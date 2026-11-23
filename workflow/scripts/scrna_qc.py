"""Shallow single-cell RNA QC on reconstructed ONT units, one FASTQ per capture."""
import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import rna_qc
import scrna_families
import scrna_regions
import scrna_resources
import scrna_species
from library_sequences import fastq_records, rc
from report_provenance import capture_analysis
from run_star_rna import align
from scrna_chemistry import PARSE, Unit, detect, expected_cdna_strand, extract, native_mgi_units


def raw_bases_beyond_r2(row, reverse):
    """Distance from the reconstructed physical Read 2 boundary to the outer ONT read end."""
    fields = ('insert_start', 'insert_end', 'raw_length')
    if any(row.get(field) in (None, '') for field in fields):
        return None
    start, end, length = (int(row[field]) for field in fields)
    if not 0 <= start <= end <= length:
        raise ValueError('Insert coordinates outside the raw read')
    return start if reverse else length - end


def load_units(raw, work):
    records = rna_qc.read_units(work / 'reads.tsv.gz', work / 'assign.tsv.gz')
    units = []
    # Dimer-only reads are also claimed by the standard reconstruction.
    with gzip.open(work / 'reads.tsv.gz', 'rt') as f:
        claimed = {row['read_id'] for row in csv.DictReader(f, delimiter='\t')}
    for name, seq, qual in rna_qc.read_fastq(work / 'inserts.fq.gz'):
        row = records.get(name)
        if row is None:
            continue
        if row['orientation'] not in ('p5_to_p7', 'p7_to_p5'):
            raise ValueError('Missing physical insert orientation')
        reverse = row['orientation'] == 'p7_to_p5'
        scope = row['scope_id'] if row['assignment_status'] == 'assigned' else 'unassigned'
        units.append(Unit(name, row['read_id'], scope, rc(seq) if reverse else seq,
                          qual[::-1] if reverse else qual, raw_bases_beyond_r2(row, reverse)))
    if len(units) != len(records):
        raise ValueError('Extracted units and adapter table disagree')
    for parent, seq, qual in fastq_records([raw]):
        if parent not in claimed:
            units.extend(native_mgi_units(parent, seq, qual))
    if len({u.name for u in units}) != len(units):
        raise ValueError('Duplicate unit IDs in FASTQ; a run requires unique read IDs')
    return units


def physical_library_summary(adapter, work, native_mgi_units):
    """Keep physical index evidence and all-unit composition at their original denominators."""
    source = adapter['index_pairing']
    pairing = {key: source[key] for key in (
        'mode', 'top_pair_fraction', 'single_library_fraction', 'informative_units',
        'minimum_informative_units', 'unassigned_units', 'ambiguous_units')}
    families = defaultdict(Counter)
    with gzip.open(work / 'reads.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['assignment_status'] == 'assigned':
                families[row['scope_id']][row['family']] += 1

    pairs, sublibraries = [], []
    for pair in source['pairs']:
        # The preparation record also holds read IDs; the report needs only aggregates.
        summary = {key: pair[key] for key in (
            'scope_id', 'name', 'i5', 'i7', 'support_units', 'informative_support', 'support_fraction')}
        observed = families[pair['scope_id']]
        if sum(observed.values()) != pair['support_units']:
            raise ValueError('Physical index counts disagree with the reconstructed unit table')
        summary['adapter_family_counts'] = dict(observed)
        pairs.append(summary)
        sublibraries.append({'scope_id': pair['scope_id'], 'library_units': pair['support_units']})
    pairing['pairs'] = pairs
    if native_mgi_units:
        # Native MGI reconstruction examines previously unclaimed reads only.
        sublibraries.append({'scope_id': 'native-unindexed', 'scope_name': 'Native unindexed',
                             'library_units': native_mgi_units})
    total = adapter['library_units'] + native_mgi_units
    assigned = sum(scope['library_units'] for scope in sublibraries)
    if assigned + pairing['ambiguous_units'] + pairing['unassigned_units'] != total:
        raise ValueError('Physical index assignments do not sum to all reconstructed units')
    return {'library_units': total, 'adapter_dimer_units': adapter['adapter_dimer_units'],
            'sub_libraries': sublibraries, 'index_pairing': pairing}


def read_well_map(path, valid):
    if path is None:
        return {}
    mapping = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            bc, well = row['round1_barcode'], row['well']
            if bc in mapping or bc not in valid or not well or not well.isascii() or not well.isalnum():
                raise ValueError('Invalid, duplicate or unknown Parse round-1 well mapping')
            mapping[bc] = well
    if set(mapping) != valid:
        raise ValueError('Provide a complete mapping for the configured Parse round-1 on-list')
    return mapping


def reference_provenance(index, gtf, assembly):
    metadata = {'contigs_sha256': scrna_resources.digest(index / 'chrNameLength.txt')}
    manifest = index / 'ontlibqc-star-index.json'
    if manifest.exists():
        index_info = json.loads(manifest.read_text())
        if index_info['reference_id'] != assembly or index_info['source_sha256']['annotation'] != scrna_resources.digest(gtf):
            raise ValueError('Assembly or GTF differs from the STAR index manifest')
        metadata.update(index_info['source_sha256'])
    return metadata


def accepted_umi(tag):
    return (set(tag.umi) <= set('ACGT') and len(set(tag.umi)) >= 2
            and tag.quality and sum(ord(q) - 33 for q in tag.quality) >= 15 * len(tag.umi))


def make_queries(units, refs, wells, cdna, single_capture=False):
    groups = defaultdict(list)
    for unit in units:
        groups[unit.scope].append(unit)
    unresolved = {'unassigned', 'native-unindexed'}
    assigned_scopes = sorted(set(groups) - unresolved)
    if single_capture and len(assigned_scopes) > 1:
        raise ValueError('--single-capture conflicts with multiple assigned physical libraries')
    capture = assigned_scopes[0] if assigned_scopes else 'single-capture'
    scopes, tagged = {}, {}
    with gzip.open(cdna, 'wt') as out:
        for scope, group in sorted(groups.items()):
            call, linkers = detect(group, refs)
            scopes[scope] = call
            counts = call['funnel'] = Counter()
            kit = call['kit']
            if kit == PARSE and linkers:
                call['linkers'] = list(linkers)
            if not kit:
                continue
            for unit in group:
                result = extract(unit, kit, refs, linkers)
                if not result:
                    continue
                tag, seq, qual = result
                counts['barcode_exact' if tag.edits == 0 else 'barcode_corrected'] += 1
                if not accepted_umi(tag):
                    continue
                counts['valid_umi'] += 1
                if len(seq) < 30:
                    continue
                counts['cdna'] += 1
                barcode = tag.barcode
                if kit == PARSE and wells:
                    first, second, rt = barcode.split('-')
                    barcode = '-'.join((first, second, wells[rt]))
                # The FASTQ written here is always in physical Read 1 orientation.
                tagged[unit.name] = {
                    'scope_id': scope, 'read_id': unit.parent, 'barcode': barcode,
                    'umi': tag.umi, 'kit': kit, 'insert_length': len(seq), 'orientation': 'p5_to_p7',
                    'assignment_status': 'assigned', 'countable': scope not in unresolved or single_capture,
                    'molecule_scope': capture if single_capture else scope,
                }
                out.write(f'@{unit.name}\n{seq}\n+\n{qual}\n')
    if single_capture and len({d['kit'] for d in scopes.values() if d['kit']}) > 1:
        raise ValueError('--single-capture conflicts with different barcode layouts')
    return scopes, tagged


def molecule_summary(rows, tagged, scopes):
    molecules = Counter()
    exonic_keys = set()
    mappings, regions, strand, rna_evidence = Counter(), Counter(), Counter(), Counter()
    for row in rows:
        unit = tagged[row['unit_id']]
        counts = scopes[unit['scope_id']]['funnel']
        status = row['mapping_status']
        mappings[status] += 1
        if status != 'unique':
            continue
        counts['unique'] += 1
        regions.update(row['region_bases'])
        if row.get('mitochondrial'):
            rna_evidence['mitochondrial'] += 1
        if row['family'] in ('rRNA', 'mt-rRNA'):
            rna_evidence['annotated_rrna'] += 1
        if row['junctions']:
            rna_evidence['spliced'] += 1
        if row['known_junctions']:
            rna_evidence['annotated_splice'] += 1
        if not row['eligible_gene']:
            continue
        strand[row['strand_evidence']] += 1
        if row['strand_evidence'] != expected_cdna_strand(unit['kit']):
            continue
        counts['gene'] += 1
        if not unit['countable']:
            continue
        exon = any(row['region_bases'].get(k, 0) for k in ('CDS', 'five_prime_UTR', 'three_prime_UTR', 'other_exon'))
        if not exon and not row['region_bases'].get('intron', 0):
            continue
        key = (unit['molecule_scope'], unit['barcode'], row['gene_id'], unit['umi'])
        molecules[key] += 1
        if exon:
            exonic_keys.add(key)
    per_barcode, gene_counts, matrix, exonic_matrix = defaultdict(Counter), Counter(), Counter(), Counter()
    for (scope, barcode, gene, umi), count in molecules.items():
        per_barcode[scope, barcode][gene] += 1
        gene_counts[gene] += count
        matrix[scope, barcode, gene] += 1
        if (scope, barcode, gene, umi) in exonic_keys:
            exonic_matrix[scope, barcode, gene] += 1
    depth = sum(molecules.values())
    draws = sorted({round(depth * n / 10) for n in range(11)})
    molecule_counts = list(molecules.values())
    genes = list(gene_counts.values())
    if not tagged:
        assessment = 'no_qualified_cdna'
    elif not any(u['countable'] for u in tagged.values()):
        assessment = 'no_assigned_capture'
    else:
        assessment = 'evaluated'
    summary = {
        'molecule_assessment': assessment,
        'mapping': dict(mappings), 'aligned_region_bases': dict(regions), 'gene_strand': dict(strand),
        'rna_evidence': dict(rna_evidence),
        'assigned_units': depth, 'molecules': len(molecules), 'genes': len(gene_counts),
        'exonic_molecules': len(exonic_keys), 'intronic_only_molecules': len(molecules) - len(exonic_keys),
        'barcodes': len(per_barcode), 'singleton_barcodes': sum(sum(v.values()) == 1 for v in per_barcode.values()),
        'barcodes_ge_2_genes': sum(len(v) >= 2 for v in per_barcode.values()),
        'ranks': sorted((sum(v.values()) for v in per_barcode.values()), reverse=True),
        'rarefaction': [{'units': d, 'molecules': round(rna_qc.expected_detected(molecule_counts, d), 2),
                         'genes': round(rna_qc.expected_detected(genes, d), 2)} for d in draws],
    }
    entries = ({'scope': s, 'barcode': b, 'gene_id': g, 'umis': n, 'exonic_umis': exonic_matrix[s, b, g]}
               for (s, b, g), n in sorted(matrix.items()))
    return summary, entries


def run(args):
    out = Path(args.out).resolve()
    if not 1 <= args.threads <= 64 or not 0 <= args.mapq <= 60:
        raise ValueError('Threads must be 1–64 and MAPQ must be 0–60')
    if out.exists():
        raise FileExistsError(f'Choose a new output directory: {out}')
    out.parent.mkdir(parents=True, exist_ok=True)
    # Publish the whole result in one rename. Private per-unit data stays in work.
    with tempfile.TemporaryDirectory(dir=out.parent, prefix='.scrna-') as temp:
        stage = Path(temp)
        work, result = stage / 'work', stage / 'result'
        work.mkdir()
        result.mkdir(mode=0o700)
        script = Path(__file__).with_name('ont_qc.py')
        files = {'inserts': 'inserts.fq.gz', 'legs': 'legs.fq.gz', 'leg1': 'leg1.fq.gz', 'leg2': 'leg2.fq.gz',
                 'read-table': 'reads.tsv.gz', 'index-assignments': 'assign.tsv.gz', 'summary': 'adapter.json',
                 'raw-lengths': 'raw.tsv', 'insert-lengths': 'insert.tsv'}
        command = [sys.executable, str(script), 'prepare', '--assay', 'rnaseq',
                   '--sample', args.sample, '--fastq', str(args.fastq)]
        for key, filename in files.items():
            command.extend(['--' + key, str(work / filename)])
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
        adapter = json.loads((work / 'adapter.json').read_text())
        units = load_units(args.fastq, work)
        observed_ends = {bc for u in units for bc in (u.seq[:16], rc(u.seq[-16:]))}
        refs = scrna_resources.load(args.barcodes, observed_ends)
        wells = read_well_map(args.parse_well_map, refs['parse'][2])
        cdna = work / 'cdna.fq.gz'
        scopes, tagged = make_queries(units, refs, wells, cdna, args.single_capture)
        rows, annotation_meta, star_version, reference, families = [], None, None, None, None
        if tagged:
            reference = reference_provenance(args.star_index, args.gtf, args.assembly)
            bam = work / 'cdna.bam'
            star_log = work / 'star.json'
            align(argparse.Namespace(fastq=cdna, index=args.star_index, output=bam, summary=star_log,
                                     star=args.star, samtools=args.samtools, threads=args.threads))
            star_version = json.loads(star_log.read_text())['version']
            alignments, chroms = rna_qc.read_alignments(bam, args.samtools)
            annotation = rna_qc.read_annotation(args.gtf, transcript_regions=not args.rna_family_gtf)
            if not set(chroms).intersection(annotation.index):
                raise ValueError('No shared chromosomes between STAR index and GTF')
            annotation_meta = {k: annotation.metadata[k] for k in ('genes', 'sha256', 'gene_families')}
            rows, _ = rna_qc.analyze_units(tagged, alignments, annotation, args.mapq, collect_profiles=False)
            if args.rna_family_gtf:
                annotation = None  # Release the STAR annotation before reading the detailed family GTF.
            families = scrna_families.analyze(args, work, rows, alignments, chroms, tagged, annotation)
            del annotation
            for row in rows:
                if row['mapping_status'] == 'unique':
                    row['mitochondrial'] = any(r['chrom'] in ('chrM', 'MT', 'chrMT', 'M')
                                               for r in alignments[row['unit_id']])
        metrics, entries = molecule_summary(rows, tagged, scopes)
        metrics['species_composition'] = scrna_species.analyze(args, work, tagged)
        if families:
            metrics['transcript_regions'] = families.pop('transcript_regions')
            metrics['transcript_metaplot'] = families.pop('transcript_metaplot')
            metrics['rna_families'] = families
        summary = {
            'schema_version': 1, 'sample': args.sample, 'assembly': args.assembly,
            'input_sha256': scrna_resources.digest(args.fastq),
            'analysis_provenance': capture_analysis(),
            'raw_reads': adapter['raw']['count'], 'non_dimer_units': len(units),
            'reconstructed_parents': len({u.parent for u in units}), 'raw_median_bp': adapter['raw']['median'],
            'native_mgi_units': sum(u.scope == 'native-unindexed' for u in units),
            'adapter_family_counts': adapter.get('adapter_families', {}),
            'physical_indexes': {pair['scope_id']: {key: pair.get(key) for key in ('i5', 'i7')}
                                 for pair in adapter.get('index_pairing', {}).get('pairs', []) if pair.get('scope_id')},
            'scopes': scopes, 'metrics': metrics, 'mapq': args.mapq, 'star_version': star_version,
            'annotation': annotation_meta, 'reference': reference, 'single_capture': args.single_capture,
            'parse_well_map_sha256': scrna_resources.digest(Path(args.parse_well_map)) if wells else None,
            'barcode_resources_sha256': scrna_resources.HASHES,
        }
        summary['physical_library'] = physical_library_summary(adapter, work, summary['native_mgi_units'])
        rna_qc.write_tsv(result / 'barcode_gene_counts.tsv.gz', entries,
                         ['scope', 'barcode', 'gene_id', 'umis', 'exonic_umis'])
        (result / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        from scrna_metaplot_plot import figures as transcript_figures
        from scrna_report import family_plot, render
        figure = family_plot(families['cohorts']) if families else None
        if figure:
            (result / 'rna_families.svg').write_bytes(figure)
        region_figure = scrna_regions.figure(metrics['transcript_regions']) if families else None
        if region_figure:
            (result / 'transcript_regions.svg').write_bytes(region_figure)
        metaplots = transcript_figures(metrics.get('transcript_metaplot') or {})
        if 'all' in metaplots:
            (result / 'transcript_metaplot.svg').write_bytes(metaplots['all'])
        (result / 'report.html').write_text(render(summary, family_figure=figure, region_figure=region_figure,
                                                  metaplot_figures=metaplots))
        if args.keep_work:
            shutil.move(work, result / 'work')
        os.rename(result, out)


def add_commands(sub):
    description = 'Shallow single-cell RNA library QC (Beta)'
    p = sub.add_parser('scrna', help=description, description=description)
    p.add_argument('--fastq', required=True, type=Path, help='One ONT FASTQ or FASTQ.gz')
    p.add_argument('--sample', required=True)
    p.add_argument('--assembly', required=True)
    p.add_argument('--barcodes', required=True, type=Path)
    p.add_argument('--star-index', required=True, type=Path)
    p.add_argument('--gtf', required=True, type=Path)
    p.add_argument('--rna-family-gtf', type=Path, help='Optional biotype-rich GTF on the same assembly for the RNA family plot')
    p.add_argument('--rna-source-panel', type=Path, help='RNA-seq rRNA reference directory with panel.fa and panel.json')
    p.add_argument('--repeatmasker', type=Path, help='UCSC rmsk.txt or rmsk.txt.gz for the STAR assembly')
    p.add_argument('--minimap2', default='minimap2', help='Align cDNA to the independent rRNA panel')
    p.add_argument('--species-panel', type=Path,
                   help='Competitive Bowtie 2 panel JSON (default: species-panel.json beside the STAR index directory)')
    p.add_argument('--bowtie2', default='bowtie2', help='Align qualified cDNA for the species screen')
    p.add_argument('--star', default='STARlong')
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--mapq', type=int, default=20)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--parse-well-map', type=Path, help='Verified round1_barcode/well TSV for this kit')
    p.add_argument('--single-capture', action='store_true',
                   help='Pool matching barcodes across index scopes only when the FASTQ is known to contain one capture')
    p.add_argument('--keep-work', action='store_true', help='Retain intermediate FASTQ, tables and BAM locally')
    p.set_defaults(func=run)
    p = sub.add_parser('prepare-scrna-barcodes', help='Fetch checksum-pinned public barcode definitions')
    p.add_argument('--out', required=True, type=Path)
    p.set_defaults(func=lambda args: scrna_resources.prepare(args.out))
