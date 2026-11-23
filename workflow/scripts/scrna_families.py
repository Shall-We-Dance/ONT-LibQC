"""Classify qualified cDNA once per fragment, with optional matched RepeatMasker evidence."""
import argparse
import bisect
import csv
from collections import Counter, defaultdict

import rna_qc
import rna_sources
import scrna_metaplot
import scrna_regions
from scrna_resources import digest


def interval_index(intervals):
    merged = rna_qc.merge_intervals(intervals)
    return [a for a, _ in merged], [b for _, b in merged]


def intersections(index, start, end):
    starts, ends = index
    lo = bisect.bisect_right(ends, start)
    hi = bisect.bisect_left(starts, end)
    for i in range(lo, hi):
        yield max(start, starts[i]), min(end, ends[i])


def repeat_loci(path, alignments, reference_chroms):
    """For each class, ≥half the M/= /X bases at a locus must fall in its union."""
    query = defaultdict(list)
    for records in alignments.values():
        for r in records:
            if not r['supplementary']:
                query[r['chrom']].extend(r['blocks'])
    windows = {chrom: interval_index(blocks) for chrom, blocks in query.items()}
    hits = defaultdict(list)
    shared_chrom = False
    with rna_qc.open_text(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 13:
                raise ValueError('Expected a UCSC rmsk table with its bin column')
            chrom = fields[5]
            shared_chrom |= chrom in reference_chroms
            if chrom in windows:
                group = 'rRNA' if fields[11] == 'rRNA' else 'other'
                hits[chrom, group].extend(intersections(windows[chrom], int(fields[6]), int(fields[7])))
    if not shared_chrom:
        raise ValueError('No shared chromosomes between RepeatMasker and STAR')
    indexed = {key: interval_index(values) for key, values in hits.items()}
    result = {}
    for name, records in alignments.items():
        loci = [r for r in records if not r['supplementary']]
        counts = Counter()
        for r in loci:
            size = rna_qc.interval_length(r['blocks'])
            for group in ('rRNA', 'other'):
                key = r['chrom'], group
                if not size or key not in indexed:
                    continue
                overlap = sum(y - x for a, b in r['blocks'] for x, y in intersections(indexed[key], a, b))
                if 2 * overlap >= size:
                    counts[group] += 1
        total = max([len(loci)] + [r['nh'] for r in loci])
        result[name] = (total, counts['rRNA'], counts['other'])
    return result


def family_counts(rows, source_rows, references, repeats, output):
    source = {row['unit_id']: row for row in source_rows}
    if len(source) != len(rows) or set(source) != {r['unit_id'] for r in rows}:
        raise ValueError('RNA panel and genome classifications have different input fragments')
    counts = {'all': Counter(), 'multimapped': Counter()}
    audit = Counter()
    fieldnames = ['unit_id', 'scope_id', 'mapping_status', 'family', 'genome_family', 'evidence',
                  'panel_reference', 'matched_fraction', 'identity', 'candidate_loci', 'rRNA_repeat_loci', 'other_repeat_loci']
    with rna_qc.open_text(output, 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            name, mapping, family = row['unit_id'], row['mapping_status'], row['family']
            hit = source[name]
            target = hit['reference_name']
            loci, rrna, other = repeats.get(name, (0, 0, 0))
            evidence = 'genome'
            if hit['status'] == 'assigned' and hit['source_class'] == 'rRNA':
                evidence = 'sequence_rRNA'
                family = 'mt-rRNA' if (family == 'mt-rRNA' or references[target].get('rna_family') == 'mt-rRNA') else 'rRNA'
            elif family in ('rRNA', 'mt-rRNA'):
                evidence = 'annotated_rRNA'
            elif loci and rrna == loci and mapping in ('unique', 'multimapped', 'low_mapq'):
                evidence, family = 'repeat_rRNA_all_loci', 'rRNA'
            counts['all'][family] += 1
            if mapping == 'multimapped':
                counts['multimapped'][family] += 1
                reason = evidence
                if evidence == 'genome':
                    if rrna:
                        reason = 'rRNA_some_loci'
                    elif loci and other == loci:
                        reason = 'other_repeat_all_loci'
                    elif other:
                        reason = 'other_repeat_some_loci'
                    else:
                        reason = 'unresolved'
                audit[reason] += 1
            writer.writerow(dict(unit_id=name, scope_id=row['scope_id'], mapping_status=mapping, family=family,
                genome_family=row['family'], evidence=evidence, panel_reference=target,
                matched_fraction=hit['matched_fraction'], identity=hit['identity'],
                candidate_loci=loci if name in repeats else '', rRNA_repeat_loci=rrna if name in repeats else '',
                other_repeat_loci=other if name in repeats else ''))
    return {k: {'units': sum(v.values()), 'counts': dict(v)} for k, v in counts.items()}, dict(audit)


def analyze(args, work, rows, alignments, chroms, tagged, annotation):
    chroms = set(chroms)
    if args.rna_family_gtf:
        annotation = rna_qc.read_annotation(args.rna_family_gtf, transcript_regions=True)
        if not set(annotation.index) <= chroms:
            raise ValueError('RNA family GTF contigs must use the STAR reference names')
        rows, _ = rna_qc.analyze_units(tagged, alignments, annotation, args.mapq, collect_profiles=False)
    regions = scrna_regions.summarize(rows, alignments, annotation)
    panel = args.rna_source_panel or rna_sources.RESOURCE
    if args.rna_source_panel and not all((panel / file).is_file() for file in ('panel.fa', 'panel.json')):
        raise ValueError('RNA source panel requires panel.fa and panel.json')
    assignments = work / 'rna_sources.tsv.gz'
    source_info = rna_sources.run(argparse.Namespace(input=work / 'cdna.fq.gz', reference_fasta=panel / 'panel.fa',
        reference_manifest=panel / 'panel.json', output=assignments, summary=work / 'rna_sources.json',
        minimap2=args.minimap2, threads=args.threads))
    references = rna_sources.load_panel(panel / 'panel.fa', panel / 'panel.json')['references'] if source_info['status'] == 'complete' else {}
    repeats = repeat_loci(args.repeatmasker, alignments, chroms) if args.repeatmasker else {}
    family_path = work / 'rna_families.tsv.gz'
    with rna_qc.open_text(assignments) as handle:
        cohorts, audit = family_counts(rows, csv.DictReader(handle, delimiter='\t'), references, repeats, family_path)
    with rna_qc.open_text(family_path) as handle:
        metaplot = scrna_metaplot.summarize(rows, alignments, annotation, tagged, csv.DictReader(handle, delimiter='\t'))
    return {'cohorts': cohorts, 'multimapping_evidence': audit, 'panel': source_info, 'transcript_regions': regions,
            'transcript_metaplot': metaplot, 'annotation_sha256': annotation.metadata['sha256'],
            'repeatmasker_sha256': digest(args.repeatmasker) if args.repeatmasker else None}
