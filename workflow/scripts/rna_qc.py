#!/usr/bin/env python3
"""RNA/CLIP QC for adapter-bounded ONT inserts from completed NGS libraries.

Input BAM and FASTQ must retain original ONT read orientation. Physical strand
is restored from the adapter table. All intervals are zero-based, half-open.
Analysis uses the standard library and samtools; scientific plots use matplotlib. SAM input works directly.
"""
from __future__ import annotations
import argparse
import bisect
import collections
import csv
import gzip
import hashlib
import html
import itertools
import json
import math
import random
import re
import statistics
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

BIOTYPES = ('rRNA', 'mt-rRNA', 'mt-mRNA', 'tRNA', 'snRNA', 'snoRNA', 'protein_coding', 'lncRNA', 'other', 'unknown', 'ambiguous_family')
REGIONS = ('CDS', 'five_prime_UTR', 'three_prime_UTR', 'other_exon', 'intron', 'genic_unannotated', 'intergenic')
LIMITATIONS = [
    'Shallow ONT QC of completed NGS libraries: counts describe observed adapter-bounded units, not original RNA molecules.',
    'Annotation completeness limits RNA-family classification. Missing rRNA/repeat annotations and mapper-suppressed alignments cannot be recovered here.',
    'Multimapping family evidence is retained separately; locus profiles use primary MAPQ-pass nonchimeric unambiguous assignments.',
    'Physical strand is the adapter-defined P5-to-P7 strand, interpreted relative to nonconflicting gene annotation.',
    'Selection chemistry, kit, cell type, RNA integrity and DNA contamination cannot be uniquely determined from these observations.',
    'Coordinate collisions and observed UMI groups are not calibrated PCR-duplicate rates or library-complexity estimates; ONT errors and UMI collisions remain possible.',
    'HOMER CLIP motifs describe enrichment against dinucleotide-preserving shuffled inserts, not IP/input specificity; identical sequences are counted once for motif discovery.',
    'RNA/repeat enrichment, intronic signal and uneven coverage can reflect target biology. No universal pass/fail thresholds are applied.',
]


def open_text(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode.replace('t', ''))


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def merge_intervals(intervals):
    output = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if output and start <= output[-1][1]:
            output[-1] = (output[-1][0], max(output[-1][1], end))
        else:
            output.append((start, end))
    return output


def intersect_intervals(left, right):
    left, right = merge_intervals(left), merge_intervals(right)
    output, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        start, end = max(left[i][0], right[j][0]), min(left[i][1], right[j][1])
        if start < end:
            output.append((start, end))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return output


def subtract_intervals(left, right):
    right, output = merge_intervals(right), []
    for start, end in merge_intervals(left):
        cursor = start
        for rstart, rend in right:
            if rend <= cursor:
                continue
            if rstart >= end:
                break
            if rstart > cursor:
                output.append((cursor, min(rstart, end)))
            cursor = max(cursor, rend)
            if cursor >= end:
                break
        if cursor < end:
            output.append((cursor, end))
    return output


def interval_length(intervals):
    return sum(end - start for start, end in intervals)


def overlap_length(left, right):
    return interval_length(intersect_intervals(left, right))


def attributes(text):
    result = {}
    for item in text.strip().split(';'):
        match = re.match(r'\s*([^\s=]+)\s+(?:"([^"]*)"|(\S+))', item)
        if match:
            result[match.group(1)] = match.group(2) if match.group(2) is not None else match.group(3)
    return result


def classify_biotype(biotype, name, chrom, mitochondrial):
    value, label = (biotype or '').lower(), (name or '').upper()
    mito = chrom in mitochondrial or value.startswith('mt_')
    if 'pseudogene' in value:
        return 'other'
    nuclear_rrna_symbol = bool(re.fullmatch(r'RNA(?:18S|28S|45S|5S|5[-_]8S)[0-9]*', label))
    mitochondrial_rrna_symbol = mito and label in ('RNR1', 'RNR2', 'MT-RNR1', 'MT-RNR2', 'MTRNR1', 'MTRNR2', 'RRNS', 'RRNL')
    if 'rrna' in value or (not value and (nuclear_rrna_symbol or mitochondrial_rrna_symbol)):
        return 'mt-rRNA' if mito else 'rRNA'
    if 'trna' in value or (not value and mito and re.fullmatch(r'(?:MT-)?TRN[A-Z][12]?', label)):
        return 'tRNA'
    sno_symbol = bool(re.fullmatch(r'(?:SNOR[AD][0-9]+[A-OQ-Z]?(?:-[0-9]+)?|SCARNA[0-9]+)', label))
    if 'snorna' in value or 'scarna' in value or (not value and sno_symbol):
        return 'snoRNA'
    if 'snrna' in value or (not value and re.fullmatch(r'RNU[0-9]+(?:ATAC)?(?:-[0-9]+)?', label)):
        return 'snRNA'
    if value in ('protein_coding', 'mrna'):
        return 'mt-mRNA' if mito else 'protein_coding'
    if 'lncrna' in value or value in ('lincrna', 'antisense', 'sense_intronic', 'sense_overlapping', 'processed_transcript'):
        return 'lncRNA'
    return 'other' if value else 'unknown'


def infer_utrs(exons, cds, strand):
    if not cds or not exons or strand not in ('+', '-'):
        return [], []
    low, high = min(s for s, _ in cds), max(e for _, e in cds)
    left = intersect_intervals(exons, [(0, low)])
    right = intersect_intervals(exons, [(high, max(e for _, e in exons))])
    return (right, left) if strand == '-' else (left, right)


@dataclass
class Gene:
    gene_id: str
    chrom: str
    strand: str
    name: str = ''
    biotype: str = ''
    family: str = 'unknown'
    family_evidence: str = 'unclassified'
    bounds: list = field(default_factory=list)
    exons: list = field(default_factory=list)
    cds: list = field(default_factory=list)
    utr5: list = field(default_factory=list)
    utr3: list = field(default_factory=list)
    transcripts: dict = field(default_factory=dict)
    representative: list = field(default_factory=list)
    segments: dict = field(default_factory=dict)

    def finish(self, mitochondrial):
        self.exons, self.cds = merge_intervals(self.exons), merge_intervals(self.cds)
        effective_type = self.biotype or ('protein_coding' if self.cds else '')
        self.family = classify_biotype(effective_type, self.name, self.chrom, mitochondrial)
        self.family_evidence = 'explicit_biotype' if self.biotype else 'CDS_annotation' if self.cds else 'gene_symbol_heuristic' if self.family != 'unknown' else 'unclassified'
        if not self.exons and self.family in ('rRNA', 'mt-rRNA', 'tRNA', 'snRNA', 'snoRNA'):
            self.exons = merge_intervals(self.bounds)
        self.bounds = [(min(s for s, _ in self.bounds), max(e for _, e in self.bounds))]
        candidates = []
        for transcript, parts in self.transcripts.items():
            exons = merge_intervals(parts['exon'])
            if exons:
                candidates.append((interval_length(exons), transcript, exons, parts))
        if candidates:
            _, _, self.representative, parts = sorted(candidates, key=lambda x: (-x[0], x[1]))[0]
            cds = merge_intervals(parts['CDS'])
        else:
            self.representative, cds = self.exons, self.cds
        five, three = infer_utrs(self.representative, cds, self.strand)
        self.segments = {'five_prime_UTR': five, 'CDS': cds, 'three_prime_UTR': three}
        five, three = infer_utrs(self.exons, self.cds, self.strand)
        self.utr5, self.utr3 = merge_intervals(self.utr5 or five), merge_intervals(self.utr3 or three)


class Annotation:
    def __init__(self, genes, known_junctions, metadata):
        self.genes, self.known_junctions, self.metadata = genes, known_junctions, metadata
        by_chrom = collections.defaultdict(list)
        for gene in genes.values():
            start, end = gene.bounds[0]
            by_chrom[gene.chrom].append((start, end, gene.gene_id))
        self.index = {}
        for chrom, rows in by_chrom.items():
            rows.sort()
            maximum, prefix = 0, []
            for _, end, _ in rows:
                maximum = max(maximum, end)
                prefix.append(maximum)
            self.index[chrom] = rows, [r[0] for r in rows], prefix

    def query(self, chrom, blocks):
        if chrom not in self.index:
            return []
        rows, starts, prefix = self.index[chrom]
        found = set()
        for start, end in blocks:
            index = bisect.bisect_left(starts, end) - 1
            while index >= 0 and prefix[index] > start:
                if rows[index][1] > start:
                    found.add(rows[index][2])
                index -= 1
        return [self.genes[key] for key in sorted(found)]


def read_annotation(path, mitochondrial=('chrM', 'MT')):
    genes, malformed, digest = {}, 0, hashlib.sha256()
    with open_text(path) as handle:
        for line in handle:
            digest.update(line.encode())
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                malformed += 1
                continue
            chrom, _, feature, start, end, _, strand, _, raw = fields
            attrs = attributes(raw)
            gene_id = attrs.get('gene_id')
            if not gene_id:
                malformed += 1
                continue
            try:
                interval = int(start) - 1, int(end)
                if interval[0] < 0 or interval[1] <= interval[0]:
                    raise ValueError()
            except ValueError:
                malformed += 1
                continue
            key = gene_id + '@' + chrom + ':' + strand
            gene = genes.setdefault(key, Gene(key, chrom, strand, attrs.get('gene_name', gene_id)))
            gene.biotype = attrs.get('gene_type', attrs.get('gene_biotype', gene.biotype))
            gene.bounds.append(interval)
            destination = {'exon': gene.exons, 'cds': gene.cds, 'five_prime_utr': gene.utr5, 'three_prime_utr': gene.utr3}.get(feature.lower())
            if destination is not None:
                destination.append(interval)
            transcript = attrs.get('transcript_id')
            if transcript and feature in ('exon', 'CDS'):
                gene.transcripts.setdefault(transcript, {'exon': [], 'CDS': []})[feature].append(interval)
    known = set()
    for gene in genes.values():
        for parts in gene.transcripts.values():
            exons = merge_intervals(parts['exon'])
            for left, right in zip(exons, exons[1:]):
                if left[1] < right[0]:
                    known.add((gene.chrom, left[1], right[0], gene.strand))
        gene.finish(set(mitochondrial))
    if not genes:
        raise ValueError('GTF contains no valid gene_id annotations')
    families = collections.Counter(g.family for g in genes.values())
    return Annotation(genes, known, {'sha256': digest.hexdigest(), 'genes': len(genes), 'known_junctions': len(known),
        'gene_families': dict(families), 'classification_evidence_counts': dict(collections.Counter(g.family_evidence for g in genes.values())),
        'genes_without_explicit_biotype': sum(not g.biotype for g in genes.values()), 'malformed_or_unassigned_lines': malformed,
        'missing_families': [f for f in ('rRNA', 'mt-rRNA', 'tRNA', 'snRNA', 'snoRNA') if not families[f]],
        'annotation_chromosomes': sorted({g.chrom for g in genes.values()})})


def parse_cigar(cigar, start):
    if cigar == '*':
        return [], [], start, 0
    tokens = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    if ''.join(n + op for n, op in tokens) != cigar:
        raise ValueError('Malformed CIGAR: ' + cigar)
    position, length, blocks, junctions = start, 0, [], []
    for number, operation in tokens:
        size = int(number)
        if operation in 'M=X':
            blocks.append((position, position + size))
            position += size
            length += size
        elif operation == 'N':
            junctions.append((position, position + size))
            position += size
        elif operation == 'D':
            position += size
        elif operation in 'IS':
            length += size
    return merge_intervals(blocks), junctions, position, length


def sam_lines(path, samtools='samtools'):
    if str(path).endswith(('.sam', '.sam.gz')):
        with open_text(path) as handle:
            yield from handle
        return
    # stderr is a temporary file to avoid pipe deadlocks on large diagnostics.
    import tempfile
    with tempfile.TemporaryFile(mode='w+t') as error:
        with subprocess.Popen([samtools, 'view', '-h', str(path)], stdout=subprocess.PIPE, stderr=error, text=True) as process:
            try:
                yield from process.stdout
            finally:
                process.stdout.close()
            status = process.wait()
            if status:
                error.seek(0)
                raise RuntimeError('samtools view failed: ' + error.read()[-4000:])


def read_alignments(path, samtools='samtools'):
    alignments, references = collections.defaultdict(list), set()
    for line in sam_lines(path, samtools):
        if line.startswith('@'):
            if line.startswith('@SQ'):
                references.update(x[3:] for x in line.rstrip().split('\t') if x.startswith('SN:'))
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 11:
            raise ValueError('Invalid SAM record')
        name, flag = fields[0], int(fields[1])
        if flag & 4:
            continue
        start = int(fields[3]) - 1
        blocks, junctions, end, length = parse_cigar(fields[5], start)
        tags = {p[0]: p[2] for value in fields[11:] if len(p := value.split(':', 2)) == 3}
        alignments[name].append({'chrom': fields[2], 'start': start, 'end': end, 'reverse': bool(flag & 16),
            'secondary': bool(flag & 256), 'supplementary': bool(flag & 2048), 'mapq': int(fields[4]),
            'blocks': blocks, 'junctions': junctions, 'query_length': length, 'cigar': fields[5],
            'nh': int(tags.get('NH', 1)), 'sa': bool(tags.get('SA'))})
    return alignments, sorted(references)


def read_units(path, index_assignments=None):
    with open_text(path) as handle:
        result, reader = {}, csv.DictReader(handle, delimiter='\t')
        if 'unit_id' not in (reader.fieldnames or []):
            raise ValueError('Adapter table must contain unit_id')
        for row in reader:
            name = row['unit_id']
            if name in result:
                raise ValueError('Duplicate adapter unit_id: ' + name)
            if str(row.get('is_adapter_dimer', '0')).lower() in ('1', 'true'):
                continue
            result[name] = row
    if index_assignments:
        with open_text(index_assignments) as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row['unit_id'] in result:
                    result[row['unit_id']].update(row)
    return result


def read_fastq(path):
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence, plus, quality = handle.readline().strip(), handle.readline(), handle.readline().strip()
            if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise ValueError('Malformed extracted FASTQ record')
            yield header[1:].split()[0], sequence.upper(), quality


def region_bases(blocks, genes):
    candidates = {key: [] for key in REGIONS[:-1]}
    for gene in genes:
        candidates['CDS'].extend(gene.cds)
        candidates['five_prime_UTR'].extend(gene.utr5)
        candidates['three_prime_UTR'].extend(gene.utr3)
        candidates['other_exon'].extend(gene.exons)
        candidates['intron' if gene.exons else 'genic_unannotated'].extend(gene.bounds)
    remaining, result = blocks, {}
    for region in REGIONS[:-1]:
        hit = intersect_intervals(remaining, candidates[region])
        result[region] = interval_length(hit)
        remaining = subtract_intervals(remaining, hit)
    result['intergenic'] = interval_length(remaining)
    return result


def classify_alignment(record, annotation):
    genes = annotation.query(record['chrom'], record['blocks'])
    exonic = [g for g in genes if overlap_length(record['blocks'], g.exons)]
    candidates = exonic or genes
    return {'genes': candidates, 'gene': candidates[0] if len(candidates) == 1 else None,
            'families': {g.family for g in candidates} or {'unknown'}, 'family_evidence': {g.family_evidence for g in candidates} or {'no_annotation'}, 'regions': region_bases(record['blocks'], genes)}


def physical_strand(record, unit):
    orientation = unit.get('orientation')
    if orientation not in ('p5_to_p7', 'p7_to_p5'):
        return None
    return '-' if record['reverse'] ^ (orientation == 'p7_to_p5') else '+'


def profile_vector(blocks, exons, strand, bins=100):
    length, vector, offset = interval_length(exons), [0.0] * bins, 0
    if not length:
        return vector
    for start, end in exons:
        for lo, hi in intersect_intervals(blocks, [(start, end)]):
            left, right = offset + lo - start, offset + hi - start
            if strand == '-':
                left, right = length - right, length - left
            first, last = int(left * bins / length), min(bins - 1, math.ceil(right * bins / length) - 1)
            for index in range(first, last + 1):
                vector[index] += max(0, min(right, (index + 1) * length / bins) - max(left, index * length / bins))
        offset += end - start
    return vector


def wilson_interval(success, total):
    if not total:
        return [None, None]
    z, p = 1.959963984540054, success / total
    scale = 1 + z * z / total
    center = (p + z * z / (2 * total)) / scale
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / scale
    return [max(0., center - half), min(1., center + half)]


def expected_detected(counts, draw):
    total, result = sum(counts), 0.
    if not total or not draw:
        return 0.
    draw = min(draw, total)
    for count in counts:
        if total - count < draw:
            result += 1.
            continue
        if min(draw, count) <= 1000:
            if count <= draw:
                logp = sum(math.log1p(-draw / (total - i)) for i in range(count))
            else:
                logp = sum(math.log1p(-count / (total - i)) for i in range(draw))
        else:
            logp = math.lgamma(total - count + 1) - math.lgamma(total - count - draw + 1) - math.lgamma(total + 1) + math.lgamma(total - draw + 1)
        result += -math.expm1(min(0., logp))
    return result


def dinucleotide_shuffle(sequence, rng):
    if len(sequence) < 3:
        return sequence
    edges = collections.defaultdict(list)
    for left, right in zip(sequence, sequence[1:]):
        edges[left].append(right)
    for choices in edges.values():
        rng.shuffle(choices)
    stack, trail = [sequence[0]], []
    while stack:
        if edges[stack[-1]]:
            stack.append(edges[stack[-1]].pop())
        else:
            trail.append(stack.pop())
    return ''.join(reversed(trail))


def present_kmers(sequence, k):
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1) if not re.search('[^ACGT]', sequence[i:i + k])}


def motif_aligned_sequence(sequence, record, gene_strand):
    """Orient aligned query sequence to RNA; omit terminal soft clips and hard-clipped records."""
    operations = re.findall(r'(\d+)([MIDNSHP=X])', record['cigar'])
    if any(op == 'H' for _, op in operations) or sum(int(n) for n, op in operations if op in 'MIS=X') != len(sequence):
        return ''
    oriented = reverse_complement(sequence) if record['reverse'] else sequence
    left = int(operations[0][0]) if operations and operations[0][1] == 'S' else 0
    right = int(operations[-1][0]) if operations and operations[-1][1] == 'S' else 0
    aligned = oriented[left:len(oriented)-right]
    return reverse_complement(aligned) if gene_strand == '-' else aligned


def discover_motifs(sequences, minimum=100, seed=1729, shuffles=8, maximum=2000):
    from rna_motifs import discover_homer
    return discover_homer(sequences, dinucleotide_shuffle, minimum, seed, shuffles, maximum)


def analyze_units(units, alignments, annotation, mapq=10):
    rows, profiles = [], {}
    for name, unit in units.items():
        records = alignments.get(name, [])
        primary = [r for r in records if not r['secondary'] and not r['supplementary']]
        record = max(primary, key=lambda r: r['mapq']) if primary else None
        scope = unit.get('scope_id', '') if unit.get('assignment_status') == 'assigned' else ''
        row = {'unit_id': name, 'read_id': unit.get('read_id', name), 'scope_id': scope,
            'insert_length': int(unit.get('insert_length') or 0), 'orientation': unit.get('orientation', ''),
            'mapping_status': 'unmapped', 'mapq': None, 'gene_id': '', 'gene_name': '', 'family': 'unknown',
            'families': [], 'gene_strand': '', 'physical_strand': '', 'strand_evidence': 'unknown',
            'region_bases': {key: 0 for key in REGIONS}, 'junctions': [], 'known_junctions': 0,
            'matched_reference_bases': 0, 'genomic_span': 0, 'umi_status': unit.get('umi_status', 'not_configured'),
            'umi_sequence': unit.get('umi_sequence', ''), 'coordinate_key': '', 'query_reverse': None,
            'eligible_gene': False, 'mapped': bool(record)}
        families, family_evidence = set(), set()
        for alternative in records:
            if not alternative['supplementary']:
                classified = classify_alignment(alternative, annotation)
                families.update(classified['families'])
                family_evidence.update(classified['family_evidence'])
        row['families'] = sorted(families or {'unknown'})
        row['family_evidence'] = sorted(family_evidence or {'no_alignment'})
        row['family'] = next(iter(families)) if len(families) == 1 else 'ambiguous_family' if families else 'unknown'
        if record:
            detail = classify_alignment(record, annotation)
            multi = record['nh'] > 1 or any(r['secondary'] for r in records)
            chimeric = record['sa'] or any(r['supplementary'] for r in records) or len(primary) > 1
            status = 'chimeric' if chimeric else 'multimapped' if multi else 'low_mapq' if record['mapq'] < mapq else 'unique'
            row.update(mapping_status=status, mapq=record['mapq'], physical_strand=physical_strand(record, unit) or '',
                region_bases=detail['regions'], matched_reference_bases=interval_length(record['blocks']),
                genomic_span=record['end'] - record['start'], query_reverse=record['reverse'])
            for start, end in record['junctions']:
                known_strands = [strand for strand in ('+', '-') if (record['chrom'], start, end, strand) in annotation.known_junctions]
                row['junctions'].append({'chrom': record['chrom'], 'start': start, 'end': end, 'known': bool(known_strands), 'known_strands': known_strands})
                row['known_junctions'] += bool(known_strands)
            gene = detail['gene']
            strands = {g.strand for g in detail['genes'] if g.strand in ('+', '-')}
            if status == 'unique':
                # Comparable molecule coordinates do not require a gene annotation.
                # Missing physical orientation cannot be mixed with recovered P5/P7 strands.
                if row['physical_strand']:
                    row['coordinate_key'] = ':'.join((record['chrom'], str(record['start']), str(record['end']), row['physical_strand'], ','.join(f'{s}-{e}' for s, e in record['junctions'])))
                if len(strands) > 1:
                    row['strand_evidence'] = 'annotation_conflict'
                elif len(strands) == 1 and row['physical_strand']:
                    row['strand_evidence'] = 'sense' if row['physical_strand'] in strands else 'antisense'
                if gene is not None:
                    row.update(gene_id=gene.gene_id, gene_name=gene.name, gene_strand=gene.strand, eligible_gene=True)
                    profiles[name] = {'gene_body': profile_vector(record['blocks'], gene.representative, gene.strand), '_lengths': {'gene_body': interval_length(gene.representative), **{key: interval_length(value) for key, value in gene.segments.items()}}}
                    for region, exons in gene.segments.items():
                        profiles[name][region] = profile_vector(record['blocks'], exons, gene.strand, 30)
        rows.append(row)
    return rows, profiles


def summarize_profiles(rows, vectors, min_gene_units=10, min_bias_genes=20):
    counts = collections.Counter(row['gene_id'] for row in rows if row['eligible_gene'] and
        row['unit_id'] in vectors and sum(vectors[row['unit_id']]['gene_body']) > 0)
    # Plot every observed gene; depth thresholds apply only to end-bias inference.
    eligible = set(counts)
    bias_eligible = {gene for gene, count in counts.items() if count >= min_gene_units}
    by_gene, lengths = collections.defaultdict(dict), {}
    for row in rows:
        if row['gene_id'] not in eligible or row['unit_id'] not in vectors:
            continue
        source = vectors[row['unit_id']]
        lengths[row['gene_id']] = source['_lengths']
        for kind, vector in source.items():
            if kind.startswith('_'):
                continue
            target = by_gene[row['gene_id']].setdefault(kind, [0.] * len(vector))
            for index, value in enumerate(vector):
                target[index] += value
    output, end_fractions, both_ends = {}, [], 0
    for kind, bins in (('gene_body', 100), ('five_prime_UTR', 30), ('CDS', 30), ('three_prime_UTR', 30)):
        normalized, missing_region = [], 0
        for gene in sorted(by_gene):
            values = by_gene[gene].get(kind, [])
            whole_total = sum(by_gene[gene].get('gene_body', []))
            whole_length, region_length = lengths[gene]['gene_body'], lengths[gene].get(kind, 0)
            if not region_length:
                missing_region += 1
                continue
            if not whole_total or not whole_length:
                continue
            # All regions share the same whole-transcript per-base denominator.
            # A present but uncovered region contributes zeros, not a missing gene.
            normalizer = whole_total / whole_length
            normalized.append([value / (region_length / bins) / normalizer for value in values])
            if kind == 'gene_body' and gene in bias_eligible:
                five, three = sum(values[:20]), sum(values[-20:])
                if five > 0 and three > 0:
                    both_ends += 1
                if five + three > 0:
                    end_fractions.append(five / (five + three))
        mean = [statistics.mean(v[i] for v in normalized) for i in range(bins)] if normalized else [None] * bins
        sem = [statistics.stdev(v[i] for v in normalized) / math.sqrt(len(normalized)) for i in range(bins)] if len(normalized) > 1 else [None] * bins
        supporting_units = sum(1 for row in rows if row['gene_id'] in eligible and row['unit_id'] in vectors
                               and sum(vectors[row['unit_id']].get(kind, [])) > 0)
        output[kind] = {'contributing_genes': len(normalized), 'genes_missing_region': missing_region,
            'covered_genes': sum(any(value > 0 for value in vector) for vector in normalized),
            'supporting_units': supporting_units, 'minimum_units_per_gene': 1,
            'mean_relative_coverage': mean, 'sem_across_genes': sem}
    available = len(end_fractions) >= min_bias_genes
    median_fraction = statistics.median(end_fractions) if available else None
    ratio = median_fraction / (1 - median_fraction) if available and median_fraction < 1 else None
    output['bias'] = {'status': 'available' if available else 'insufficient_terminal_supported_genes',
        'eligible_terminal_supported_genes': len(end_fractions), 'eligible_genes_with_both_ends': both_ends,
        'minimum_genes': min_bias_genes, 'minimum_assigned_units_per_gene': min_gene_units,
        'median_five_prime_fraction': median_fraction,
        'median_end_balance': 2 * median_fraction - 1 if available else None,
        'five_to_three_ratio': ratio,
        'five_to_three_ratio_status': 'five_prime_only' if available and median_fraction == 1 else 'available' if available else 'insufficient_support',
        'definition': 'Per-gene first/last 20% matched coverage of concatenated representative-transcript exons. Genes need coverage at either end. Median fraction = five/(five+three); ratio = median_fraction/(1-median_fraction); balance = 2*median_fraction-1. A 5-prime-only median has no finite ratio.',
        'interpretation': 'Conditional coverage evidence; not an RNA-integrity or selection-chemistry diagnosis'}

    return output


def summarize_scope(scope_id, rows, vectors, annotation, sequences, assay, min_gene_units=10, min_bias_genes=20, min_motif_units=100):
    total = len(rows)
    mapping = collections.Counter(row['mapping_status'] for row in rows)
    mapped = sum(row['mapped'] for row in rows)
    family = collections.Counter(row['family'] for row in rows)
    family_sets = collections.Counter(';'.join(row['families']) for row in rows if len(row['families']) > 1)
    region = collections.Counter()
    for row in rows:
        if row['mapping_status'] == 'unique':
            region.update(row['region_bases'])
    strands = collections.Counter(row['strand_evidence'] for row in rows)
    informed = strands['sense'] + strands['antisense']
    gene_counts = collections.Counter(row['gene_id'] for row in rows if row['eligible_gene'])
    assigned = sum(gene_counts.values())
    splice_rows = [row for row in rows if row['mapping_status'] == 'unique']
    junction_total = sum(len(row['junctions']) for row in splice_rows)
    known = sum(row['known_junctions'] for row in splice_rows)
    lengths = [row['insert_length'] for row in rows if row['insert_length'] > 0]
    coordinates = collections.Counter(row['coordinate_key'] for row in rows if row['coordinate_key'])
    umi_rows = [row for row in rows if row['coordinate_key'] and row['umi_status'] == 'complete' and row['umi_sequence']]
    exact = collections.Counter((row['coordinate_key'], row['umi_sequence']) for row in umi_rows)
    rarefaction = [{'assigned_units': draw, 'expected_detected_genes': expected_detected(list(gene_counts.values()), draw)} for draw in sorted({0, assigned, *[int(assigned * i / 10) for i in range(1, 10)]})]
    profiles = summarize_profiles(rows, vectors, min_gene_units, min_bias_genes)
    motifs = discover_motifs([(row['unit_id'], sequences[row['unit_id']]) for row in rows if row['unit_id'] in sequences], min_motif_units) if assay == 'clip' else {'status': 'not_requested', 'motifs': []}
    interval = wilson_interval(strands['sense'], informed)
    observations = []
    if informed >= 100:
        if interval[0] > .8:
            observations.append('The physical P5-to-P7 strand predominantly agrees with annotated RNA orientation.')
        elif interval[1] < .2:
            observations.append('The physical P5-to-P7 strand predominantly opposes annotated RNA orientation.')
        elif interval[0] > .35 and interval[1] < .65:
            observations.append('Both physical strand orientations are represented; strong strand specificity is not supported.')
        else:
            observations.append('Physical strand evidence is mixed; annotation overlap, orientation recovery and protocol may contribute.')
    else:
        observations.append('Too few informative units to characterize strand specificity confidently.')
    if mapped or any(row.get('source_status') == 'assigned' for row in rows):
        observations.append(f"Annotation or independent sequence-panel ribosomal evidence accounts for {100 * (family['rRNA'] + family['mt-rRNA']) / total:.1f}% of extracted non-dimer units; incomplete references can underestimate this fraction.")
    if lengths:
        observations.append(f"Median insert length is {statistics.median(lengths):.0f} bp; {100 * sum(x < 100 for x in lengths) / len(lengths):.1f}% are shorter than 100 bp.")
    balance = profiles['bias']['median_end_balance']
    if balance is not None:
        if balance < -1 / 3:
            observations.append('Supported coverage is enriched toward 3′ transcript ends. This is compatible with 3′ capture, incomplete transcripts or other library/target biases; chemistry is not uniquely identified.')
        elif balance > 1 / 3:
            observations.append('Supported coverage is enriched toward 5′ transcript ends. This is compatible with end-biased capture or target distribution; chemistry is not uniquely identified.')
        else:
            observations.append('Supported transcript coverage does not show a greater than twofold median 5′/3′ end imbalance.')
    observations.append('Library chemistry and kit remain undetermined; these are descriptive features, not required metadata or protocol identification.')
    return {'scope_id': scope_id, 'units': total, 'mapped_units': mapped, 'mapping_counts': dict(mapping),
        'mapped_fraction': mapped / total if total else None,
        'family_counts': {key: family[key] for key in BIOTYPES},
        'classification_evidence_unit_counts': dict(collections.Counter(';'.join(row['family_evidence']) for row in rows)), 'ambiguous_family_sets': dict(family_sets),
        'family_denominator': 'All extracted non-dimer units; each unit counted once; multifamily mappings remain ambiguous',
        'region_bases': {key: region[key] for key in REGIONS}, 'region_denominator_bases': sum(region.values()),
        'region_denominator': 'Primary unique MAPQ-pass nonchimeric matched reference bases; CIGAR N and D excluded',
        'strand_counts': {key: strands[key] for key in ('sense', 'antisense', 'unknown', 'annotation_conflict')},
        'strand_informative_units': informed, 'physical_sense_fraction': strands['sense'] / informed if informed else None,
        'physical_sense_wilson95': interval, 'spliced_units': sum(bool(row['junctions']) for row in splice_rows),
        'splice_denominator_units': len(splice_rows), 'junction_observations': junction_total,
        'annotated_junction_observations': known, 'annotated_junction_fraction': known / junction_total if junction_total else None,
        'gene_assigned_units': assigned, 'detected_genes': len(gene_counts), 'gene_counts': dict(gene_counts),
        'top10_gene_fraction': sum(count for _, count in gene_counts.most_common(10)) / assigned if assigned else None,
        'rarefaction': rarefaction, 'profiles': profiles, 'motifs': motifs,
        'insert_median_bp': statistics.median(lengths) if lengths else None,
        'short_insert_fraction_lt100': sum(value < 100 for value in lengths) / len(lengths) if lengths else None,
        'umi': {'status_counts': dict(collections.Counter(row['umi_status'] for row in rows)),
                'eligible_coordinate_umi_units': len(umi_rows), 'observed_exact_coordinate_umi_groups': len(exact),
                'eligible_coordinate_definition': 'Unique MAPQ-pass nonchimeric genome alignments with recovered physical P5/P7 strand and complete observed UMI; annotation assignment is not required. Coordinates use reference start/end and splice junction chain, not source-panel alignments.',
                'exact_repeat_units': len(umi_rows) - len(exact), 'coordinate_collision_units': sum(coordinates.values()) - len(coordinates),
                'interpretation': 'Observed groups only; no deduplication or complexity extrapolation. Complete means the expected tag bases were recovered, not error-free or uniquely molecular. ONT errors may split UMI groups and collisions may merge them.'},
        'inferred_characteristics': observations}


def write_tsv(path, rows, fields):
    with open_text(path, 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False, separators=(',', ':')) if isinstance(value, (dict, list)) else '' if value is None else value for key, value in row.items()})


def _svg_header(width, height, title):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{html.escape(title, quote=True)}">',
        '<rect width="100%" height="100%" fill="#ffffff"/><style>text{font-family:Arial,sans-serif;fill:#203344;font-size:12px}.title{font-size:18px;font-weight:600}.note{font-size:10px;fill:#526577}</style>', f'<text class="title" x="24" y="30">{html.escape(title)}</text>']


def bars_svg(values, title, footnote='', percentage=True):
    rows = [(name, value) for name, value in values.items() if value]
    height = max(180, 100 + 29 * len(rows))
    svg, maximum, total = _svg_header(720, height, title), max((v for _, v in rows), default=1), sum(v for _, v in rows)
    for index, (name, value) in enumerate(rows):
        y, width = 55 + index * 29, 420 * value / maximum
        svg.append(f'<text x="20" y="{y + 15}">{html.escape(str(name))}</text><rect x="185" y="{y}" width="{width:.2f}" height="20" rx="2" fill="#267f9e"/>')
        label = f'{100 * value / total:.1f}%' if percentage and total else f'{value:g}'
        svg.append(f'<text x="{min(650, 190 + width):.2f}" y="{y + 15}">{label}</text>')
    if not rows:
        svg.append('<text x="24" y="85">No eligible observations</text>')
    svg.append(f'<text class="note" x="24" y="{height - 18}">{html.escape(footnote[:115])}</text></svg>')
    return ''.join(svg)


def profile_svg(profiles, title):
    """Display observed 5′ UTR/CDS/3′ UTR coverage, including single-insert genes."""
    import io
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    color = '#3276b5'
    with plt.rc_context({'svg.fonttype': 'none', 'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.spines.top': False, 'axes.spines.right': False}):
        fig, ax = plt.subplots(figsize=(7.6, 3.1))
        labels, positions, means, errors = [], [], [], []
        for index, (kind, label) in enumerate((('five_prime_UTR', '5′ UTR'), ('CDS', 'CDS'), ('three_prime_UTR', '3′ UTR'))):
            profile = profiles[kind]
            mean = np.array([np.nan if value is None else value for value in profile['mean_relative_coverage']], dtype=float)
            sem = np.array([np.nan if value is None else value for value in profile['sem_across_genes']], dtype=float)
            x = index + (np.arange(len(mean)) + .5) / len(mean)
            positions.append(x)
            means.append(mean)
            errors.append(sem if profile['contributing_genes'] >= 2 else np.full_like(mean, np.nan))
            if not np.any(np.isfinite(mean)):
                empty = 'No annotated region' if profile.get('genes_missing_region', 0) else 'No usable alignments'
                ax.text(index + .5, .5, empty, ha='center', va='center', fontsize=8, color='#64748b', transform=ax.get_xaxis_transform())
            support = profile.get('supporting_units')
            support_text = '' if support is None else ' · ' + str(support) + ' inserts'
            labels.append(label + '\n' + str(profile['contributing_genes']) + ' genes' + support_text)
        x, mean, sem = map(np.concatenate, (positions, means, errors))
        ax.plot(x, mean, color=color, linewidth=1.7)
        ax.fill_between(x, np.maximum(0, mean - sem), mean + sem, color=color, alpha=.18, linewidth=0)
        ax.set(xlim=(0, 3), xticks=[.5, 1.5, 2.5], xticklabels=labels, ylabel='Relative coverage')
        for boundary in (1, 2):
            ax.axvline(boundary, color='#aebac4', linewidth=.7, linestyle='--')
        ax.set_ylim(0, max(1.05, ax.get_ylim()[1]))
        ax.axhline(1, color='#b8c0c8', linestyle=':', linewidth=.8)
        ax.grid(axis='y', alpha=.12)
        ax.set_title(title, loc='left', fontsize=10, pad=12)
        fig.text(.1, .025, 'Observed coverage · equal gene weighting · shading: SEM across genes', fontsize=7, color='#526577')
        fig.subplots_adjust(left=.1, right=.98, bottom=.24, top=.85)
        buffer = io.StringIO()
        fig.savefig(buffer, format='svg', metadata={'Title': title, 'Description': '5′ UTR/CDS/3′ UTR metaplot. All observed genes, including single-insert genes; whole-transcript per-base normalization; intron gaps excluded.'})
        plt.close(fig)
        return buffer.getvalue()


def rarefaction_svg(points, title):
    """Observed gene-detection rarefaction; intentionally no extrapolated curve."""
    import io
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with plt.rc_context({'svg.fonttype': 'none', 'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.spines.top': False, 'axes.spines.right': False}):
        fig, ax = plt.subplots(figsize=(7.2, 3.5))
        xs = [point['assigned_units'] for point in points]
        ys = [point['expected_detected_genes'] for point in points]
        ax.plot(xs, ys, color='#244a64', marker='o', markersize=3, linewidth=1.7)
        ax.set(xlabel='Observed gene-assigned library units sampled', ylabel='Expected detected genes', title=title)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(alpha=.15)
        fig.text(.11, .025, 'Subsampling observed units only. This does not estimate unseen molecules or PCR library complexity.', fontsize=7, color='#526577')
        fig.subplots_adjust(left=.11, right=.97, bottom=.21, top=.86)
        buffer = io.StringIO()
        fig.savefig(buffer, format='svg', metadata={'Title': title})
        plt.close(fig)
        return buffer.getvalue()


def reference_source_comparison(rows, assembly):
    """Compare known assembly labels with candidate references, never diagnose contamination."""
    host = ({'hg19': 'Homo sapiens', 'hg38': 'Homo sapiens', 'hs1': 'Homo sapiens',
             'mm9': 'Mus musculus', 'mm10': 'Mus musculus', 'mm39': 'Mus musculus'}).get(str(assembly).lower())
    known_taxa = {'Homo sapiens', 'Mus musculus', 'Escherichia coli', 'Bacillus subtilis',
                  'Pseudomonas aeruginosa', 'Staphylococcus aureus'}
    assigned = [row for row in rows if row.get('source_status') == 'assigned' and row.get('source_class') == 'rRNA']
    unavailable = bool(rows) and all(row.get('source_status') in (None, 'not_available') for row in rows)
    result = {'assembly': assembly, 'host_reference_taxon': host,
        'status': 'source_screen_not_available' if unavailable else 'assembly_taxon_unknown' if host is None else 'available' if assigned else 'no_assigned_source_evidence',
        'assigned_rRNA_units': len(assigned), 'denominator_units': len(rows),
        'host_only_candidate_units': None, 'host_and_other_candidate_units': None,
        'non_host_only_candidate_units': None, 'unresolved_candidate_units': None,
        'non_host_only_fraction_of_all_units': None,
        'interpretation': 'Counts compare the selected assembly taxon with near-best rRNA reference candidates. A candidate set containing the host is never counted as non-host-only. Limited, conserved reference matches do not establish species identity, contamination, or the cause of low genome mapping.'}
    if host is None or unavailable:
        return result
    categories = collections.Counter()
    for row in assigned:
        taxa = {part.strip() for part in row.get('source_taxon', '').split(',') if part.strip()}
        if not taxa or not taxa <= known_taxa:
            categories['unresolved_candidate_units'] += 1
        elif taxa == {host}:
            categories['host_only_candidate_units'] += 1
        elif host in taxa:
            categories['host_and_other_candidate_units'] += 1
        else:
            categories['non_host_only_candidate_units'] += 1
    for key in ('host_only_candidate_units', 'host_and_other_candidate_units', 'non_host_only_candidate_units', 'unresolved_candidate_units'):
        result[key] = categories[key]
    result['non_host_only_fraction_of_all_units'] = categories['non_host_only_candidate_units'] / len(rows) if rows else None
    return result


def run_analysis(bam, gtf, read_table, fastq, sample, assay, out_dir, samtools='samtools', mapq=10,
                 mitochondrial=('chrM', 'MT'), index_assignments=None, min_gene_units=10, min_bias_genes=20, min_motif_units=100, assembly='', preparation_json=None, preparation_tsv=None, source_assignments=None, source_summary=None):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    preparation = json.loads(Path(preparation_json).read_text()) if preparation_json else None
    annotation = read_annotation(gtf, mitochondrial)
    units = read_units(read_table, index_assignments)
    preparation_units = {}
    if preparation_tsv:
        with open_text(preparation_tsv) as handle:
            for record in csv.DictReader(handle, delimiter='\t'):
                name = record['unit_id']
                if name in preparation_units:
                    raise ValueError('Duplicate preparation unit: ' + name)
                if name not in units:
                    raise ValueError('Preparation unit absent from non-dimer adapter table: ' + name)
                for key in ('input_length', 'analysis_length', 'p5_trim_bp', 'p7_trim_bp'):
                    record[key] = int(record[key])
                if record['input_length'] != int(units[name]['insert_length']):
                    raise ValueError('Preparation input length differs from adapter unit: ' + name)
                if record['analysis_length'] != max(0, record['input_length'] - record['p5_trim_bp'] - record['p7_trim_bp']):
                    raise ValueError('Inconsistent preparation trim coordinates: ' + name)
                preparation_units[name] = record
                units[name]['umi_status'] = record['umi_status']
                units[name]['umi_sequence'] = record['umi_sequence']
        if set(preparation_units) != set(units):
            raise ValueError('Preparation sidecar does not cover every non-dimer adapter unit')
        if preparation and preparation.get('units_tsv_sha256') and hashlib.sha256(Path(preparation_tsv).read_bytes()).hexdigest() != preparation['units_tsv_sha256']:
            raise ValueError('Preparation sidecar checksum mismatch')
    elif preparation and preparation.get('schema_version', 1) >= 2:
        raise ValueError('Version 2 RNA preparation requires --preparation-tsv')
    alignments, references = read_alignments(bam, samtools)
    rows, vectors = analyze_units(units, alignments, annotation, mapq)
    by_id = {row['unit_id']: row for row in rows}
    for name, record in preparation_units.items():
        by_id[name].update(preparation_status=record['status'], p5_random_tag=record.get('p5_random_tag', ''),
            fixed_barcode=record.get('fixed_barcode', ''), p7_umi=record.get('p7_umi', ''),
            p5_trim_bp=record['p5_trim_bp'], p7_trim_bp=record['p7_trim_bp'])
    source_units = {}
    source_metadata = json.loads(Path(source_summary).read_text()) if source_summary else None
    if source_assignments:
        with open_text(source_assignments) as handle:
            for record in csv.DictReader(handle, delimiter='\t'):
                name = record['unit_id']
                if name in source_units or name not in by_id:
                    raise ValueError('Duplicate or unknown source-screen unit: ' + name)
                source_units[name] = record
                row = by_id[name]
                row.update(source_class=record.get('source_class', ''), source_status=record['status'],
                    source_reference_name=record.get('reference_name', ''), source_taxon=record.get('taxon', ''),
                    source_taxon_ambiguous=record.get('taxon_ambiguous', '') == '1',
                    source_matched_fraction=record.get('matched_fraction', ''), source_identity=record.get('identity', ''),
                    genome_family=row['family'], genome_families=list(row['families']))
                if record['status'] == 'assigned' and record.get('source_class') == 'rRNA':
                    row.update(family='rRNA', families=['rRNA'], eligible_gene=False,
                        strand_evidence='unknown', excluded_from_gene_qc_reason='rRNA_source_panel')
                    row['family_evidence'] = sorted(set(row['family_evidence']) | {'independent_rRNA_panel'})
                    vectors.pop(name, None)
    # Annotated ribosomal genes also do not define ordinary transcript-body or CLIP motif profiles.
    for row in rows:
        if row['family'] in ('rRNA', 'mt-rRNA'):
            vectors.pop(row['unit_id'], None)
    sequences, seen_fastq, terminal_a = {}, set(), {}
    for name, sequence, quality in read_fastq(fastq):
        if name in seen_fastq:
            raise ValueError('Duplicate extracted FASTQ unit_id: ' + name)
        seen_fastq.add(name)
        row = by_id.get(name)
        if row is None:
            raise ValueError('Analysis FASTQ unit is absent from non-dimer adapter table: ' + name)
        if row['insert_length'] and row['insert_length'] < len(sequence):
            raise ValueError('Analysis FASTQ exceeds physical adapter-bounded insert length: ' + name)
        if preparation_units and len(sequence) != preparation_units[name]['analysis_length']:
            raise ValueError('Analysis length disagrees with per-unit preparation: ' + name)
        if preparation and not preparation_units and len(sequence) != max(0, row['insert_length'] - preparation['p5_trim_bp'] - preparation['p7_trim_bp']):
            raise ValueError('Analysis length disagrees with configured preparation: ' + name)
        row['analysis_length'] = len(sequence)
        row['technical_bases_removed'] = max(0, row['insert_length'] - len(sequence))
        if not row['insert_length']:
            row['insert_length'] = len(sequence)
        if row['eligible_gene'] and row['gene_strand'] in ('+', '-') and row['family'] not in ('rRNA', 'mt-rRNA'):
            oriented = reverse_complement(sequence) if row['query_reverse'] != (row['gene_strand'] == '-') else sequence
            if assay == 'clip':
                primary = [r for r in alignments[name] if not r['secondary'] and not r['supplementary']]
                sequences[name] = motif_aligned_sequence(sequence, max(primary, key=lambda r: r['mapq']), row['gene_strand'])
            terminal_a[name] = int(oriented.endswith('A' * 8))
    missing = sorted(set(units) - seen_fastq)
    if missing:
        trim = (preparation['p5_trim_bp'] + preparation['p7_trim_bp']) if preparation else 0
        allowed = (all(preparation_units[name]['analysis_length'] == 0 and not by_id[name]['mapped'] for name in missing)
                   if preparation_units else bool(preparation) and all(0 < by_id[name]['insert_length'] <= trim and not by_id[name]['mapped'] for name in missing))
        if not allowed or (preparation is not None and len(missing) != preparation.get('empty_after_trim', -1)):
            raise ValueError(f'{len(missing)} unexpected adapter units missing from analysis FASTQ (first: {missing[0]})')
        for name in missing:
            by_id[name].update(mapping_status='empty_after_technical_trim', analysis_length=0, technical_bases_removed=by_id[name]['insert_length'])
    if preparation and (preparation.get('input_units') != len(units) or preparation.get('retained_units') != len(seen_fastq) or preparation.get('empty_after_trim') != len(missing)):
        raise ValueError('Preparation summary counts disagree with adapter units or analysis FASTQ')
    if source_assignments and set(source_units) != seen_fastq:
        raise ValueError('Source-screen sidecar does not cover exactly the analysis FASTQ units')
    summaries = {}
    for scope in ['total'] + sorted({row['scope_id'] for row in rows if row['scope_id']}):
        selected = rows if scope == 'total' else [row for row in rows if row['scope_id'] == scope]
        summaries[scope] = summarize_scope(scope, selected, vectors, annotation, sequences, assay, min_gene_units, min_bias_genes, min_motif_units)
        tails = [terminal_a[row['unit_id']] for row in selected if row['unit_id'] in terminal_a]
        summaries[scope]['source_evidence'] = {
            'status': source_metadata.get('status', 'complete') if source_metadata else 'complete' if source_assignments else 'not_screened',
            'status_counts': dict(collections.Counter(row.get('source_status', 'not_screened') for row in selected)),
            'class_counts': dict(collections.Counter(row.get('source_class') for row in selected if row.get('source_status') == 'assigned')),
            'assigned_units': sum(row.get('source_status') == 'assigned' for row in selected),
            'ambiguous_taxon_units': sum(bool(row.get('source_taxon_ambiguous')) for row in selected),
            'candidate_taxa_counts': dict(collections.Counter(row.get('source_taxon', '') for row in selected if row.get('source_status') == 'assigned')),
            'denominator_units': len(selected),
            'interpretation': 'Independent panel sequence evidence; does not change genome mapping counts. Candidate reference taxa are not confirmed organism identity or contamination.'}
        summaries[scope]['source_evidence']['reference_comparison'] = reference_source_comparison(selected, assembly)
        summaries[scope]['umi']['auto_grip_p7_complete_units'] = sum(bool(row.get('p7_umi')) and row['umi_status'] == 'complete' for row in selected)
        summaries[scope]['terminal_poly_a_evidence'] = {'oriented_gene_assigned_units': len(tails), 'terminal_a_run_ge8_units': sum(tails), 'interpretation': 'Observed terminal homopolymer, not a direct poly(A)-selection or tail-length assay'}
    summary = {'schema_version': 1, 'sample': sample, 'assay': assay, 'assembly': assembly, 'status': 'beta', 'annotation': annotation.metadata,
        'alignment_reference_contigs': references,
        'annotation_contigs_missing_from_bam_header': sorted(set(annotation.metadata['annotation_chromosomes']) - set(references)),
        'alignment_contigs_missing_from_annotation': sorted(set(references) - set(annotation.metadata['annotation_chromosomes'])),
        'settings': {'mapq': mapq, 'min_gene_units': min_gene_units, 'metaplot_min_gene_units': 1, 'min_bias_genes': min_bias_genes, 'min_motif_units': min_motif_units, 'motif_engine': 'HOMER', 'mitochondrial_contigs': list(mitochondrial)},
        'preparation': preparation, 'source_screen': source_metadata,
        'definitions': {'gene_assignment': 'Single overlapping exonic gene, or single overlapping gene if no exon is hit, on a unique MAPQ-pass nonchimeric primary alignment',
            'region_priority': list(REGIONS), 'junctions': 'CIGAR N gaps; exact annotation matches are support, not independent validation',
            'gene_counts': 'Observed assigned library units, not transcript abundance or deduplicated RNA molecules',
            'family_classification': 'Explicit biotype first; if missing, CDS features support coding classification and conservative conventional symbols support selected ncRNA families. Other untyped genes remain unknown; evidence is recorded per gene.',
            'metaplot': 'All genes with observed coverage on the representative transcript, including single-insert genes. Equal gene weighting; each region shares the whole-transcript per-base denominator. Existing unique MAPQ-pass and RNA-source filters retained. End-bias inference has separate support thresholds.',
            'rarefaction': 'Expected genes detected in subsamples of observed gene-assigned units; no extrapolation'},
        'limitations': LIMITATIONS, 'scopes': summaries,
        'files': ['summary.json', 'units.tsv.gz', 'gene_counts.tsv', 'junctions.tsv', 'rarefaction.tsv', 'profiles.tsv', 'motifs.tsv']}
    unit_fields = ['unit_id', 'read_id', 'scope_id', 'insert_length', 'analysis_length', 'technical_bases_removed', 'mapping_status', 'mapq', 'family', 'families', 'family_evidence', 'gene_id', 'gene_name', 'gene_strand', 'orientation', 'physical_strand', 'strand_evidence', 'matched_reference_bases', 'genomic_span', 'region_bases', 'junctions', 'known_junctions', 'umi_status', 'umi_sequence', 'preparation_status', 'p5_random_tag', 'fixed_barcode', 'p7_umi', 'p5_trim_bp', 'p7_trim_bp', 'genome_family', 'genome_families', 'source_status', 'source_class', 'source_reference_name', 'source_taxon', 'source_taxon_ambiguous', 'source_matched_fraction', 'source_identity', 'excluded_from_gene_qc_reason']
    write_tsv(out / 'units.tsv.gz', rows, unit_fields)
    if preparation_tsv:
        if str(preparation_tsv).endswith('.gz'):
            (out / 'preparation_units.tsv.gz').write_bytes(Path(preparation_tsv).read_bytes())
        else:
            with open_text(preparation_tsv) as handle, gzip.open(out / 'preparation_units.tsv.gz', 'wt') as copied:
                copied.write(handle.read())
        summary['files'].append('preparation_units.tsv.gz')
    if source_assignments:
        with open_text(source_assignments) as handle, gzip.open(out / 'source_assignments.tsv.gz', 'wt') as copied:
            copied.write(handle.read())
        summary['files'].append('source_assignments.tsv.gz')
    if source_metadata is not None:
        (out / 'source_summary.json').write_text(json.dumps(source_metadata, indent=2) + '\n')
        summary['files'].append('source_summary.json')
    genes_out, junctions_out, rare_out, profile_out, motif_out = [], [], [], [], []
    for scope, data in summaries.items():
        for gene_id, count in sorted(data['gene_counts'].items(), key=lambda item: (-item[1], item[0])):
            gene = annotation.genes[gene_id]
            genes_out.append({'scope_id': scope, 'gene_id': gene_id, 'gene_name': gene.name, 'biotype': gene.biotype, 'family': gene.family, 'family_evidence': gene.family_evidence, 'units': count})
        rare_out.extend(dict(scope_id=scope, **record) for record in data['rarefaction'])
        for kind in ('gene_body', 'five_prime_UTR', 'CDS', 'three_prime_UTR'):
            profile = data['profiles'][kind]
            for index, (mean, sem) in enumerate(zip(profile['mean_relative_coverage'], profile['sem_across_genes'])):
                profile_out.append({'scope_id': scope, 'region': kind, 'bin': index, 'mean_relative_coverage': mean, 'sem_across_genes': sem, 'contributing_genes': profile['contributing_genes'], 'supporting_units': profile['supporting_units'], 'covered_genes': profile['covered_genes']})
        motif_out.extend(dict(scope_id=scope, **record) for record in data['motifs']['motifs'])
        jcounts = collections.Counter()
        for row in rows:
            if row['mapping_status'] == 'unique' and (scope == 'total' or row['scope_id'] == scope):
                for junction in row['junctions']:
                    jcounts[(junction['chrom'], junction['start'], junction['end'], junction['known'])] += 1
        for (chrom, start, end, known), count in sorted(jcounts.items()):
            junctions_out.append({'scope_id': scope, 'chrom': chrom, 'start': start, 'end': end, 'annotated': known, 'supporting_units': count})
        slug = 'scope_' + hashlib.sha256(scope.encode()).hexdigest()[:12]
        plots = {'composition': slug + '.composition.svg', 'regions': slug + '.regions.svg', 'profiles': slug + '.profiles.svg', 'rarefaction': slug + '.rarefaction.svg'}
        (out / plots['composition']).write_text(bars_svg(data['family_counts'], sample + ' · ' + scope + ': RNA families', 'All extracted non-dimer units; ambiguous family mappings retained.'))
        (out / plots['regions']).write_text(bars_svg(data['region_bases'], sample + ' · ' + scope + ': annotated regions', 'Unique MAPQ-pass matched bases; intron skips do not contribute coverage.'))
        (out / plots['profiles']).write_text(profile_svg(data['profiles'], sample + ' · ' + scope + ': transcript coverage'))
        (out / plots['rarefaction']).write_text(rarefaction_svg(data['rarefaction'], sample + ' · ' + scope + ': observed gene detection'))
        data['plots'] = plots
        summary['files'].extend(plots.values())
    write_tsv(out / 'gene_counts.tsv', genes_out, ['scope_id', 'gene_id', 'gene_name', 'biotype', 'family', 'family_evidence', 'units'])
    write_tsv(out / 'junctions.tsv', junctions_out, ['scope_id', 'chrom', 'start', 'end', 'annotated', 'supporting_units'])
    write_tsv(out / 'rarefaction.tsv', rare_out, ['scope_id', 'assigned_units', 'expected_detected_genes'])
    write_tsv(out / 'profiles.tsv', profile_out, ['scope_id', 'region', 'bin', 'mean_relative_coverage', 'sem_across_genes', 'contributing_genes', 'supporting_units', 'covered_genes'])
    write_tsv(out / 'motifs.tsv', motif_out, ['scope_id', 'motif', 'dna_kmer', 'observed_units', 'expected_units', 'target_fraction', 'background_fraction', 'enrichment', 'homer_log10_p', 'adapter_sequence_match', 'pwm'])
    (out / 'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    return summary


def cohort_analysis(summary_paths, out_dir, min_assigned_units=200, min_shared_genes=20):
    items = []
    for path in summary_paths:
        data = json.loads(Path(path).read_text())
        for scope, metrics in data['scopes'].items():
            items.append({'sample': data['sample'], 'scope': scope, 'label': data['sample'] + ' / ' + scope, 'assay': data['assay'], 'assembly': data.get('assembly', ''), 'annotation': data['annotation']['sha256'], 'metrics': metrics})
    pairs = []
    for left, right in itertools.combinations(items, 2):
        base = {'sample_a': left['sample'], 'scope_a': left['scope'], 'sample_b': right['sample'], 'scope_b': right['scope'], 'status': 'available', 'log_cpm_pearson': None, 'gene_jaccard': None, 'related_sample_scopes': left['sample'] == right['sample'], 'same_assay': left['assay'] == right['assay']}
        if left['assay'] != right['assay']:
            base['status'] = 'different_assay'
        elif not left['assembly'] or not right['assembly'] or left['assembly'] != right['assembly']:
            base['status'] = 'different_or_unknown_assembly'
        elif left['annotation'] != right['annotation']:
            base['status'] = 'different_annotation'
        elif min(left['metrics']['gene_assigned_units'], right['metrics']['gene_assigned_units']) < min_assigned_units:
            base['status'] = 'insufficient_assigned_units'
        else:
            a, b = left['metrics']['gene_counts'], right['metrics']['gene_counts']
            genes = sorted(set(a) | set(b))
            base['gene_jaccard'] = len(set(a) & set(b)) / len(genes) if genes else None
            total_a, total_b = sum(a.values()), sum(b.values())
            xs, ys = [math.log1p(a.get(g, 0) * 1e6 / total_a) for g in genes], [math.log1p(b.get(g, 0) * 1e6 / total_b) for g in genes]
            base['shared_detected_genes'] = len(set(a) & set(b))
            if base['shared_detected_genes'] >= min_shared_genes and len(genes) >= 3 and len(set(xs)) > 1 and len(set(ys)) > 1:
                mx, my = statistics.mean(xs), statistics.mean(ys)
                covariance = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
                norm = math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
                base['log_cpm_pearson'] = max(-1., min(1., covariance / norm)) if norm else None
            if base['log_cpm_pearson'] is None:
                base['status'] = 'insufficient_variable_genes'
        pairs.append(base)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    result = {'samples_scopes': [{k: item[k] for k in ('sample', 'scope', 'label', 'assay', 'assembly', 'annotation')} for item in items], 'pairs': pairs, 'minimum_assigned_units': min_assigned_units, 'minimum_shared_genes': min_shared_genes,
        'interpretation': 'Descriptive gene-profile similarity. No replicate grouping or sample identity is inferred. Total and child scopes are dependent. Only matching assay, assembly and annotation are compared; biology and library depth affect similarity.'}
    (out / 'cohort.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    write_tsv(out / 'cohort_similarity.tsv', pairs, ['sample_a', 'scope_a', 'sample_b', 'scope_b', 'status', 'log_cpm_pearson', 'gene_jaccard', 'related_sample_scopes', 'same_assay'])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('bam', 'gtf', 'read-table', 'fastq', 'sample', 'index-assignments', 'assembly', 'preparation-json', 'preparation-tsv', 'source-assignments', 'source-summary'):
        parser.add_argument('--' + name)
    parser.add_argument('--assay', choices=('rnaseq', 'clip'), default='rnaseq')
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--samtools', default='samtools')
    parser.add_argument('--mapq', type=int, default=10)
    parser.add_argument('--mitochondrial-contigs', default='chrM,MT')
    parser.add_argument('--min-gene-units', type=int, default=10, help='Minimum units per gene for end-bias inference; metaplots include every observed gene')
    parser.add_argument('--min-bias-genes', type=int, default=20)
    parser.add_argument('--min-motif-units', type=int, default=100)
    parser.add_argument('--cohort', nargs='+')
    args = parser.parse_args()
    if args.cohort:
        cohort_analysis(args.cohort, args.out_dir)
        return
    for name in ('bam', 'gtf', 'read_table', 'fastq', 'sample'):
        if not getattr(args, name):
            parser.error('--' + name.replace('_', '-') + ' is required outside cohort mode')
    if min(args.min_gene_units, args.min_bias_genes, args.min_motif_units) < 1 or args.mapq < 0:
        parser.error('Support thresholds must be positive and MAPQ nonnegative')
    run_analysis(args.bam, args.gtf, args.read_table, args.fastq, args.sample, args.assay, args.out_dir,
                 args.samtools, args.mapq, tuple(filter(None, args.mitochondrial_contigs.split(','))),
                 args.index_assignments, args.min_gene_units, args.min_bias_genes, args.min_motif_units, args.assembly or '', args.preparation_json, args.preparation_tsv, args.source_assignments, args.source_summary)


if __name__ == '__main__':
    main()
