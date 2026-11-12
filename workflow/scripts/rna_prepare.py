#!/usr/bin/env python3
"""Preserve ONT direction while removing supported, scope-specific RNA library tags."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

BARCODES = ('CCTATAT', 'TGCTATT')
SOURCES = [
    'https://www.nature.com/articles/s41557-022-01038-4',
    'https://static-content.springer.com/esm/art%3A10.1038%2Fs41557-022-01038-4/MediaObjects/41557_2022_1038_MOESM1_ESM.pdf',
    'https://www.encodeproject.org/documents/842f7424-5396-424a-a1a3-3f18707c3222/@@download/attachment/eCLIP_SOP_v1.P_110915.pdf',
]
FIELDS = ['unit_id', 'scope_id', 'orientation', 'input_length', 'analysis_length', 'status',
          'p5_trim_bp', 'p7_trim_bp', 'physical_analysis_start', 'physical_analysis_end',
          'ont_analysis_start', 'ont_analysis_end', 'p5_random_tag', 'fixed_barcode',
          'barcode_expected', 'barcode_start', 'barcode_end', 'barcode_edits',
          'p7_umi', 'p7_umi_quality', 'umi_status', 'umi_sequence']


def open_text(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def fastq_records(path):
    with open_text(path) as handle:
        while True:
            name = handle.readline().rstrip()
            if not name:
                return
            seq, plus, qual = handle.readline().rstrip(), handle.readline().rstrip(), handle.readline().rstrip()
            if not name.startswith('@') or not plus.startswith('+') or len(seq) != len(qual):
                raise ValueError('Malformed RNA preparation FASTQ: ' + name[:80])
            yield name[1:].split()[0], seq.upper(), qual


def trim_insert(sequence, quality, orientation, p5, p7):
    if orientation not in {'p5_to_p7', 'p7_to_p5'} and (p5 or p7):
        raise ValueError('Physical orientation is required for end trimming')
    front, tail = (p7, p5) if orientation == 'p7_to_p5' else (p5, p7)
    stop = max(front, len(sequence) - tail)
    return sequence[front:stop], quality[front:stop]


def scope_of(row):
    return row.get('scope_id', '') if row.get('assignment_status') == 'assigned' and row.get('scope_id') else 'unassigned'


def edit_distance_one(left, right):
    """Small exact Levenshtein distance, used only for seven-base adapter tags."""
    previous = list(range(len(right) + 1))
    for i, a in enumerate(left, 1):
        current = [i]
        for j, b in enumerate(right, 1):
            current.append(min(current[-1] + 1, previous[j] + 1, previous[j - 1] + (a != b)))
        previous = current
    return previous[-1]


def find_barcode(sequence):
    """Find a bounded tag; indels in the ONT tag or its prefix can shift the cut."""
    hits = []
    for start in range(2, 9):
        for length in (6, 7, 8):
            end = start + length
            if end + 10 > len(sequence):
                continue
            observed = sequence[start:end]
            for expected in BARCODES:
                edits = edit_distance_one(observed, expected)
                if edits <= 1:
                    hits.append(((edits, abs(start - 5), abs(length - 7)), start, end, expected, observed, edits))
    if not hits:
        return None
    best_score = min(hit[0] for hit in hits)
    best = [hit for hit in hits if hit[0] == best_score]
    if len({hit[1:4] for hit in best}) != 1:
        return None
    _, start, end, expected, observed, edits = best[0]
    return {'barcode_start': start, 'barcode_end': end, 'barcode_expected': expected,
            'fixed_barcode': observed, 'barcode_edits': edits}


def infer_scopes(records, rows, minimum_units=20, minimum_exact=10, minimum_fraction=.5,
                 minimum_distinct_tags=4, maximum_joint_tag_fraction=.5):
    evidence = collections.defaultdict(lambda: {'input_units': 0, 'eligible_units': 0,
        'exact_barcode_units': 0, 'valid_random_tags_units': 0, 'barcodes': collections.Counter(),
        'p5_tags': collections.Counter(), 'p7_tags': collections.Counter(), 'joint_tags': collections.Counter()})
    for name, sequence, _ in records:
        row = rows[name]
        ev = evidence[scope_of(row)]
        ev['input_units'] += 1
        if row.get('orientation') not in ('p5_to_p7', 'p7_to_p5') or len(sequence) < 22:
            continue
        seq = reverse_complement(sequence) if row['orientation'] == 'p7_to_p5' else sequence
        ev['eligible_units'] += 1
        if seq[5:12] not in BARCODES:
            continue
        ev['exact_barcode_units'] += 1
        ev['barcodes'][seq[5:12]] += 1
        p5, p7 = seq[:5], seq[-10:]
        if set(p5 + p7) <= set('ACGT'):
            ev['valid_random_tags_units'] += 1
            ev['p5_tags'][p5] += 1
            ev['p7_tags'][p7] += 1
            ev['joint_tags'][(p5, p7)] += 1
    result = {}
    for scope, ev in evidence.items():
        eligible, exact = ev['eligible_units'], ev['exact_barcode_units']
        fraction = exact / eligible if eligible else None
        dominant = max(ev['joint_tags'].values(), default=0) / ev['valid_random_tags_units'] if ev['valid_random_tags_units'] else None
        diverse = len(ev['p5_tags']) >= minimum_distinct_tags and len(ev['p7_tags']) >= minimum_distinct_tags and dominant is not None and dominant <= maximum_joint_tag_fraction
        decision = ('insufficient_eligible_units' if eligible < minimum_units else
                    'insufficient_exact_barcode_support' if exact < minimum_exact or fraction < minimum_fraction else
                    'insufficient_random_tag_diversity' if not diverse else 'supported')
        result[scope] = {key: value for key, value in ev.items() if not key.endswith('_tags')}
        result[scope].update(exact_fraction=fraction, p5_random_distinct=len(ev['p5_tags']),
            p7_random_distinct=len(ev['p7_tags']), max_joint_tag_fraction=dominant,
            supported=decision == 'supported', decision=decision)
    return result


def prepare(fastq, table, output, summary, p5=0, p7=0, units_tsv=None, auto_grip=True,
            structure_json=None, minimum_units=20, minimum_exact=10, minimum_fraction=.5):
    if min(p5, p7) < 0 or max(p5, p7) > 1000:
        raise ValueError('Technical end trimming must be within 0..1000 bases')
    if minimum_units < 1 or minimum_exact < 1 or not 0 < minimum_fraction <= 1:
        raise ValueError('Invalid GRIP scope evidence thresholds')
    rows = {}
    with open_text(table) as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['unit_id'] in rows:
                raise ValueError('Duplicate adapter unit: ' + row['unit_id'])
            rows[row['unit_id']] = row
    records, seen = [], set()
    for record in fastq_records(fastq):
        name, seq, qual = record
        if name in seen:
            raise ValueError('Duplicate preparation FASTQ unit: ' + name)
        seen.add(name)
        if name not in rows:
            raise ValueError('Preparation FASTQ unit missing from adapter table: ' + name)
        if int(rows[name].get('insert_length') or len(seq)) != len(seq):
            raise ValueError('Preparation FASTQ length differs from adapter table: ' + name)
        records.append(record)
    structure = (json.loads(str(structure_json)) if str(structure_json).lstrip().startswith('{') else json.loads(Path(structure_json).read_text())) if structure_json else {}
    if 'library_structure' in structure:
        structure = structure['library_structure']
    has_custom_umi = bool(structure.get('umis')) or any(row.get('umi_status', '') not in ('', 'not_configured') for row in rows.values())
    bypass = ('manual_end_trimming' if p5 or p7 else 'configured_umi_extraction' if has_custom_umi else 'auto_disabled' if not auto_grip else None)
    evidence = infer_scopes(records, rows, minimum_units, minimum_exact, minimum_fraction)
    if bypass:
        for ev in evidence.values():
            ev['supported'] = False
            ev['decision'] = 'bypassed_' + bypass
    units_tsv = Path(units_tsv) if units_tsv else Path(summary).with_suffix('.units.tsv.gz')
    for path in (output, summary, units_tsv):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    counts = {'schema_version': 2, 'input_units': len(records), 'retained_units': 0, 'empty_after_trim': 0,
              'p5_trim_bp': p5, 'p7_trim_bp': p7, 'mode': 'manual' if p5 or p7 else 'passthrough' if bypass else 'auto_grip',
              'units_tsv_file': units_tsv.name, 'umi_complete_units': 0,
              'method': 'Scope-specific published GRIP/eCLIP X1 structure; observed P7 N10 UMI, P5 random tag and fixed barcode recorded separately',
              'auto_grip': {'enabled': auto_grip, 'bypass_reason': bypass,
                  'thresholds': {'minimum_eligible_units': minimum_units, 'minimum_exact_barcode_units': minimum_exact,
                      'minimum_exact_fraction': minimum_fraction, 'minimum_distinct_tags_each_end': 4, 'maximum_joint_tag_fraction': .5},
                  'sources': SOURCES, 'scope_evidence': evidence,
                  'expected_structure': 'P5-N5-[CCTATAT|TGCTATT]-RNA-derived insert-N10-P7',
                  'limitations': ['N5 is recorded as a random RNA-adapter tag, not a second validated UMI.',
                      'Observed N10 tags are not error-corrected ONT molecule counts.',
                      'Evidence is inferred independently within each assigned or unassigned scope; unsupported units remain untrimmed.',
                      'Barcode support does not identify RNA species, protein specificity, or a unique library kit.']}}
    statuses = collections.Counter()
    with gzip.open(output, 'wt') as handle, open_text(units_tsv, 'wt') as sidecar:
        writer = csv.DictWriter(sidecar, delimiter='\t', fieldnames=FIELDS, extrasaction='ignore')
        writer.writeheader()
        for name, sequence, quality in records:
            source = rows[name]
            scope, orientation = scope_of(source), source.get('orientation', '')
            length = len(sequence)
            normalized = reverse_complement(sequence) if orientation == 'p7_to_p5' else sequence
            normalized_quality = quality[::-1] if orientation == 'p7_to_p5' else quality
            row = {key: '' for key in FIELDS}
            row.update(unit_id=name, scope_id=scope, orientation=orientation, input_length=length,
                       p5_trim_bp=0, p7_trim_bp=0, umi_status=source.get('umi_status') or 'not_configured',
                       umi_sequence=source.get('umi_sequence') or '')
            if p5 or p7:
                row.update(status='manual_trim', p5_trim_bp=p5, p7_trim_bp=p7)
            elif bypass:
                row['status'] = 'bypassed_' + bypass
            elif not evidence[scope]['supported']:
                row['status'] = 'structure_not_supported'
            else:
                hit = find_barcode(normalized) if orientation in ('p5_to_p7', 'p7_to_p5') else None
                if hit is None:
                    row.update(status='uncertain_barcode', umi_status='structure_uncertain', umi_sequence='')
                elif set(normalized[:hit['barcode_start']] + normalized[-10:]) - set('ACGT'):
                    row.update(status='uncertain_random_tag', umi_status='structure_uncertain', umi_sequence='')
                else:
                    row.update(hit)
                    row.update(status='grip_trimmed', p5_trim_bp=hit['barcode_end'], p7_trim_bp=10,
                               p5_random_tag=normalized[:hit['barcode_start']], p7_umi=normalized[-10:],
                               p7_umi_quality=normalized_quality[-10:], umi_status='complete', umi_sequence=normalized[-10:])
                    counts['umi_complete_units'] += 1
            front, tail = (row['p7_trim_bp'], row['p5_trim_bp']) if orientation == 'p7_to_p5' else (row['p5_trim_bp'], row['p7_trim_bp'])
            row.update(physical_analysis_start=min(length, row['p5_trim_bp']),
                       physical_analysis_end=max(min(length, row['p5_trim_bp']), length - row['p7_trim_bp']),
                       ont_analysis_start=min(length, front), ont_analysis_end=max(min(length, front), length - tail))
            seq, qual = trim_insert(sequence, quality, orientation, row['p5_trim_bp'], row['p7_trim_bp'])
            row['analysis_length'] = len(seq)
            if not seq:
                counts['empty_after_trim'] += 1
                row['status'] += '_empty'
            else:
                counts['retained_units'] += 1
                handle.write(f'@{name}\n{seq}\n+\n{qual}\n')
            statuses[row['status']] += 1
            writer.writerow(row)
    counts['status_counts'] = dict(statuses)
    counts['units_tsv_sha256'] = hashlib.sha256(units_tsv.read_bytes()).hexdigest()
    Path(summary).write_text(json.dumps(counts, indent=2) + '\n')
    return counts


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('fastq', 'read-table', 'output', 'summary'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--units-tsv')
    p.add_argument('--structure-json')
    p.add_argument('--p5', type=int, default=0)
    p.add_argument('--p7', type=int, default=0)
    p.add_argument('--auto-grip', action=argparse.BooleanOptionalAction, default=True)
    p.add_argument('--minimum-scope-units', type=int, default=20)
    p.add_argument('--minimum-exact-barcodes', type=int, default=10)
    p.add_argument('--minimum-barcode-fraction', type=float, default=.5)
    a = p.parse_args()
    prepare(a.fastq, a.read_table, a.output, a.summary, a.p5, a.p7, a.units_tsv, a.auto_grip,
            a.structure_json, a.minimum_scope_units, a.minimum_exact_barcodes, a.minimum_barcode_fraction)
