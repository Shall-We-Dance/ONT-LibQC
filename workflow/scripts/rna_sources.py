#!/usr/bin/env python3
"""Conservative per-insert rRNA evidence from a small public reference panel.

Reference-organism labels are candidates, never a taxonomic census or a
contamination diagnosis. Panel and host-annotation counts overlap.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

RESOURCE = Path(__file__).resolve().parents[1] / 'resources' / 'rna_sources'
FIELDS = ['unit_id', 'source_class', 'reference_name', 'matched_fraction',
          'identity', 'aligned_bases', 'taxon', 'taxon_ambiguous', 'status']
THRESHOLDS = {'minimum_query_coverage': 0.8, 'minimum_identity': 0.85,
              'minimum_aligned_matching_bases': 30, 'near_best_score_fraction': 0.95}
LIMITATIONS = [
    'Small representative rRNA panel; unassigned does not mean non-ribosomal RNA.',
    'Conserved rRNA may match several organisms. Taxon labels identify reference candidates only.',
    'Panel evidence and host-genome annotation overlap and their percentages must not be added.',
]


def open_text(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def input_ids(path):
    ids = []
    with open_text(path) as handle:
        while True:
            name = handle.readline().rstrip()
            if not name:
                break
            sequence = handle.readline().rstrip()
            plus = handle.readline().rstrip()
            quality = handle.readline().rstrip()
            if not name.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise ValueError('Malformed RNA FASTQ')
            ids.append(name[1:].split()[0])
    if len(set(ids)) != len(ids):
        raise ValueError('RNA source evidence requires unique unit identifiers')
    return ids


def load_panel(fasta, manifest):
    metadata = json.loads(Path(manifest).read_text())
    digest = hashlib.sha256(Path(fasta).read_bytes()).hexdigest()
    if digest != metadata.get('fasta_sha256'):
        raise ValueError('RNA source panel FASTA does not match its provenance manifest')
    references = metadata.get('references', {})
    if not references or any(r.get('source_class') != 'rRNA' or not r.get('taxon') for r in references.values()):
        raise ValueError('RNA source panel lacks RNA class or reference-organism metadata')
    return metadata


def parse_paf(lines, references, thresholds=None):
    """Require long, accurate evidence before comparing competing references."""
    thresholds = thresholds or THRESHOLDS
    hits = {}
    for line in lines:
        if not line.strip():
            continue
        fields = line.rstrip().split('\t')
        if len(fields) < 12:
            raise ValueError('Malformed minimap2 PAF')
        unit, target = fields[0], fields[5]
        if target not in references:
            raise ValueError('PAF target is absent from RNA source provenance')
        query_length = int(fields[1])
        matches, block_length = int(fields[9]), int(fields[10])
        coverage = (int(fields[3]) - int(fields[2])) / query_length if query_length else 0
        identity = matches / block_length if block_length else 0
        if (coverage < thresholds['minimum_query_coverage'] or
                identity < thresholds['minimum_identity'] or
                matches < thresholds['minimum_aligned_matching_bases']):
            continue
        tags = {tag.split(':', 2)[0]: tag.split(':', 2)[2] for tag in fields[12:] if tag.count(':') >= 2}
        score = int(tags.get('AS', matches))
        candidate = {'reference_name': target, 'matched_fraction': coverage,
                     'identity': identity, 'aligned_bases': matches, 'score': score}
        rank = (-score, -coverage, -identity, target)
        state = hits.setdefault(unit, {'best': candidate, 'rank': rank, 'taxon_scores': {}})
        if rank < state['rank']:
            state['best'], state['rank'] = candidate, rank
        taxon = references[target]['taxon']
        state['taxon_scores'][taxon] = max(score, state['taxon_scores'].get(taxon, score))
    assignments = {}
    for unit, state in hits.items():
        best = state['best']
        taxa = sorted(taxon for taxon, score in state['taxon_scores'].items()
                      if score >= best['score'] * thresholds['near_best_score_fraction'])
        assignments[unit] = {
            'unit_id': unit, 'source_class': references[best['reference_name']]['source_class'],
            'reference_name': best['reference_name'],
            'matched_fraction': round(best['matched_fraction'], 6),
            'identity': round(best['identity'], 6), 'aligned_bases': best['aligned_bases'],
            'taxon': ','.join(taxa), 'taxon_ambiguous': int(len(taxa) > 1), 'status': 'assigned',
        }
    return assignments


def write_results(ids, assignments, output, summary_path, metadata=None, unavailable_reason=None):
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    counts, taxa = collections.Counter(), collections.Counter()
    ambiguous = 0
    with open_text(output, 'wt') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        for unit in ids:
            row = assignments.get(unit)
            if row is None:
                row = {'unit_id': unit, 'status': 'not_available' if unavailable_reason else 'unassigned',
                       'taxon_ambiguous': 0}
            else:
                counts[row['source_class']] += 1
                taxa[row['taxon']] += 1
                ambiguous += int(row['taxon_ambiguous'])
            writer.writerow(row)
    assigned = sum(counts.values())
    summary = {
        'schema_version': 1, 'status': 'not_available' if unavailable_reason else 'complete',
        'unavailable_reason': unavailable_reason, 'denominator': 'analysis_input_inserts',
        'input_units': len(ids), 'assigned_units': assigned,
        'unassigned_units': len(ids) - assigned if not unavailable_reason else 0,
        'not_available_units': len(ids) if unavailable_reason else 0,
        'source_counts': dict(counts), 'candidate_taxon_counts': dict(taxa),
        'taxon_ambiguous_units': ambiguous, 'thresholds': THRESHOLDS,
        'panel_version': (metadata or {}).get('panel_version'),
        'panel_fasta_sha256': (metadata or {}).get('fasta_sha256'),
        'panel_reference_count': len((metadata or {}).get('references', {})),
        'reference_source_urls': sorted({r.get('source_url') for r in (metadata or {}).get('references', {}).values() if r.get('source_url')}),
        'mapper_version': (metadata or {}).get('_mapper_version'),
        'method': 'minimap2 nucleotide alignment; query coverage and identity from base-aligned PAF',
        'limitations': LIMITATIONS,
    }
    Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
    Path(summary_path).write_text(json.dumps(summary, indent=2) + '\n')
    return summary


def run(args):
    ids = input_ids(args.input)
    metadata, reason, assignments = None, None, {}
    if not Path(args.reference_fasta).is_file() or not Path(args.reference_manifest).is_file():
        reason = 'reference_panel_unavailable'
    elif not shutil.which(args.minimap2):
        reason = 'minimap2_unavailable'
    else:
        metadata = load_panel(args.reference_fasta, args.reference_manifest)
        metadata['_mapper_version'] = subprocess.check_output([args.minimap2, '--version'], text=True).strip()
        if ids:
            command = [args.minimap2, '-x', 'map-ont', '-k', '9', '-w', '5', '-m', '20',
                       '-s', '20', '-n', '2', '-N', '100', '--secondary=yes', '-c', '--eqx',
                       '-t', str(max(1, min(args.threads, 4))), args.reference_fasta, args.input]
            with tempfile.TemporaryFile(mode='w+t') as paf:
                subprocess.run(command, stdout=paf, check=True)
                paf.seek(0)
                assignments = parse_paf(paf, metadata['references'])
            if not set(assignments).issubset(ids):
                raise ValueError('RNA source mapper returned unknown input identifiers')
    return write_results(ids, assignments, args.output, args.summary, metadata, reason)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--reference-fasta', default=str(RESOURCE / 'panel.fa'))
    parser.add_argument('--reference-manifest', default=str(RESOURCE / 'panel.json'))
    parser.add_argument('--output', required=True)
    parser.add_argument('--summary', required=True)
    parser.add_argument('--minimap2', default='minimap2')
    parser.add_argument('--threads', type=int, default=1)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
