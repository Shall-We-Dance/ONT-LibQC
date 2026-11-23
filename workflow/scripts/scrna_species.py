"""Run the whole-library species screen on barcode- and UMI-qualified cDNA."""
import csv
import gzip
import json
import subprocess
from pathlib import Path

from library_sequences import fastq_records
from scrna_chemistry import scope_label
from scrna_resources import digest
from species_composition import AMBIG, NONE, RAW, TOTAL, align_hits, classify


def read_panel(path, assembly):
    config = json.loads(path.read_text())
    indexes = config['indexes']
    aliases = config.get('assembly_aliases', {})
    reference = assembly if assembly in indexes else aliases.get(assembly)
    if reference not in indexes or len(indexes) < 2:
        raise ValueError('Species panel must contain the expected species and a competitor')
    prefixes = {}
    for name, value in indexes.items():
        if not name or name in {AMBIG, NONE, TOTAL, RAW} or not isinstance(value, str) or not value.strip():
            raise ValueError('Invalid reference in the species panel')
        prefix = (path.parent / value).resolve()
        if not any(all(Path(f'{prefix}.{part}.{suffix}').is_file() for part in ('1', '2', '3', '4', 'rev.1', 'rev.2'))
                   for suffix in ('bt2', 'bt2l')):
            raise ValueError(f'Incomplete Bowtie 2 species index for {name}')
        prefixes[name] = str(prefix)
    if len(set(prefixes.values())) != len(prefixes):
        raise ValueError('Species panel has duplicate index prefixes')
    return reference, prefixes


def make_queries(cdna, tagged, output):
    names, scopes, seen = [], [], set()
    with output.open('w') as handle:
        for name, sequence, _ in fastq_records([cdna]):
            if name in seen or name not in tagged:
                raise ValueError('Species FASTQ and qualified unit IDs disagree')
            seen.add(name)
            scope = tagged[name]['scope_id']
            scopes.append(scope_label(scope))
            handle.write(f'>{len(names)}\n{sequence}\n')
            names.append(name)
    if seen != tagged.keys():
        raise ValueError('Species FASTQ is missing qualified cDNA units')
    return names, scopes


def analyze(args, work, tagged):
    panel = args.species_panel or args.star_index.parent / 'species-panel.json'
    result = {'method': 'competitive-cdna-v1', 'host_assembly': args.assembly,
              'input_units': len(tagged), 'denominator': 'barcode_and_UMI_qualified_cDNA'}
    if not args.species_panel and not panel.exists():
        return dict(result, status='not_configured')
    reference, indexes = read_panel(panel, args.assembly)
    result.update(reference=reference, panel=list(indexes), panel_config_sha256=digest(panel))
    if not tagged:
        return dict(result, status='no_qualified_cdna')
    fasta = work / 'species_queries.fa'
    try:
        names, scopes = make_queries(work / 'cdna.fq.gz', tagged, fasta)
        hits = align_hits(fasta, len(names), indexes, args.bowtie2, args.threads)
    finally:
        fasta.unlink(missing_ok=True)
    rows = classify(scopes, hits, indexes, reference)['rows']
    # The qualified FASTQ contains only recovered units. RAW is a synthetic
    # zero row from the whole-library classifier, with no matching cohort here.
    result.update(rows=[r for r in rows if r['scope'] != RAW], status='complete',
                  mapper_version=subprocess.check_output([args.bowtie2, '--version'], text=True).splitlines()[0].split()[-1])
    with gzip.open(work / 'rna_families.tsv.gz', 'rt') as handle:
        families = {r['unit_id']: r['family'] for r in csv.DictReader(handle, delimiter='\t')}
    if families.keys() != tagged.keys():
        raise ValueError('RNA family and species screen inputs disagree')
    retained = [i for i, name in enumerate(names) if families[name] not in ('rRNA', 'mt-rRNA')]
    row = classify([scopes[i] for i in retained], [hits[i] for i in retained], indexes, reference)['rows'][0]
    result['non_rrna_row'] = dict(row, scope='After excluding identified rRNA')
    return result
