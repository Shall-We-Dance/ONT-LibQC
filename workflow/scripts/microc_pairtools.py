#!/usr/bin/env python3
"""Pairtools categories from overlapping, fastp-treated 150-bp insert ends."""
import argparse
from collections import Counter
import csv
import gzip
import html
import json
from pathlib import Path
import subprocess
import tempfile

CATEGORIES = ('UU', 'UR/RU', 'DD', 'WW', 'Other')
COLORS = ('#466e8a', '#71a8a0', '#c4866b', '#a197bb', '#c9ced3')
SETTINGS = dict(end_length_bp=150, minimum_insert_bp=50, overlap_allowed=True,
                fastp_detect_adapter_for_pe=True, aligner='bwa-mem2 mem -SP',
                min_mapq=1, walks_policy='mask', max_inter_align_gap=20,
                max_molecule_size=2000, duplicate_tolerance_bp=1,
                duplicate_method='max', duplicate_backend='scipy',
                duplicate_scope='sample', denominator='all parsed pairs, including unmapped and masked walks')


def fastq(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as handle:
        while header := handle.readline():
            seq, plus, qual = (handle.readline().rstrip() for _ in range(3))
            if not header.startswith('@') or not plus.startswith('+') or len(seq) != len(qual):
                raise ValueError('Invalid four-line insert FASTQ')
            yield header[1:].split()[0], seq, qual


def emit_ends(inserts, r1, r2):
    eligible, short = set(), set()
    with gzip.open(r1, 'wt') as left, gzip.open(r2, 'wt') as right:
        for unit, seq, qual in fastq(inserts):
            if unit in eligible or unit in short:
                raise ValueError('Repeated library unit identifier')
            if len(seq) < SETTINGS['minimum_insert_bp']:
                short.add(unit)
                continue
            eligible.add(unit)
            n = min(len(seq), SETTINGS['end_length_bp'])
            reverse = seq[-n:].translate(str.maketrans('ACGTRYMKBDHVacgtrymkbdhv', 'TGCAYRKMVHDBtgcayrkmvhdb'))[::-1]
            left.write(f'@{unit}\n{seq[:n]}\n+\n{qual[:n]}\n')
            right.write(f'@{unit}\n{reverse}\n+\n{qual[-n:][::-1]}\n')
    return eligible, short


def pair_types(path):
    result, columns = {}, None
    with open(path) as handle:
        for line in handle:
            if line.startswith('#columns:'):
                columns = line.split(':', 1)[1].split()
            elif line.startswith('#') or not line.strip():
                continue
            else:
                if columns is None:
                    raise ValueError('Pairs file lacks a columns header')
                row = dict(zip(columns, line.rstrip().split('\t')))
                unit = row['readID']
                if unit in result:
                    raise ValueError('Pairtools emitted more than one pair per unit')
                result[unit] = row['pair_type']
    if columns is None:
        raise ValueError('Pairs file lacks a columns header')
    return result


def category(pair_type):
    return 'UR/RU' if pair_type in {'UR', 'RU'} else pair_type if pair_type in CATEGORIES else 'Other'


def summarize(sample, before, after, eligible, short, assignments, names=None):
    if set(before) != set(after) or not set(before) <= set(eligible):
        raise ValueError('Pairtools changed the unit population or emitted an unknown unit')
    if not (set(eligible) | set(short)) <= set(assignments):
        raise ValueError('Insert units are missing from the reconstruction table')
    names = names or {}
    scope_ids = sorted({assignments[u]['scope_id'] for u in eligible | short
                        if assignments[u].get('assignment_status') == 'assigned'} | set(names))
    scopes = []
    for scope_id in ['total'] + scope_ids:
        def selected(u):
            return scope_id == 'total' or (assignments[u].get('assignment_status') == 'assigned'
                                          and assignments[u]['scope_id'] == scope_id)
        ids = {u for u in before if selected(u)}
        input_count = sum(selected(u) for u in eligible)
        scope = dict(scope_id=scope_id, scope_name='Total library' if scope_id == 'total' else names.get(scope_id, scope_id),
                     input_units=input_count, short_inserts_excluded=sum(selected(u) for u in short),
                     parsed_pairs=len(ids), fastp_excluded_units=input_count - len(ids))
        counts = Counter(category(after[u]) for u in ids)
        scope['after'] = dict(counts={c: counts[c] for c in CATEGORIES},
                              fractions={c: counts[c] / len(ids) if ids else None for c in CATEGORIES},
                              pair_types=dict(sorted(Counter(after[u] for u in ids).items())),
                              uu_to_ur_ru=counts['UU'] / counts['UR/RU'] if counts['UR/RU'] else None)
        scopes.append(scope)
    return dict(schema_version=2, sample=sample, settings=dict(SETTINGS), scopes=scopes)


def after_only(result):
    """Normalize saved reports to the final duplicate-marked results."""
    return dict(result, schema_version=2,
                scopes=[{key: value for key, value in scope.items() if key != 'before'}
                        for scope in result['scopes']])


def save_results(result, after, eligible, short, assignments, output):
    result = after_only(result)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    # Append suffixes: sample names may themselves contain periods.
    Path(str(output)+'.json').write_text(json.dumps(result, indent=2) + '\n')
    with Path(str(output)+'.summary.tsv').open('w') as handle:
        fields = ['sample', 'scope_id', 'stage', 'category', 'count', 'denominator', 'fraction',
                  'input_units', 'fastp_excluded_units', 'short_inserts_excluded']
        writer = csv.DictWriter(handle, fields, delimiter='\t'); writer.writeheader()
        for scope in result['scopes']:
            for cat in CATEGORIES:
                writer.writerow(dict(sample=result['sample'], scope_id=scope['scope_id'], stage='after',
                                     category=cat, count=scope['after']['counts'][cat], denominator=scope['parsed_pairs'],
                                     fraction=scope['after']['fractions'][cat], input_units=scope['input_units'],
                                     fastp_excluded_units=scope['fastp_excluded_units'],
                                     short_inserts_excluded=scope['short_inserts_excluded']))
    with gzip.open(str(output)+'.units.tsv.gz', 'wt') as handle:
        fields = ['unit_id', 'read_id', 'scope_id', 'assignment_status', 'status',
                  'after_pair_type', 'after_category']
        writer = csv.DictWriter(handle, fields, delimiter='\t'); writer.writeheader()
        for unit in sorted(eligible | short):
            row = {key: assignments[unit].get(key, '') for key in fields[:4]}
            row['unit_id'] = unit
            row['status'] = 'short_insert' if unit in short else 'parsed' if unit in after else 'fastp_excluded'
            row.update(after_pair_type=after.get(unit, ''),
                       after_category=category(after[unit]) if unit in after else '')
            writer.writerow(row)


def run_analysis(args):
    output = Path(args.output_prefix)
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.read_table, 'rt') as handle:
        assignments = {row['unit_id']: row for row in csv.DictReader(handle, delimiter='\t')}
    adapters = json.loads(Path(args.adapter_summary).read_text())
    names = {p['scope_id']: p.get('name', p['scope_id']) for p in adapters.get('index_pairing', {}).get('pairs', [])}
    with tempfile.TemporaryDirectory(prefix='pairtools-', dir=output.parent) as temp:
        work = Path(temp)
        r1, r2 = work/'R1.fastq.gz', work/'R2.fastq.gz'
        eligible, short = emit_ends(args.inserts, r1, r2)
        before, after, fastp_summary = {}, {}, {}
        def run(command, stdout=None):
            subprocess.run(list(map(str, command)), stdout=stdout, check=True)
        if eligible:
            clean1, clean2 = work/'clean_R1.fastq.gz', work/'clean_R2.fastq.gz'
            run([args.fastp, '-i', r1, '-I', r2, '-o', clean1, '-O', clean2,
                 '--thread', min(args.threads, 2), '--detect_adapter_for_pe',
                 '--html', work/'fastp.html', '--json', work/'fastp.json'])
            fastp_summary = json.loads((work/'fastp.json').read_text())['summary']
            passing = {u for u, _, _ in fastq(clean1)}
            if passing != {u for u, _, _ in fastq(clean2)} or not passing <= eligible:
                raise ValueError('fastp paired outputs do not match eligible units')
            if passing:
                sam = work/'aligned.sam'
                with sam.open('w') as handle:
                    run([args.bwa_mem2, 'mem', '-SP', '-t', args.threads, args.index, clean1, clean2], handle)
                run([args.pairtools, 'parse', '--chroms-path', args.chrom_sizes, '--drop-sam', '--add-columns', 'mapq',
                     '--walks-policy', 'mask', '--min-mapq', '1', '--max-inter-align-gap', '20',
                     '--max-molecule-size', '2000', '--nproc-in', '1', '--nproc-out', '1',
                     sam, '-o', work/'parsed.pairs'])
                run([args.pairtools, 'sort', '--nproc', min(args.threads, 2), '--tmpdir', work,
                     work/'parsed.pairs', '-o', work/'sorted.pairs'])
                marked = work/'all_marked.pairs'
                run([args.pairtools, 'dedup', '--mark-dups', '--max-mismatch', '1', '--method', 'max',
                     '--backend', 'scipy', '--output', marked, '--output-dups', marked, '--output-unmapped', marked,
                     '--nproc-in', '1', '--nproc-out', '1', work/'sorted.pairs'])
                before, after = pair_types(work/'parsed.pairs'), pair_types(marked)
                if set(before) != passing:
                    raise ValueError('Pairtools did not retain every fastp-passing pair')
        result = summarize(args.sample, before, after, eligible, short, assignments, names)
        result['fastp_summary'] = fastp_summary
        result['versions'] = {}
        for label, executable in [('pairtools', args.pairtools), ('fastp', args.fastp)]:
            proc = subprocess.run([executable, '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, check=True)
            result['versions'][label] = proc.stdout.strip()
        save_results(result, after, eligible, short, assignments, output)
    print(f"{args.sample}: {len(before):,} parsed pairs from {len(eligible):,} eligible inserts", flush=True)


def svg_categories(scope):
    esc = html.escape
    result = scope['after']
    parts = ["<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 840 164' role='img' "
             "style='font-family:DejaVu Sans,Arial,sans-serif' aria-label='Pairtools category fractions after duplicate marking'>",
             "<rect width='840' height='164' fill='white'/>",
             f"<text x='22' y='30' font-size='17' fill='#243d4b'>{esc(scope['scope_name'])}: Pairtools categories</text>",
             "<text x='22' y='80' font-size='13' fill='#445966'>After marking</text>"]
    x = 159
    if not scope['parsed_pairs']:
        parts.append(f"<text x='{x}' y='80' font-size='14'>No parsed pairs</text>")
    for cat, color in zip(CATEGORIES, COLORS):
        frac = result['fractions'][cat] or 0
        width = 650 * frac
        parts.append(f"<rect x='{x:.3f}' y='57' width='{width:.3f}' height='34' fill='{color}'><title>{esc(cat)}: "
                     f"{result['counts'][cat]:,} ({frac:.2%})</title></rect>")
        if frac >= .065:
            parts.append(f"<text x='{x+width/2:.3f}' y='79' text-anchor='middle' font-size='12' "
                         f"fill='{'white' if cat == 'UU' else '#172e3b'}'>{frac:.1%}</text>")
        x += width
    for i, (cat, color) in enumerate(zip(CATEGORIES, COLORS)):
        x = 159 + i*126
        parts.append(f"<rect x='{x}' y='112' width='13' height='13' fill='{color}'/>"
                     f"<text x='{x+19}' y='123' font-size='13'>{esc(cat)}</text>")
    parts.append(f"<text x='159' y='151' font-size='12' fill='#526675'>Denominator: all {scope['parsed_pairs']:,} parsed pairs</text></svg>")
    return ''.join(parts)


def report_html(result, assets, scopes, data_url='../qc/pairtools'):
    """Shared renderer for new analyses and existing portable examples."""
    esc = html.escape
    sample = result['sample']
    panels = []
    data_scopes = {s['scope_id']: s for s in result['scopes']}
    for requested in scopes:
        sid = requested.get('scope_id', 'total')
        if sid not in data_scopes:
            raise ValueError(f'Missing Pairtools scope: {sid}')
        scope = dict(data_scopes[sid], scope_name=requested.get('scope_name', data_scopes[sid]['scope_name']))
        name = f'{sample}.{sid}.pairtools.svg'
        Path(assets).mkdir(parents=True, exist_ok=True)
        (Path(assets)/name).write_text(svg_categories(scope))
        rows = []
        for cat in CATEGORIES:
            frac = scope['after']['fractions'][cat]
            cells = f"<td>{scope['after']['counts'][cat]:,}</td><td>{frac:.2%}</td>" if frac is not None else '<td>0</td><td>—</td>'
            rows.append(f'<tr><th>{esc(cat)}</th>{cells}</tr>')
        ratio = f"{scope['after']['uu_to_ur_ru']:.2f}" if scope['after']['uu_to_ur_ru'] is not None else '—'
        panels.append(f"<div class='pairtools-scope'><h3>{esc(scope['scope_name'])}</h3>"
                      f"<p>{scope['input_units']:,} eligible inserts → {scope['parsed_pairs']:,} parsed pairs; "
                      f"{scope['fastp_excluded_units']:,} excluded by fastp. {scope['short_inserts_excluded']:,} inserts &lt;50 bp excluded before end extraction.</p>"
                      f"<figure class='plot-card'><a href='assets/{name}' target='_blank' rel='noopener'>"
                      f"<img src='assets/{name}' alt='Pairtools UU, UR/RU, DD, WW and Other fractions' style='width:100%'></a></figure>"
                      "<div class='table-scroll'><table><tr><th>Category</th><th>Count</th><th>Fraction</th></tr>"
                      + ''.join(rows) + "</table></div>"
                      f"<p class='muted'>UU / (UR + RU): {ratio}.</p></div>")
    links = ' · '.join(f"<a href='{esc(data_url, quote=True)}/{sample}.{suffix}'>{label}</a>" for suffix, label in
                       [('summary.tsv','Category counts and fractions'),('units.tsv.gz','Per-unit pair types'),('json','Settings and metrics')])
    return ("<section class='report-section' id='pairtools'><span class='eyebrow'>PAIRED-END PARSING</span>"
            "<h2>Pairtools pair types</h2><p>Reconstructed ONT inserts are viewed as overlapping terminal reads of at most 150 bp, "
            "filtered with fastp, aligned with BWA-MEM2 and classified by Pairtools. Results are shown after duplicate marking. Percentages use every parsed pair, including unmapped pairs and masked walks.</p>"
            + ''.join(panels) + "<details><summary>Pair types and analysis settings</summary><p><strong>UU</strong>: unique–unique; "
            "<strong>UR/RU</strong>: rescued single-ligation pairs; <strong>DD</strong>: coordinate duplicates marked within the sample; "
            "<strong>WW</strong>: masked complex walks; <strong>Other</strong>: remaining types, including unmapped or non-unique ends.</p>"
            "<p>fastp paired-end adapter detection; BWA-MEM2 mem -SP; Pairtools parse with MAPQ ≥1, walks-policy mask, "
            "maximum inter-alignment gap 20 bp and maximum molecule size 2,000 bp; duplicate marking with a 1-bp tolerance, max method and scipy backend. "
            "Duplicates stay in the denominator. Sublibraries inherit marks from the full sample.</p>"
            "<p>These are parsing outcomes, not direct ligation counts or a calibrated library-quality score. DD depends on sampling depth and does not estimate PCR duplication or library complexity. "
            "The native valid-contact metrics above use their own end geometry and filters. "
            "<a href='https://pairtools.readthedocs.io/en/latest/parsing.html' target='_blank' rel='noopener'>Pairtools definitions</a>.</p></details>"
            f"<p>{links}</p></section>")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('sample', 'inserts', 'read-table', 'adapter-summary', 'index', 'chrom-sizes', 'output-prefix'):
        parser.add_argument('--'+flag, required=True)
    parser.add_argument('--bwa-mem2', default='bwa-mem2')
    parser.add_argument('--fastp', default='fastp')
    parser.add_argument('--pairtools', default='pairtools')
    parser.add_argument('--threads', type=int, default=2)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    run_analysis(args)


if __name__ == '__main__':
    main()
