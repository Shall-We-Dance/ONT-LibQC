"""Sequence-supported single-cell chemistry calls and read-oriented insert schematics."""
from base64 import b64encode
from collections import Counter
from html import escape

from report_sublibrary import percent
from scrna_chemistry import C4, C4_COMPACT, LABELS, PARSE, TENX, TENX_TSO, parse_offsets, scope_label


CSS = """
.sc-chem-group+.sc-chem-group {margin-top:1.3rem;border-top:1px solid var(--line);padding-top:1rem}
.sc-chem-group h3 {color:var(--ink);font-size:1.15rem;margin:0 0 .6rem}
.sc-chem-matches {max-width:42rem}
.sc-chem-matches td {white-space:nowrap}
.sc-chem-charts {display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1rem;margin-bottom:.65rem}
.sc-protocol-figure,.sc-partial-figure {display:flex;flex-direction:column;margin:0 0 1.4rem;max-width:42rem;break-inside:avoid}
.sc-chem-charts>figure {margin:0;max-width:none;padding:1rem 1.15rem;border:1px solid var(--line);border-radius:6px}
.sc-protocol-figure h3,.sc-partial-figure h3 {color:var(--ink);font-size:1.03rem;margin:0}
.sc-chem-charts figcaption {min-height:4em}
.sc-protocol-figure figcaption p,.sc-partial-figure figcaption p {color:var(--muted);margin:.25rem 0 .85rem;line-height:1.45}
.sc-protocol-figure svg,.sc-partial-figure svg {display:block;width:100%;max-width:100%}
.sc-protocol-figure .figure-link,.sc-partial-figure .figure-link {align-self:flex-end;margin-top:auto;padding-top:.5rem;font-size:.82rem}
.sc-site-detail {margin-bottom:1.3rem}
.sc-site-detail .sc-wrap {max-width:52rem}
.sc-site-detail .sc-table th {vertical-align:bottom}
.sc-molecule {margin:1rem 0;padding:1rem;border:1px solid var(--line);border-radius:6px;break-inside:avoid}
.sc-molecule figcaption {font-weight:600;color:var(--ink)}
.sc-molecule h4 {margin:.7rem 0 .15rem;color:var(--muted);font-size:.88rem;font-weight:500}
.sc-molecule .library-diagram {grid-template-columns:max-content repeat(var(--blocks),minmax(78px,108px) max-content)}
.sc-molecule .library-diagram-note {margin-bottom:0}
.sc-chem-details .sc-table :is(td,th):first-child {width:auto;white-space:normal}
.sc-chem-details h4 {margin-bottom:.25rem}
@media screen and (max-width:680px) {
.sc-molecule {padding:.7rem}
.sc-molecule .library-diagram {justify-content:start}
}
@media (max-width:840px) {
.sc-chem-charts {grid-template-columns:minmax(0,1fr)}
.sc-chem-charts figcaption {min-height:0}
.sc-chem-charts>figure {padding:.75rem .85rem}
}
@media screen and (max-width:400px) {
.sc-chem-charts svg {font-size:12px}
}
@media print {
.sc-chem-charts {display:block}
.sc-chem-charts>figure {margin-bottom:.6rem}
.sc-chem-charts figcaption {min-height:0}
.sc-protocol-figure .figure-link,.sc-partial-figure .figure-link {display:none}
.sc-molecule {-webkit-print-color-adjust:exact;print-color-adjust:exact}
.sc-molecule .library-diagram-scroll {overflow:visible}
.sc-molecule .library-diagram {min-width:0;grid-template-columns:max-content repeat(var(--blocks),minmax(0,1fr) max-content);gap:2px}
.sc-molecule .library-block {min-width:0;padding:3px 2px}
}
"""


def _segments(kit):
    # All schematics start at the barcode-bearing sequencing end. The Parse
    # primer side lies between its round-1 barcode and the transcript.
    if kit == PARSE:
        blocks = [('UMI', 'umi', 10)]
        last = 10
        for round_no, (start, end) in zip((3, 2, 1), parse_offsets()):
            if start > last:
                blocks.append((f'Linker {3-round_no}', 'anchor', start - last))
            blocks.append((f'BC{round_no}', 'barcode', end - start))
            last = end
        note = 'Read 2 carries the three cell-barcode rounds. Read 1 starts at the cDNA end. A poly(T) tract may lie between BC1 and cDNA.'
    elif kit in TENX:
        _, umi, anchor = TENX[kit]
        blocks = [('Cell barcode', 'barcode', 16), ('UMI', 'umi', umi)]
        blocks.append(('Poly(T)', 'anchor', '≥6 bp') if anchor == 'polyT' else
                      ('TSO', 'anchor', len(TENX_TSO)))
        note = 'Read 1 carries the cell barcode and UMI. Read 2 starts at the cDNA end.'
    elif kit in (C4, C4_COMPACT):
        if kit == C4:
            blocks = [('Linker 1', 'anchor', 6), ('Barcode 1', 'barcode', 10), ('Linker 2', 'anchor', 6),
                      ('Barcode 2', 'barcode', 10), ('Linker 3', 'anchor', 5)]
        else:
            blocks = [('Barcode 1', 'barcode', 10), ('Barcode 2', 'barcode', 10)]
        blocks.extend([('UMI', 'umi', 10), ('Poly(T)', 'anchor', '≥6 bp')])
        note = 'Read 1 carries both cell-barcode segments and the UMI. Read 2 starts at the cDNA end.'
    else:
        raise ValueError('Unknown single-cell chemistry')
    return blocks + [('cDNA', 'cdna', 'Variable length')], note


def _diagram(title, blocks, read_span):
    """Use the same component and arrow styling as Chromatin, in P5-to-P7 orientation."""
    first, last = read_span
    parts = [f'<h4>{escape(title)}</h4><div class="library-diagram-scroll" tabindex="0" role="region" '
             f'aria-label="{escape(title)}"><div class="library-diagram" style="--blocks:{len(blocks)}">',
             '<span class="library-end">5′</span>']
    for i, (label, kind, content) in enumerate(blocks):
        if i:
            parts.append('<span class="library-connector" aria-hidden="true">—</span>')
        parts.append(f'<span class="library-block {kind}"><strong>{label}</strong>{content}</span>')
    parts.extend(['<span class="library-end">3′</span>',
                  f'<div class="library-read-directions" style="grid-column:{2 + 2 * first}/{3 + 2 * last}">'
                  '<span class="library-read forward" aria-label="Read 1, forward">R1<i aria-hidden="true"></i></span>'
                  '<span class="library-read reverse" aria-label="Read 2, reverse">R2<i aria-hidden="true"></i></span>'
                  '</div></div></div>'])
    return ''.join(parts)


def _index_text(seq):
    if seq is None:
        return 'Unresolved'
    return seq or 'No index'


def _index_content(pairs, key):
    values = {pair.get(key) for pair in pairs.values()}
    if len(values) == 1:
        seq = values.pop()
        if seq:
            return f'<code class="library-sequence">{escape(seq)}</code><small>{len(seq)} bp</small>'
        return f'<small>{_index_text(seq)}</small>'
    return '<small>Varies by library</small>'


def structure(kit, pairs=None, adapter_family=None):
    segments, note = _segments(kit)
    # _segments follows the barcode read for extraction; Parse reads that end from R2.
    if kit == PARSE:
        segments.reverse()
    parts = ['<figure class="sc-molecule"><figcaption>Library structure</figcaption>']
    if pairs:
        family = f'<small>{escape(adapter_family)}</small>' if adapter_family else ''
        outer = [('P5', 'handle', '<small>Expected handle</small>'),
                 ('i5', 'index', _index_content(pairs, 'i5')),
                 ('Adapter', 'adapter', '<small>i5 side</small>' + family),
                 ('Insert', 'insert', '<small>Barcoded cDNA</small>'),
                 ('Adapter', 'adapter', '<small>i7 side</small>' + family),
                 ('i7', 'index', _index_content(pairs, 'i7')),
                 ('P7', 'handle', '<small>Expected handle</small>')]
        scope = next(iter(pairs)) if len(pairs) == 1 else f'{len(pairs)} libraries'
        parts.append(_diagram(f'Physical library · {scope}', outer, (3, 3)))
    classes = {'cdna': 'insert', 'umi': 'umi', 'anchor': 'adapter', 'barcode': 'index'}
    blocks = [(label, classes[kind], f'<small>{str(size) + " bp" if isinstance(size, int) else size}</small>')
              for label, kind, size in segments]
    parts.extend([_diagram('Reconstructed insert', blocks, (0, len(blocks) - 1)),
                  f'<p class="library-diagram-note">{note}</p></figure>'])
    return ''.join(parts)


_UNRESOLVED = {'unresolved': 'No supported layout', 'insufficient_evidence': 'Too few supporting reads',
               'mixed_or_ambiguous': 'Multiple supported layouts'}

_PROTOCOLS = (('parse', 'Parse Evercode WT', '#226e63'), ('10x', '10x Chromium', '#4787a7'),
              ('mgi', 'MGI DNBelab C4', '#8364a2'), ('multiple', 'Multiple families', '#946b3f'),
              ('none', 'No tested layout matched', '#93a3b1'))

_PARTIAL_PATTERNS = (('barcode', 'Barcode mismatch', '#447897'), ('linker', 'Linker mismatch', '#947349'),
                     ('shifted', 'Shifted cassette or spacing', '#648692'),
                     ('short_end', 'Possible Read 2 truncation', '#8a72a1'),
                     ('mixed', 'Other partial structure', '#639187'),
                     ('unresolved', 'No clear Parse pattern', '#93a3b1'))
_SITES = (('bc3', 'BC3'), ('linker1', 'Linker 1'), ('bc2', 'BC2'), ('linker2', 'Linker 2'), ('bc1', 'BC1'))


def _count_figure(kinds, counts, title, subtitle, filename, css_class):
    total = sum(counts.values())
    figure_id = filename.removesuffix('.svg')
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="640" height="{36 * len(kinds)}" role="img" '
           f'aria-labelledby="{figure_id}-title {figure_id}-desc" font-family="Arial, sans-serif" font-size="14" fill="#24384c">',
           f'<title id="{figure_id}-title">{escape(title)}</title>',
           f'<desc id="{figure_id}-desc">{escape(subtitle)} ' +
           '; '.join(f'{label}: {counts[key]:,} ({escape(percent(counts[key], total))})' for key, label, _ in kinds) + '.</desc>']
    for i, (key, label, color) in enumerate(kinds):
        count = counts[key]
        width = 100 * count / total
        y = 14 + 36 * i
        svg.extend([f'<text y="{y}">{label}</text>',
                    f'<text x="100%" y="{y}" text-anchor="end">{count:,} ({escape(percent(count, total))})</text>',
                    f'<rect y="{y+6}" width="100%" height="6" rx="2" fill="#edf1f5"/>',
                    f'<rect y="{y+6}" width="{width:.8f}%" height="6" rx="2" fill="{color}"/>'])
        if 0 < width < 1:
            svg.append(f'<line x1="{width:.8f}%" x2="{width:.8f}%" y1="{y+6}" y2="{y+12}" stroke="{color}" stroke-width="2"/>')
    graphic = ''.join(svg) + '</svg>'
    uri = 'data:image/svg+xml;base64,' + b64encode(graphic.encode()).decode('ascii')
    return (f'<figure class="{css_class}"><figcaption><h3>{escape(title)}</h3>'
            f'<p>{escape(subtitle)}</p></figcaption>{graphic}'
            f'<a class="figure-link" href="{uri}" download="{filename}">Download SVG</a></figure>')


def _protocol_figure(scopes):
    # Legacy per-kit counts can count the same unit more than once.
    if any('protocol_support' not in scope for scope in scopes.values()):
        return ''
    counts = {key: 0 for key, _, _ in _PROTOCOLS}
    for scope in scopes.values():
        support = scope['protocol_support']
        if (support.keys() != counts.keys() or any(n < 0 for n in support.values())
                or sum(support.values()) != scope['units']):
            raise ValueError('Protocol support counts must partition all tested library units')
        for key, n in support.items():
            counts[key] += n
    total = sum(counts.values())
    if not total:
        return ''
    return _count_figure(_PROTOCOLS, counts, 'Support by library preparation',
                         f'Structure matches across {total:,} reconstructed non-dimer units.',
                         'chemistry_support.svg', 'sc-protocol-figure')


def _parse_evidence(scopes, table):
    diagnostics = []
    for scope in scopes.values():
        if not scope.get('protocol_support', {}).get('none'):
            continue
        result = scope.get('unmatched_parse')
        if result is None:
            return '', ''
        if (result['units'] != scope['protocol_support']['none']
                or sum(result['categories'].values()) != result['assessed']
                or result['partial'] != result['assessed'] - result['categories']['unresolved']):
            raise ValueError('Parse diagnostics must account for the unmatched units in each scope')
        diagnostics.append(result)
    assessed = sum(d['assessed'] for d in diagnostics)
    if not assessed:
        return '', ''
    units = sum(d['units'] for d in diagnostics)
    partial = sum(d['partial'] for d in diagnostics)
    counts = Counter()
    for d in diagnostics:
        counts.update(d['categories'])
    patterns = _PARTIAL_PATTERNS
    if units > assessed:
        patterns += (('untested', 'No linker consensus', '#cad2d9'),)
        counts['untested'] = units - assessed
    subtitle = (f'{partial:,} of {units:,} unmatched units ({percent(partial, units)}) '
                'retain recognizable cassette features.')
    figure = _count_figure(patterns, counts, 'Parse structure in unmatched units', subtitle,
                           'unmatched_parse.svg', 'sc-partial-figure')
    rows = []
    for site, name in _SITES:
        rows.append([name, *(sum(d['sites'][site][key] for d in diagnostics)
                             for key in ('exact', 'near', 'mismatch', 'incomplete'))])
    short = sum(d['shorter_than_58'] for d in diagnostics)
    short_supported = sum(d['short_end_patterns']['short_insert'] for d in diagnostics)
    early = sum(d['short_end_patterns']['early_anchors'] for d in diagnostics)
    opposite = sum(d['opposite_end_partial'] for d in diagnostics)
    detail = [
        '<details class="sc-site-detail"><summary>Barcode, linker and end evidence</summary>',
        f'<p>Each row covers {assessed:,} units with a known linker pair. Sites are checked at their expected Read 2 positions. '
        'A one-base barcode difference must point to a single listed barcode.</p>',
        table(['Read 2 site', 'Exact', 'One base away', 'Mismatch', 'Incomplete'], rows),
        '<p>Matching linkers can pinpoint a barcode mismatch. Matching barcodes can pinpoint a linker mismatch. '
        'A change in spacing can move several sites out of position at once. Other partial structures have recognizable '
        'features but no single clear explanation.</p>',
        f'<p><strong>End evidence.</strong> {short:,} of {units:,} inserts were shorter than the 58-base cassette. '
        f'Among the {counts["short_end"]:,} possible truncations, {short_supported:,} showed cassette landmarks in a short insert '
        f'and {early:,} showed landmarks closer to the Read 2 boundary than expected. These patterns can also arise during reconstruction.</p>',
    ]
    boundaries = [d['raw_r2_boundary'] for d in diagnostics if d.get('raw_r2_boundary', {}).get('measured')]
    if boundaries:
        measured = sum(b['measured'] for b in boundaries)
        continued = sum(b['at_least_1'] for b in boundaries)
        min_bp = min(b['min_bp'] for b in boundaries)
        extent = (f' The shortest extension was {min_bp:,} bp.' if min_bp else
                  f' {sum(b["at_least_12"] for b in boundaries):,} extended at least 12 bp.')
        detail.append(f'<p><strong>Raw-read boundary.</strong> In {continued:,} of {measured:,} units checked, the raw ONT '
                      f'read continued beyond the reconstructed Read 2 boundary.{extent}</p>')
    detail.append(
        f'<p><strong>Opposite-end check.</strong> The same search found partial structure in {opposite:,} of '
        f'{assessed:,} units ({escape(percent(opposite, assessed))}) at the other end.</p></details>')
    return figure, ''.join(detail)


def _scope_results(scopes, table, physical, family, shown_kits, adapter_note):
    groups, detail, unresolved = {}, [], []
    for scope, data in scopes.items():
        kit = data.get('kit')
        if kit in LABELS:
            groups.setdefault(kit, []).append((scope, data))
        else:
            unresolved.append([scope_label(scope), _UNRESOLVED.get(data.get('status'), 'Unresolved'), data['units']])
        rows = []
        for k, ev in sorted(data['evidence'].items(), key=lambda item: (-item[1]['supported_units'], item[0])):
            rows.append([LABELS.get(k, k), ev['exact_barcode_units'], ev['opposite_end_exact'], ev['supported_units'],
                         ev['parents'], ev['distinct_barcodes']])
        detail.extend([f'<h4 data-sublibrary="{escape(scope, quote=True)}">{escape(scope_label(scope))}</h4>',
                       table(['Layout', 'Exact barcode', 'Opposite end', 'Barcode + anchors', 'Distinct ONT reads',
                              'Distinct barcodes'], rows)])
    content = []
    for kit, members in groups.items():
        parts = [f'<div class="sc-chem-group"><h3>{escape(LABELS[kit])}</h3>']
        rule = ('A match requires all three cell barcodes and both linkers at their expected positions in Read 2.'
                if kit == PARSE else 'A match requires the expected cell barcode and flanking sequence.')
        parts.append(f'<p>{rule} Each row counts non-dimer units in that sample-index group.</p>')
        matches = []
        for scope, data in members:
            n, total = data['evidence'][kit]['supported_units'], data['units']
            percentage = f'{100 * n / total:.1f}%' if total else '—'
            matches.append([scope_label(scope), f'{n:,} / {total:,}', percentage])
        parts.extend(['<div class="sc-chem-matches">', table(['Index group', 'Matched / tested', 'Match rate'], matches), '</div>'])
        if kit not in shown_kits:
            pairs = {scope: physical[scope] for scope, _ in members if scope in physical}
            parts.append(structure(kit, pairs, family))
        parts.append('</div>')
        content.extend(parts)
    if unresolved:
        content.append(table(['Library', 'Detection result', 'Units screened'], unresolved))
    content += ['<details class="sc-chem-details"><summary>Detection evidence</summary>']
    if adapter_note:
        content.append(adapter_note)
    if any(d.get('linkers') and d['kit'] == PARSE for d in scopes.values()):
        content.append(table(['Physical library', 'Read 2 linker 1', 'Read 2 linker 2'],
                             [[scope_label(scope), *d['linkers']] for scope, d in scopes.items() if d.get('linkers') and d['kit'] == PARSE]))
    return ''.join(content + detail + ['</details>'])


def render_chemistry(summary, table):
    scopes = summary['scopes']
    if not scopes:
        return '<p>No complete non-dimer units were recovered.</p>'
    content = _protocol_figure(scopes)
    if content:
        partial, site_detail = _parse_evidence(scopes, table)
        if partial:
            content = '<div class="sc-chem-charts" id="chemistry-evidence">' + content + partial + '</div>' + site_detail
    physical = summary.get('physical_indexes', {})
    adapters = {k: n for k, n in summary.get('adapter_family_counts', {}).items() if n}
    labels = {'truseq': 'TruSeq', 'nextera': 'Nextera'}
    family = labels.get(next(iter(adapters))) if len(adapters) == 1 else None
    adapter_note = '<p>Whole-library outer adapters · ' + ', '.join(
        f'{escape(labels.get(k, k))} ({n:,} units)' for k, n in sorted(adapters.items())) + '</p>' if adapters else ''
    primary = {name: group for name, group in scopes.items() if name != 'unassigned'}
    if primary:
        content += _scope_results(primary, table, physical, family, (), adapter_note)
    if 'unassigned' in scopes:
        shown = {group.get('kit') for group in primary.values()}
        content += '<details class="sc-unassigned"><summary>Unassigned index</summary>'
        content += _scope_results({'unassigned': scopes['unassigned']}, table, physical, family, shown,
                                  '' if primary else adapter_note) + '</details>'
    return content
