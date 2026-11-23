"""Sublibrary comparisons from the retained Micro-C analysis outputs."""
import csv
import gzip
import html
import json
import math
import statistics
from pathlib import Path

from microc_ligation import (CATEGORIES as BLOCK_CATEGORIES, CATEGORY_COLORS as BLOCK_COLORS,
                             summarize as summarize_ligation)
from microc_pairtools import CATEGORIES as PAIR_CATEGORIES, COLORS as PAIR_COLORS
from qc_plot_metrics import size_trace_data
from report_palette import bar_ink, scope_color
from report_svg import svg_start, text


BLOCK_LABELS = ('1 block', '2', '3', '4', '≥5', 'Unresolved')


def _scope_label(scope_id):
    if scope_id == 'total':
        return 'Total library'
    return scope_id if len(scope_id) <= 12 else scope_id[:6] + '…' + scope_id[-5:]


def _scope_group(scope):
    sid = html.escape(scope['scope_id'], quote=True)
    label = html.escape(scope.get('scope_name') or scope['scope_id'])
    if len(scope['scope_id']) > 12 and scope.get('scope_name') not in (None, '', scope['scope_id']):
        label += f' ({sid})'
    return f"<g data-scope='{sid}'><title>{label}</title>"


def _stacked_segments(sid, categories, labels, colors, counts, denominator, left, right, y, height=22,
                      label_size=12):
    svg = []
    x = left
    for key, label, color in zip(categories, labels, colors):
        count = counts[key]
        if not count:
            continue
        fraction = count / denominator
        width = fraction * (right - left)
        tooltip = html.escape(f'{sid} · {label} · {count:,} / {denominator:,} ({fraction:.2%})')
        svg.append(f"<rect x='{x:.3f}' y='{y-height/2:g}' width='{width:.3f}' height='{height}' "
                   f"fill='{color}' data-category='{html.escape(key, quote=True)}' data-count='{count}'><title>{tooltip}</title></rect>")
        if width >= 49:
            svg.append(text(x + width / 2, y + 4, f'{fraction:.1%}', size=label_size, color=bar_ink(color)))
        x += width
    return svg


def _stacked_chart(scopes, title, description, categories, labels, colors, rows,
                   denominator_label, statistic_label):
    """A missing measurement differs from a measured zero."""
    left = 222 if any(len(scope['scope_id']) > 5 for scope in scopes) else 110
    right = 588
    height = 131 + 35 * len(scopes)
    svg = svg_start(840, height, title)
    svg += [text(64, 49, description, 'start', 13),
            text(702, 75, denominator_label, 'end', 12),
            text(819, 75, statistic_label, 'end', 12)]
    for i, scope in enumerate(scopes):
        y = 97 + 35 * i
        sid = scope['scope_id']
        svg += [_scope_group(scope),
                f"<circle cx='28' cy='{y}' r='4' fill='{scope_color(sid)}'/>",
                text(40, y + 4, _scope_label(sid), 'start', 13, '#263541')]
        row = rows.get(sid)
        if row is None:
            svg += [text(left, y + 4, 'Unavailable', 'start'), '</g>']
            continue
        denominator, counts, statistic = row
        if sum(counts[key] for key in categories) != denominator:
            raise ValueError(f'{title}: category counts do not match the denominator for {sid}')
        svg += [text(702, y + 4, f'{denominator:,}', 'end', 14, '#263541'),
                text(819, y + 4, statistic, 'end', 14, '#263541')]
        if not denominator:
            svg.append(text(left, y + 4, 'No observations', 'start'))
        svg.extend(_stacked_segments(sid, categories, labels, colors, counts, denominator, left, right, y))
        svg.append('</g>')
    for i, (label, color) in enumerate(zip(labels, colors)):
        x, y = left + i * (right - left) / len(labels), height - 20
        svg += [f"<rect x='{x:.1f}' y='{y-10}' width='11' height='11' fill='{color}'/>",
                text(x + 16, y, label, 'start', 12, '#263541')]
    return ''.join(svg) + '</svg>'


def _pairtools_chart(scopes, result):
    rows = {}
    for stats in result['scopes']:
        counts = stats['after']['counts']
        ratio = f"{counts['UU'] / counts['UR/RU']:.2f}" if counts['UR/RU'] else '—'
        rows[stats['scope_id']] = (stats['parsed_pairs'], counts, ratio)
    return _stacked_chart(scopes, 'Pairtools pair types',
                          'Share of all parsed pairs after sample-level duplicate marking',
                          PAIR_CATEGORIES, PAIR_CATEGORIES, PAIR_COLORS, rows,
                          'Parsed pairs', 'UU / (UR + RU)')


def _length_box(values):
    """Matplotlib’s default inclusive quartiles and 1.5-IQR whiskers, without fliers."""
    if len(values) == 1:
        q1 = median = q3 = values[0]
    else:
        q1, median, q3 = statistics.quantiles(values, n=4, method='inclusive')
    extent = 1.5 * (q3 - q1)
    low = min(v for v in values if v >= q1 - extent)
    high = max(v for v in values if v <= q3 + extent)
    return dict(count=len(values), whislo=min(q1, low), q1=q1, med=median,
                q3=q3, whishi=max(q3, high))


def _length_boxes(scopes, units):
    grouped = {s['scope_id']: {category: [] for category in BLOCK_CATEGORIES} for s in scopes}
    with gzip.open(units, 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            sid, category, length = row['scope_id'], row['category'], int(row['insert_bp'])
            if 'total' in grouped:
                grouped['total'][category].append(length)
            if sid != 'total' and sid in grouped:
                grouped[sid][category].append(length)
    return {(sid, category): _length_box(values)
            for sid, categories in grouped.items() for category, values in categories.items() if values}


def _length_box_svg(scope, category, stats, x, width, y, color):
    sid = html.escape(scope['scope_id'], quote=True)
    label = html.escape(scope.get('scope_name') or scope['scope_id'])
    attributes = ' '.join(f"data-{key}='{value}'" for key, value in stats.items())
    tip = (f"{label}, {category} blocks · n = {stats['count']:,} · median {stats['med']:g} bp · "
           f"quartiles {stats['q1']:g}–{stats['q3']:g} bp · whiskers {stats['whislo']:g}–{stats['whishi']:g} bp")
    half = width / 2
    lo, q1, med, q3, hi = (y(stats[key]) for key in ('whislo', 'q1', 'med', 'q3', 'whishi'))
    return (f"<g data-boxplot='insert-length' data-scope='{sid}' data-category='{category}' {attributes}>"
            f"<title>{tip}</title>"
            f"<path d='M{x:.2f},{lo:.2f} V{q1:.2f} M{x:.2f},{q3:.2f} V{hi:.2f} M{x-half/2:.2f},{lo:.2f} h{half:.2f} "
            f"M{x-half/2:.2f},{hi:.2f} h{half:.2f}' fill='none' stroke='#000000' stroke-width='1' stroke-linecap='square'/>"
            f"<rect x='{x-half:.2f}' y='{q3:.2f}' width='{width:.2f}' height='{q1-q3:.2f}' "
            f"fill='{color}' stroke='#000000' stroke-width='1'/>"
            f"<path d='M{x-half:.2f},{med:.2f} h{width:.2f}' stroke='#ff7f0e' stroke-width='1.45'/></g>")


def _ligation_lengths_chart(scopes, units):
    width, margin = 500, 18
    title = 'Insert length by block count'
    subtitle = text(margin, 49, 'Adapter-bounded insert length (bp)', 'start', 14)
    boxes = _length_boxes(scopes, units)
    if not boxes:
        return ''.join(svg_start(width, 90, title, title_x=margin) + [
            subtitle, text(margin, 75, 'No units in these libraries', 'start', 13), '</svg>'])
    largest = max(1, max(b['whishi'] for b in boxes.values()))
    magnitude = 10 ** math.floor(math.log10(largest / 4))
    step = max(1, next(n * magnitude for n in (1, 2, 2.5, 5, 10) if n * magnitude >= largest / 4))
    ticks = math.ceil(largest / step)
    left, right = 52, width - margin
    pitch = (right - left) / len(BLOCK_CATEGORIES)
    frame_height = 122
    offset = 49
    svg = ["<g data-panel='insert-length-by-blocks'>"]
    # Split large pools into groups of four on a common bp scale.
    for start in range(0, len(scopes), 4):
        chunk = scopes[start:start + 4]
        long_labels = any(len(_scope_label(s['scope_id'])) > 7 or not s['scope_id'].startswith('lib') for s in chunk)
        columns, column_width = (2, 172) if long_labels else (4, 62)
        legend_rows = math.ceil(len(chunk) / columns)
        longest_row = min(len(chunk), columns) * column_width
        # Short keys share the unit-label line. Longer keys wrap below it.
        legend_y = offset + (20 if start == 0 and longest_row > right - 254 else 0)
        for i, scope in enumerate(chunk):
            row = i // columns
            row_count = min(columns, len(chunk) - row * columns)
            lx = right - (row_count - i % columns) * column_width
            ly = legend_y + row * 19
            sid = scope['scope_id']
            svg += [f"<rect data-library-key='{html.escape(sid, quote=True)}' x='{lx}' y='{ly-10}' "
                    f"width='11' height='9' fill='{scope_color(sid)}' stroke='#000000' stroke-width='.8'/>",
                    text(lx + 17, ly, _scope_label(sid), 'start', 13, '#263541')]
        top = legend_y + 19 * legend_rows
        bottom = top + frame_height

        def y(bp):
            return bottom - bp / (ticks * step) * (bottom - top)

        svg.append(f"<path d='M{left},{top} V{bottom} H{right}' fill='none' stroke='#000000' stroke-width='.8'/>")
        for index in range(ticks + 1):
            value = index * step
            yy = y(value)
            svg += [f"<path d='M{left-4},{yy:.2f} h4' stroke='#000000' stroke-width='.8'/>",
                    text(left - 8, yy + 4, f'{value:,g}', 'end', 13, '#263541')]
        spacing = min(25, pitch * .8 / len(chunk))
        for j, (category, label) in enumerate(zip(BLOCK_CATEGORIES, BLOCK_LABELS)):
            center = left + (j + .5) * pitch
            svg += [f"<path d='M{center:.2f},{bottom} v4' stroke='#000000' stroke-width='.8'/>",
                    text(center, bottom + 18, label, size=13, color='#263541')]
            for i, scope in enumerate(chunk):
                sid = scope['scope_id']
                if (sid, category) in boxes:
                    x = center + (i - (len(chunk) - 1) / 2) * spacing
                    svg.append(_length_box_svg(scope, category, boxes[sid, category], x, spacing * .73, y,
                                               scope_color(sid)))
        offset = bottom + 31
    svg.append(text(left, offset + 4, 'Whiskers within 1.5 × IQR. Outliers hidden.', 'start', 13))
    svg.append('</g></svg>')
    return ''.join(svg_start(width, offset + 20, title, title_x=margin) + [subtitle] + svg)


def _ligation_chart(scopes, result):
    stats_by_id = {'total': result, **result['scopes']}
    width, left, right, count_x = 500, 18, 482, 342
    first_row_y, row_pitch = 101, 56
    legend_top = first_row_y + row_pitch * len(scopes)
    svg = svg_start(width, legend_top + 33, 'Genomic blocks per library unit', title_x=left)
    svg += [text(left, 49, 'Screened non-dimer units, including unresolved units', 'start', 14),
            text(left, 75, 'Library', 'start', 13), text(count_x, 75, 'Screened units', 'end', 13),
            text(right, 75, '≥2 blocks', 'end', 13)]
    for i, scope in enumerate(scopes):
        sid = scope['scope_id']
        stats = stats_by_id.get(sid)
        # The ligation JSON lists observed scopes. QC metrics also include configured empty libraries.
        if stats is None and scope.get('ligation_screened_units') == 0:
            stats = summarize_ligation([])
        y = first_row_y + row_pitch * i
        svg += [_scope_group(scope),
                f"<circle cx='{left+3}' cy='{y-4}' r='3.2' fill='{scope_color(sid)}'/>",
                text(left + 13, y, _scope_label(sid), 'start', 14, '#263541')]
        if stats is None:
            svg += [text(right, y, 'Unavailable', 'end'), '</g>']
            continue
        counts = stats['category_counts']
        n = stats['ligation_screened_units']
        if sum(counts[key] for key in BLOCK_CATEGORIES) != n:
            raise ValueError(f'Ligation structure: category counts do not match the denominator for {sid}')
        multi = sum(counts[key] for key in ('2', '3', '4', '5+'))
        svg += [text(count_x, y, f'{n:,}', 'end', 15, '#263541'),
                text(right, y, f'{multi / n:.1%}' if n else '—', 'end', 15, '#263541')]
        if n:
            svg.extend(_stacked_segments(sid, BLOCK_CATEGORIES, BLOCK_LABELS, BLOCK_COLORS, counts,
                                         n, left, right, y + 17, label_size=13))
        else:
            svg.append(text(left, y + 21, 'No observations', 'start', 12))
        svg.append('</g>')
    for i, (label, color) in enumerate(zip(BLOCK_LABELS, BLOCK_COLORS)):
        x, y = left + (i % 3) * (right - left) / 3, legend_top + (i // 3) * 18
        svg += [f"<rect x='{x}' y='{y-9}' width='10' height='10' fill='{color}'/>",
                text(x + 17, y, label, 'start', 13, '#263541')]
    return ''.join(svg) + '</svg>'


def _mnase_chart(scopes, fragment_file):
    lengths = {scope['scope_id']: [] for scope in scopes}
    with gzip.open(fragment_file, 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            sid = row['scope_id']
            value = int(row['reference_bp'])
            if 'total' in lengths:
                lengths['total'].append(value)
            if sid != 'total' and sid in lengths:
                lengths[sid].append(value)

    # Match the individual panels' grid, display range, bandwidth and count weighting.
    maximum = 1000
    left, right, top, bottom = 72, 806, 82, 270
    height = 361 + 24 * len(scopes)
    svg = svg_start(840, height, 'Inferred pre-ligation fragment size')
    svg += [text(64, 49, 'Reference spans of blocks in resolved multi-block units', 'start', 13),
            text(left, top - 11, 'Count-weighted signal (peak = 100 per curve)', 'start', 12)]
    for tick in (0, 25, 50, 75, 100):
        y = bottom - tick / 100 * (bottom - top)
        svg += [f"<line x1='{left}' x2='{right}' y1='{y}' y2='{y}' stroke='#e1e7eb'/>",
                text(left - 10, y + 4, tick, 'end')]
    for bp in (0, 50, 100, 150, 200, 300, 500, 700, 1000):
        x = left + math.log1p(bp / 50) / math.log1p(maximum / 50) * (right - left)
        svg += [f"<line x1='{x:.2f}' x2='{x:.2f}' y1='{top}' y2='{bottom}' stroke='#edf0f2'/>",
                text(x, bottom + 18, f'{bp:,}')]
    svg += [text((left + right) / 2, bottom + 38, 'Fragment size (bp; nonlinear scale)', size=13),
            text(130, 340, 'Library', 'start', 12),
            text(330, 340, 'Fragments', 'end', 12),
            text(463, 340, 'Median (bp)', 'end', 12),
            text(610, 340, '100–<200 bp', 'end', 12),
            text(751, 340, '>1,000 bp (n)', 'end', 12)]
    visible = False
    for i, scope in enumerate(scopes):
        sid = scope['scope_id']
        color = scope_color(sid)
        values = lengths[sid]
        y = 361 + i * 24
        svg.append(_scope_group(scope))
        trace = size_trace_data(values, maximum, .045)
        if trace['shown']:
            visible = True
            signal = trace['signal']
            points = ' '.join(f'{left+j*(right-left)/(len(signal)-1):.2f},{bottom-s*(bottom-top)/100:.2f}'
                              for j, s in enumerate(signal))
            svg.append(f"<polyline points='{points}' fill='none' stroke='{color}' "
                       "stroke-width='2.2' stroke-linejoin='round'/>")
        svg += [f"<line x1='104' x2='123' y1='{y-4}' y2='{y-4}' stroke='{color}' stroke-width='2.5'/>",
                text(130, y, _scope_label(sid), 'start', 13, '#263541'),
                text(330, y, f'{len(values):,}', 'end', 13, '#263541'),
                text(463, y, f'{statistics.median(values):,g}' if values else '—', 'end', 13, '#263541'),
                text(610, y, f'{sum(100 <= bp < 200 for bp in values) / len(values):.1%}' if values else '—',
                     'end', 13, '#263541'),
                text(751, y, f'{trace["above_range"]:,}', 'end', 13, '#263541'), '</g>']
    if not visible:
        svg.append(text((left + right) / 2, (bottom + top) / 2,
                        'No resolved fragments in the displayed range', size=14))
    return ''.join(svg) + '</svg>'


def write_microc_comparisons(sample, scopes, data_dir, assets):
    """Write requested scope comparisons for a working run (qc/) or a bundle (data/)."""
    if not scopes:
        return {}
    data_dir, assets = Path(data_dir), Path(assets)
    assets.mkdir(parents=True, exist_ok=True)
    figures = {}
    pairtools = data_dir / 'pairtools' / f'{sample}.json'
    if pairtools.is_file():
        figures['pairtools'] = _pairtools_chart(scopes, json.loads(pairtools.read_text()))
    ligation = data_dir / 'ligation' / f'{sample}.json'
    if ligation.is_file():
        figures['ligation'] = _ligation_chart(scopes, json.loads(ligation.read_text()))
    units = data_dir / 'ligation' / f'{sample}.units.tsv.gz'
    if units.is_file():
        figures['ligation_lengths'] = _ligation_lengths_chart(scopes, units)
    fragments = data_dir / 'ligation' / f'{sample}.mnase_fragments.tsv.gz'
    if fragments.is_file():
        figures['mnase'] = _mnase_chart(scopes, fragments)
    result = {}
    for key, figure in figures.items():
        result[key] = f'{sample}.microc.{key}.svg'
        (assets / result[key]).write_text(figure)
    return result

def native_contact_summary(scopes):
    """Native endpoint QC, with the counts behind every saved fraction."""
    rates = [('Valid contact fraction', 'valid_contact_rate', 'unique_valid_contacts', 'mappable_library_units'),
             ('Cis fraction', 'cis_ratio', 'cis_contacts', 'unique_valid_contacts'),
             ('Long-range cis fraction', 'long_range_cis_fraction', 'long_range_cis_contacts', 'unique_valid_contacts')]
    if not any(s.get('unique_valid_contacts') is not None or
               any(s.get(key) is not None for _, key, _, _ in rates) for s in scopes):
        return None

    def count(value):
        return '—' if value is None else f'{value:,}'

    cutoffs = {s.get('long_range_cis_min_bp') for s in scopes}
    shared_cutoff = next(iter(cutoffs)) if len(cutoffs) == 1 else None
    subtitle = 'Valid contact fraction uses mappable units. Cis fractions use all unique valid contacts.'
    if shared_cutoff is not None:
        subtitle += f' Long-range cis ≥{shared_cutoff:,} bp.'
    row_height = 58 if len(cutoffs) > 1 else 45
    svg = svg_start(840, 90 + row_height * len(scopes), 'Native Micro-C contacts')
    svg += [text(64, 49, subtitle, 'start', 13),
            text(26, 78, 'Library', 'start', 12),
            text(250, 78, 'Unique valid contacts', 'end', 12)]
    columns = (445, 600, 812)
    for x, (label, _, _, _) in zip(columns, rates):
        svg.append(text(x, 78, label, 'end', 12))
    for i, scope in enumerate(scopes):
        sid = scope['scope_id']
        y = 103 + i * row_height
        svg += [_scope_group(scope),
                f"<line x1='24' x2='815' y1='{y-17}' y2='{y-17}' stroke='#e1e7eb'/>",
                f"<circle cx='28' cy='{y-3}' r='4' fill='{scope_color(sid)}'/>",
                text(40, y, _scope_label(sid), 'start', 13, '#263541'),
                text(250, y, count(scope.get('unique_valid_contacts')), 'end', 14, '#263541')]
        for x, (_, key, numerator_key, denominator_key) in zip(columns, rates):
            fraction = scope.get(key)
            numerator = scope.get(numerator_key)
            denominator = scope.get(denominator_key)
            display = f'{fraction:.1%}' if fraction is not None and denominator != 0 else '—'
            svg += [text(x, y, display, 'end', 14, '#263541'),
                    text(x, y + 16, f'{count(numerator)} / {count(denominator)}', 'end', 12)]
        if len(cutoffs) > 1 and scope.get('long_range_cis_min_bp') is not None:
            svg.append(text(columns[-1], y + 31, f'≥{scope["long_range_cis_min_bp"]:,} bp', 'end', 12))
        svg.append('</g>')
    return ''.join(svg) + '</svg>'
