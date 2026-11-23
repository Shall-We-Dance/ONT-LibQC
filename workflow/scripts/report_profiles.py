"""Compare recovered sublibraries using the retained numeric plot data."""
import csv
import json
import xml.etree.ElementTree as ET
import html
import math
from pathlib import Path
from urllib.parse import quote
import re

from report_palette import scope_color
from report_svg import svg_start, text
from report_regions import REGION_ROWS, REGION_LEGEND, region_track, region_percent
from report_profile_band import read_profile_band
from profile_data import read_gene_profile, read_numeric_profile, read_data
from qc_plot_metrics import refresh_scope_size_traces

CSS = Path(__file__).with_name("report_profiles.css").read_text()


INTERACTION_SCRIPT = "<script>" + Path(__file__).with_name("report_profiles.js").read_text() + "</script>"


def size_figure(svg, url, legend=""):
    """Inline the generated vector plot so highlighting also works in offline ZIPs."""
    escaped = html.escape(url, quote=True)
    # Keep the standalone SVG intact; use a responsive HTML heading in the report.
    svg = re.sub(r"<text\b[^>]*>Full-library size profiles</text>", '', svg)
    svg = re.sub(r"<text\b[^>]*>ONT-derived · length-weighted</text>", '', svg)
    svg = re.sub(r"viewBox='0 0 780 (\d+)'",
                 lambda match: f"viewBox='0 54 780 {int(match[1]) - 54}'", svg, count=1)
    heading = ("<div class='profile-size-heading'><div><h3>Full-library size profiles</h3>"
               "<p>ONT-derived · length-weighted</p></div>"+legend+"</div>")
    return ("<figure class='plot-card comparison-size'>"+heading+"<div class='profile-size-chart'>"+svg+
            "</div><figcaption class='figure-link'><a href='"+escaped+
            "' target='_blank' rel='noopener'>Open figure ↗</a></figcaption></figure>")


def legend_item(scope, label, color):
    esc = html.escape
    return (f"<li><button type='button' data-scope='{esc(scope,quote=True)}' aria-label='Highlight {esc(label,quote=True)}' aria-pressed='false' style='--legend-tint:{color}0f;--legend-active:{color}24;--legend-border:{color}30'>"
            f"<i style='background:{color}' aria-hidden='true'></i>{esc(label)}</button></li>")


def polyline(xs, ys, color, label):
    # Missing bins stay missing; never draw across an unavailable interval.
    paths, current = [], []
    for x, y in zip(xs, ys):
        if math.isfinite(y): current.append(f'{x:.3f},{y:.3f}')
        elif current: paths.append(current); current = []
    if current: paths.append(current)
    return ''.join(f"<polyline data-scope='{html.escape(label, quote=True)}' points='{' '.join(points)}' fill='none' stroke='{color}' stroke-width='2' stroke-linejoin='round'/>" for points in paths)


def read_size_profile(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    xs = [float(row['size_bp']) for row in rows]
    ys = [float(row['length_weighted_relative_signal']) for row in rows]
    if len(xs) < 2 or any(not math.isfinite(x) or x < 0 for x in xs) or any(b <= a for a, b in zip(xs, xs[1:])):
        raise ValueError('Invalid size profile axis')
    if any(not math.isfinite(y) or not 0 <= y <= 100.000001 for y in ys):
        raise ValueError('Invalid size profile signal')
    return xs, ys


def size_comparison(series, *, compact=False):
    # Match the individual electropherogram: the same grids, type and baseline.
    left, right, top = 64, 744, 64
    bottom, height = (250, 308) if compact else (314, 390)
    mid = (top + bottom) / 2
    signal_label = 'Relative signal (peak = 100)' if compact else 'Relative signal (a.u.; peak = 100)'
    maximum = max(row['sizes'][0][-1] for row in series)
    upper = math.log1p(maximum / 50)
    x = lambda bp: left + math.log1p(bp / 50) / upper * (right-left)
    svg = svg_start(780, height, 'Full-library size profiles')
    svg[3] = "<text x='64' y='28' font-size='19' font-weight='600'>Full-library size profiles</text>"
    svg.append(text(left, 46, 'ONT-derived · length-weighted', 'start', 13, '#576574'))
    for value in (0, 25, 50, 75, 100):
        y = bottom - value / 100 * (bottom-top)
        svg += [f"<line x1='{left}' x2='{right}' y1='{y}' y2='{y}' stroke='#e1e7eb'/>", text(left-12, y+4, value, 'end', 13, '#485460')]
    for value in (0, 50, 100, 150, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000, 10000):
        if value > maximum: continue
        svg += [f"<line x1='{x(value):.3f}' x2='{x(value):.3f}' y1='{top}' y2='{bottom}' stroke='#edf0f2'/>", text(x(value), bottom+21, f'{value:,}', size=13, color='#111')]
    for row in series:
        xs, ys = row['sizes']
        points = ' '.join(f'{x(bp):.3f},{bottom-v/100*(bottom-top):.3f}' for bp,v in zip(xs,ys))
        svg.append(f"<polygon data-scope='{html.escape(row['id'],quote=True)}' points='{x(xs[0]):.3f},{bottom} {points} {x(xs[-1]):.3f},{bottom}' fill='{row['color']}' fill-opacity='.06'/>")
    for row in series:
        xs, ys = row['sizes']
        svg.append(polyline([x(v) for v in xs], [bottom-v/100*(bottom-top) for v in ys], row['color'], row['id']))
    svg += [f"<line x1='{left}' x2='{right}' y1='{bottom}' y2='{bottom}' stroke='#586069'/>",
            text((left+right)/2, bottom+45, 'Fragment size (bp; nonlinear size axis)', size=14, color='#111'),
            f"<text x='23' y='{mid:g}' transform='rotate(-90 23 {mid:g})' text-anchor='middle' font-size='14'>{signal_label}</text>", '</svg>']
    return ''.join(svg)


def virtual_lanes(series, *, compact=False):
    width = max(190, 100 + 44*len(series))
    top = 64
    bottom, height = (250, 308) if compact else (314, 390)
    maximum = max(row['sizes'][0][-1] for row in series);upper = math.log1p(maximum/50)
    y = lambda bp: bottom - math.log1p(bp/50)/upper*(bottom-top)
    svg = svg_start(width, height, 'Virtual lanes')
    svg[3] = "<text x='68' y='28' font-size='19' font-weight='600'>Virtual lanes</text>"
    for value in sorted({0, 100, 300, 1000, maximum}):
        if value <= maximum:
            svg += [text(54, y(value)+4, f'{value:,.0f}', 'end', 13, '#485460'),
                    f"<line x1='59' x2='64' y1='{y(value):.3f}' y2='{y(value):.3f}' stroke='#586069'/>"]
    for i, row in enumerate(series):
        xs, ys = row['sizes'];lane_x=68+i*44;lane_w=32
        svg.append(f"<g data-scope='{html.escape(row['id'], quote=True)}'><title>{html.escape(row['label'])}</title>")
        for j, signal in enumerate(ys):
            y0 = y(xs[j]) if j == len(xs)-1 else (y(xs[j])+y(xs[j+1]))/2
            y1 = y(xs[j]) if j == 0 else (y(xs[j])+y(xs[j-1]))/2
            shade=round(248-205*signal/100)
            svg.append(f"<rect x='{lane_x}' y='{y0:.3f}' width='{lane_w}' height='{y1-y0:.3f}' fill='rgb({shade},{shade},{shade})'/>")
        svg += [f"<rect x='{lane_x}' y='{top}' width='{lane_w}' height='{bottom-top}' fill='none' stroke='#aab4bc' stroke-width='.8'/>",
                f"<rect x='{lane_x}' y='{top-7}' width='{lane_w}' height='4' fill='{row['color']}'/>",
                text(lane_x+lane_w/2, bottom+23, row['id'], size=13, color='#34495e'), '</g>']
    svg += [text(width/2, bottom+45, 'Fragment size (bp)', size=14, color='#111'), '</svg>']
    return ''.join(svg)


def _gene_panel(series):
    rows = [row for row in series if row.get('gene')]
    if not rows: return None
    labels, bins, _ = rows[0]['gene']
    if any(row['gene'][:2] != (labels, bins) for row in rows):
        raise ValueError('Sublibrary gene profiles must share the same bin coordinates')
    left, right, top, bottom = 80, 792, 52, 350
    maximum = max((v for row in rows for v in (row['gene'][2] + (row['band'][1] if row.get('band') else [])) if math.isfinite(v)), default=0) or 1
    # A little headroom and round tick values keep the axes readable at report size.
    magnitude = 10 ** math.floor(math.log10(maximum / 4))
    step = min((v * magnitude for v in (1, 2, 2.5, 5, 10)), key=lambda v: abs(math.log(v / (maximum / 4))))
    maximum = math.ceil(maximum * 1.06 / step) * step
    x = lambda v: left+(v-bins[0])/(bins[-1]-bins[0])*(right-left)
    y = lambda v: bottom-v/maximum*(bottom-top)
    legend_rows = math.ceil(len(rows) / 8)
    svg = svg_start(820, 440 + 24*(legend_rows-1), 'TSS–TES metaplot')
    svg.append(f"<defs><clipPath id='profile-frame'><rect x='{left}' y='{top}' width='{right-left}' height='{bottom-top}'/></clipPath></defs>")
    for tick in range(round(maximum / step) + 1):
        value = tick * step
        svg += [f"<line x1='{left-5}' x2='{left}' y1='{y(value):.3f}' y2='{y(value):.3f}' stroke='#202b36'/>", text(left-11, y(value)+5, f'{value:g}', 'end', 14, '#202b36')]
    for label, b in zip(labels, bins):
        if label:
            label = re.sub(r'([+-]?\d+)(?:\.0)?Kb', r'\1 kb', label).replace('-', '−')
            svg += [f"<line x1='{x(b):.3f}' x2='{x(b):.3f}' y1='{bottom}' y2='{bottom+5}' stroke='#202b36'/>", text(x(b), bottom+28, label, size=15, color='#202b36')]
    # Bands are recovered from the original deepTools SE plot, never synthesized.
    for row in rows:
        if not row.get('band'): continue
        lower, upper = row['band']
        segment = []
        for i in range(len(bins)+1):
            if i < len(bins) and math.isfinite(lower[i]) and math.isfinite(upper[i]):
                segment.append(i)
            elif segment:
                points = [(x(bins[j]), y(upper[j])) for j in segment]+[(x(bins[j]), y(lower[j])) for j in reversed(segment)]
                coords = ' '.join(f'{px:.3f},{py:.3f}' for px,py in points)
                svg.append(f"<polygon data-band='{html.escape(row['id'],quote=True)}' points='{coords}' fill='{row['color']}' fill-opacity='.16' clip-path='url(#profile-frame)'/>")
                segment = []
    for row in rows:
        svg.append(polyline([x(b) for b in bins], [y(v) for v in row['gene'][2]], row['color'], row['id']).replace("stroke-width='2'", "stroke-width='2.4' clip-path='url(#profile-frame)'"))
    svg.append(f"<rect x='{left}' y='{top}' width='{right-left}' height='{bottom-top}' fill='none' stroke='#202b36' stroke-width='1.2'/>")
    svg.append("<text x='22' y='201' transform='rotate(-90 22 201)' text-anchor='middle' font-size='15' fill='#202b36'>CPM-normalized signal</text>")
    for i, row in enumerate(rows):
        count = min(8, len(rows) - (i//8)*8);gap = 86
        lx = (820-gap*count)/2 + gap*(i%8);ly = 417 + 24*(i//8)
        svg += [f"<line x1='{lx}' x2='{lx+20}' y1='{ly}' y2='{ly}' stroke='{row['color']}' stroke-width='2.4'/>", text(lx+27, ly+4, row['id'], 'start', 13, '#202b36')]
    svg.append('</svg>')
    return ''.join(svg)



def gene_comparison(series):
    """Small multiples retain CPM values and use an independent y axis per panel."""
    rows = [row for row in series if row.get('gene')]
    if not rows:
        return None
    labels, bins, _ = rows[0]['gene']
    if any(row['gene'][:2] != (labels, bins) for row in rows):
        raise ValueError('Sublibrary gene profiles must share the same bin coordinates')
    columns = min(2, len(rows))
    width, height = 820 * columns, 440 * math.ceil(len(rows) / columns)
    result = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='Sublibrary TSS–TES profiles with independent y axes'>",
              '<title>Sublibrary TSS–TES profiles · independent CPM axes</title>']
    for i, row in enumerate(rows):
        panel = _gene_panel([row])
        panel = panel.replace('profile-frame', f'profile-frame-{i}')
        panel = panel.replace('TSS–TES metaplot', html.escape(str(row['label']), quote=True))
        panel = panel.replace('<svg ', f"<svg x='{820*(i%columns)}' y='{440*(i//columns)}' width='820' height='440' data-panel='{html.escape(row['id'],quote=True)}' ", 1)
        result.append(panel)
    result.append('</svg>')
    return ''.join(result)


def read_region_metrics(path):
    """Use the same saved values as each individual genomic-region figure."""
    path = Path(path)
    numeric = path.with_suffix('.json')
    if numeric.is_file(): return read_data(numeric, 'genomic_regions')
    if not path.is_file():
        return None
    node = ET.fromstring(path.read_text()).find("{*}metadata[@id='genomic-region-metrics']")
    return json.loads(node.text) if node is not None else None


def region_comparison(series):
    if not any(row.get('regions') is not None for row in series):
        return None
    height = 98 + 3 * (32 + 25 * len(series))
    svg = svg_start(820, height, 'Sublibrary genomic regions')
    svg.append("<metadata id='sublibrary-region-metrics'>" + html.escape(json.dumps([
        dict(scope_id=row['id'], color=row['color'], metrics=row.get('regions')) for row in series])) + '</metadata>')
    for x, label in [(535, 'Sample'), (658, 'Background'), (796, 'Enrichment')]:
        svg.append(text(x, 60, label, 'end', 13))
    for group_index, (label, key) in enumerate(REGION_ROWS):
        top = 84 + group_index * (32 + 25 * len(series))
        svg.append(text(24, top, label, 'start', 16, '#203040'))
        for i, row in enumerate(series):
            y = top + 23 + i * 25
            metric = row.get('regions') or {}
            rate, bg, fold = (metric.get(key + suffix) for suffix in ('_rate', '_background_rate', '_enrichment'))
            tip = f"{row['label']} · {label}: {region_percent(rate)}, background {region_percent(bg)}"
            svg.append(f"<g data-scope='{html.escape(row['id'],quote=True)}' data-region='{key}'><title>{html.escape(tip)}</title>")
            svg.append(text(28, y + 4, row['id'], 'start', 13, '#34495e'))
            svg.append(region_track(rate, bg, row['color'], 108, y, 322))
            for x, value in [(535,region_percent(rate)),(658,region_percent(bg)),(796,'—' if fold is None else f'{fold:,.2f}×')]:
                svg.append(text(x,y+4,value,'end',14,'#34495e'))
            svg.append('</g>')
    svg.append(text(24,height-13,REGION_LEGEND,'start',12))
    svg.append('</svg>')
    return ''.join(svg)

def rna_region_comparison(series, metrics):
    from rna_report import REGION_LABELS
    group_height = 30 + 24 * len(series)
    height = 92 + len(REGION_LABELS) * group_height
    svg = svg_start(760, height, 'RNA region distribution')
    svg.extend([text(590, 55, 'Matched bases', 'end', 13), text(737, 55, 'Fraction', 'end', 13)])
    for i, (key, label) in enumerate(REGION_LABELS.items()):
        top = 80 + i * group_height
        svg.append(text(24, top, label, 'start', 15, '#203040'))
        for j, row in enumerate(series):
            y = top + 20 + 24 * j
            scope = metrics.get(row['id'], {})
            counts = scope.get('region_bases')
            total = scope.get('region_denominator_bases')
            count = counts.get(key, 0) if counts is not None else None
            fraction = count / total if total and count is not None else None
            svg.append(text(28, y + 4, row['id'], 'start', 13, '#34495e'))
            svg.append(region_track(fraction, None, row['color'], 110, y, 290))
            svg.append(text(590, y + 4, '—' if count is None else f'{count:,}', 'end', 14))
            svg.append(text(737, y + 4, region_percent(fraction), 'end', 14))
    svg.append(text(24, height - 12, 'Fraction of matched bases within each library', 'start', 12))
    return ''.join(svg) + '</svg>'


def transcript_panels(series, metrics, prefix):
    panels = ['<section class="rna-profile-panels plot-card"><h3>Transcript metaplots</h3><div>']
    for row in series:
        label = html.escape(row['label'])
        plot = metrics.get(row['id'], {}).get('plots', {}).get('profiles')
        panels.append(f'<figure><figcaption><i style="background:{row["color"]}"></i>{label}</figcaption>')
        if plot:
            url = html.escape(prefix + '/' + quote(plot, safe=''), quote=True)
            panels.append(f'<a href="{url}" target="_blank" rel="noopener"><img src="{url}" alt="{label} transcript coverage"></a>')
        else:
            panels.append('<p class="muted">Transcript profile unavailable</p>')
        panels.append('</figure>')
    return ''.join(panels) + '</div></section>'


def render_profiles_card(sample, scopes, scope_html, asset_dir, gene_dir, asset_prefix='assets/', *, rna_summary=None, rna_prefix=None, microc=None):
    """Works during report creation and on a retained, report-only bundle."""
    asset_dir, gene_dir = Path(asset_dir), Path(gene_dir)
    esc = html.escape
    def figure(name, caption, cls=''):
        url=esc(asset_prefix+name,quote=True)
        return f"<figure class='plot-card {cls}'><a href='{url}' target='_blank' rel='noopener'><img src='{url}' alt='{esc(caption,quote=True)}'></a><figcaption class='figure-link'><a href='{url}' target='_blank' rel='noopener'>Open figure ↗</a></figcaption></figure>"
    title = "Size and RNA profiles" if rna_summary is not None else "Size and genomic profiles"
    heading = f'<h2>{title}</h2>'
    if microc:
        link = esc(asset_prefix + f'{sample}.microc_qc.svg', quote=True)
        heading = ("<div class='microc-heading'><h2>Size and Micro-C QC</h2>"
                   f"<a href='{link}' download>Download overview (SVG)</a></div>")
    section_class = 'report-section microc-profiles' if microc else 'report-section'
    opening = f"<section class='{section_class}' id='profiles'>"
    sub=scopes[1:]
    refresh_scope_size_traces(asset_dir, sample, sub)
    if not sub and rna_summary is None and microc is None:
        return opening + heading + scope_html[0] + '</section>'
    parts = [opening + f'<style>{CSS}</style>' + heading]
    series=[]
    for scope in sub or scopes[:1]:
        sid=scope['scope_id'];label=scope.get('scope_name',sid)
        series.append(dict(id=sid,label=label,color=scope_color(sid),
                           sizes=read_size_profile(asset_dir/f'{sample}.{sid}.full_library_profiles.tsv'),
                           gene=None))
        if rna_summary is not None or microc:
            continue
        series[-1]['gene'] = read_gene_profile(gene_dir/f'{sample}.{sid}.profile.tsv')
        series[-1]['regions'] = read_region_metrics(asset_dir/f'{sample}.{sid}.regions.svg')
        numeric = gene_dir/f'{sample}.{sid}.profile.json'
        if numeric.is_file():
            series[-1]['gene'], series[-1]['band'] = read_numeric_profile(numeric, series[-1]['gene'])
        else:
            series[-1]['band'] = read_profile_band(asset_dir/f'{sample}.{sid}.TSS_TES.svg', series[-1]['gene'])
    parts.append("<div class='profile-comparison'>")
    legend = ["<ul class='profile-comparison-legend' aria-label='Size profile legend'>"]
    for row in series:
        label=row['label'] if row['id'] == 'total' or row['label']==row['id'] else row['id']+' · '+row['label']
        legend.append(legend_item(row['id'],label,row['color']))
    legend.append('</ul>')
    size_name=f'{sample}.sublibraries.size_profiles.svg';gel_name=f'{sample}.sublibraries.virtual_gel.svg'
    size_svg = size_comparison(series, compact=bool(microc))
    (asset_dir/size_name).write_text(size_svg)
    (asset_dir/gel_name).write_text(virtual_lanes(series, compact=bool(microc)))
    parts.append("<div class='profile-comparison-grid'>"+size_figure(size_svg,asset_prefix+size_name,''.join(legend))+figure(gel_name,'Sublibrary virtual lanes','comparison-gel')+'</div>')
    parts.append("<p class='comparison-note'>ONT-derived full-library sizes, length-weighted; each trace and lane is scaled to its own peak of 100.</p>")
    parts.append(INTERACTION_SCRIPT)
    if microc:
        from report_microc_layout import overview_html, write_overview_figure
        write_overview_figure(sample, asset_dir, size_name, gel_name, microc, sub or scopes)
        parts.append(overview_html(asset_prefix, microc))
    if rna_summary is not None:
        from rna_report import rna_family_comparison
        metrics = rna_summary['scopes']
        family = rna_family_comparison(metrics, [row['id'] for row in series],
                                       {row['id']: row['label'] for row in series},
                                       {row['id']: row['color'] for row in series})
        parts.append('<style>' + Path(__file__).with_name('rna_comparison.css').read_text() + '</style>')
        name = f'{sample}.rna_region_distribution.svg'
        (asset_dir/name).write_text(rna_region_comparison(series, metrics))
        prefix = rna_prefix or '../qc/rna/' + quote(sample, safe='')
        parts.extend(["<div class='profile-genomic-grid rna-distributions'>", figure(name, 'RNA region distribution', 'comparison-regions'),
                      family, '</div>', transcript_panels(series, metrics, prefix)])
    meta=gene_comparison(series)
    regions=region_comparison(series)
    if regions:
        name=f'{sample}.sublibraries.regions.svg';(asset_dir/name).write_text(regions)
        parts.append("<div class='profile-genomic-grid'>"+figure(name,'Sublibrary genomic region proportions, backgrounds and enrichment','comparison-regions'))
    if meta:
        name=f'{sample}.sublibraries.TSS_TES.svg';(asset_dir/name).write_text(meta)
        parts.append(figure(name,'Sublibrary TSS–TES profiles with independent CPM axes','comparison-metaplot'))
        parts.append("<p class='comparison-note'>Each metaplot has its own y-axis range; CPM values and SEM bands are unchanged.</p>")
        missing=[r['label'] for r in series if not r['gene']]
        if missing:parts.append("<p class='comparison-note'>No mapped profile: "+esc(', '.join(missing))+'.</p>')
    if regions:parts.append('</div>')
    parts.append('</div>')
    for i, scope in enumerate(sub):
        label=scope.get('scope_name',scope['scope_id']);color=series[i]['color']
        content=re.sub(r'<h3>.*?</h3>','',scope_html[i+1],count=1,flags=re.S)
        if microc:
            content += "<div class='microc-scope'>" + microc.details[scope['scope_id']] + '</div>'
        count=int(scope.get('library_units',0))
        parts.append(f"<details class='profile-detail' data-scope='{esc(scope['scope_id'],quote=True)}' style='--scope-tint:{color}0f;--scope-border:{color}38'><summary><i style='background:{color}'></i>{esc(label)}<small>{count:,} units</small></summary>{content}</details>")
    total_label = 'Total library · pooled profiles' if sub else 'Total library details'
    total_content = scope_html[0]
    if microc:
        total_label = 'Total library' if sub else 'Total library details'
        total_content += "<div class='microc-scope'>" + microc.details['total'] + '</div>'
    expanded = ' open' if microc and not sub else ''
    parts.append(f"<details class='profile-detail pooled-profile'{expanded}><summary>{total_label}</summary>{total_content}</details>")
    return ''.join(parts)+'</section>'
