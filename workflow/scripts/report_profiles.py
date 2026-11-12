"""Compare recovered sublibraries using the retained numeric plot data."""
import csv
import html
import math
from pathlib import Path
import re

from report_sublibrary import COLORS
from report_profile_band import read_profile_band

CSS = '''
.profile-comparison-legend{display:flex;flex-wrap:wrap;gap:8px 22px;list-style:none;margin:8px 0 14px;padding:0;font-size:13px}
.profile-comparison-legend i,.profile-detail summary i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:7px;print-color-adjust:exact;-webkit-print-color-adjust:exact}
.profile-comparison-grid{display:grid;grid-template-columns:minmax(0,1fr) minmax(220px,.43fr);gap:12px;align-items:stretch}
.profile-comparison .plot-card{margin:0;min-width:0;padding:12px;border-radius:9px}.profile-comparison img{display:block;width:100%;height:auto;min-height:0}
.profile-comparison-grid .plot-card>a{display:flex;align-items:center;justify-content:center;height:clamp(240px,25vw,350px)}
.profile-comparison-grid img{width:100%;height:100%;object-fit:contain}
.profile-comparison .comparison-metaplot{max-width:850px;margin:12px auto 0}
.profile-comparison .figure-link{margin-top:6px;font-size:11px}
.profile-comparison .comparison-note{font-size:12px;color:var(--muted);margin:8px 0 16px}
.profile-detail{margin:12px 0 0;border:1px solid var(--line);border-radius:8px;padding:0;background:#fff}
.profile-detail>summary{padding:13px 16px;color:var(--deep);font-weight:550}.profile-detail>summary small{font-weight:400;color:var(--muted);margin-left:12px}
.profile-detail>.profile-scope{margin:0;padding:0 16px 16px;border:0}.profile-detail .primary-size-plots{margin-top:8px}
#profiles{scroll-margin-top:80px}#profiles>.profile-scope{border:0;margin-top:0;padding-top:0}
@media(max-width:700px){.profile-comparison-grid{grid-template-columns:minmax(0,1fr)}.profile-comparison-grid .plot-card>a{height:auto}.profile-comparison-grid .comparison-gel>a{height:290px}.profile-detail>summary small{display:block;margin:5px 0 0 20px}}
@media print{.profile-comparison-grid{display:block}.profile-comparison .plot-card{max-width:100%;break-inside:avoid}.profile-detail{border:0}.profile-detail>summary{padding:2mm 0}.profile-detail>.profile-scope{padding:0}}
'''


def svg_start(width, height, title):
    return [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='{html.escape(title, quote=True)}' style='font-family:Arial,Helvetica,sans-serif'>",
            f'<title>{html.escape(title)}</title>', "<rect width='100%' height='100%' fill='white'/>",
            f"<text x='64' y='27' font-size='17' fill='#293f53' font-weight='600'>{html.escape(title)}</text>"]


def text(x, y, label, anchor='middle', size=12, color='#526579'):
    return f"<text x='{x:.3f}' y='{y:.3f}' text-anchor='{anchor}' font-size='{size}' fill='{color}'>{html.escape(str(label))}</text>"


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


def read_gene_profile(path):
    if not Path(path).is_file(): return None
    with Path(path).open() as handle: rows = list(csv.reader(handle, delimiter='\t'))
    if len(rows) == 1 and rows[0][0].startswith('No '): return None
    if len(rows) != 3 or rows[0][:1] != ['bin labels'] or rows[1][:1] != ['bins']:
        raise ValueError('Unsupported gene profile table')
    labels, xs, ys = rows[0][2:], [float(v) for v in rows[1][2:]], [float(v) for v in rows[2][2:]]
    if len(xs) < 2 or len(xs) != len(ys) or len(xs) != len(labels) or any(not math.isfinite(x) for x in xs) or any(b <= a for a, b in zip(xs, xs[1:])):
        raise ValueError('Invalid gene profile axis')
    if any(math.isinf(y) or y < 0 for y in ys): raise ValueError('Invalid gene coverage')
    return labels, xs, ys


def size_comparison(series):
    # Match the individual electropherogram: the same grids, type and baseline.
    left, right, top, bottom = 64, 744, 64, 314
    maximum = max(row['sizes'][0][-1] for row in series)
    upper = math.log1p(maximum / 50)
    x = lambda bp: left + math.log1p(bp / 50) / upper * (right-left)
    svg = svg_start(780, 390, 'Full-library size profiles')
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
        svg.append(f"<polygon points='{x(xs[0]):.3f},{bottom} {points} {x(xs[-1]):.3f},{bottom}' fill='{row['color']}' opacity='.06'/>")
    for row in series:
        xs, ys = row['sizes']
        svg.append(polyline([x(v) for v in xs], [bottom-v/100*(bottom-top) for v in ys], row['color'], row['id']))
    svg += [f"<line x1='{left}' x2='{right}' y1='{bottom}' y2='{bottom}' stroke='#586069'/>",
            text((left+right)/2, bottom+45, 'Fragment size (bp; nonlinear size axis)', size=14, color='#111'),
            "<text x='23' y='189' transform='rotate(-90 23 189)' text-anchor='middle' font-size='14'>Relative signal (a.u.; peak = 100)</text>", '</svg>']
    return ''.join(svg)


def virtual_lanes(series):
    width = 100 + 44*len(series);top, bottom = 64, 314
    maximum = max(row['sizes'][0][-1] for row in series);upper = math.log1p(maximum/50)
    y = lambda bp: bottom - math.log1p(bp/50)/upper*(bottom-top)
    svg = svg_start(width, 390, 'Virtual gel')
    svg[3] = "<text x='68' y='28' font-size='19' font-weight='600'>Virtual gel</text>"
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
                text(lane_x+lane_w/2, bottom+23, row['id'], size=13, color=row['color']), '</g>']
    svg += [text(width/2, bottom+45, 'Fragment size (bp)', size=14, color='#111'), '</svg>']
    return ''.join(svg)


def gene_comparison(series):
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


def render_profiles_card(sample, scopes, scope_html, asset_dir, gene_dir, asset_prefix='assets/'):
    """Works during report creation and on a retained, report-only bundle."""
    asset_dir, gene_dir = Path(asset_dir), Path(gene_dir)
    esc = html.escape
    def figure(name, caption, cls=''):
        url=esc(asset_prefix+name,quote=True)
        return f"<figure class='plot-card {cls}'><a href='{url}' target='_blank' rel='noopener'><img src='{url}' alt='{esc(caption,quote=True)}'></a><figcaption class='figure-link'><a href='{url}' target='_blank' rel='noopener'>Open figure ↗</a></figcaption></figure>"
    parts=["<section class='report-section' id='profiles'><style>"+CSS+"</style><h2>Size and genomic profiles</h2>"]
    sub=scopes[1:]
    if not sub: return ''.join(parts)+scope_html[0]+'</section>'
    series=[]
    for i, scope in enumerate(sub):
        sid=scope['scope_id'];label=scope.get('scope_name',sid)
        series.append(dict(id=sid,label=label,color=COLORS[i%len(COLORS)],
                           sizes=read_size_profile(asset_dir/f'{sample}.{sid}.full_library_profiles.tsv'),
                           gene=read_gene_profile(gene_dir/f'{sample}.{sid}.profile.tsv')))
        series[-1]['band'] = read_profile_band(asset_dir/f'{sample}.{sid}.TSS_TES.svg', series[-1]['gene'])
    parts += ["<div class='profile-comparison'><ul class='profile-comparison-legend'>"]
    for row in series:
        label=row['label'] if row['label']==row['id'] else row['id']+' · '+row['label']
        parts.append(f"<li><i style='background:{row['color']}'></i>{esc(label)}</li>")
    parts.append('</ul>')
    size_name=f'{sample}.sublibraries.size_profiles.svg';gel_name=f'{sample}.sublibraries.virtual_gel.svg'
    (asset_dir/size_name).write_text(size_comparison(series));(asset_dir/gel_name).write_text(virtual_lanes(series))
    parts.append("<div class='profile-comparison-grid'>"+figure(size_name,'Sublibrary full-library size traces')+figure(gel_name,'Sublibrary virtual gel lanes','comparison-gel')+'</div>')
    parts.append("<p class='comparison-note'>ONT-derived full-library sizes, length-weighted; each trace and lane is scaled to its own peak of 100.</p>")
    meta=gene_comparison(series)
    if meta:
        name=f'{sample}.sublibraries.TSS_TES.svg';(asset_dir/name).write_text(meta)
        parts.append(figure(name,'Sublibrary TSS–TES profiles on a shared CPM axis','comparison-metaplot'))
        missing=[r['label'] for r in series if not r['gene']]
        if missing:parts.append("<p class='comparison-note'>No mapped profile: "+esc(', '.join(missing))+'.</p>')
    parts.append('</div>')
    for i, scope in enumerate(sub):
        label=scope.get('scope_name',scope['scope_id']);color=series[i]['color']
        content=re.sub(r'<h3>.*?</h3>','',scope_html[i+1],count=1,flags=re.S)
        count=int(scope.get('library_units',0))
        parts.append(f"<details class='profile-detail' data-scope='{esc(scope['scope_id'],quote=True)}'><summary><i style='background:{color}'></i>{esc(label)}<small>{count:,} units · individual profiles</small></summary>{content}</details>")
    parts.append("<details class='profile-detail pooled-profile'><summary>Total library · pooled profiles</summary>"+scope_html[0]+'</details>')
    return ''.join(parts)+'</section>'
