"""Shared, alignment-independent mixture proportions for every assay report."""
import html

COLORS = ('#2878c8', '#df7a14', '#219357', '#9356c0', '#ca4d87', '#1196aa', '#819610', '#d34d45')
TOTAL_COLOR = '#496576'
AMBIGUOUS_COLOR = '#a4a49c'
UNASSIGNED_COLOR = '#dddcd4'


def bar_ink(fill):
    """Choose a label color with at least 4.5:1 contrast against its fill."""
    rgb = (int(fill[i:i + 2], 16) / 255 for i in (1, 3, 5))
    linear = (v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb)
    luminance = sum(weight * v for weight, v in zip((0.2126, 0.7152, 0.0722), linear))
    return '#ffffff' if 1.05 / (luminance + 0.05) >= 4.5 else '#000000'


def scope_color(scope, index):
    # Standard libNN identifiers retain their color when a scope has no RNA reads.
    if scope == 'total':
        return TOTAL_COLOR
    suffix = scope.removeprefix('lib')
    if suffix != scope and suffix.isascii() and suffix.isdecimal() and int(suffix) > 0:
        index = int(suffix) - 1
    return COLORS[index % len(COLORS)]


CSS = '''
.sublibrary-mixture{break-inside:avoid}
.mixture-bar{display:flex;height:38px;width:100%;overflow:hidden;border-radius:6px;background:#edf2f7;print-color-adjust:exact;-webkit-print-color-adjust:exact}
.mixture-segment{display:flex;align-items:center;justify-content:center;min-width:0;flex-shrink:0;font-size:13px;font-weight:600;overflow:hidden}
.mixture-axis{display:flex;justify-content:space-between;font-size:11px;color:var(--muted,#5b6d82);margin-top:4px}
.mixture-legend{list-style:none;padding:0;margin:16px 0 0;display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px 20px}
.mixture-legend li{display:grid;grid-template-columns:10px minmax(0,1fr);gap:8px;align-items:baseline;min-width:0}
.mixture-swatch{width:10px;height:10px;border-radius:2px;print-color-adjust:exact;-webkit-print-color-adjust:exact}
.mixture-label{overflow-wrap:anywhere}.mixture-label strong{font-weight:600}.mixture-count{display:block;font-size:12px;color:var(--muted,#5b6d82)}
@media(max-width:600px){.mixture-segment-text{display:none}}
'''


def mixture_counts(metric):
    """Never renormalize assigned units or infer mixing shares from raw reads."""
    total = int(metric.get('library_units', 0))
    rows = [(str(s.get('scope_name') or s['scope_id']), int(s['library_units']), scope_color(s.get('scope_id', ''), i))
            for i, s in enumerate(metric.get('sub_libraries', []))]
    pairing = metric.get('index_pairing', {})
    if any(n < 0 for _, n, _ in rows) or total < 0:
        raise ValueError('Negative library-unit count')
    remaining = total - sum(n for _, n, _ in rows)
    if remaining < 0:
        raise ValueError('Sublibrary counts exceed total library units')
    if 'unassigned_units' in pairing and 'ambiguous_units' in pairing:
        unassigned, ambiguous = int(pairing['unassigned_units']), int(pairing['ambiguous_units'])
        if min(unassigned, ambiguous) < 0 or unassigned + ambiguous != remaining:
            raise ValueError('Index-assignment counts do not sum to total library units')
        rows += [(label, n, color) for label, n, color in
                 [('Ambiguous', ambiguous, AMBIGUOUS_COLOR), ('Unassigned', unassigned, UNASSIGNED_COLOR)] if n]
    elif remaining:
        rows.append(('Unassigned / ambiguous', remaining, '#bcbab1'))
    return total, rows


def percent(count, total):
    value = 100 * count / total if total else 0
    return '<0.1%' if 0 < value < .05 else f'{value:.1f}%'


def render_sublibrary_composition(metric):
    if len(metric.get('sub_libraries', [])) < 2:
        return ''
    total, rows = mixture_counts(metric)
    if not total:
        return ''
    esc = html.escape
    description = '; '.join(f'{label}: {percent(n, total)}, {n:,} units' for label, n, _ in rows)
    parts = ['<section class="report-section sublibrary-mixture" id="sublibrary-composition">',
             '<style>' + CSS + '</style><h2>Sublibrary composition</h2>',
             f'<p class="muted">Share of all {total:,} reconstructed library units, including adapter dimers and unresolved index assignments. Counts are library units, not raw ONT reads.</p>',
             '<div class="mixture-bar" role="img" aria-label="' + esc(description, quote=True) + '">']
    for label, n, color in rows:
        if not n:
            continue
        width = 100 * n / total
        title = esc(f'{label}: {percent(n, total)} ({n:,} units)', quote=True)
        text = esc(percent(n, total)) if width >= 8 else ''
        parts.append(f'<span class="mixture-segment" style="width:{width:.10f}%;background:{color};color:{bar_ink(color)}" title="{title}" aria-hidden="true"><span class="mixture-segment-text">{text}</span></span>')
    parts.append('</div><div class="mixture-axis" aria-hidden="true"><span>0%</span><span>100%</span></div><ul class="mixture-legend">')
    for label, n, color in rows:
        parts.append(f'<li><span class="mixture-swatch" style="background:{color}" aria-hidden="true"></span><span class="mixture-label">{esc(label)} <strong>{esc(percent(n, total))}</strong><span class="mixture-count">{n:,} units</span></span></li>')
    return ''.join(parts) + '</ul></section>'
