"""Shared, alignment-independent mixture proportions for every assay report."""
import html

from report_palette import NPG, bar_ink, scope_color


CSS = '''
.sublibrary-mixture{break-inside:avoid}
.mixture-bar{display:flex;height:38px;width:100%;overflow:hidden;border-radius:6px;background:#edf2f7;print-color-adjust:exact;-webkit-print-color-adjust:exact}
.mixture-segment{container:mixture-segment / inline-size;display:flex;align-items:center;justify-content:center;min-width:0;flex-shrink:0;font-size:13px;font-weight:600;overflow:hidden}
.mixture-axis{display:flex;justify-content:space-between;font-size:11px;color:var(--muted,#5b6d82);margin-top:4px}
.mixture-legend{list-style:none;padding:0;margin:16px 0 0;display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px 20px}
.mixture-legend li{display:grid;grid-template-columns:10px minmax(0,1fr);gap:8px;align-items:baseline;min-width:0}
.mixture-swatch{width:10px;height:10px;border-radius:2px;print-color-adjust:exact;-webkit-print-color-adjust:exact}
.mixture-label{overflow-wrap:anywhere}.mixture-label strong{font-weight:600}.mixture-count{display:block;font-size:12px;color:var(--muted,#5b6d82)}
.mixture-segment-text{display:none;white-space:nowrap;font-variant-numeric:tabular-nums}
.mixture-segment-text.outlined{-webkit-text-stroke:1px rgb(0 0 0 / 42%);paint-order:stroke fill}
.sublibrary-mixture:has(.mixture-other:not([open])) .mixture-deferred{visibility:hidden}
@container mixture-segment (min-width:6.5ch){.mixture-segment-text{display:inline}}
'''


def mixture_counts(metric, color_for_scope=None):
    """Never renormalize assigned units or infer mixing shares from raw reads."""
    total = int(metric.get('library_units', 0))
    color_for_scope = color_for_scope or scope_color
    rows = []
    for sublib in metric.get('sub_libraries', []):
        scope = sublib.get('scope_id') or sublib['scope_name']
        color = color_for_scope(scope)
        rows.append((str(sublib.get('scope_name') or scope), int(sublib['library_units']), color))
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
                 [('Ambiguous', ambiguous, '#717d92'), ('Unassigned', unassigned, '#aab5c4')] if n]
    elif remaining:
        rows.append(('Unassigned / ambiguous', remaining, '#8794a7'))
    return total, rows


def percent(count, total):
    value = 100 * count / total if total else 0
    return '<0.1%' if 0 < value < .05 else f'{value:.1f}%'


def unit_count(count):
    return f'{count:,} unit' if count == 1 else f'{count:,} units'


def mixture_chart(total, rows, *, folded=()):
    esc = html.escape
    description = '; '.join(f'{label}: {percent(n, total)}, {unit_count(n)}' for label, n, _ in [*rows, *folded])
    parts = ['<style>' + CSS + '</style>',
             '<div class="mixture-bar" role="img" aria-label="' + esc(description, quote=True) + '">']
    for i, (label, n, color) in enumerate([*rows, *folded]):
        if not n:
            continue
        width = 100 * n / total
        title = esc(f'{label}: {percent(n, total)} ({unit_count(n)})', quote=True)
        ink = bar_ink(color)
        label_class = 'mixture-segment-text'
        if color in NPG[:4] and ink != '#ffffff':
            ink = '#ffffff'
            label_class += ' outlined'
        segment_class = 'mixture-segment' if i < len(rows) else 'mixture-segment mixture-deferred'
        parts.append(f'<span class="{segment_class}" style="width:{width:.10f}%;background:{color};color:{ink}" title="{title}" aria-hidden="true"><span class="{label_class}">{esc(percent(n, total))}</span></span>')
    parts.append('</div><div class="mixture-axis" aria-hidden="true"><span>0%</span><span>100%</span></div>')

    def legend(items):
        entries = ''.join(
            f'<li><span class="mixture-swatch" style="background:{color}" aria-hidden="true"></span>'
            f'<span class="mixture-label">{esc(label)} <strong>{esc(percent(n, total))}</strong>'
            f'<span class="mixture-count">{unit_count(n)}</span></span></li>' for label, n, color in items)
        return '<ul class="mixture-legend">' + entries + '</ul>' if entries else ''

    parts.append(legend(rows))
    if folded:
        parts.append('<details class="mixture-other"><summary>Unresolved indexes</summary>' + legend(folded) + '</details>')
    return ''.join(parts)


def render_sublibrary_composition(metric):
    if len(metric.get('sub_libraries', [])) < 2:
        return ''
    total, rows = mixture_counts(metric)
    if not total:
        return ''
    return ('<section class="report-section sublibrary-mixture" id="sublibrary-composition">'
            '<h2>Sublibrary composition</h2>'
            f'<p class="muted">Share of all {total:,} reconstructed library units, including adapter dimers and unresolved index assignments. Counts are library units, not raw ONT reads.</p>'
            + mixture_chart(total, rows) + '</section>')
