"""Legacy mitochondrial profile metadata and saved-report fallback."""
import html
import json
import math


def coverage_html(records, sizes, mitochondrial, mapq, metric=None, primary=None):
    """Bin alignment midpoints; None records means positions were not retained.

    Keep mitochondrial tracks independent of the nuclear region/collision metrics.
    Saved reports without positions show only their measured mitochondrial fraction.
    """
    names = sorted(name for name in mitochondrial if name in sizes and sizes[name] > 0)
    if not names:
        return ''
    metric = metric or {}
    escape = lambda value: html.escape(str(value), quote=True)
    parts = ["<div class='mitochondrial-coverage' style='margin:16px 0 4px;padding:14px;border-top:1px solid #dce5ef'>"]
    if records is None:
        fraction = metric.get('mitochondrial_fraction')
        label = f'{fraction:.1%}' if fraction is not None else 'Unavailable'
        parts.append(f"<strong>{escape(' / '.join(names))} · Mitochondrial fraction: {label}</strong>"
                     "<p>Position data unavailable in this saved report.</p></div>")
        return ''.join(parts)
    eligible = [r for r in records.values() if mapq <= r['mapq'] < 255]
    allowed = set(primary if primary is not None else sizes) | set(names)
    denominator = sum(r['chrom'] in allowed for r in eligible)
    for name in names:
        size = sizes[name]
        step = max(1, math.ceil(size / 80))
        bins = [0] * math.ceil(size / step)
        for record in eligible:
            if record['chrom'] != name:
                continue
            position = min(size - 1, max(0, (record['start'] + record['end']) // 2))
            bins[position // step] += 1
        count = sum(bins)
        fraction = f'{count / denominator:.1%}' if denominator else '—'
        metadata = dict(chromosome=name, length_bp=size, bin_size_bp=step, midpoint_counts=bins,
                        mapped_alignments=count, denominator=denominator, mapq=mapq)
        parts.append(f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 900 175' role='img' aria-label='{escape(name)} mitochondrial coverage' style='width:100%;font-family:Arial,Helvetica,sans-serif'>"
                     f"<title>{escape(name)} mitochondrial coverage</title>"
                     f"<metadata>{escape(json.dumps(metadata))}</metadata>"
                     "<rect width='900' height='175' fill='white'/>"
                     f"<text x='50' y='24' font-size='16' fill='#203040'>{escape(name)} · {count:,} alignments · {fraction}</text>"
                     "<path d='M50 40V120H875' stroke='#536475' fill='none'/>")
        peak = max(1, max(bins))
        for i, value in enumerate(bins):
            x = 50 + 825 * i * step / size
            width = 825 * min(step, size - i * step) / size
            height = 70 * value / peak
            parts.append(f"<rect x='{x:.3f}' y='{120-height:.3f}' width='{width:.3f}' height='{height:.3f}' fill='#a82030' fill-opacity='.65'/>")
        for y, value in [(120, 0), (50, peak)]:
            parts.append(f"<text x='40' y='{y+4}' text-anchor='end' font-size='12' fill='#536475'>{value}</text>")
        for fraction_x in [0, .25, .5, .75, 1]:
            x = 50 + 825 * fraction_x
            parts.append(f"<text x='{x}' y='139' text-anchor='middle' font-size='12' fill='#536475'>{size*fraction_x/1000:g}</text>")
        parts.append(f"<text x='50' y='164' font-size='12' fill='#536475'>Position (kb) · alignment midpoints per {step:,} bp bin · separate coordinate scale</text></svg>")
    parts.append('</div>')
    return ''.join(parts)
