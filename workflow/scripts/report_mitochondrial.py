"""Mitochondrial coverage on an independent scale within the chromosome figure."""
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


def integrate_karyotype(svg, coverage):
    """Put mitochondrial rows inside the existing SVG, before its footer.

    Preserve the nuclear drawing and embedded count metadata. Integration is
    idempotent, including when previously generated plot assets are reused.
    """
    import xml.etree.ElementTree as ET
    ns = 'http://www.w3.org/2000/svg'
    ET.register_namespace('', ns)
    ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
    root = ET.fromstring(svg)
    if not coverage:
        return svg
    source = ET.fromstring(coverage)
    tracks = [json.loads(node.text) for node in source.iter(f'{{{ns}}}metadata')]
    previous = root.find(f"{{{ns}}}g[@id='mitochondrial-tracks']")
    if previous is not None:
        root.remove(previous)
    view = [float(value) for value in root.attrib['viewBox'].split()]
    width = view[2]
    height = float(root.get('data-nuclear-height', view[3]))
    root.set('data-nuclear-height', f'{height:g}')
    extra = 54 * max(1, len(tracks))
    root.set('height', f'{height + extra:g}')
    root.set('viewBox', f'{view[0]:g} {view[1]:g} {width:g} {height + extra:g}')
    # Cairo/R writes the outer caption as a top-level group of glyph uses.
    # Move only that caption; chromosome geometry and densities stay untouched.
    for child in root:
        uses = list(child.iter(f'{{{ns}}}use'))
        if uses and all(float(u.get('y', 0)) > height - 20 for u in uses):
            child.set('transform', f'translate(0 {extra:g})')
    group = ET.SubElement(root, f'{{{ns}}}g', {'id': 'mitochondrial-tracks',
                         'font-family': 'Arial,Helvetica,sans-serif'})
    def element(tag, attrs, text=None):
        node = ET.SubElement(group, f'{{{ns}}}{tag}', {k: str(v) for k, v in attrs.items()})
        node.text = text
        return node
    left, right = width * .11, width * .955
    span = right - left
    if not tracks:
        label = ' '.join(source.itertext()).replace('Position data unavailable in this saved report.', '').strip()
        element('text', {'x': left, 'y': height, 'font-size': 7, 'fill': '#536475'},
                label + ' · Position data unavailable')
    for index, track in enumerate(tracks):
        baseline = height - 10 + index * 54
        bins = track['midpoint_counts']; size = track['length_bp']; step = track['bin_size_bp']
        peak = max(1, max(bins, default=0))
        element('metadata', {}, json.dumps(track))
        element('text', {'x': width * .055, 'y': baseline + 5, 'text-anchor': 'middle',
                         'font-size': 7, 'fill': '#34495e'}, track['chromosome'])
        element('rect', {'x': left, 'y': baseline - 20, 'width': span, 'height': 20, 'fill': '#f5f8fa'})
        points = [(left, baseline)] + [(left + span * min(size, (i + .5) * step) / size,
                   baseline - 19 * value / peak) for i, value in enumerate(bins)] + [(right, baseline)]
        path = 'M ' + ' L '.join(f'{x:.3f} {y:.3f}' for x, y in points)
        element('path', {'d': path, 'fill': '#e45756', 'fill-opacity': '.4', 'stroke': '#a82030', 'stroke-width': '.65'})
        element('rect', {'x': left, 'y': baseline + 3, 'width': span, 'height': 6,
                         'fill': '#f5f8fa', 'stroke': '#34495e', 'stroke-width': '.6'})
        for fraction in (0, .25, .5, .75, 1):
            x = left + span * fraction
            element('path', {'d': f'M{x:.3f} {baseline+9:g}v2', 'stroke': '#607080', 'stroke-width': '.4'})
            element('text', {'x': x, 'y': baseline + 16, 'text-anchor': 'middle',
                             'font-size': 5, 'fill': '#536475'}, f'{size*fraction/1000:.1f} kb')
        fraction = f"{track['mapped_alignments']/track['denominator']:.1%}" if track['denominator'] else '—'
        element('text', {'x': left, 'y': baseline + 28, 'font-size': 6.5, 'fill': '#536475'},
                f"{track['chromosome']}: independent coordinate and density scales · {track['mapped_alignments']:,} alignments ({fraction}) · {step:,} bp bins")
    return ET.tostring(root, encoding='unicode')
