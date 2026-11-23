"""Micro-C overview and supporting figures for live and portable reports."""
from dataclasses import dataclass
import html
import json
import re
import shutil
import xml.etree.ElementTree as ET
from xml.parsers import expat

from microc_pairtools import report_parts
from microc_ligation import CATEGORY_COLORS
from report_microc import native_contact_summary, write_microc_comparisons
from report_metrics import fmt, METRIC_META, scoped_table_html
from report_palette import scope_color
from report_svg import text


COMPARISON_GROUPS = {
    'pairtools': ('Pairtools pair types', {'pairtools': 'Pairtools pair types'}),
    'ligation': ('Micro-C ligation structure', {
        'ligation': 'Distribution of inferred genomic block counts',
        'ligation_lengths': 'Adapter-bounded insert length by inferred block count',
    }),
    'mnase': ('Inferred MNase fragment sizes', {'mnase': 'Inferred MNase fragment sizes'}),
}


@dataclass
class MicrocProfiles:
    charts: dict
    notes: dict
    details: dict

    def figures(self, section):
        return [(self.charts[key], label) for key, label in COMPARISON_GROUPS[section][1].items()
                if key in self.charts]

    @property
    def section_ids(self):
        return tuple(section for section in COMPARISON_GROUPS if self.figures(section))


def collect_ligation_assets(out, sample, scopes):
    source = out / 'plots/ligation'
    assets = out / 'reports/assets'
    summary = out/'qc/ligation'/f'{sample}.json'
    observed = json.loads(summary.read_text()).get('scopes') if summary.is_file() else None

    def sync(src, dest, produced=True):
        target = assets / dest
        if produced and (source / src).is_file():
            shutil.copyfile(source / src, target)
        elif observed is not None:
            # This is a report copy from an earlier run. Source plots and
            # archival assets lacking an authoritative scope list stay intact.
            target.unlink(missing_ok=True)

    for scope in scopes:
        sid = scope['scope_id']
        produced = sid == 'total' or observed is None or sid in observed
        prefix = f'{sample}.{sid}'
        for src, dest in [('.svg', '.ligation.svg'), ('.mnase.svg', '.mnase.svg'),
                          ('.single_block.svg', '.single_block.svg')]:
            sync(prefix + src, prefix + dest, produced)
    sync(f'{sample}.examples.svg', f'{sample}.ligation_examples.svg')


def plot(url, label, caption=''):
    url = html.escape(url, quote=True)
    return (f"<figure class='plot-card'><a href='{url}' target='_blank' rel='noopener'>"
            f"<img src='{url}' alt='{html.escape(label, quote=True)}'></a>{caption}</figure>")


# Expat identifies the element. These tokens locate its quoted source attribute without serializing the SVG.
_XML_TAG_NAME = re.compile(rb'<[^\s/>]+')
_XML_ATTRIBUTE = re.compile(rb"""\s+([^\s=/>]+)\s*=\s*(["'])(.*?)\2""", re.DOTALL)
_FILL = re.compile(r'(?:^|;)\s*fill\s*:\s*(#[0-9a-fA-F]{6})(?=\s*(?:;|$))')
# Older refreshed reports retain Matplotlib's patch_N ids with a former scope color.
# A box in this exporter is a closed rectangle, including one with zero-height IQR.
_BOX_RECTANGLE = re.compile(r'\s*M\s+(\S+)\s+(\S+)\s+L\s+(\S+)\s+\2\s+L\s+\3\s+(\S+)'
                            r'\s+L\s+\1\s+\4\s+L\s+\1\s+\2\s+[zZ]\s*')


def color_retained_ligation(path, color):
    """Preserve the original SVG bytes outside the box face's style attribute."""
    if not path.is_file():
        return
    source = path.read_bytes()
    parser = expat.ParserCreate(namespace_separator='}')
    groups, edits = [], []

    def start(tag, attributes):
        name = tag.rsplit('}', 1)[-1]
        parent = groups[-1] if groups else ''
        in_box_axes = len(groups) >= 2 and groups[-2] == 'axes_2'
        groups.append(attributes.get('id', '') if name == 'g' else '')
        if name != 'path' or not in_box_axes:
            return
        marked = re.fullmatch(r'ligation-box-\d+', parent)
        if not marked and not re.fullmatch(r'patch_\d+', parent):
            return
        style = attributes.get('style', '')
        match = _FILL.search(style)
        if not match or match[1].lower() == color.lower():
            return
        shape = attributes.get('d', '')
        if not marked and (not attributes.get('clip-path') or not shape.rstrip().endswith(('z', 'Z'))
                           or (match[1].lower() not in CATEGORY_COLORS and not _BOX_RECTANGLE.fullmatch(shape))
                           or not re.search(r'(?:^|;)\s*stroke\s*:\s*#000000(?=\s*(?:;|$))', style)):
            return
        # The old exporter uses UTF-8. Leave unfamiliar encodings intact.
        tag_name = _XML_TAG_NAME.match(source, parser.CurrentByteIndex)
        if tag_name is None:
            return
        pos = tag_name.end()
        while attr := _XML_ATTRIBUTE.match(source, pos):
            if attr[1] == b'style':
                replacement = style[:match.start(1)] + color + style[match.end(1):]
                escaped = html.escape(replacement, quote=True).encode('utf-8')
                edits.append((*attr.span(3), escaped))
                break
            pos = attr.end()

    parser.StartElementHandler = start
    parser.EndElementHandler = lambda tag: groups.pop()
    parser.Parse(source, True)
    if edits:
        out = bytearray(source)
        for first, last, replacement in reversed(edits):
            out[first:last] = replacement
        path.write_bytes(out)


def mnase_table(scopes, result):
    rows = []
    for scope in scopes:
        sid = scope['scope_id']
        stats = result if sid == 'total' else result['scopes'].get(sid, {})
        cells = [html.escape(scope.get('scope_name', sid)), fmt(stats.get('mnase_resolved_fragments')),
                 fmt(stats.get('mnase_contributing_molecules')), fmt(stats.get('mnase_median_bp'))]
        for key in ('short_fraction', 'mono_range_fraction', 'intermediate_fraction', 'long_fraction'):
            value = stats.get('mnase_' + key)
            cells.append(f'{value:.1%}' if value is not None else 'NA')
        rows.append('<tr>' + ''.join(f'<td>{v}</td>' for v in cells) + '</tr>')
    return ("<div class='table-scroll'><table><tr><th>Library</th><th>Resolved fragments</th>"
            "<th>Raw molecules</th><th>Median bp</th><th>&lt;100 bp</th><th>100–200 bp</th>"
            "<th>200–300 bp</th><th>≥300 bp</th></tr>" + ''.join(rows) + '</table></div>')


def prepare_microc_profiles(sample, scopes, data_dir, assets, asset_url, data_url):
    comparisons = write_microc_comparisons(sample, scopes[1:] or scopes, data_dir, assets)
    if not comparisons:
        return None
    notes, details = {}, {s['scope_id']: [] for s in scopes}
    if 'pairtools' in comparisons:
        pairtools = json.loads((data_dir / 'pairtools' / f'{sample}.json').read_text())
        known = {s['scope_id'] for s in pairtools['scopes']}
        available = [s for s in scopes if s['scope_id'] in known]
        intro, panels, methods = report_parts(pairtools, assets, available, data_url + '/pairtools', asset_url,
                                              data_dir / 'pairtools')
        notes['pairtools'] = intro + methods
        for sid, panel in panels.items():
            details[sid].append('<h4>Pairtools pair types</h4>' + panel)
    def links(files):
        anchors = []
        for suffix, label in files:
            name = f'{sample}.{suffix}'
            if (data_dir / 'ligation' / name).is_file():
                url = html.escape(f'{data_url}/ligation/{name}', quote=True)
                anchors.append(f"<a href='{url}'>{label}</a>")
        return '<p>' + ' · '.join(anchors) + '</p>' if anchors else ''

    has_ligation = 'ligation' in comparisons or 'ligation_lengths' in comparisons
    if has_ligation:
        for scope in scopes[1:]:
            sid = scope['scope_id']
            color_retained_ligation(assets / f'{sample}.{sid}.ligation.svg', scope_color(sid))
        table = ''
        if any(s.get('ligation_screened_units') is not None for s in scopes):
            table = scoped_table_html(scopes, [k for k in METRIC_META if k.startswith('ligation_')], {})
        notes['ligation'] = (
            '<p>Whole-insert split alignments estimate the minimum number of resolved genomic blocks per '
            'non-dimer library unit. The denominator includes unresolved units. Adjacent and short fragments '
            'can appear as a single block. These counts cannot establish biochemical ligation efficiency.</p>' + table +
            links([('units.tsv.gz', 'Per-unit insert lengths and inferred block counts'), ('segments.tsv.gz', 'Genomic blocks and query coordinates')]))
    if 'mnase' in comparisons:
        summary = data_dir / 'ligation' / f'{sample}.mnase.json'
        notes['mnase'] = ''
        if summary.is_file():
            mnase = json.loads(summary.read_text())
            notes['mnase'] = (mnase_table(scopes, mnase) + '<p>' + html.escape(mnase.get('limitations', '')) +
                              ' Size intervals include the lower bound and exclude the upper bound.</p>')
        notes['mnase'] += links([('mnase_fragments.tsv.gz', 'Fragment lengths and parent molecules'),
                                 ('mnase.json', 'MNase size metrics')])

    for sid, parts in details.items():
        def figure(suffix, label, caption=''):
            name = f'{sample}.{sid}.{suffix}.svg'
            return plot(asset_url + name, label, caption) if (assets / name).is_file() else ''

        if has_ligation:
            lig = figure('ligation', 'Genomic blocks per library unit and insert sizes')
            examples = f'{sample}.ligation_examples.svg'
            if sid == 'total' and (assets / examples).is_file():
                lig += plot(asset_url + examples, 'Representative query block structures')
            if lig:
                parts.append('<h4>Micro-C ligation structure</h4>' + lig)
        if 'mnase' in comparisons:
            traces = (figure('mnase', 'Constituent fragment lengths',
                             '<figcaption>Reference spans of individual blocks in resolved multi-block units.</figcaption>') +
                      figure('single_block', 'Whole insert lengths in single-block units',
                             '<figcaption>Whole inserts resolved as one block, including indistinguishable junctions.</figcaption>'))
            if traces:
                parts.append("<h4>MNase fragment and single-block sizes</h4><div class='plots size-plots'>" + traces + '</div>')
        contacts = (figure('contact_funnel', 'Valid-contact funnel') +
                    figure('contact_distances', 'Contact counts by genomic distance'))
        if contacts:
            parts.append("<h4>Native contact filtering and distance</h4><div class='plots'>" + contacts + '</div>')
    return MicrocProfiles(comparisons, notes, {k: ''.join(v) for k, v in details.items()})


def write_overview_figure(sample, assets, size_name, gel_name, profiles, scopes):
    """Embed each vector panel so the downloaded figure works by itself."""
    width = 1060
    panels = [text(24, 37, sample + ' · Micro-C QC', 'start', 24, '#263944')]

    def add_row(sources, y, gap=0):
        sizes = [tuple(map(float, ET.fromstring(source).attrib['viewBox'].split()[2:])) for source in sources]
        scale = (width - gap * (len(sizes) - 1)) / sum(w for w, h in sizes)
        x = 0
        for source, (w, h) in zip(sources, sizes):
            position = f"<svg x='{x:.2f}' y='{y:.2f}' width='{w * scale:.2f}' height='{h * scale:.2f}' "
            panels.append(source.replace('<svg ', position, 1))
            x += w * scale + gap
        return y + max(h for w, h in sizes) * scale + 12

    y = add_row([(assets / size_name).read_text(), (assets / gel_name).read_text()], 56)
    for section in profiles.section_ids:
        sources = [(assets / name).read_text() for name, _ in profiles.figures(section)]
        y = add_row(sources, y, gap=12)
    native = native_contact_summary(scopes)
    if native:
        y = add_row([native], y)
    name = f'{sample}.microc_qc.svg'
    with (assets / name).open('w') as output:
        output.write(f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {y:.2f}' "
                     "font-family='Arial,Helvetica,sans-serif' role='img'><title>Micro-C QC overview</title>"
                     "<rect width='100%' height='100%' fill='white'/>" + ''.join(panels) + '</svg>')
    return name


def overview_html(asset_url, profiles):
    sections = []
    for key in profiles.section_ids:
        title = COMPARISON_GROUPS[key][0]
        figures = ''.join(plot(asset_url + name, label) for name, label in profiles.figures(key))
        if key == 'ligation':
            figures = f"<div class='ligation-comparisons'>{figures}</div>"
        sections.append(f"<div class='microc-comparison' id='{key}'>" + figures +
                        f"<details><summary>{title} · counts and methods</summary>{profiles.notes[key]}</details></div>")
    return ''.join(sections)
