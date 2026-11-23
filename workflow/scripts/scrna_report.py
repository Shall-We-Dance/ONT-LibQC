"""Offline aggregate report for shallow single-cell libraries."""
import base64
import html
import io
from pathlib import Path

from report_brand import MARK, brand_bar, favicon
from report_provenance import analysis_label
from rna_report import FAMILIES
import scrna_metaplot_plot
import scrna_regions
from scrna_composition import sublibrary_composition
from scrna_library_report import library_identity, naming_assets
from species_composition import render as render_species
from scrna_chemistry import PARSE, scope_label
from scrna_chemistry_report import CSS as CHEMISTRY_CSS, render_chemistry

CSS = """
.sc-report .report-hero {display:flex;align-items:center;justify-content:space-between;gap:1rem 2rem;flex-wrap:wrap}
.sc-report .report-hero h1 {margin-bottom:.15rem}
.sc-report .report-hero small {display:block;margin-top:.4rem;color:var(--muted);font-size:.77rem}
.sc-overview {display:flex;flex-wrap:wrap;gap:.75rem 1.5rem;margin:0;font-variant-numeric:tabular-nums}
.sc-overview>div {padding-left:1.1rem;border-left:1px solid var(--line)}
.sc-overview dt {color:var(--muted);font-size:.8rem}
.sc-overview dd {margin:0;color:var(--ink);font-size:1.55rem;font-weight:600;letter-spacing:-.035em}
.sc-report .report-nav {flex-wrap:nowrap;overflow-x:auto;scrollbar-width:thin}
.sc-report .report-nav a {white-space:nowrap}
.sc-grid {display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:.5rem}
.sc-stat {padding:.8rem 1rem;background:var(--wash);border-radius:5px}
.sc-stat strong {display:block;font-size:1.55rem;font-variant-numeric:tabular-nums}
.sc-stat span {color:var(--muted)}
.sc-plots {display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1rem}
.sc-plots figure {min-width:0;margin:0;border:1px solid var(--line);border-radius:6px;padding:.85rem}
.sc-plots h3 {margin:0 0 .4rem;color:var(--ink)}
.sc-plots img {display:block;width:100%;height:auto}
.sc-plots figcaption {color:var(--muted);font-size:.92rem}
.sc-family-figure {margin:0;max-width:52rem}
.sc-family-figure img {width:100%;height:auto;display:block}
.sc-region-plot {margin:1rem 0;max-width:46rem;break-inside:avoid;font-variant-numeric:tabular-nums}
.sc-region-track {display:flex;height:1rem;border-radius:3px;overflow:hidden}
.sc-region-labels {display:grid;grid-template-columns:repeat(3,1fr);gap:.5rem;margin-top:.65rem}
.sc-region-labels b,.sc-region-labels small {display:block}
.sc-region-labels b::before {content:"";display:inline-block;width:.6em;height:.6em;border-radius:2px;margin-right:.4em;background:var(--region-color)}
.sc-region-labels small {color:var(--muted)}
.sc-metaplot {margin:0;break-inside:avoid}
.sc-metaplot img {width:100%;height:auto;display:block}
.sc-metaplot figcaption {display:flex;align-items:baseline;justify-content:space-between;gap:.5rem 1rem;flex-wrap:wrap;color:var(--muted);font-size:.85rem}
@media screen and (min-width:960px) {
  .sc-section-pair {display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1rem;margin:18px 0}
  .sc-section-pair>.report-section {min-width:0;margin:0}
}
@media print {
  .sc-section-pair>.report-section {break-inside:avoid}
  .sc-metaplot a[download] {display:none}
}
.sc-table {width:100%;border-collapse:collapse;font-variant-numeric:tabular-nums;margin:.6rem 0}
.sc-table :is(th,td) {padding:.5rem .75rem;border-bottom:1px solid var(--line);text-align:right}
.sc-table :is(th,td):first-child {text-align:left}
.sc-table thead {background:var(--wash)}
.sc-wrap {max-width:100%;overflow:auto}
.sc-report #species-composition .table-scroll {max-width:100%;overflow-x:auto}
@media screen and (max-width:680px) {
  .sc-report .report-hero {display:block}
  .sc-overview {display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:.5rem;border-top:1px solid var(--line);padding-top:.65rem;margin-top:.75rem}
  .sc-overview>div {padding-left:.5rem}
  .sc-overview>div:first-child {padding-left:0;border-left:0}
  .sc-overview dt {font-size:.72rem}
  .sc-overview dd {font-size:1.27rem}
  .sc-report #species-composition .species-counts th:first-child {min-width:6.5em}
  .sc-report #species-composition .species-counts th[scope="col"] {min-width:5.5em}
}
.sc-report details {margin-top:1rem}
.sc-report summary {color:var(--blue);cursor:pointer}
@media(max-width:680px){.sc-plots{grid-template-columns:1fr}.sc-table :is(th,td){padding:.35rem}.sc-grid{grid-template-columns:repeat(2,minmax(0,1fr))}}
@media print {body{max-width:none;margin:0;padding:0}.report-section,.report-hero{box-shadow:none}figure,tr,.sc-stat{break-inside:avoid}.sc-plots img{max-height:70mm}.sc-table{font-size:9pt}.sc-wrap{overflow:visible}.report-section{break-inside:auto}}
"""


def image_uri(content):
    return 'data:image/svg+xml;base64,' + base64.b64encode(content).decode('ascii')


def percent(n, d):
    return f'{100 * n / d:.1f}%' if d else '—'


def table(headers, rows):
    head = ''.join('<th scope="col">' + html.escape(h) + '</th>' for h in headers)
    body = []
    for fields in rows:
        first, *rest = [html.escape(f'{v:,}' if isinstance(v, int) else str(v)) for v in fields]
        body.append('<tr><th scope="row">' + first + '</th>' + ''.join('<td>' + v + '</td>' for v in rest) + '</tr>')
    return '<div class="sc-wrap"><table class="sc-table"><thead><tr>' + head + '</tr></thead><tbody>' + ''.join(body) + '</tbody></table></div>'


def plot(metrics, rank):
    from matplotlib.backends.backend_svg import FigureCanvasSVG
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator, ScalarFormatter
    fig = Figure(figsize=(5.4, 3.0), layout='constrained')
    ax = fig.subplots()
    if rank:
        ranks = metrics['ranks']
        ax.step(range(1, len(ranks) + 1), ranks, where='post', color='#266c9e', linewidth=2, marker='.' if len(ranks) == 1 else None)
        ax.set(xlabel='Barcode combination rank', ylabel='Observed UMIs', xscale='log', yscale='log')
        if max(ranks) < 15:
            ax.set_yticks(sorted(set(ranks)))
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.minorticks_off()
    else:
        points = metrics['rarefaction']
        for key, label, color in (('molecules', 'UMIs', '#266c9e'), ('genes', 'Genes', '#168472')):
            ax.plot([p['units'] for p in points], [p[key] for p in points], label=label, color=color, marker='o', markersize=3)
        ax.set(xlabel='Gene-assigned units sampled', ylabel='Expected number observed', xlim=(0, None), ylim=(0, None))
        ax.legend(frameon=False)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    ax.spines[['right', 'top']].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.grid(axis='y', alpha=.17)
    buf = io.BytesIO()
    FigureCanvasSVG(fig).print_svg(buf, metadata={'Date': None})
    return image_uri(buf.getvalue())


FAMILY_LABELS = {'protein_coding': 'Protein-coding', 'other': 'Other annotated RNA',
                 'unknown': 'Unclassified', 'ambiguous_family': 'Ambiguous family'}


def family_plot(cohorts):
    """Each series uses its own input denominator, including unclassified fragments."""
    from matplotlib.backends.backend_svg import FigureCanvasSVG
    from matplotlib.figure import Figure
    from matplotlib.ticker import PercentFormatter
    groups = [(cohorts[k], label, color) for k, label, color in
              [('all', 'All cDNA', '#316a9a'), ('multimapped', 'Multimapped', '#bc662d')] if cohorts[k]['units']]
    if not groups:
        return None
    keys = [key for key in FAMILIES if cohorts['all']['counts'].get(key)]
    fig = Figure(figsize=(7.8, 1.15 + 0.53 * len(keys)), layout='constrained')
    ax = fig.subplots()
    for i, (group, label, color) in enumerate(groups):
        den = group['units']
        rates = [100 * group['counts'].get(k, 0) / den for k in keys]
        positions = [j + (i - (len(groups) - 1) / 2) * .34 for j in range(len(keys))]
        ax.barh(positions, rates, height=.29, color=color, label=f'{label} (n = {den:,})')
        for pos, rate in zip(positions, rates):
            if rate:
                value = '<0.1%' if rate < .1 else f'{rate:.1f}%'
                ax.text(rate + 1, pos, value, fontsize=8, va='center', color=color)
    ax.set_yticks(range(len(keys)), [FAMILY_LABELS.get(k, k) for k in keys])
    ax.invert_yaxis()
    ax.set_xlim(0, 112)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_axisbelow(True)
    ax.grid(axis='x', alpha=.17)
    ax.tick_params(length=0, labelsize=9)
    ax.spines[:].set_visible(False)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1), ncol=len(groups), frameon=False, fontsize=9)
    buf = io.BytesIO()
    FigureCanvasSVG(fig).print_svg(buf, metadata={'Title': 'RNA family distribution', 'Date': None})
    return buf.getvalue()


def family_section(data, figure, audit=None):
    cohorts = data['cohorts']
    def value(key, group):
        n = group['counts'].get(key, 0)
        return f'{n:,} ({percent(n, group["units"])})' if group['units'] else '—'
    rows = [(FAMILY_LABELS.get(k, k), value(k, cohorts['all']), value(k, cohorts['multimapped'])) for k in FAMILIES]
    uri = image_uri(figure)
    content = '<figure class="sc-family-figure"><img alt="RNA family percentages for all qualified cDNA and the multimapped subset" src="' + uri + '"></figure>'
    content += '<p>Each fragment is counted once in its group. Conflicting gene families across alignments are grouped as ambiguous.</p>'
    content += '<details><summary>Counts and classification</summary>' + table(['RNA family', 'All qualified cDNA', 'Multimapped cDNA'], rows)
    if data['panel']['status'] == 'complete':
        content += '<p>Ribosomal sequence matches cover at least 80% of the fragment at 85% identity or higher, with at least 30 matching bases.</p>'
    else:
        content += '<p>The independent rRNA sequence screen was unavailable. Ribosomal assignments use the reference annotations.</p>'
    if data['repeatmasker_sha256']:
        content += '<p>A fragment can also be assigned to rRNA when at least half of its aligned bases overlap annotated rRNA repeats at every reported locus.</p>'
    content += '<p><a download="rna_families.svg" href="' + uri + '">Download figure (SVG)</a></p></details>'
    if cohorts['multimapped']['units']:
        labels = [('sequence_rRNA', 'Mature rRNA sequence match'), ('annotated_rRNA', 'rRNA gene annotation at all loci'),
                  ('repeat_rRNA_all_loci', 'rRNA repeat at all loci'), ('rRNA_some_loci', 'rRNA repeat at some loci'),
                  ('other_repeat_all_loci', 'Other repeat at all loci'), ('other_repeat_some_loci', 'Other repeat at some loci'),
                  ('unresolved', 'Origin unresolved')]
        counts, total = data['multimapping_evidence'], cohorts['multimapped']['units']
        content += '<details><summary>Evidence for multimapping</summary>' + table(['Strongest evidence', 'Fragments', '% of multimapped cDNA'],
            [(label, counts[k], percent(counts[k], total)) for k, label in labels if counts.get(k)])
        if audit and audit.get('panel_subunits'):
            subunits = ', '.join(f'{n:,} {html.escape(name)}' for name, n in sorted(audit['panel_subunits'].items()))
            content += f'<p>Mature rRNA matches included {subunits}. '
            content += f'{audit["host_panel_matches"]:,} also match the host-only rRNA panel.</p>'
        content += '</details>'
    if audit:
        n = audit['input']
        if n != cohorts['multimapped']['units']:
            raise ValueError('Terminal audit and RNA family plot have different multimapped inputs')
        content += '<details><summary>Terminal sequence and remapping checks</summary>'
        content += (f'<p>Before cDNA extraction, {audit["barcode_before"]:,} of {n:,} fragments matched their own barcode and adjacent sequence. '
                    f'After extraction, {audit["any_signature"]:,} retained a terminal adapter, template-switch or linked-barcode match.</p>')
        content += table(['Additional trimming', 'Input', 'Still multimapped', 'Unique', 'Other mapping outcome', 'Below 30 bp'],
                         [[v['label'], v['input'], v['multimapped'], v['unique'], v['other'], v['too_short']] for v in audit['experiments']])
        content += '<p>Both checks used the original STARlong settings and genome index.</p></details>'
    return content


def region_section(metrics, figure, mapq):
    den = metrics['coding_region_bases']
    total = metrics['aligned_bases']
    if not total:
        return '<p>No uniquely mapped cDNA bases passed the alignment filters.</p>'
    content = ''
    if den:
        bars, labels = [], []
        for key, label, color in scrna_regions.REGIONS:
            n = metrics['counts'][key]
            bars.append(f'<span style="width:{n / den:.6%};background:{color}"></span>')
            labels.append(f'<span style="--region-color:{color}"><b>{label}</b>{percent(n, den)} <small>{n:,} bp</small></span>')
        content = ('<figure class="sc-region-plot" aria-label="Base composition in coding transcript regions">'
                   '<div class="sc-region-track" aria-hidden="true">' + ''.join(bars) + '</div>'
                   '<figcaption class="sc-region-labels">' + ''.join(labels) + '</figcaption></figure>')
        content += (f'<p>{den:,} of {total:,} aligned bases ({percent(den, total)}) fall in these three regions, '
                    f'across {metrics["unique_units"]:,} uniquely mapped cDNA fragments (MAPQ ≥{mapq}).</p>')
    else:
        content += f'<p>None of the {total:,} uniquely aligned bases overlap an annotated UTR or CDS.</p>'
    fields = [(k, label) for k, label, _ in scrna_regions.REGIONS] + list(scrna_regions.OTHERS)
    content += '<details><summary>All aligned bases</summary>'
    content += table(['Region', 'Aligned bases', '% of all aligned bases'],
                     [(label, metrics['counts'][k], percent(metrics['counts'][k], total)) for k, label in fields])
    content += '<p>Each aligned base is counted once. UTRs follow the transcript strand. CDS includes annotated start and stop codons and takes precedence where isoforms overlap.</p>'
    if len(metrics['scopes']) > 1:
        content += table(['Library', '5′ UTR', 'CDS', '3′ UTR', 'Bases in these regions'],
                         [[scope, *(percent(group['counts'][k], group['coding_region_bases']) for k, _, _ in scrna_regions.REGIONS),
                           group['coding_region_bases']] for scope, group in metrics['scopes'].items()])
    if figure:
        content += '<p><a download="transcript_regions.svg" href="' + image_uri(figure) + '">Download figure (SVG)</a></p>'
    return content + '</details>'


def metaplot_section(metrics, figures, mapq):
    if 'all' not in figures:
        return '<p>No eligible cDNA fragments overlap representative transcript exons.</p>'
    scopes = metrics['scopes']

    def chart(scope):
        profile = scopes[scope]
        body = profile['gene_body']
        label = html.escape(scope_label(scope))
        axis = '5′ UTR, CDS and 3′ UTR' if profile['CDS']['contributing_genes'] else 'spliced transcript, 5′ to 3′'
        alt = f'{label} mean relative coverage by transcript position ({axis})'
        svg = image_uri(figures[scope])
        basename = 'transcript_metaplot' if scope == 'all' else 'transcript_metaplot_' + ''.join(
            ch if ch.isascii() and ch.isalnum() else '_' for ch in scope)
        content = (f'<figure class="sc-metaplot"><img src="{svg}" alt="{alt}"><figcaption>'
                   f'<span>{" · ".join(scrna_metaplot_plot.support_counts(body))}</span>'
                   f'<a href="{svg}" download="{basename}.svg" aria-label="Download {label} metaplot as SVG">Download SVG</a>'
                   '</figcaption></figure>')
        if profile['CDS']['contributing_genes']:
            rows = [(region, profile[key]['covered_genes'], profile[key]['contributing_genes'],
                     profile[key]['supporting_units']) for key, region, _ in scrna_regions.REGIONS]
            content += '<details><summary>Coverage by region</summary>' + table(
                ['Region', 'Genes with coverage', 'Genes with region', 'Fragments'], rows)
            content += ('<p>An annotated region with no coverage contributes zero to the mean. '
                        'A fragment can overlap more than one region.</p></details>')
        return content

    content = chart('all')
    others = [scope for scope in figures if scope != 'all']
    if len(others) == 1 and scopes[others[0]] == scopes['all']:
        others = []
    for scope in others:
        name = html.escape(scope_label(scope))
        content += '<details><summary>' + name + '</summary>' + chart(scope) + '</details>'
    content += (f'<details><summary>How coverage is calculated</summary><p>Each gene contributes equally. '
                f'Coverage uses the longest spliced transcript per gene and uniquely mapped cDNA (MAPQ ≥{mapq}) '
                'assigned to one gene on the expected strand. Ribosomal RNA is excluded. A value of 1 is the '
                'gene’s average exonic coverage. Shading shows one standard error across genes.</p>')
    bias = scopes['all'].get('bias', {})
    if bias.get('status') == 'insufficient_terminal_supported_genes':
        content += (f'<p>5′/3′ bias needs at least {bias["minimum_genes"]} genes with '
                    f'{bias["minimum_assigned_units_per_gene"]} or more fragments and coverage at either end. '
                    f'Genes meeting these criteria: {bias["eligible_terminal_supported_genes"]}.</p>')
    return content + '</details>'


def render(s, family_figure=None, region_figure=None, metaplot_figures=None):
    esc = html.escape
    m, scopes = s['metrics'], s['scopes']
    mark = image_uri(MARK.read_bytes())
    css = Path(__file__).with_name('report_modern.css').read_text() + CHEMISTRY_CSS
    sample = esc(s['sample'])
    provenance = s['analysis_provenance']
    version = esc(analysis_label(provenance, provenance['pipeline_version'] + ' · ' + provenance['pipeline_build']))
    eligible = any(d['kit'] for d in scopes.values())
    overview = [('ONT reads', s['raw_reads']), ('Non-dimer units', s['non_dimer_units']),
                ('Qualified cDNA', sum(m['mapping'].values()) if eligible else None)]
    overview_html = ''.join(f'<div><dt>{label}</dt><dd>{format(n, ",") if n is not None else "—"}</dd></div>'
                            for label, n in overview)
    pieces = [f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{sample} · Single-cell RNA QC (Beta)</title>{favicon(mark)}<style>{css}{CSS}</style></head><body class="sc-report">',
              brand_bar(mark, 'Single-cell RNA QC'),
              f'<header class="report-hero"><div><h1>{sample}</h1><p>Shallow scRNA-seq <span class="rna-beta">Beta</span> · {esc(s["assembly"])}</p>'
              f'<small>{version}</small></div><dl class="sc-overview">{overview_html}</dl></header>']
    navigation, sections = [], []
    short_labels = {'sublibrary-composition': 'Sublibraries', 'chemistry': 'Chemistry', 'gene-yield': 'Gene yield',
                    'unit-to-gene': 'Unit yield', 'sequencing-depth': 'Depth', 'rna-families': 'RNA families',
                    'transcript-regions': 'Transcript regions', 'transcript-metaplot': 'Metaplot'}
    def add_panel(markup, beside_previous=False):
        if beside_previous:
            markup = '<div class="sc-section-pair">' + sections.pop() + markup + '</div>'
        sections.append(markup)

    def section(title, content, anchor=None, *, beside_previous=False):
        target = f' id="{anchor}"' if anchor else ''
        if anchor:
            navigation.append((anchor, short_labels.get(anchor, title)))
        add_panel(f'<section class="report-section"{target}><h2>' + title + '</h2>' + content + '</section>', beside_previous)

    section('Library quality', library_identity(s), 'library-quality')
    composition = sublibrary_composition(s)
    if composition:
        section('Sublibrary composition', composition, 'sublibrary-composition')
    if m.get('species_composition'):
        navigation.append(('species-composition', 'Species'))
        add_panel(render_species(m['species_composition']), beside_previous=bool(composition))
    section('Chemistry detection', render_chemistry(s, table), 'chemistry')
    if not eligible:
        section('Barcode and gene QC', '<p>Barcode-resolved metrics require a supported library layout.</p>')
    else:
        stats = [('Detected genes', m['genes']), ('Distinct UMIs', m['molecules']),
                 ('Barcode combinations with a UMI', m['barcodes']), ('Combinations with ≥2 genes', m['barcodes_ge_2_genes'])]
        cards = ''.join(f'<div class="sc-stat"><strong>{v:,}</strong><span>{label}</span></div>' for label, v in stats)
        notes = []
        if any(d['kit'] == PARSE for d in scopes.values()) and not s['parse_well_map_sha256']:
            notes.append('Parse uses both random and poly(T) primers. Without a kit-specific well map, their round-1 barcodes are counted separately.')
        notes.append(f'{m["singleton_barcodes"]:,} combinations have one observed UMI. Cell recovery was not estimated.')
        if m['assigned_units']:
            repeated = m['assigned_units'] - m['molecules']
            notes.append(f'{repeated:,} of {m["assigned_units"]:,} countable gene-assigned units repeated a barcode–gene–UMI key ({percent(repeated, m["assigned_units"])}).')
        if not s['single_capture'] and any(scope in ('unassigned', 'native-unindexed') for scope in scopes):
            notes.append('Molecule counts include only units with an assigned sample index.')
        if m['molecule_assessment'] == 'no_assigned_capture':
            section('Barcode and gene yield', '<p>Only index-unassigned units passed QC. Molecule counting requires a confirmed sample or capture.</p>')
        elif m['molecule_assessment'] == 'no_qualified_cdna':
            section('Barcode and gene yield', '<p>No units passed the barcode, UMI and cDNA length filters.</p>')
        else:
            detail = '<details><summary>Counting details</summary><ul>' + ''.join(f'<li>{note}</li>' for note in notes) + '</ul></details>'
            section('Barcode and gene yield', '<div class="sc-grid">' + cards + '</div>' + detail, 'gene-yield')
        stages = [('Barcode decoded', 'decoded'), ('UMI passes sequence QC', 'valid_umi'),
                  ('cDNA ≥30 bp sent to STARlong', 'cdna'), (f'Unique alignment, MAPQ ≥{s["mapq"]}', 'unique'),
                  ('One gene on the expected strand', 'gene')]
        stage_rows, unassigned_rows = [], []
        for label, key in stages:
            ns = {}
            for scope, d in scopes.items():
                counts = d['funnel']
                n = sum(counts.get(k, 0) for k in ('barcode_exact', 'barcode_corrected')) if key == 'decoded' else counts.get(key, 0)
                ns[scope] = n
            total = sum(ns.values())
            if 'unassigned' in ns:
                unassigned_rows.append([label, ns.pop('unassigned')])
            stage_rows.append([label, f'{total:,} ({percent(total, s["non_dimer_units"])})', *ns.values()])
        info = f'<p>{s["reconstructed_parents"]:,} of {s["raw_reads"]:,} ONT reads yielded {s["non_dimer_units"]:,} non-dimer units. Percentages below use those units.</p>'
        info += table(['Stage', 'All units', *(scope_label(scope) for scope in scopes if scope != 'unassigned')], stage_rows)
        if unassigned_rows:
            info += ('<details><summary>Unassigned index</summary>' +
                     table(['Stage', 'Units'], unassigned_rows) + '</details>')
        exact = sum(d['funnel'].get('barcode_exact', 0) for d in scopes.values())
        corrected = sum(d['funnel'].get('barcode_corrected', 0) for d in scopes.values())
        info += f'<p>Decoded barcodes include {exact:,} exact matches and {corrected:,} unique one-base corrections.</p>'
        section('From library unit to gene', info, 'unit-to-gene', beside_previous=True)
        if m['molecules']:
            figures = []
            for rank, title, caption in ((True, 'Barcode rank', 'Each combination is ordered by its observed UMIs.'),
                                         (False, 'Observed-depth rarefaction', 'Subsampling uses the gene-assigned units that enter the molecule counts.')):
                figures.append(f'<figure><h3>{title}</h3><img src="{plot(m, rank)}" alt="{title}"><figcaption>{caption}</figcaption></figure>')
            section('Sequencing depth', '<div class="sc-plots">' + ''.join(figures) + '</div>', 'sequencing-depth')
        mapping = m['mapping']
        cdna = sum(mapping.values())
        mapped = cdna - mapping.get('unmapped', 0)
        text = f'<p>{mapped:,} of {cdna:,} barcode- and UMI-qualified cDNA inserts mapped to {esc(s["assembly"])} ({percent(mapped, cdna)}).'
        if m['molecule_assessment'] == 'evaluated':
            text += f' Of {m["molecules"]:,} counted UMIs, {m["exonic_molecules"]:,} have exonic evidence and {m["intronic_only_molecules"]:,} have intronic evidence only.'
        text += '</p>'
        text += table(['cDNA alignment', 'Units', '% of cDNA input'],
                      [(status, mapping.get(key, 0), percent(mapping.get(key, 0), cdna))
                       for key, status in (('unique', 'Unique, MAPQ-pass'), ('multimapped', 'Multiple loci'),
                                           ('low_mapq', 'Low MAPQ'), ('chimeric', 'Chimeric'), ('unmapped', 'Unmapped'))])
        ev, unique = m['rna_evidence'], mapping.get('unique', 0)
        families = (s['annotation'] or {}).get('gene_families', {})
        evidence_rows = [('Mitochondrial chromosome', ev.get('mitochondrial', 0)),
                         ('Spliced', ev.get('spliced', 0)),
                         ('Contains a known splice junction', ev.get('annotated_splice', 0))]
        if families.get('rRNA') or families.get('mt-rRNA'):
            evidence_rows.append(('Annotated rRNA', ev.get('annotated_rrna', 0)))
        rna = table(['Unique cDNA alignments', 'Units', '% of unique'],
                    [(label, n, percent(n, unique)) for label, n in evidence_rows])
        bases = m['aligned_region_bases']
        regions = [('Exonic', sum(bases.get(k, 0) for k in ('CDS', 'five_prime_UTR', 'three_prime_UTR', 'other_exon'))),
                   ('Intronic', bases.get('intron', 0)), ('Other genic', bases.get('genic_unannotated', 0)),
                   ('Intergenic', bases.get('intergenic', 0))]
        detail_title = 'Splicing and mitochondrial alignments'
        if not m.get('transcript_regions'):
            detail_title = 'Splicing and annotated regions'
            rna += table(['Annotated region', 'Aligned bases', '% of unique aligned bases'],
                         [(k, n, percent(n, sum(bases.values()))) for k, n in regions])
        section('RNA alignment', text + '<details><summary>' + detail_title + '</summary>' + rna + '</details>', 'rna-alignment')
        families = m.get('rna_families')
        if families and families['cohorts']['all']['units']:
            figure = family_figure or family_plot(families['cohorts'])
            section('RNA family distribution', family_section(families, figure, s.get('multimapping_audit')), 'rna-families', beside_previous=True)
        regions = m.get('transcript_regions')
        metaplot = m.get('transcript_metaplot')
        if regions:
            section('Coding transcript regions', region_section(regions, region_figure or scrna_regions.figure(regions), s['mapq']), 'transcript-regions')
        if metaplot:
            charts = metaplot_figures if metaplot_figures is not None else scrna_metaplot_plot.figures(metaplot)
            section('Transcript metaplot', metaplot_section(metaplot, charts, s['mapq']), 'transcript-metaplot', beside_previous=bool(regions))
    if len(navigation) > 1:
        links = ''.join(f'<a href="#{anchor}">{esc(label)}</a>' for anchor, label in navigation)
        pieces.append('<nav class="report-nav" aria-label="Report sections">' + links + '</nav>')
    return ''.join(pieces) + '<main>' + ''.join(sections) + '</main>' + naming_assets(s) + '</body></html>'
