"""Self-contained RNA/CLIP report fragments, with explicit support denominators."""
import html
import csv
import gzip
from collections import Counter, defaultdict
import json
from pathlib import Path
from urllib.parse import quote, urlsplit


def esc(value):
    return html.escape(str(value), quote=True)


def number(value, digits=0):
    return 'Not evaluable' if value is None else f'{value:,.{digits}f}'


def percent(value):
    return 'Not evaluable' if value is None else '&lt;0.1%' if 0 < value < .001 else f'{100 * value:.1f}%'


def table(headers, rows):
    return '<div class="table-scroll rna-table" tabindex="0"><table><thead><tr>' + ''.join('<th>' + esc(x) + '</th>' for x in headers) + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join('<td>' + esc(x) + '</td>' for x in row) + '</tr>' for row in rows) + '</tbody></table></div>'

def render_structure_evidence(preparation, scope, scope_labels, prefix):
    """Keep observed scope evidence distinct from the protocol model it tests."""
    if not preparation:
        return ''
    auto = preparation.get('auto_grip') or {}
    if not auto:
        return ''
    mode = preparation.get('mode', 'passthrough')
    parts = ['<section class="rna-evidence-block rna-structure-evidence"><h3>RNA-adapter structure and UMI evidence</h3>']
    bypass = auto.get('bypass_reason')
    if bypass:
        descriptions = {'configured_umi_extraction': 'The configured UMI structure takes precedence over automatic GRIP inference.',
                        'manual_end_trimming': 'The configured end trimming takes precedence over automatic GRIP inference.',
                        'auto_disabled': 'Automatic GRIP structure inference is disabled.'}
        parts.append('<p class="muted">' + esc(descriptions.get(bypass, str(bypass).replace('_', ' '))) + '</p>')
        if mode == 'manual':
            parts.append('<p>Configured trimming: ' + number(preparation.get('p5_trim_bp', 0)) + ' bp at P5 and ' + number(preparation.get('p7_trim_bp', 0)) + ' bp at P7. These lengths are settings, not inferred structure evidence.</p>')
    else:
        parts.append('<p class="rna-structure-model"><strong>Protocol model tested</strong><code>' + esc(auto.get('expected_structure', 'P5–N5–fixed barcode–RNA insert–N10–P7')) + '</code></p>')
        parts.append('<p class="muted">The 7-nt fixed barcode is a structure landmark; P5 N5 is a random adapter tag. Only P7 N10 is interpreted as the protocol UMI. Structure support is assessed independently within each index scope; it is not a sample-wide assumption.</p>')
    evidence = auto.get('scope_evidence', {})
    selected = evidence.items() if scope == 'total' else [(scope, evidence[scope])] if scope in evidence else []
    decisions = {'supported': 'Supported', 'insufficient_eligible_units': 'Too few eligible units',
                 'insufficient_exact_barcode_support': 'Exact-barcode support insufficient',
                 'insufficient_random_tag_diversity': 'Tag diversity insufficient',
                 'bypassed_configured_umi_extraction': 'Configured UMI structure used',
                 'bypassed_manual_end_trimming': 'Configured trimming used',
                 'bypassed_auto_disabled': 'Automatic inference disabled'}
    rows = []
    for key, values in selected:
        label = scope_labels.get(key) or ('Unassigned units' if key == 'unassigned' else 'Sublibrary ' + key)
        exact, eligible = values.get('exact_barcode_units', 0), values.get('eligible_units', 0)
        rows.append((label, number(exact) + ' / ' + number(eligible) + ' · ' + percent(values.get('exact_fraction')),
                     number(values.get('p5_random_distinct', 0)) + ' / ' + number(values.get('p7_random_distinct', 0)),
                     percent(values.get('max_joint_tag_fraction')),
                     decisions.get(values.get('decision'), str(values.get('decision', 'Not evaluable')).replace('_', ' '))))
    if rows:
        parts.append(table(('Scope', 'Exact barcode / eligible units', 'Distinct P5 / P7 tags', 'Largest tag-pair fraction', 'Structure evidence'), rows))
    else:
        parts.append('<p class="muted">No per-scope structure evidence is available for these units.</p>')
    if not bypass:
        threshold = auto.get('thresholds', {})
        parts.append('<p class="muted">Support requires ≥' + number(threshold.get('minimum_eligible_units', 20)) + ' eligible units, ≥' + number(threshold.get('minimum_exact_barcode_units', 10)) + ' exact barcode observations and ≥' + percent(threshold.get('minimum_exact_fraction', .5)) + ' exact support, with ≥' + number(threshold.get('minimum_distinct_tags_each_end', 4)) + ' distinct tags at each end and no tag pair above ' + percent(threshold.get('maximum_joint_tag_fraction', .5)) + '.</p>')
        parts.append('<p class="muted">Trimming requires both a supported scope and a recoverable barcode and tags in that unit. Unsupported or uncertain units remain untrimmed. Observed P7 tags are not error-corrected molecule counts; ONT errors and tag collisions can alter apparent UMI diversity.</p>')
    if preparation.get('units_tsv_file'):
        parts.append('<p class="rna-evidence-download"><a href="' + esc(prefix + '/preparation_units.tsv.gz') + '">Per-unit structure, trimming and UMI evidence</a></p>')
    parts.append('</section>')
    return ''.join(parts)

def render_source_evidence(data, metrics, prefix):
    """Show the independent RNA panel beside genome evidence without a diagnosis."""
    panel = data.get('source_screen') or {}
    source = metrics.get('source_evidence') or {}
    if not panel:
        return ''
    parts = ['<section class="rna-evidence-block rna-source-evidence"><h3>Selected reference and RNA-source evidence</h3>']
    if panel.get('status') != 'complete' or source.get('status') in {'not_available', 'not_screened'}:
        parts.append('<p class="muted">The RNA-source panel screen was unavailable. A zero source count must not be interpreted as absence of ribosomal RNA.</p></section>')
        return ''.join(parts)
    denominator = source.get('denominator_units', metrics['units'])
    assigned = source.get('assigned_units', 0)
    ambiguous = source.get('ambiguous_taxon_units', 0)
    status = source.get('status_counts', {})
    rows = [
        ('Selected-genome mapping', number(metrics['mapped_units']) + ' / ' + number(metrics['units']), percent(metrics['mapped_fraction']) + ' · ' + str(data.get('assembly') or 'Unspecified assembly')),
        ('RNA-panel assigned units', number(assigned) + ' / ' + number(denominator), percent(assigned / denominator if denominator else None) + ' of all recovered units'),
        ('Several candidate reference organisms', number(ambiguous) + ' / ' + number(assigned), percent(ambiguous / assigned if assigned else None) + ' of panel-assigned units'),
        ('No qualifying panel match', number(status.get('unassigned', 0)), 'A limited panel cannot establish that these units are non-ribosomal'),
        ('Not screened / unavailable', number(sum(count for state, count in status.items() if state not in {'assigned', 'unassigned'})), 'Retained in the all-recovered-unit denominator'),
    ]
    parts.append(table(('Evidence', 'Observed units', 'Denominator / meaning'), rows))
    classes = source.get('class_counts', {})
    if classes:
        parts.append(table(('RNA class assigned by the panel', 'Units', 'Fraction of recovered units'),
                           [(name, number(count), percent(count / denominator if denominator else None)) for name, count in sorted(classes.items(), key=lambda item: -item[1])]))
    candidates = source.get('candidate_taxa_counts', {})
    if candidates:
        parts.append(table(('Candidate reference organism(s)', 'Unit count', 'Fraction of recovered units'),
                           [(name.replace(',', ' / ') or 'Unresolved candidate', number(count), percent(count / denominator if denominator else None)) for name, count in sorted(candidates.items(), key=lambda item: -item[1])]))
    comparison = source.get('reference_comparison') or {}
    if comparison:
        parts.append('<h4>Reference comparison</h4>')
        if comparison.get('status') == 'source_screen_not_available':
            parts.append('<p class="muted">RNA-source evidence is unavailable for these units; host versus non-host candidate fractions cannot be estimated.</p>')
        elif comparison.get('status') == 'assembly_taxon_unknown':
            parts.append('<p class="muted">No supported organism mapping is available for this assembly; host versus non-host reference candidates are not evaluated.</p>')
        else:
            comparison_rows = []
            for key, label in [('host_only_candidate_units', 'Host reference candidates only'), ('host_and_other_candidate_units', 'Host and other reference candidates'), ('non_host_only_candidate_units', 'Non-host reference candidates only'), ('unresolved_candidate_units', 'Unresolved reference candidates')]:
                count = comparison.get(key)
                comparison_rows.append((label, number(count), percent(count / denominator if count is not None and denominator else None)))
            parts.append('<p class="muted">Selected assembly: ' + esc(comparison.get('assembly') or data.get('assembly') or 'Unspecified') + ' · expected reference organism: ' + esc(comparison.get('host_reference_taxon') or 'Unspecified') + '.</p>')
            parts.append(table(('RNA-panel candidate relationship', 'Unit count', 'Fraction of recovered units'), comparison_rows))
        parts.append('<p class="muted">' + esc(comparison.get('interpretation', 'Candidate relationships describe the limited reference panel; they do not establish a sample organism, contamination, or an incorrect reference.')) + '</p>')
    parts.append('<p class="muted">Use these two evidence sources when reviewing the selected reference. Low genome mapping also depends on sequence quality and alignment filters. Conserved rRNA may support several reference organisms; candidate labels do not confirm a species or diagnose contamination.</p>')
    parts.append('<p class="muted">Panel assignments and genome-annotation evidence overlap. Their percentages must not be added. Each row listing several candidate organisms represents one unresolved candidate set, not separate assignments to every organism.</p>')
    thresholds = panel.get('thresholds', {})
    if thresholds:
        parts.append('<p class="muted">Panel matches require ≥' + percent(thresholds.get('minimum_query_coverage')) + ' query coverage, ≥' + percent(thresholds.get('minimum_identity')) + ' identity and ≥' + number(thresholds.get('minimum_aligned_matching_bases')) + ' matching bases. Competing candidates retain ≥' + percent(thresholds.get('near_best_score_fraction')) + ' of the best alignment score.</p>')
    files = data.get('files', [])
    downloads = [('source_assignments.tsv.gz', 'Per-unit RNA-source evidence'), ('source_summary.json', 'Panel settings and reference provenance')]
    links = ['<a href="' + esc(prefix + '/' + name) + '">' + label + '</a>' for name, label in downloads if name in files]
    if links:
        parts.append('<p class="rna-evidence-download">' + ' · '.join(links) + '</p>')
    parts.append('</section>')
    return ''.join(parts)


# Source evidence has precedence over alignment to a selected genome.
SOURCE_COLORS = {'Target genome': '#3276b5', 'E. coli candidate': '#df8b36',
                 'Other genome candidate': '#8b68aa', 'Uncertain / unassigned': '#a4afba'}
FAMILIES = ('rRNA', 'mt-rRNA', 'mt-mRNA', 'tRNA', 'protein_coding', 'snoRNA', 'snRNA', 'lncRNA', 'other', 'unknown', 'ambiguous_family')
REGION_LABELS = {'CDS': 'CDS', 'five_prime_UTR': '5′ UTR', 'three_prime_UTR': '3′ UTR',
                 'other_exon': 'Other exon', 'intron': 'Intron', 'genic_unannotated': 'Other genic', 'intergenic': 'Intergenic'}


def source_category(row, host):
    if row.get('source_status') == 'assigned':
        taxa = {x.strip() for x in row.get('source_taxon', '').split(',') if x.strip()}
        if host and taxa == {host}:
            return 'Target genome'
        if host and taxa and host not in taxa:
            return 'E. coli candidate' if taxa == {'Escherichia coli'} else 'Other genome candidate'
        return 'Uncertain / unassigned'
    if row.get('mapping_status') == 'unique':
        return 'Target genome'
    return 'Uncertain / unassigned'


def source_breakdown(rows, host):
    families, regions = defaultdict(Counter), defaultdict(Counter)
    for row in rows:
        origin = source_category(row, host)
        families[row['family']][origin] += 1
        if row['mapping_status'] == 'unique':
            for region, count in json.loads(row.get('region_bases') or '{}').items():
                regions[region][origin] += count
    return families, regions


def evidence_bars(counts, denominator, breakdown=None, labels=None):
    """Zero categories remain visible; unavailable denominators never become 0%."""
    parts = ['<div class="rna-bars">']
    for key, value in counts.items():
        value = value or 0
        pieces = (breakdown or {}).get(key) or {'Target genome': value}
        segments = ''.join('<span style="width:' + str(100 * n / denominator if denominator else 0) + '%;background:' + SOURCE_COLORS[origin] + '" title="' + esc(origin + ': ' + number(n)) + '"></span>' for origin, n in pieces.items())
        parts.append('<div class="rna-bar-row"><span>' + esc((labels or {}).get(key, key)) + '</span><div class="rna-bar-track">' + segments + '</div><span class="rna-bar-value">' + percent(value / denominator if denominator else None) + ' <small>(' + number(value) + ')</small></span></div>')
    return ''.join(parts) + '</div>'


def core_figures(metrics, rows, prefix, label):
    host = ((metrics.get('source_evidence') or {}).get('reference_comparison') or {}).get('host_reference_taxon')
    family, regions = source_breakdown(rows, host)
    # Legacy reports without unit records retain counts but cannot infer source colors.
    if not rows:
        family = {k: {'Uncertain / unassigned': v} for k, v in metrics['family_counts'].items()}
        regions = {k: {'Uncertain / unassigned': v} for k, v in metrics['region_bases'].items()}
    legend = '<div class="rna-source-legend">' + ''.join('<span><i style="background:' + color + '"></i>' + esc(name) + '</span>' for name, color in SOURCE_COLORS.items()) + '</div>'
    parts = ['<div class="rna-core"><h4>RNA family</h4>', legend,
        evidence_bars({k: metrics['family_counts'].get(k, 0) for k in FAMILIES}, metrics['units'], family),
        '<p class="muted">All recovered non-dimer units. Colors indicate reference evidence, not confirmed species identity. Zero means no observed assignment.</p>',
        '<h4>RNA regions</h4>', evidence_bars(metrics['region_bases'], metrics['region_denominator_bases'], regions, REGION_LABELS),
        '<p class="muted">' + number(metrics['region_denominator_bases']) + ' matched bases on the selected genome; skipped introns and deletions excluded. Exonic bases are split into CDS, UTRs and other exon.</p>',
        '<h4>Strand orientation</h4>', evidence_bars(metrics['strand_counts'], metrics['units'], labels={'sense':'Sense', 'antisense':'Antisense', 'unknown':'Unknown', 'annotation_conflict':'Annotation conflict'}),
        '<p class="muted">Physical P5→P7 direction relative to annotated RNA. Sense among ' + number(metrics['strand_informative_units']) + ' informative units: <strong>' + percent(metrics['physical_sense_fraction']) + '</strong>.</p>',
        '<h4>5′ / 3′ coverage</h4>']
    bias = metrics['profiles']['bias']
    fraction = bias.get('median_five_prime_fraction')
    if fraction is None:
        parts.append('<p class="rna-empty">Insufficient end coverage: ' + number(bias.get('eligible_terminal_supported_genes', bias.get('eligible_genes_with_both_ends', 0))) + ' supported genes; ' + number(bias['minimum_genes']) + ' required.</p>')
    else:
        parts.append('<div class="rna-end-bar"><span style="width:' + str(100*fraction) + '%"></span></div><p>5′ <strong>' + percent(fraction) + '</strong> · 3′ <strong>' + percent(1-fraction) + '</strong></p>')
    parts.append('<p class="muted">Median per-gene share of coverage in the first / last 20% of transcript exons; does not identify library chemistry.</p>')
    target = esc(prefix + '/' + metrics['plots']['profiles'])
    parts.append('<h4>Transcript metaplot</h4><figure class="plot-card"><a href="' + target + '"><img loading="lazy" src="' + target + '" alt="' + esc(label) + ': transcript metaplot"></a></figure>')
    parts.append('<h4>Splice evidence</h4><p class="muted">' + number(metrics['splice_denominator_units']) + ' unique MAPQ-pass units</p>')
    parts.append(evidence_bars({'Spliced':metrics['spliced_units'], 'No splice observed':metrics['splice_denominator_units']-metrics['spliced_units']}, metrics['splice_denominator_units']))
    parts.append('<p class="muted">Junction observations</p>' + evidence_bars({'Annotated':metrics['annotated_junction_observations'], 'Not in annotation':metrics['junction_observations']-metrics['annotated_junction_observations']}, metrics['junction_observations']))
    return ''.join(parts) + '</div>'


def render_motifs(motif):
    parts = ['<h4>CLIP motif calling</h4>']
    if motif.get('engine') != 'HOMER':
        parts.append('<div class="rna-motif-notice" role="note">HOMER analysis has not been run for this report.</div>')
    elif motif['status'] == 'insufficient_units':
        parts.append('<div class="rna-motif-notice" role="note">Not enough distinct sequences for HOMER (' + number(motif.get('sampled_units', 0)) + ' / ' + number(motif['minimum_units']) + ' required).</div>')
    elif motif['status'] != 'complete':
        label = {'engine_unavailable':'HOMER is unavailable.', 'timed_out':'HOMER reached its time limit.', 'engine_failed':'HOMER analysis did not complete.'}.get(motif['status'], 'HOMER results are unavailable.')
        parts.append('<div class="rna-motif-notice" role="note">' + label + '</div>')
    else:
        from rna_motifs import logo_svg
        parts.append('<p class="muted">HOMER · ' + number(motif['sampled_units']) + ' distinct RNA sequences</p>')
        if motif['motifs']:
            parts.append('<div class="table-scroll rna-table"><table class="rna-motif-table"><thead><tr><th>Motif</th><th>Sequences with motif</th><th>Shuffled background</th><th>Enrichment</th></tr></thead><tbody>')
            for row in motif['motifs'][:3]:
                adapter_note = '<br><small style="color:#a65f18" title="Consensus also occurs in a standard library adapter">Adapter sequence match</small>' if row.get('adapter_sequence_match') else ''
                parts.append('<tr><td><div style="width:190px;max-width:100%" role="img" aria-label="' + esc(row['motif']) + '">' + logo_svg(row['pwm']) + '</div><code>' + esc(row['motif']) + '</code>' + adapter_note + '</td><td>' + percent(row['target_fraction']) + ' (' + number(row['observed_units']) + ')</td><td>' + percent(row['background_fraction']) + '</td><td>' + number(row['enrichment'], 2) + '×</td></tr>')
            parts.append('</tbody></table></div>')
        else:
            parts.append('<p class="muted">No recurrent enriched motif detected.</p>')
    return ''.join(parts)


def render_rna_section(summary_path, relative_prefix=None, scope_labels=None):
    path = Path(summary_path)
    if not path.exists():
        return '<section class="report-section" id="rna-qc"><h2>RNA QC · Beta</h2><p>RNA QC results are unavailable.</p></section>'
    data = json.loads(path.read_text())
    prefix = (relative_prefix or '../qc/rna/' + quote(data['sample'], safe='')).rstrip('/')
    title = 'CLIP-seq QC' if data['assay'] == 'clip' else 'RNA-seq QC'
    scope_labels = scope_labels or {}
    parts = ['<section class="report-section rna-report" id="rna-qc"><h2>' + title + ' <span class="rna-beta">Beta</span></h2>',
        '<p class="muted">Beta version — still being tested.</p>']
    unit_path = path.parent / 'units.tsv.gz'
    with gzip.open(unit_path, 'rt') if unit_path.exists() else __import__('io').StringIO('') as handle:
        unit_rows = list(csv.DictReader(handle, delimiter='\t'))
    annotation = data['annotation']
    missing = annotation.get('missing_families', [])
    parts.append('<p class="rna-metadata"><span>Assembly <strong>' + esc(data.get('assembly') or 'Unspecified') + '</strong></span><span>' + number(annotation.get('genes', 0)) + ' annotated genes</span><span>' + number(data['scopes'].get('total', {}).get('units', 0)) + ' recovered units</span></p>')
    annotation_notes = []
    if annotation.get('genes_without_explicit_biotype'):
        evidence = annotation.get('classification_evidence_counts', {})
        annotation_notes.append('<p class="note">' + number(annotation['genes_without_explicit_biotype']) + ' genes lack an explicit RNA biotype. Where supported, CDS annotations (' + number(evidence.get('CDS_annotation', 0)) + ' genes) or RNA gene-symbol heuristics (' + number(evidence.get('gene_symbol_heuristic', 0)) + ' genes) supply descriptive family evidence; others remain unknown.</p>')
    if missing:
        annotation_notes.append('<p class="note"><strong>Incomplete RNA-family annotation:</strong> ' + esc(', '.join(missing)) + '. Zero observed assignments do not establish absence of these RNAs. Multi-family mappings remain visible.</p>')
    def scope_label(scope):
        return scope_labels.get(scope) or ('Whole sample' if scope == 'total' else 'Sublibrary ' + scope)

    if len(data['scopes']) > 1:
        parts.append('<nav class="rna-scope-nav" aria-label="RNA QC scopes">' + ''.join('<a href="#rna-scope-' + esc(quote(scope, safe='')) + '">' + esc(scope_label(scope)) + '</a>' for scope in data['scopes']) + '</nav>')
    for scope, metrics in data['scopes'].items():
        label = scope_label(scope)
        parts.append('<section class="profile-scope rna-scope' + (' rna-primary' if scope == 'total' else '') + '" id="rna-scope-' + esc(quote(scope, safe='')) + '" aria-label="' + esc(label) + '"><header class="rna-scope-heading"><h3>' + esc(label) + '</h3><span>' + number(metrics['units']) + ' recovered units</span></header>')
        strand = metrics['physical_sense_wilson95']
        confidence = 'Not evaluable' if strand[0] is None else f'{100 * strand[0]:.1f}–{100 * strand[1]:.1f}%'
        bias = metrics['profiles']['bias']
        bases = metrics['region_denominator_bases']
        regions = metrics['region_bases']
        exon = sum(regions.get(key, 0) for key in ('CDS', 'five_prime_UTR', 'three_prime_UTR', 'other_exon'))
        tails = metrics.get('terminal_poly_a_evidence', {})
        tail_n, tail_total = tails.get('terminal_a_run_ge8_units', 0), tails.get('oriented_gene_assigned_units', 0)
        ratio_display = '5′ only (no finite ratio)' if bias.get('five_to_three_ratio_status') == 'five_prime_only' else number(bias['five_to_three_ratio'], 2)
        terminal_supported = bias.get('eligible_terminal_supported_genes', bias.get('eligible_genes_with_both_ends', 0))
        umi = metrics['umi']
        preparation = data.get('preparation') or {}
        automatic_grip = preparation.get('mode') == 'auto_grip' and not (preparation.get('auto_grip') or {}).get('bypass_reason')
        umi_status = umi.get('status_counts', {})
        umi_total = sum(umi_status.values())
        umi_configured = any(count for status, count in umi_status.items() if status != 'not_configured')
        umi_label = 'Complete UMI recovery'
        umi_support = 'Recovered units with a complete configured UMI' if umi_configured else ('No UMI structure configured for these units' if umi_total else 'No recovered units with UMI evidence')
        umi_recovery = (number(umi_status.get('complete', 0)) + ' / ' + number(umi_total) if umi_configured else 'Not configured') if umi_total else 'Not evaluable'
        if automatic_grip:
            valid_p7 = umi.get('auto_grip_p7_complete_units')
            umi_label = 'Valid P7 N10 UMI recovery'
            umi_recovery = number(valid_p7) + ' / ' + number(metrics['units'])
            umi_support = percent(valid_p7 / metrics['units'] if valid_p7 is not None and metrics['units'] else None) + '; all recovered units; structure-supported, observed A/C/G/T tags only'
            umi_configured = valid_p7 is not None and valid_p7 > 0
        ribosomal = sum(metrics.get('family_counts', {}).get(family, 0) for family in ('rRNA', 'mt-rRNA'))
        rows = [
            ('Mapped units', number(metrics['mapped_units']) + ' / ' + number(metrics['units']), percent(metrics['mapped_fraction'])),
            ('Usable gene assignments', number(metrics['gene_assigned_units']), number(metrics['detected_genes']) + ' detected genes'),
            ('Ribosomal-family evidence', number(ribosomal) + ' / ' + number(metrics['units']), percent(ribosomal / metrics['units'] if metrics['units'] else None) + '; annotation or RNA-panel evidence, each unit counted once'),
            ('Physical P5→P7 sense fraction', percent(metrics['physical_sense_fraction']), number(metrics['strand_informative_units']) + ' informative units; 95% CI ' + confidence),
            ('Spliced units', number(metrics['spliced_units']) + ' / ' + number(metrics['splice_denominator_units']), 'Primary unique MAPQ-pass units'),
            ('Annotated junctions', number(metrics['annotated_junction_observations']) + ' / ' + number(metrics['junction_observations']), percent(metrics['annotated_junction_fraction'])),
            ('Exonic / intronic / intergenic matched bases', ' / '.join(percent(value / bases if bases else None) for value in (exon, regions['intron'], regions['intergenic'])), number(bases) + ' matched bases; untyped genic bases reported separately'),
            ('Terminal A-run ≥8 bp', number(tail_n) + ' / ' + number(tail_total), percent(tail_n / tail_total if tail_total else None) + '; observed homopolymer, not poly(A) selection or tail-length identification'),
            ('Median physical insert', number(metrics['insert_median_bp']) + ' bp', percent(metrics['short_insert_fraction_lt100']) + ' shorter than 100 bp'),
            ('Top 10 genes', percent(metrics['top10_gene_fraction']), 'Of assigned units; expression/target concentration is not PCR duplication'),
            (umi_label, umi_recovery, umi_support),
            ('Observed coordinate + UMI groups', number(umi['observed_exact_coordinate_umi_groups']) if umi_configured else 'Not evaluable', number(umi['eligible_coordinate_umi_units']) + ' unique MAPQ-pass alignments with physical orientation and complete UMI; no deduplication' if umi_configured else 'Requires a complete recovered UMI and usable genome coordinates'),
            ('5′ / 3′ end coverage ratio', ratio_display, number(terminal_supported) + ' terminal-supported genes; minimum ' + number(bias['minimum_genes'])),
        ]
        comparison = (metrics.get('source_evidence') or {}).get('reference_comparison') or {}
        if comparison.get('status') in {'available', 'no_assigned_source_evidence'}:
            rows.insert(3, ('Non-host-only RNA-panel candidates', number(comparison.get('non_host_only_candidate_units')) + ' / ' + number(comparison.get('denominator_units')), percent(comparison.get('non_host_only_fraction_of_all_units')) + '; candidate reference evidence, not a contamination or mismatch diagnosis'))
        selected_rows = unit_rows if scope == 'total' else [row for row in unit_rows if row['scope_id'] == scope]
        parts.append(core_figures(metrics, selected_rows, prefix, label))
        motif = metrics['motifs']
        if data['assay'] == 'clip':
            parts.append(render_motifs(motif))
        parts.append('<details class="rna-auxiliary"><summary>Counts, UMI and methods</summary>')
        if data['assay'] == 'clip':
            parts.append('<h4>Motif method</h4><p>' + esc(motif['interpretation']) + '. ' + esc(motif.get('background', '')) + '.</p>')
            if motif.get('engine') == 'HOMER' and motif.get('motifs'):
                parts.append(table(('Motif', 'HOMER log10 P (discovery score)'), [(row['motif'], number(row['homer_log10_p'], 2)) for row in motif['motifs']]))
        parts.append(table(('Evidence', 'Observed', 'Denominator / interpretation'), rows))
        parts.append(render_structure_evidence(preparation, scope, scope_labels, prefix))
        parts.append(render_source_evidence(data, metrics, prefix))
        if scope == 'total':
            parts.extend(annotation_notes)
        parts.append('</details>')
        if metrics['ambiguous_family_sets']:
            parts.append('<details><summary>Ambiguous RNA families</summary>' + table(('Candidate family set', 'Units'), sorted(metrics['ambiguous_family_sets'].items(), key=lambda item: -item[1])) + '</details>')
        parts.append('</section>')
    parts.append('<details class="metric-group"><summary>Methods, limitations and source data</summary><ul>' + ''.join('<li>' + esc(text) + '</li>' for text in data['limitations']) + '</ul>')
    model_sources = []
    for index, url in enumerate(((data.get('preparation') or {}).get('auto_grip') or {}).get('sources', []), 1):
        if isinstance(url, str) and urlsplit(url).scheme in {'https', 'http'} and urlsplit(url).netloc:
            model_sources.append('<a href="' + esc(url) + '" target="_blank" rel="noopener">Protocol source ' + str(index) + '</a>')
    if model_sources:
        parts.append('<p>Structure-model references: ' + ' · '.join(model_sources) + '</p>')
    parts.append('<p>Annotation SHA-256: <code>' + esc(annotation['sha256']) + '</code>. Assembly: ' + esc(data.get('assembly') or 'unspecified') + '.</p><p>')
    parts.append(' · '.join('<a href="' + esc(prefix + '/' + file) + '">' + esc(file) + '</a>' for file in data['files'] if not file.endswith('.svg')))
    parts.append('</p></details></section>')
    return ''.join(parts)


def render_cohort_qc(summary_paths, out_dir, relative_prefix='../qc/rna/cohort'):
    from collections import Counter
    from rna_qc import cohort_analysis
    paths = list(summary_paths)
    data = cohort_analysis(paths, out_dir)
    parts = ['<section class="report-section rna-report" id="rna-cohort"><h2>RNA sample overview</h2>',
             '<p class="muted">Whole-sample RNA evidence. Open a sample for its expanded RNA results and independently assessed index scopes.</p>',
             '<div class="table-scroll rna-table" tabindex="0"><table><thead><tr>' + ''.join('<th>' + esc(label) + '</th>' for label in ('Sample', 'Recovered units', 'Ribosomal-family evidence', 'Genome mapped', 'Valid P7 N10 UMI recovery')) + '</tr></thead><tbody>']
    for path in paths:
        summary = json.loads(Path(path).read_text())
        metrics = summary['scopes'].get('total', {})
        units = metrics.get('units', 0)
        ribosomal = sum(metrics.get('family_counts', {}).get(family, 0) for family in ('rRNA', 'mt-rRNA'))
        preparation = summary.get('preparation') or {}
        automatic_grip = preparation.get('mode') == 'auto_grip' and not (preparation.get('auto_grip') or {}).get('bypass_reason')
        valid_p7 = metrics.get('umi', {}).get('auto_grip_p7_complete_units') if automatic_grip else None
        umi_display = (number(valid_p7) + ' / ' + number(units) + ' (' + percent(valid_p7 / units if units else None) + ')') if valid_p7 is not None else 'Not evaluable'
        target = '../reports/' + quote(summary['sample'], safe='') + '.html#rna-qc'
        parts.append('<tr><td><a href="' + esc(target) + '">' + esc(summary['sample']) + '</a></td>' + ''.join('<td>' + esc(value) + '</td>' for value in (number(units), percent(ribosomal / units if units else None), percent(metrics.get('mapped_fraction')), umi_display)) + '</tr>')
    parts.append('</tbody></table></div><p class="muted">Fractions use all recovered units. Ribosomal-family evidence combines annotation and the independent RNA panel, counting each unit once. P7 recovery reflects observed tags from supported index scopes, not a sample-wide structure assumption or deduplicated molecule count.</p>')
    parts.append('<h3>RNA cohort similarity</h3><p class="note">' + esc(data['interpretation']) + '</p>')
    available = [pair for pair in data['pairs'] if pair['status'] == 'available' and pair.get('log_cpm_pearson') is not None]
    excluded = Counter(pair['status'] if pair['status'] != 'available' else 'not_evaluable' for pair in data['pairs'] if pair not in available)
    labels = {'different_assay': 'Different assays', 'different_or_unknown_assembly': 'Different or unknown assembly',
              'different_annotation': 'Different annotation', 'insufficient_assigned_units': 'Too few gene-assigned units',
              'insufficient_variable_genes': 'Too few shared or variable genes', 'not_evaluable': 'Not evaluable'}
    parts.append('<p>' + number(len(available)) + ' / ' + number(len(data['pairs'])) + ' scope pairs meet the comparison criteria: at least ' + number(data['minimum_assigned_units']) + ' gene-assigned units per scope and ' + number(data['minimum_shared_genes']) + ' shared genes, with matching assay, assembly and annotation.</p>')
    if available:
        rows = []
        for pair in available:
            label_a = pair['sample_a'] + ' / ' + ('whole sample' if pair['scope_a'] == 'total' else pair['scope_a'])
            label_b = pair['sample_b'] + ' / ' + ('whole sample' if pair['scope_b'] == 'total' else pair['scope_b'])
            rows.append((label_a, label_b, number(pair['log_cpm_pearson'], 3), percent(pair['gene_jaccard']),
                         'Related scopes within one sample' if pair.get('related_sample_scopes') else 'Different samples'))
        parts.append(table(('Library A', 'Library B', 'log1p CPM Pearson', 'Detected gene Jaccard', 'Scope relationship'), rows))
    else:
        parts.append('<p class="muted">No scope pair has sufficient comparable gene-profile evidence for a correlation estimate.</p>')
    if excluded:
        parts.append('<p class="muted">Other pairs: ' + esc('; '.join(labels.get(status, status.replace('_', ' ')) + ': ' + number(count) for status, count in sorted(excluded.items()))) + '. Complete pairwise statuses remain available in the downloads.</p>')
    parts.append('<p><a href="' + esc(relative_prefix + '/cohort.json') + '">Cohort JSON</a> · <a href="' + esc(relative_prefix + '/cohort_similarity.tsv') + '">All pairwise results</a></p></section>')
    return ''.join(parts)
