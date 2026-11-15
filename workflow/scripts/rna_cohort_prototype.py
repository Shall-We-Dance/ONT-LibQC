"""Inactive RNA cohort presentation prototype; not imported by report generation."""
from rna_report import esc, number, percent, table
from urllib.parse import quote
from pathlib import Path
import json

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
