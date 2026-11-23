"""Physical sample-index composition in the single-cell report."""
from report_palette import scope_color
from report_sublibrary import mixture_chart, mixture_counts


def sublibrary_composition(summary):
    physical = summary.get('physical_library')
    if not physical:
        return ''
    total, rows = mixture_counts(physical, scope_color)
    if not total:
        return ''
    assigned = len(physical.get('sub_libraries', []))
    return (f'<p>Physical sample-index assignments across all {total:,} reconstructed library units, '
            'including adapter dimers.</p>'
            '<div class="sublibrary-mixture">' + mixture_chart(total, rows[:assigned], folded=rows[assigned:]) + '</div>')
