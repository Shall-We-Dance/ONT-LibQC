"""Candidate fractions within each resolved sublibrary."""
import html

SHORT_NAMES = {
    'Escherichia coli': 'E. coli',
    'Saccharomyces cerevisiae': 'S. cerevisiae',
    'Schizosaccharomyces pombe': 'S. pombe',
    'Drosophila melanogaster': 'D. melanogaster',
}



MATRIX_CSS = """
#spike-screen .spike-matrix-frame{border:1px solid #dce5ef;border-radius:8px;overflow:auto;background:#fff;scrollbar-width:thin;scrollbar-color:#bac7d6 #f6f8fb}
#spike-screen .spike-matrix{margin:0;border:0;border-collapse:separate;border-spacing:0;border-radius:0;box-shadow:none}
#spike-screen .spike-matrix th,#spike-screen .spike-matrix td{border:0;border-bottom:1px solid #e8edf3;padding:12px 14px}
#spike-screen .spike-matrix thead th{background:#f5f7fa;border-bottom:1px solid #dce5ef}
#spike-screen .spike-matrix tbody th{background:transparent}
#spike-screen .spike-matrix tr:last-child>*{border-bottom:0}
#spike-screen .spike-matrix th:first-child{border-right:1px solid #e8edf3}
#spike-screen .spike-matrix tbody tr:hover>*{background:#f8fafc}
#spike-screen .spike-matrix-frame:focus-visible{outline:2px solid #5176a5;outline-offset:3px}

@media screen {
  #spike-screen .spike-matrix{width:100%;min-width:1000px;table-layout:fixed}
  #spike-screen .spike-matrix thead th:first-child{width:112px!important}
  #spike-screen .spike-matrix thead th:nth-child(2){width:132px}
}
@media print {
  #spike-screen .spike-matrix-frame{width:100%;overflow:visible;border:0;border-radius:0}
  #spike-screen .spike-matrix{width:100%;min-width:0;table-layout:fixed;border-collapse:collapse}
  #spike-screen .spike-matrix th,#spike-screen .spike-matrix td{padding:2mm 1mm;font-size:8pt;line-height:1.35;vertical-align:middle;overflow-wrap:normal;word-break:normal}
  #spike-screen .spike-matrix th:first-child{width:12%!important;min-width:0!important;max-width:none!important}
  #spike-screen .spike-matrix thead th:first-child{white-space:nowrap}
  #spike-screen .spike-matrix thead th:nth-child(2){width:14%}
  #spike-screen .spike-matrix tbody th{white-space:normal!important;overflow-wrap:anywhere}
}
"""


def spike_matrix_html(screen, sublibraries, sample, note=''):
    candidates = screen.get('candidates', {})
    headers = ''.join('<th scope="col" title="'+html.escape(c.get('label', key), quote=True)+'">'+html.escape(SHORT_NAMES.get(c.get('label', key), c.get('label', key)))+'</th>' for key,c in candidates.items())
    rows = []
    for lib in sublibraries:
        n = lib.get('spike_in_screened_reads')
        counts = lib.get('spike_in_species_counts')
        cells = []
        for key,candidate in candidates.items():
            count = (counts or {}).get(key, 0)
            if candidate.get('status') == 'HOST_NOT_DISTINGUISHABLE':
                value, color, tip = 'N/A', 'gray', 'Same species as the host; spike-in cannot be distinguished.'
            elif not n or counts is None:
                value, color, tip = 'N/A', 'gray', 'No screened inserts or per-candidate evidence available.'
            else:
                percent = count / n * 100
                value = ('<0.01%' if percent < .01 else f'{percent:.2f}%') if count else '0%'
                color = 'amber' if count else 'gray'
                tip = f'{count:,} / {n:,} screened inserts in this sublibrary.'
            cells.append('<td style="white-space:nowrap;font-variant-numeric:tabular-nums" title="'+html.escape(tip,quote=True)+'"><span class="light '+color+'" aria-hidden="true" style="display:inline-block;margin-right:7px"></span>'+html.escape(value)+'</td>')
        label = html.escape(str(lib.get('scope_name',lib.get('scope_id',''))))
        rows.append('<tr><th scope="row" style="white-space:nowrap">'+label+'</th><td>'+ (f'{n:,}' if n is not None else 'N/A')+'</td>'+''.join(cells)+'</tr>')
    return ("<section class='report-section' id='spike-screen'><style>"+MATRIX_CSS+"</style><h2>Spike-in candidate screen</h2>"
            "<p class='muted'>Percentages are within each sublibrary’s screened inserts. Yellow: signal detected; gray: no signal. N/A: not assessable.</p>"
            '<div class="table-scroll spike-matrix-frame" tabindex="0" aria-label="Spike-in candidates by sublibrary"><table class="spike-matrix"><thead><tr><th scope="col">Sublibrary</th><th scope="col">Screened inserts</th>'+headers+'</tr></thead><tbody>'+''.join(rows)+'</tbody></table></div>'
            +note
            +"</section>")
