"""Self-contained PDF/print controls for generated and existing sample reports."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent


def add_pdf_export(document):
    """Add once to sample reports; assets are inline so downloaded reports work offline."""
    if 'id="report-export-script"' in document or "id='report-export-script'" in document:
        return document
    # Summary/download pages have no sample hero and should keep their existing layout.
    hero = re.search(r'<header\b[^>]*\bclass=[\'"][^\'"]*\breport-hero\b[^\'"]*[\'"][^>]*>', document, re.I)
    if not hero:
        return document
    controls = ('<div class="report-export-controls">'
                '<button type="button" id="report-export-pdf" aria-describedby="report-export-status">'
                '<svg width="16" height="16" viewBox="0 0 20 20" fill="none" aria-hidden="true">'
                '<path d="M10 2v10m-4-4 4 4 4-4M3 13v4h14v-4" stroke="currentColor" '
                'stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/></svg>'
                '<span>Export PDF</span></button></div>')
    document = document[:hero.end()] + controls + document[hero.end():]
    feedback = '<p id="report-export-status" class="report-export-status" role="status" aria-live="polite"></p>'
    document = re.sub(r'</header\s*>', lambda match: match[0] + feedback, document, count=1, flags=re.I)
    screen_css = '''
.report-export-controls{float:right;margin:0 0 12px 18px}
.report-export-controls button{display:inline-flex;align-items:center;gap:7px;background:#fff;color:var(--blue,#2a6cbd);border:1px solid var(--line,#dce5f0);border-radius:6px;padding:8px 12px;font-family:inherit;font-size:13px;font-weight:500;line-height:1.4;white-space:nowrap}
.report-export-controls button:hover{background:var(--blue-soft,#edf4fc);border-color:var(--blue,#2a6cbd)}
.report-export-controls button:focus-visible{outline:2px solid var(--blue,#2a6cbd);outline-offset:3px}
.report-export-controls button:disabled{cursor:wait;opacity:.65}
.report-export-status{font-size:12px;color:var(--muted,#5b6d82);margin:8px 0 0}
.report-export-status:empty{display:none}
.report-print-definitions,.report-print-note{display:none}
.report-hero::after{content:"";display:table;clear:both}
@media(max-width:520px){.report-export-controls{float:none;display:flex;justify-content:flex-end;margin:0 0 12px}}
'''
    styles = screen_css + '\n' + (ROOT / 'report_print.css').read_text()
    document = document.replace('</head>', '<style id="report-export-style">' + styles + '</style></head>', 1)
    script = (ROOT / 'report_export.js').read_text()
    return document.replace('</body>', '<script id="report-export-script">' + script + '</script></body>', 1)
