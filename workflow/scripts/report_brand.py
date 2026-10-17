"""Shared, offline presentation for ONT-LibQC reports."""
import html
from pathlib import Path
import shutil

MARK = Path(__file__).with_name("report_mark.svg")
MARK_NAME = "ont-libqc-mark.svg"


def copy_mark(assets):
    assets = Path(assets)
    assets.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(MARK, assets / MARK_NAME)


def favicon(mark_url):
    return '<link rel="icon" type="image/svg+xml" href="' + html.escape(mark_url, quote=True) + '">'


def brand_bar(mark_url, context="Analysis report"):
    return ('<div class="report-brand"><div class="brand-name">'
            '<img src="' + html.escape(mark_url, quote=True) + '" width="36" height="36" alt="">'
            '<span>ONT<span class="brand-light">-LibQC</span></span></div>'
            '<span class="brand-context">' + html.escape(context) + '</span></div>')


def downloads_page(items, *, title="Report data", back_label="All samples"):
    """Wrap escaped download-list items; keep navigation local to the bundle."""
    css = Path(__file__).with_name("report_modern.css").read_text()
    mark_url = "assets/" + MARK_NAME
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{html.escape(title)} · ONT-LibQC</title>{favicon(mark_url)}<style>{css}</style>'
            f'</head><body class="report-directory">{brand_bar(mark_url, "Report data")}'
            f'<header class="report-hero">'
            f'<h1>{html.escape(title)}</h1><p>Metrics, source tables and analysis settings.</p></header>'
            f'<nav class="report-nav" aria-label="Report navigation"><a href="index.html">{html.escape(back_label)}</a></nav>'
            '<main class="report-section"><h2>Files</h2>'
            '<p class="muted">Open a file below, or choose “Save link as” to download it.</p>'
            f'<ul class="download-list">{items}</ul></main>'
            '<footer class="report-footer"><strong>ONT-LibQC</strong> · Library quality control'
            '</footer></body></html>')
