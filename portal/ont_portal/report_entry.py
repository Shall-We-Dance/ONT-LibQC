"""Remove the web job's redundant all-sample overview, retaining scientific data."""
import html, json, re
from pathlib import Path


def direct_report(root, sample):
    root = Path(root)
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]{0,47}', sample):
        raise ValueError('Invalid report sample')
    page = root / 'samples' / (sample + '.html')
    if not page.is_file() or page.is_symlink():
        raise ValueError('Sample report is missing')
    text = page.read_text()
    text = re.sub(r"<a\s+href=['\"]\.\./index\.html['\"]>All samples(?: ↗)?</a>",
                  '<a href="../downloads.html">Report data</a>', text)
    page.write_text(text)
    entry = 'samples/' + sample + '.html'
    (root / 'index.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta http-equiv="refresh" content="0;url=' + html.escape(entry, quote=True) + '">'
        '<title>' + html.escape(sample) + ' report</title></head><body><a href="' + html.escape(entry, quote=True)
        + '">Open ' + html.escape(sample) + ' report</a></body></html>')
    (root / 'assets/all_samples_qc.svg').unlink(missing_ok=True)
    downloads = root / 'downloads.html'
    if downloads.is_file():
        text = downloads.read_text().replace('<a href="index.html">All samples</a>',
            '<a href="' + entry + '">Open report</a>')
        if (root / 'data/submission.json').is_file() and 'href="data/submission.json"' not in text:
            text = text.replace('</ul>', '<li><a href="data/submission.json">Submitted settings</a></li></ul>')
        downloads.write_text(text)
    (root / 'README.md').write_text('# ONT-LibQC report: ' + sample + '\n\nOpen `index.html` or `' + entry
        + '` in a browser. The entry point opens the sample report directly.\n\n'
        'Keep this folder together: `assets/` contains figures; `data/` and `downloads.html` contain report evidence and settings.\n')
    return entry
