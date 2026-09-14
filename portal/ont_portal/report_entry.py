"""Remove the web job's redundant all-sample overview, retaining scientific data."""
import html, json, re
from pathlib import Path



def simplify_download_footer(text):
    def footer(match):
        return re.sub(r"<p>(?:(?!</p>).)*Adapter table(?:(?!</p>).)*</p>",
            '<p class="report-download-links"><a href="../downloads.html">Browse report files</a></p>',
            match.group(0), flags=re.S)
    return re.sub(r"<footer\b[^>]*class=['\"]report-footer['\"][^>]*>.*?</footer>",footer,text,flags=re.S)


def online_report_links(text, file):
    """Add the server's archive link only online; extracted reports stay portable."""
    prefix='../'*(len(Path(file).parts)-1)
    archive=html.escape(prefix+'report.zip',quote=True)
    browse=html.escape(prefix+'downloads.html',quote=True)
    button=('<a href="'+archive+'" download style="display:inline-block;padding:10px 16px;'
            'border-radius:7px;background:#176e87;color:#fff;text-decoration:none;font-weight:600">'
            'Download complete report (ZIP)</a>')
    block=('<div class="report-download-actions" style="margin:20px 0">'
           '<p style="display:flex;flex-wrap:wrap;align-items:center;gap:12px">'+button)
    if file!='downloads.html':block+='<a href="'+browse+'">Browse report files</a>'
    block+='</p><p class="muted" style="margin-top:8px;font-size:13px">Saved HTML report, figures, tables and settings in one ZIP.</p></div>'
    if file=='downloads.html':
        return re.sub(r'(<main\b[^>]*>)',lambda match:match.group(1)+block,text,count=1)
    text=simplify_download_footer(text)
    return re.sub(r'<p class="report-download-links">.*?</p>',lambda _:block,text,count=1,flags=re.S)


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
    page.write_text(simplify_download_footer(text))
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
