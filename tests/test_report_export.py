"""PDF controls stay self-contained without disturbing the scientific report."""
from html.parser import HTMLParser
from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from report_export import add_pdf_export


class ReportMarkup(HTMLParser):
    """Collect elements and their context without adding a parser dependency."""
    VOID_TAGS = {'area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input',
                 'link', 'meta', 'param', 'source', 'track', 'wbr'}

    def __init__(self, document):
        super().__init__()
        self.elements = []
        self.stack = []
        self.feed(document)
        self.close()

    def handle_starttag(self, tag, attrs):
        node = {'tag': tag, 'attrs': dict(attrs), 'parents': tuple(self.stack), 'text': ''}
        self.elements.append(node)
        if tag not in self.VOID_TAGS:
            self.stack.append(node)

    def handle_endtag(self, tag):
        for index in range(len(self.stack) - 1, -1, -1):
            if self.stack[index]['tag'] == tag:
                del self.stack[index:]
                break

    def handle_data(self, text):
        for node in self.stack:
            node['text'] += text

    def identified(self, identifier):
        return [node for node in self.elements if node['attrs'].get('id') == identifier]

    def resource_references(self):
        return [(node['tag'], key, value) for node in self.elements
                for key, value in node['attrs'].items() if key in {'src', 'href'}]


class ReportExportTests(unittest.TestCase):
    SAMPLE = '<h1>Sample A &amp; B</h1><p class="hero-meta">Chromatin · dm6</p>'
    METRICS = '''<section class="report-section" id="mapping">
<h2>Mapping</h2><table class="qc-metrics"><tbody><tr>
<th><span data-tip="Mapped inserts divided by recovered inserts">Mapped inserts</span></th>
<td>80 / 100</td></tr></tbody></table>
<details><summary>Chromosome coverage</summary>
<img loading="lazy" src="../assets/chromosome.svg" alt="Chromosome coverage">
<a href="../data/mapping.tsv?view=example-grant" download>Download table</a></details>
</section>'''
    METRIC_SCRIPTS = '''<script id="metric-data" type="application/json">{"sample":"A & B","mapped":80}</script>
<script id="metric-ui">document.querySelector('#mapping').dataset.ready = 'true';</script>'''

    def document(self, hero='<header class="report-hero">'):
        return ('<!doctype html><html><head><meta charset="utf-8"><title>Sample report</title>'
                '<style>.report-section{margin:20px 0}</style></head><body>'
                + hero + self.SAMPLE + '</header>' + self.METRICS
                + '<footer><a href="../downloads.html">Browse report files</a></footer>'
                + self.METRIC_SCRIPTS + '</body></html>')

    def test_sample_report_receives_inline_accessible_pdf_controls(self):
        original = self.document()
        markup = ReportMarkup(add_pdf_export(original))
        button, = markup.identified('report-export-pdf')
        status, = markup.identified(button['attrs']['aria-describedby'])
        script, = markup.identified('report-export-script')
        style, = markup.identified('report-export-style')

        self.assertEqual(button['tag'], 'button')
        self.assertEqual(button['attrs']['type'], 'button')
        self.assertIn('Export PDF', button['text'])
        self.assertTrue(any(parent['tag'] == 'header' and
                            'report-hero' in parent['attrs'].get('class', '').split()
                            for parent in button['parents']))
        self.assertEqual(status['attrs']['role'], 'status')
        self.assertEqual(status['attrs']['aria-live'], 'polite')
        self.assertEqual(script['tag'], 'script')
        self.assertNotIn('src', script['attrs'])
        self.assertIn('window.print()', script['text'])
        self.assertEqual(style['tag'], 'style')
        self.assertIn('@media print', style['text'])
        self.assertIn('@page', style['text'])
        self.assertTrue(any(parent['tag'] == 'head' for parent in style['parents']))
        self.assertNotRegex(style['text'], r'@import\b|url\(\s*[\'"]?https?://')
        self.assertEqual(markup.resource_references(), ReportMarkup(original).resource_references())

    def test_repeated_upgrade_is_byte_for_byte_idempotent(self):
        upgraded = add_pdf_export(self.document())
        for quote in ('"', "'"):
            existing = upgraded.replace('id="report-export-script"',
                                        f'id={quote}report-export-script{quote}')
            with self.subTest(script_id_quote=quote):
                self.assertEqual(add_pdf_export(existing), existing)
                self.assertEqual(len(ReportMarkup(existing).identified('report-export-pdf')), 1)

    def test_sample_content_links_and_existing_metric_scripts_are_preserved(self):
        original = self.document()
        upgraded = add_pdf_export(original)
        for fragment in (self.SAMPLE, self.METRICS, self.METRIC_SCRIPTS):
            self.assertIn(fragment, upgraded)
        before = ReportMarkup(original)
        after = ReportMarkup(upgraded)
        self.assertEqual(after.resource_references(), before.resource_references())
        for identifier in ('metric-data', 'metric-ui'):
            self.assertEqual(after.identified(identifier)[0]['text'],
                             before.identified(identifier)[0]['text'])

    def test_hero_class_accepts_both_quote_styles_and_multiple_classes(self):
        for quote in ('"', "'"):
            hero = f'<header data-sample="example" class={quote}wide report-hero compact{quote}>'
            with self.subTest(class_quote=quote):
                markup = ReportMarkup(add_pdf_export(self.document(hero)))
                self.assertEqual(len(markup.identified('report-export-pdf')), 1)
                self.assertEqual(len(markup.identified('report-export-script')), 1)
                self.assertEqual(len(markup.identified('report-export-style')), 1)

    def test_pages_without_a_sample_hero_are_unchanged(self):
        documents = (
            '<html><head><title>Downloads</title></head><body><a href="data.tsv">Data</a></body></html>',
            self.document('<header class="summary-header">'),
            '<html><head></head><body><div class="report-hero">Summary</div></body></html>',
        )
        for document in documents:
            with self.subTest(document=document[:80]):
                self.assertEqual(add_pdf_export(document), document)


if __name__ == '__main__':
    unittest.main()
