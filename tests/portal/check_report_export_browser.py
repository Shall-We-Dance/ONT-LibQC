"""Exercise PDF export on public reports through a temporary local portal.

Usage: PYTHONPATH=portal python tests/portal/check_report_export_browser.py OUTDIR [chromium webkit]
No submissions, worker, mail, external requests, or private reports are used.
ONT_WEBKIT_EXECUTABLE optionally selects a locally configured WebKit executable.
"""
import json
import os
from pathlib import Path
import socket
import shutil
import sys
import tempfile
import threading
import time

import uvicorn
from playwright.sync_api import sync_playwright
from ont_portal.app import create_app
from ont_portal.settings import Settings


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'workflow/scripts'))
from report_refresh import refresh
from report_spike_matrix import spike_matrix_html
from report_sublibrary import render_sublibrary_composition
from report_metrics import REPORT_CSS, scoped_table_html
from species_composition import classify, render as render_species
from report_export import add_pdf_export


OUTPUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('test-results/report-export')
OUTPUT.mkdir(parents=True, exist_ok=True)
GRANT = 'synthetic-print-check'


def wait_js(page, expression, timeout=30000):
    deadline = time.monotonic() + timeout / 1000
    while time.monotonic() < deadline:
        if page.evaluate(expression):
            return
        page.wait_for_timeout(50)
    raise AssertionError('Timed out: ' + expression)


def snapshot(page):
    return page.evaluate("""() => ({
        title: document.title,
        details: Array.from(document.querySelectorAll('details'), d => [d.open, d.getAttribute('name')]),
        links: Array.from(document.querySelectorAll('a'), a => a.getAttribute('href')),
        images: Array.from(document.images, i => i.getAttribute('loading')),
        tables: Array.from(document.querySelectorAll('table'), t => t.innerHTML)
    })""")


def assert_prepared(page):
    assert page.evaluate("Array.from(document.querySelectorAll('details')).every(d => d.open)")
    assert page.locator('.report-print-definitions').count() == 1
    assert page.locator('.report-print-definitions dt').count() > 0
    assert page.evaluate("""Array.from(document.querySelectorAll('a[href]')).every(a =>
        !a.href.includes('view=') && !a.href.includes('synthetic-print-check'))"""), 'Private grant leaked into a printable link'
    assert page.locator('a[href^="#"]').count() == 0, 'Fragments inherit the current page query'
    assert page.evaluate("""Array.from(document.querySelectorAll('table')).every(t =>
        !t.rows.length || !Array.from(t.rows[0].cells).every(c => c.tagName === 'TH') || !!t.tHead)""")


def restore_and_check(page, before):
    page.evaluate('window.dispatchEvent(new Event("afterprint"))')
    page.emulate_media(media='screen')
    assert snapshot(page) == before, 'Print preparation changed the report after printing/cancellation'
    assert page.locator('#report-export-pdf').is_enabled()
    assert page.locator('.report-print-definitions').count() == 0
    assert not page.locator('#report-export-pdf').get_attribute('aria-busy')
    assert page.locator('#report-export-status').inner_text() == ''


def draw_custom_bed(page):
    bed = page.evaluate(r"""() => {
        const p = JSON.parse(document.getElementById('custom-payload').textContent);
        const row = p.scopes[document.getElementById('custom-scope').value].records[0];
        return row[0] + '\t' + row[1] + '\t' + row[2] + '\n';
    }""")
    page.locator('#custom-bed').set_input_files({
        'name': 'export-example.bed', 'mimeType': 'text/plain', 'buffer': bed.encode(),
    })
    wait_js(page, '!document.getElementById("custom-run").disabled')
    page.locator('#custom-assembly').check()
    page.locator('#custom-run').click()
    page.locator('#custom-output svg').wait_for(state='visible')
    canonical_svg = page.locator('#custom-output svg').evaluate('(node) => node.outerHTML')
    chrom, rest = bed.split('\t', 1)
    alternate = chrom[3:] if chrom.startswith('chr') else 'chr' + chrom
    page.locator('#custom-bed').set_input_files({
        'name': 'export-example.bed', 'mimeType': 'text/plain',
        'buffer': (alternate + '\t' + rest).encode(),
    })
    wait_js(page, '!document.getElementById("custom-run").disabled')
    page.locator('#custom-assembly').check()
    page.locator('#custom-run').click()
    page.locator('#custom-output svg').wait_for(state='visible')
    assert page.locator('#custom-output svg').evaluate('(node) => node.outerHTML') == canonical_svg, 'Chromosome prefix changed the custom profile'
    return page.locator('#custom-output').inner_html()


def check_wide_table(page):
    # Six sublibraries exceed the paper width: every original value must occur
    # exactly once across the printed column groups, with labels repeated.
    page.evaluate("""() => {
        const table = document.createElement('table');
        table.id = 'wide-print-fixture'; table.className = 'qc-metrics';
        for (let r = 0; r < 3; r++) {
            const row = table.insertRow();
            for (let c = 0; c < 9; c++) {
                const cell = document.createElement(r ? 'td' : 'th');
                cell.textContent = r + ':' + c; row.append(cell);
            }
        }
        document.body.append(table);
    }""")
    before = snapshot(page)
    page.evaluate('window.ONTReportExport.prepare()')
    assert page.evaluate("""() => {
        const source = document.getElementById('wide-print-fixture');
        const copies = [...source.nextElementSibling.querySelectorAll('table')];
        return copies.length === 2 && [...source.rows].every((row, index) => {
            const parts = copies.map(copy => [...copy.rows[index].cells].map(c => c.textContent));
            return parts.every(p => p[0] === row.cells[0].textContent && p.at(-1) === row.cells[8].textContent)
                && parts.flatMap(p => p.slice(1,-1)).join('|') === [...row.cells].slice(1,-1).map(c => c.textContent).join('|');
        });
    }""")
    restore_and_check(page, before)
    page.locator('#wide-print-fixture').evaluate('(node) => node.remove()')


def check_report(browser, origin, sample, engine):
    context = browser.new_context(viewport={'width': 1440, 'height': 1000})
    page = context.new_page()
    errors, failed_requests, external, posts, held = [], [], [], [], []
    page.on('pageerror', lambda error: errors.append(str(error)))
    context.on('response', lambda response: failed_requests.append(response.status) if response.status >= 400 else None)
    context.on('request', lambda request: posts.append(request.url) if request.method == 'POST' else None)

    def intercept(route):
        url = route.request.url
        if not url.startswith(origin + '/'):
            external.append(url)
            route.abort()
        elif not held and '/assets/' in url and '.svg' in url and 'ont-libqc-mark' not in url:
            # Hold one real plot until the button has been pressed. This checks
            # preparation waits for assets rather than exporting a partial page.
            held.append(route)
        else:
            route.continue_()

    context.route('**/*', intercept)
    try:
        response = page.goto(origin + f'/example/samples/{sample}.html?view={GRANT}', wait_until='domcontentloaded')
        assert response.status == 200
        csp = response.headers['content-security-policy']
        assert 'allow-modals' in csp and 'allow-same-origin' not in csp
        assert "connect-src 'none'" in csp
        wait_js(page, 'Array.from(document.images).some(i => !i.complete)')
        assert held, 'Expected one delayed report plot'
        assert page.locator('details:not([open])').count() > 0
        assert page.locator('a[href^="#"]').count() > 0
        assert page.evaluate("Array.from(document.querySelectorAll('a[href^=\"#\"]')).some(a => a.href.includes('view='))")
        assert page.locator('#report-export-content').input_value() == 'overview'
        page.select_option('#report-export-content', 'full')
        before = snapshot(page)
        page.evaluate("""() => {
            window.printCalls = 0;
            window.print = () => {window.printCalls += 1;};
        }""")
        button = page.locator('#report-export-pdf')
        button.click()
        wait_js(page, 'document.getElementById("report-export-pdf").disabled')
        assert button.get_attribute('aria-busy') == 'true'
        assert page.evaluate('window.printCalls') == 0
        held[0].continue_()
        wait_js(page, 'window.printCalls === 1')
        assert page.locator('#report-export-status').inner_text() == ''
        assert page.evaluate('Array.from(document.images).every(i => i.complete && i.naturalWidth > 0)')
        assert_prepared(page)
        restore_and_check(page, before)

        # The default overview retains the user's expanded sections, excludes
        # unopened details, and does not append a long metric glossary.
        page.select_option('#report-export-content', 'overview')
        before = snapshot(page)
        page.evaluate('window.ONTReportExport.prepare()')
        assert page.evaluate("Array.from(document.querySelectorAll('details:not([open])')).every(d => d.classList.contains('report-print-omit'))")
        assert page.locator('.report-print-definitions').count() == 0
        assert page.locator('#custom-regions.report-print-omit').count() == (1 if sample == 'HL240' else 0)
        page.emulate_media(media='print')
        assert page.evaluate("Array.from(document.querySelectorAll('.report-print-omit')).every(d => getComputedStyle(d).display === 'none')")
        if engine == 'chromium':
            page.pdf(path=str(OUTPUT / f'{sample}-overview.pdf'), prefer_css_page_size=True,
                     print_background=True, display_header_footer=False)
        restore_and_check(page, before)
        page.select_option('#report-export-content', 'full')

        custom = draw_custom_bed(page) if sample == 'HL240' else None
        before = snapshot(page)
        # Keyboard/OS printing must prepare the same full report, idempotently.
        page.evaluate('window.dispatchEvent(new Event("beforeprint")); window.dispatchEvent(new Event("beforeprint"))')
        assert_prepared(page)
        page.emulate_media(media='print')
        for width in (703, 704):
            page.set_viewport_size({'width': width, 'height': 1000})
            assert page.locator('.report-print-definitions').is_visible()
            assert button.is_hidden() and page.locator('#report-export-status').is_hidden()
            assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
            edges = page.evaluate("""() => Array.from(document.querySelectorAll(
                'table,.library-diagram,.library-end,.identity-card,.library-sublibrary header>span'
            )).filter(e => e.getBoundingClientRect().width).map(e => {
                const r = e.getBoundingClientRect();
                return {class: e.className, left: r.left, right: innerWidth-r.right,
                        overflow: e.scrollWidth-e.clientWidth};
            })""")
            assert edges and all(e['left'] >= 3 and e['right'] >= 3 and e['overflow'] <= 2 for e in edges), (engine, sample, width, edges)
            if custom is not None:
                assert page.locator('#custom-output svg').is_visible()
                assert page.locator('#custom-output').inner_html() == custom
                assert 'export-example.bed' in page.locator('#custom-current-source').inner_text()
            page.screenshot(path=str(OUTPUT / f'{sample}-{engine}-print-{width}.png'))
        if engine == 'chromium':
            page.pdf(path=str(OUTPUT / f'{sample}.pdf'), prefer_css_page_size=True,
                     print_background=True, display_header_footer=False)
        restore_and_check(page, before)
        if custom is not None:
            assert page.locator('#custom-output svg').is_visible()
            assert page.locator('#custom-output').inner_html() == custom
        check_wide_table(page)
        assert not errors and not failed_requests, (errors, failed_requests)
        assert not external and not posts, (external, posts)
        assert not context.cookies()
        return {'browser': engine, 'sample': sample, 'waits_for_plot': True,
                'all_details_and_metric_notes': True, 'grant_and_fragment_links_removed': True,
                'state_and_links_restored': True, 'print_widths': [703, 704],
                'minimum_inset_px': 3, 'custom_bed_preserved': custom is not None}
    finally:
        context.close()


def check_refreshed_layout(browser, original, refreshed, sample):
    """Retained reports must keep export controls and print geometry after refresh."""
    pages = [browser.new_page() for _ in range(2)]
    selectors = ['.report-hero', '.report-hero h1', '.hero-meta',
                 '.report-export-controls', '#report-export-content', '#report-export-pdf',
                 '.report-nav', '.scope-card', '.qc-metrics', '.metric-group',
                 '.identity-card', '.profile-scope', '.plot-card']
    properties = ['display', 'position', 'float', 'grid-template-columns', 'flex-direction',
                  'gap', 'margin', 'padding', 'font-size', 'line-height', 'color',
                  'background-color', 'border', 'width', 'height', 'max-width',
                  'overflow', 'break-inside', 'break-before', 'break-after']
    try:
        for page, root in zip(pages, (original, refreshed)):
            page.goto((root/'samples'/f'{sample}.html').as_uri())
            page.evaluate('document.fonts.ready')
        for media, width in [('screen',1440),('screen',900),('screen',390),('screen',320),('print',794)]:
            snapshots=[]
            for page in pages:
                page.set_viewport_size({'width':width,'height':1000})
                page.emulate_media(media=media)
                if media == 'print':
                    # The print media change event may run after Playwright reads the layout.
                    page.evaluate('window.ONTReportExport.prepare()')
                snapshots.append(page.evaluate("""({selectors,properties}) => selectors.map(selector =>
                    [...document.querySelectorAll(selector)].map(node => {
                        const style=getComputedStyle(node), rect=node.getBoundingClientRect();
                        return {style:properties.map(p=>style.getPropertyValue(p)),
                                box:[rect.x,rect.width,rect.height],
                                top:node.closest('.report-hero') ? rect.y : null};
                    }))""", {'selectors':selectors,'properties':properties}))
            differences=[(selectors[i],a,b) for i,(a,b) in enumerate(zip(*snapshots)) if a!=b]
            assert not differences, (sample,media,width,'Refresh changed report layout',differences[:1])
            button=pages[1].locator('#report-export-pdf')
            if media=='print':
                assert button.is_hidden()
                assert pages[1].locator('.plot-card').first.evaluate('(e)=>getComputedStyle(e).breakInside')=='avoid'
            else:
                assert button.is_visible()
                assert pages[1].evaluate('document.documentElement.scrollWidth <= innerWidth+1')
                if width>=900:
                    controls=pages[1].locator('.report-export-controls').bounding_box()
                    hero=pages[1].locator('.report-hero').bounding_box()
                    assert controls['x'] > hero['x']+hero['width']/2
    finally:
        for page in pages: page.close()


def check_composition_tables(browser):
    ids = [f'lib{i:02}' for i in range(1, 7)]
    panel = ['hg38', 'dm6', 'sacCer3']
    hits = [[(species, 60, True, ('chr1', 0))] for species in panel]
    data = classify([sid for sid in ids for _ in panel], hits * len(ids), panel, 'dm6')
    scopes = [dict(scope_id=sid, collision_observed_pairs=i, collision_oe=2.5,
                   mycoplasma_status='NOT_DETECTED', mycoplasma_screened_reads=1)
              for i, sid in enumerate(['total', *ids])]
    metrics = scoped_table_html(scopes, ['collision_oe', 'mycoplasma_status'], {})
    document = add_pdf_export("<html><head><style>" + REPORT_CSS + "</style></head><body>"
                              "<header class='report-hero'><h1>Sample</h1></header>"
                              + render_species(data) + metrics + "</body></html>")
    page = browser.new_page(viewport={'width': 1250, 'height': 900})
    def check_bar_labels():
        segments = page.locator('.species-bar').first.locator('.species-segment')
        assert segments.locator('b').all_inner_texts() == ['33.3%'] * len(panel)
        assert all(label.is_visible() for label in segments.locator('b').all())
        assert segments.evaluate_all('''bars => bars.map(bar => [
            getComputedStyle(bar).backgroundColor, getComputedStyle(bar.querySelector("b")).color
        ])''') == [['rgb(205, 32, 44)', 'rgb(255, 255, 255)'],
                  ['rgb(0, 178, 169)', 'rgb(0, 0, 0)'],
                  ['rgb(116, 118, 120)', 'rgb(255, 255, 255)']]

    try:
        page.set_content(document)
        check_bar_labels()
        page.locator('details').evaluate('(node) => node.open = true')
        assert page.locator('.species-counts').first.evaluate('(t) => t.clientWidth <= t.parentNode.clientWidth + 1')
        assert page.locator('.metric-note').first.inner_text() == 'Low depth'
        assert page.locator('.metric-reading').filter(has_text='Not detected').first.evaluate("""node => {
            const label = node.lastElementChild, range = document.createRange();
            range.selectNodeContents(label);
            return range.getClientRects().length === 1 &&
                label.getBoundingClientRect().left - node.firstElementChild.getBoundingClientRect().right < 10;
        }""")
        page.select_option('#report-export-content', 'full')
        page.evaluate('window.ONTReportExport.prepare()')
        page.emulate_media(media='print')
        page.set_viewport_size({'width': 703, 'height': 1000})
        check_bar_labels()
        for kind, count in [('species-counts', 2), ('qc-metrics', 2)]:
            copies = page.locator('.report-print-tables .' + kind)
            assert copies.count() == count
            assert copies.evaluate_all('(tables) => tables.every(t => t.scrollWidth <= t.parentNode.clientWidth + 2)')
        assert page.locator('.report-print-tables .species-counts thead th').all_text_contents() == [
            'Species', *ids[:4], 'Species', *ids[4:]]
        assert page.locator('.report-print-tables .metric-note').count() == 5
        page.evaluate('window.ONTReportExport.restore()')
        assert page.locator('.report-print-tables').count() == 0
    finally:
        page.close()


def check_sublibrary_percentages(browser, engine):
    scopes = [dict(scope_id=sid, library_units=count)
              for sid, count in [('lib01', 400), ('lib02', 240), ('lib03', 180), ('lib04', 120)]]
    section = render_sublibrary_composition(dict(library_units=1000, sub_libraries=scopes,
                                                 index_pairing=dict(ambiguous_units=0, unassigned_units=60)))
    doc = add_pdf_export(f"<html><head><style>{REPORT_CSS}</style></head><body>"
                         f"<header class='report-hero'><h1>Sample</h1></header>{section}</body></html>")
    page = browser.new_page()
    try:
        page.set_content(doc)
        labels = page.locator('.mixture-segment-text')
        percentages = ['40.0%', '24.0%', '18.0%', '12.0%', '6.0%']
        assert labels.all_text_contents() == percentages
        for media, width, visible in [('screen', 1250, [True, True, True, True, True]),
                                      ('screen', 375, [True, True, True, False, False]),
                                      ('print', 703, [True, True, True, True, False])]:
            page.set_viewport_size({'width': width, 'height': 900})
            if media == 'print':
                page.evaluate('window.ONTReportExport.prepare()')
            page.emulate_media(media=media)
            assert [label.is_visible() for label in labels.all()] == visible
            assert labels.evaluate_all("""nodes => nodes.map(label => [
                getComputedStyle(label.parentElement).backgroundColor, getComputedStyle(label).color
            ])""") == [['rgb(230, 75, 53)', 'rgb(255, 255, 255)'],
                       ['rgb(77, 187, 213)', 'rgb(255, 255, 255)'],
                       ['rgb(0, 160, 135)', 'rgb(255, 255, 255)'],
                       ['rgb(60, 84, 136)', 'rgb(255, 255, 255)'],
                       ['rgb(170, 181, 196)', 'rgb(0, 0, 0)']]
            assert page.locator('.mixture-legend strong').all_inner_texts() == percentages
        if engine == 'chromium':
            page.pdf(path=str(OUTPUT/'sublibrary-composition.pdf'), format='A4', print_background=True)
    finally:
        page.close()


def check_spike_matrix_print(browser, engine):
    labels = ['Saccharomyces cerevisiae', 'Schizosaccharomyces pombe',
              'Lambda phage', 'PhiX174', 'Escherichia coli', 'Drosophila melanogaster']
    screen = {'candidates': {str(i): {'label': name} for i, name in enumerate(labels)}}
    libs = [{'scope_id': f'lib{i+1:02}', 'spike_in_screened_reads': 2300-i*100,
             'spike_in_species_counts': {}} for i in range(6)]
    table = spike_matrix_html(screen, libs, 'synthetic')
    # Retained reports still have the old inline 1% width: print must repair them too.
    table = table.replace('<th scope="col">Sublibrary',
                          '<th scope="col" style="width:1%;min-width:76px;max-width:120px">Sublibrary')
    document = add_pdf_export('<html><head><meta name="viewport" content="width=device-width,initial-scale=1">'
                             '<style>'+REPORT_CSS+'</style></head><body>'
                             '<header class="report-hero"><h1>Print layout check</h1></header>'+table+'</body></html>')
    page = browser.new_page()
    try:
        page.set_content(document)
        for media, width in [('screen',1440), ('screen',390), ('print',703), ('print',794)]:
            page.set_viewport_size({'width':width,'height':1000})
            page.emulate_media(media=media)
            widths = page.locator('.spike-matrix thead th').evaluate_all(
                '(es)=>es.slice(2).map(e=>e.getBoundingClientRect().width)')
            assert max(widths)-min(widths)<1
            assert page.evaluate('document.documentElement.scrollWidth<=innerWidth+1')
            if media == 'print':
                first = page.locator('.spike-matrix thead th').first
                assert first.evaluate('(e)=>e.clientWidth>=e.scrollWidth'), 'Sublibrary heading clipped'
                assert page.locator('.spike-matrix thead').bounding_box()['height']<50, 'Header wraps letter by letter'
                assert page.locator('.spike-matrix th,.spike-matrix td').evaluate_all(
                    '(es)=>es.every(e=>e.scrollWidth<=e.clientWidth+1)'), 'Printed matrix cell overflows'
                page.screenshot(path=str(OUTPUT/f'spike-matrix-{engine}-{width}.png'))
        if engine == 'chromium':
            for paper in ('A4', 'Letter'):
                page.pdf(path=str(OUTPUT/f'spike-matrix-{paper}.pdf'), format=paper,
                         print_background=True, prefer_css_page_size=False)
    finally:
        page.close()


results = []
with tempfile.TemporaryDirectory() as directory:
    sock = socket.socket()
    sock.bind(('127.0.0.1', 0))
    origin = f'http://127.0.0.1:{sock.getsockname()[1]}'
    examples=Path(directory)/'examples'
    shutil.copytree(ROOT/'example_report',examples)
    refresh(examples, apply=True)
    settings = Settings(examples=examples, state=Path(directory) / 'state', public_url=origin,
                        allowed_hosts=['127.0.0.1'], secret='local-export-test-' * 4,
                        submission_key='', public_catalog={'references': []},
                        min_free_bytes=0, require_https=False, mail_enabled=False)
    server = uvicorn.Server(uvicorn.Config(create_app(settings), log_level='error',
                                          access_log=False, limit_concurrency=256))
    thread = threading.Thread(target=server.run, kwargs={'sockets': [sock]}, daemon=True)
    thread.start()
    deadline = time.monotonic() + 15
    while not server.started:
        if not thread.is_alive() or time.monotonic() > deadline:
            raise RuntimeError('Temporary portal failed to start')
        time.sleep(.02)
    try:
        with sync_playwright() as playwright:
            for engine in sys.argv[2:] or ['chromium', 'webkit']:
                options = {'executable_path': os.environ['ONT_WEBKIT_EXECUTABLE']} if engine == 'webkit' and os.environ.get('ONT_WEBKIT_EXECUTABLE') else {}
                browser = getattr(playwright, engine).launch(**options)
                try:
                    check_composition_tables(browser)
                    check_sublibrary_percentages(browser, engine)
                    check_spike_matrix_print(browser, engine)
                    for sample in ('HL120', 'HL240'):
                        check_refreshed_layout(browser, ROOT/'example_report', examples, sample)
                        results.append(check_report(browser, origin, sample, engine))
                finally:
                    browser.close()
    finally:
        server.should_exit = True
        thread.join(timeout=10)
        sock.close()
(OUTPUT / 'validation.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
