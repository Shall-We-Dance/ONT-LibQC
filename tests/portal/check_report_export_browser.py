"""Exercise PDF export on public reports through a temporary local portal.

Usage: PYTHONPATH=portal python tests/portal/check_report_export_browser.py OUTDIR [chromium webkit]
No submissions, worker, mail, external requests, or private reports are used.
ONT_WEBKIT_EXECUTABLE optionally selects a locally configured WebKit executable.
"""
import json
import os
from pathlib import Path
import socket
import sys
import tempfile
import threading
import time

import uvicorn
from playwright.sync_api import sync_playwright
from ont_portal.app import create_app
from ont_portal.settings import Settings


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


results = []
with tempfile.TemporaryDirectory() as directory:
    sock = socket.socket()
    sock.bind(('127.0.0.1', 0))
    origin = f'http://127.0.0.1:{sock.getsockname()[1]}'
    settings = Settings(state=Path(directory) / 'state', public_url=origin,
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
                    for sample in ('HL120', 'HL240'):
                        results.append(check_report(browser, origin, sample, engine))
                finally:
                    browser.close()
    finally:
        server.should_exit = True
        thread.join(timeout=10)
        sock.close()
(OUTPUT / 'validation.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
