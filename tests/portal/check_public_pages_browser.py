"""Run real portal routes and complete example reports in Chromium and WebKit.

No external data, live submissions or emails are sent. A temporary local server
accepts the example uploads; no analysis worker is started.
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

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(sys.argv[1]) if len(sys.argv)>1 else Path('test-results/public-pages')
OUTPUT.mkdir(parents=True, exist_ok=True)
catalog={'references':[{'id':name,'assembly':name,'label':label,'ready':True,
    'mappers':['bowtie2','bwa-mem2'],'microc_mapper':'bwa-mem2',
    'primary_chromosomes':['chr1'] if name=='hg38' else ['chr2L'],
    'mycoplasma':False,'spike_in':False} for name,label in [('hg38','Human · hg38'),('dm6','Fruit fly · dm6')]]}
results=[]
def wait_js(page, expression, timeout=30000):
    # Poll through the browser inspector, without injecting eval under the site's CSP.
    deadline=time.monotonic()+timeout/1000
    while time.monotonic()<deadline:
        if page.evaluate(expression): return
        page.wait_for_timeout(50)
    raise AssertionError('Timed out: '+expression)

with tempfile.TemporaryDirectory() as directory:
    sock=socket.socket();sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
    origin=f'http://127.0.0.1:{port}'
    config=Settings(state=Path(directory)/'state',public_url=origin,allowed_hosts=['127.0.0.1'],
        secret='local-browser-test-'*4,submission_key='',public_catalog=catalog,
        min_free_bytes=0,require_https=False,mail_enabled=False)
    app=create_app(config);app.state.store.heartbeat()
    server=uvicorn.Server(uvicorn.Config(app,log_level='error',access_log=False,limit_concurrency=256))
    thread=threading.Thread(target=server.run,kwargs={'sockets':[sock]},daemon=True);thread.start()
    while not server.started:time.sleep(.02)
    try:
        with sync_playwright() as p:
            for engine in sys.argv[2:] or ['chromium','webkit']:
                browser=getattr(p,engine).launch(**({'executable_path':os.environ['ONT_WEBKIT_EXECUTABLE']} if engine=='webkit' and os.environ.get('ONT_WEBKIT_EXECUTABLE') else {}))
                try:
                    for width in [1440,390]:
                        context=browser.new_context(viewport={'width':width,'height':960})
                        page=context.new_page();failures=[];posts=[]
                        context.on('response',lambda response: failures.append({'url':response.url,'status':response.status}) if response.status>=400 else None)
                        context.on('weberror',lambda error: failures.append({'error':str(error.error)}))
                        context.on('request',lambda request: posts.append(request.url) if request.method=='POST' else None)
                        page.goto(origin+'/?example=HL240')
                        wait_js(page, "files.length === 1 && !loadingExample")
                        assert page.evaluate("files[0].sample.value === 'HL240' && files[0].ref.value === 'hg38' && files[0].assay.value === 'chromatin' && files[0].structure.value().mode === 'auto'")
                        assert page.evaluate("location.search === ''")
                        assert not page.locator('#example-inputs').count()
                        assert not page.get_by_role('button',name='Load HL240 example').count()
                        assert not page.locator('#consent').is_checked() and not page.locator('#data-policy').is_checked()
                        assert page.locator('.topbar nav [aria-current=page]').inner_text()=='New analysis'
                        page.locator('.topbar nav').get_by_role('link',name='Example reports',exact=True).click()
                        page.wait_for_url(origin+'/example/index.html')
                        examples=page
                        assert examples.url==origin+'/example/index.html'
                        assert len(context.pages)==1
                        assert examples.evaluate('document.documentElement.scrollWidth <= innerWidth')
                        accent=examples.evaluate("getComputedStyle(document.documentElement).getPropertyValue('--blue').trim()")
                        accent_rgb='rgb('+', '.join(str(int(accent[i:i+2],16)) for i in [1,3,5])+')'
                        for sample in ['HL240','HL120']:
                            primary=examples.get_by_role('link',name='Open '+sample+' report (opens in a new tab)',exact=True)
                            assert primary.evaluate('(link) => getComputedStyle(link).backgroundColor')==accent_rgb
                            secondary=examples.get_by_role('link',name='Load '+sample+' example',exact=True)
                            assert secondary.evaluate('(link) => getComputedStyle(link).backgroundColor')=='rgba(0, 0, 0, 0)'
                        examples.screenshot(path=str(OUTPUT/f'example-actions-{engine}-{width}.png'),full_page=True)
                        for sample in ['HL240','HL120']:
                            with examples.expect_popup() as opened:
                                examples.get_by_role('link',name='Open '+sample+' report (opens in a new tab)',exact=True).click()
                            report=opened.value;report.wait_for_load_state('load')
                            wait_js(report, "Array.from(document.images).every(img => img.complete && img.naturalWidth > 0)")
                            assert report.url==origin+f'/example/samples/{sample}.html'
                            assert report.evaluate('window.opener === null')
                            assert examples.url==origin+'/example/index.html'
                            report.close()
                        # Keep the portable index's absolute public links local in this test.
                        examples.locator('.card-actions a[href^="https://ontlibqc.org/"]').evaluate_all(
                            '(links, origin) => links.forEach(link => link.href = origin + new URL(link.href).pathname + new URL(link.href).search)',origin)
                        with examples.expect_download() as pending:
                            examples.get_by_role('link',name='Download HL120 FASTQ',exact=True).click()
                        assert pending.value.failure() is None
                        assert pending.value.suggested_filename=='HL120.fastq.gz'
                        examples.get_by_role('link',name='Load HL120 example',exact=True).click()
                        examples.wait_for_load_state('load')
                        wait_js(examples, "files.length === 1 && !loadingExample")
                        assert examples.evaluate("files[0].ref.value === 'dm6' && files[0].assay.value === 'microc'")
                        assert not examples.locator('#consent').is_checked() and not examples.locator('#data-policy').is_checked()
                        assert page.url==origin+'/' and page.evaluate("files.length === 1 && files[0].sample.value === 'HL120'")
                        assert page.locator('.info-panel a').count()==1
                        page.locator('.info-panel a').click();page.wait_for_url(origin+'/example/index.html')
                        assert len(context.pages)==1
                        page.goto(origin+'/?example=HL240')
                        wait_js(page, "files.length === 1 && !loadingExample")
                        # Exercise duplicate protection and a mixed batch without restoring inline controls.
                        page.evaluate("loadExample('HL240')")
                        assert page.evaluate('files.length')==1
                        page.evaluate("loadExample('HL120')")
                        wait_js(page, "files.length === 2 && !loadingExample")
                        assert page.evaluate("files[1].ref.value === 'dm6' && files[1].assay.value === 'microc'")
                        assert not posts,'Loading an example must not submit a job.'
                        assert not page.locator('#consent').is_checked() and not page.locator('#data-policy').is_checked()
                        assert page.evaluate('document.documentElement.scrollWidth <= innerWidth'), (engine,width,'upload overflow')
                        page.screenshot(path=str(OUTPUT/f'upload-{engine}-{width}.png'))
                        # One real local submission verifies that loaded files use the normal upload route.
                        if engine=='chromium' and width==1440:
                            page.locator('#email').fill('test@example.org');page.locator('#consent').check();page.locator('#data-policy').check()
                            page.get_by_role('button',name='Start analysis').click()
                            wait_js(page, "tasks.length === 2 && !submitting",timeout=30000)
                            assert page.evaluate("tasks.every(task => task.job.status === 'queued' && !task.failure)")
                            assert page.evaluate("tasks.map(task => task.job.reference + '/' + task.job.assay).join(',')")=='hg38/chromatin,dm6/microc'
                            assert page.locator('.topbar nav [aria-current]').count()==0
                        report_counts={}
                        for sample in ['HL120','HL240']:
                            page.goto(origin+f'/example/samples/{sample}.html',wait_until='load')
                            page.locator('details').evaluate_all('(nodes) => nodes.forEach(node => node.open = true)')
                            wait_js(page, "Array.from(document.images).every(img => img.complete && img.naturalWidth > 0)",timeout=30000)
                            report_counts[sample]=page.locator('img').count()
                        for route,active in [('/','New analysis'),('/help','Help'),('/license',None),('/example/index.html','Example reports')]:
                            page.goto(origin+route,wait_until='load')
                            assert page.evaluate('document.documentElement.scrollWidth <= innerWidth'),(engine,width,route,'overflow')
                            selected=page.locator('.topbar nav [aria-current=page]')
                            assert selected.count()==int(active is not None)
                            if active:
                                assert selected.inner_text()==active
                                assert selected.evaluate('(link) => getComputedStyle(link).color')==accent_rgb
                                assert selected.evaluate('(link) => getComputedStyle(link).fontWeight')=='600'
                            page.screenshot(path=str(OUTPUT/(route.strip('/').replace('/','-')+f'-{engine}-{width}.png')))
                        assert not failures, failures
                        results.append({'browser':engine,'width':width,'example_autofill':True,'primary_open_report':True,
                            'same_tab_example_index':True,'new_tab_reports':True,'consistent_active_navigation':True,
                            'fastq_download':True,'report_images':report_counts,'failed_requests':failures})
                        context.close()
                finally:browser.close()
    finally:
        server.should_exit=True;thread.join(timeout=10);sock.close()
(OUTPUT/'validation.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
