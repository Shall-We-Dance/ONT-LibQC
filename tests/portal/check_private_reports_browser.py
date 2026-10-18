"""Exercise complete private reports without cookies in Chromium and WebKit.

Uses public example reports in temporary local jobs; no worker or mail is started.
All browser requests have Cookie removed, and every response must omit Set-Cookie.
"""
import json
import os
from pathlib import Path
import socket
import sys
import tempfile
import threading
import time
import zipfile

import uvicorn
from playwright.sync_api import sync_playwright
from ont_portal.app import create_app
from ont_portal.settings import Settings
from ont_portal.store import token_for
from ont_portal.worker import publish_report

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(sys.argv[1]) if len(sys.argv)>1 else Path('test-results/private-reports')
OUTPUT.mkdir(parents=True, exist_ok=True)


def wait_js(page, expression, timeout=30000):
    deadline=time.monotonic()+timeout/1000
    while time.monotonic()<deadline:
        if page.evaluate(expression): return
        page.wait_for_timeout(50)
    raise AssertionError('Timed out: '+expression)


results=[]
with tempfile.TemporaryDirectory() as directory:
    root=Path(directory)
    sock=socket.socket();sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
    origin=f'http://127.0.0.1:{port}'
    config=Settings(state=root/'state',public_url=origin,allowed_hosts=['127.0.0.1'],
        secret='local-browser-test-'*4,submission_key='',public_catalog={'references':[]},
        min_free_bytes=0,require_https=False,mail_enabled=False)
    app=create_app(config);now=time.time();jobs={}
    for index,sample in enumerate(['HL120','HL240']):
        job_id=str(index+1)*32;token=token_for(config.secret,job_id)
        with app.state.store.connect() as db:
            db.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires) VALUES (?,?,?,?,?,?,?,?,?,?)',
                (job_id,now,now,'succeeded','',sample+'.fastq',10,json.dumps({'sample':sample}),'test',now+7*86400))
        publish_report(ROOT/'example_report',config.state/'reports'/job_id)
        jobs[sample]=(job_id,token)
    server=uvicorn.Server(uvicorn.Config(app,log_level='error',access_log=False,limit_concurrency=256))
    thread=threading.Thread(target=server.run,kwargs={'sockets':[sock]},daemon=True);thread.start()
    while not server.started:time.sleep(.02)
    try:
        with sync_playwright() as p:
            for engine in sys.argv[2:] or ['chromium','webkit']:
                browser=getattr(p,engine).launch(**({'executable_path':os.environ['ONT_WEBKIT_EXECUTABLE']} if engine=='webkit' and os.environ.get('ONT_WEBKIT_EXECUTABLE') else {}))
                try:
                    context=browser.new_context(accept_downloads=True)
                    context.route('**/*',lambda route:route.continue_(headers={k:v for k,v in route.request.all_headers().items() if k.lower()!='cookie'}))
                    failures=[];cookie_responses=[];requests=[]
                    context.on('response',lambda response:failures.append({'status':response.status,'path':response.url.split('?')[0]}) if response.status>=400 else None)
                    context.on('response',lambda response:cookie_responses.append(response.url.split('?')[0]) if 'set-cookie' in response.all_headers() else None)
                    context.on('request',lambda request:requests.append(request.url))
                    page=context.new_page();images={}
                    for sample,(job_id,token) in jobs.items():
                        private=origin+f'/jobs/{job_id}#'+token
                        page.goto(private)
                        report_link=page.get_by_role('link',name='Open '+sample+' report')
                        report_link.wait_for()
                        assert '?view=' in report_link.get_attribute('href')
                        with page.expect_popup() as opened:report_link.click()
                        report=opened.value;report.wait_for_load_state('load')
                        wait_js(report,"Array.from(document.images).every(img => img.complete && img.naturalWidth > 0)")
                        images[sample]=report.locator('img').count()
                        assert images[sample]>0
                        assert token not in report.url
                        assert report.evaluate('window.opener === null')
                        report.close()
                        with page.expect_download() as pending:page.get_by_role('link',name='Download report ZIP').click()
                        download=pending.value;assert download.failure() is None
                        with zipfile.ZipFile(download.path()) as archive:
                            assert 'report/samples/'+sample+'.html' in archive.namelist()
                            assert '?view=' not in archive.read('report/samples/'+sample+'.html').decode()
                        # Simulate a tab resuming after the temporary report links expired.
                        page.evaluate('tasks[0].job.report_access_expires = 0')
                        with page.expect_popup() as renewed:report_link.click()
                        report=renewed.value
                        report.wait_for_url('**/reports/**');report.wait_for_load_state('load')
                        wait_js(report,"Array.from(document.images).every(img => img.complete && img.naturalWidth > 0)")
                        assert report.evaluate('window.opener === null')
                        report.close()
                        page.evaluate('tasks[0].job.report_access_expires = 0')
                        with page.expect_download() as renewed_download:page.get_by_role('link',name='Download report ZIP').click()
                        assert renewed_download.value.failure() is None
                        # A newly opened private link must work with empty browser storage.
                        page.reload();page.get_by_role('link',name='Open '+sample+' report').wait_for()
                        assert page.evaluate('localStorage.length === 0 && sessionStorage.length === 0')
                    assert not context.cookies(), 'The portal created a browser cookie.'
                    assert not cookie_responses, cookie_responses
                    assert not failures, failures
                    assert all(token not in url for _,token in jobs.values() for url in requests), 'Owner key leaked into an HTTP URL.'
                    assert not any('/session' in url for url in requests)
                    results.append({'browser':engine,'report_images':images,'zip_downloads':True,
                        'expired_link_renewal':True,'private_link_reload':True,'cookies':0,'set_cookie_responses':0,'http_errors':failures})
                    context.close()
                finally:browser.close()
    finally:
        server.should_exit=True;thread.join(timeout=10);sock.close()
(OUTPUT/'validation.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
