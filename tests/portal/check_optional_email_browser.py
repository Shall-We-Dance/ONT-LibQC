"""Submit without email, save private links, then open completed example reports.

Runs a temporary local portal with mail enabled and a mocked sender. Completion is
simulated by publishing public example reports; no analysis or real email is run.
"""
import json
import os
from pathlib import Path
import socket
import sys
import tempfile
import threading
import time
from unittest.mock import patch

import uvicorn
from playwright.sync_api import sync_playwright
from ont_portal.app import create_app
from ont_portal.notifications import notify
from ont_portal.settings import Settings
from ont_portal.worker import publish_report

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(sys.argv[1]) if len(sys.argv)>1 else Path('test-results/optional-email')
OUTPUT.mkdir(parents=True,exist_ok=True)
CATALOG={'references':[{'id':name,'assembly':name,'label':name,'ready':True,
    'mappers':['bowtie2','bwa-mem2'],'microc_mapper':'bwa-mem2','primary_chromosomes':['chr1'],
    'mycoplasma':False,'spike_in':False} for name in ['dm6','hg38']]}


def wait_js(page,expression,timeout=30000):
    deadline=time.monotonic()+timeout/1000
    while time.monotonic()<deadline:
        if page.evaluate(expression):return
        page.wait_for_timeout(50)
    raise AssertionError('Timed out: '+expression)


results=[]
with tempfile.TemporaryDirectory() as directory:
    root=Path(directory);sock=socket.socket();sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
    origin=f'http://127.0.0.1:{port}'
    settings=Settings(state=root/'state',public_url=origin,allowed_hosts=['127.0.0.1'],
        secret='local-browser-test-'*4,submission_key='',public_catalog=CATALOG,
        min_free_bytes=0,require_https=False,mail_enabled=True)
    app=create_app(settings);store=app.state.store;store.heartbeat()
    (root/'smtp.yaml').write_text('enabled: true\nfrom_address: reports@example.org\n')
    (root/'secrets.json').write_text(json.dumps({'signing_key':settings.secret}))
    mail_config={'state':str(settings.state),'mail_enabled':True,'smtp_file':str(root/'smtp.yaml'),
                 'secrets_file':str(root/'secrets.json'),'public_url':origin}
    server=uvicorn.Server(uvicorn.Config(app,log_level='error',access_log=False,limit_concurrency=256))
    thread=threading.Thread(target=server.run,kwargs={'sockets':[sock]},daemon=True);thread.start()
    while not server.started:time.sleep(.02)
    try:
        with sync_playwright() as p, patch('ont_portal.notifications.send') as sender:
            for engine in sys.argv[2:] or ['chromium','webkit']:
                browser=getattr(p,engine).launch(**({'executable_path':os.environ['ONT_WEBKIT_EXECUTABLE']} if engine=='webkit' and os.environ.get('ONT_WEBKIT_EXECUTABLE') else {}))
                try:
                    for amount,width in [(1,390),(2,1440)]:
                        store.heartbeat();context=browser.new_context(viewport={'width':width,'height':960},accept_downloads=True)
                        page=context.new_page();bad=[];posts=[]
                        context.on('response',lambda response:bad.append(response.status) if response.status>=400 else None)
                        page.on('request',lambda request:posts.append(request.url) if request.method=='POST' else None)
                        page.goto(origin+'/?example=HL120')
                        wait_js(page,'files.length === 1 && !loadingExample')
                        if amount==2:
                            # Seed a mixed batch; the public example chooser lives on its own page.
                            page.evaluate("loadExample('HL240')")
                            wait_js(page,'files.length === 2 && !loadingExample')
                        email=page.get_by_label('Email address (optional)',exact=True)
                        assert not email.evaluate('(input) => input.required')
                        assert email.input_value()==''
                        page.locator('#consent').check();page.locator('#data-policy').check()
                        email.fill('invalid-address')
                        page.get_by_role('button',name='Start analysis').click()
                        assert email.evaluate('(input) => input.validity.typeMismatch')
                        assert not posts,'Invalid supplied email must not create a submission.'
                        email.fill('')
                        assert page.locator('#submission-form').evaluate('(form) => form.checkValidity()')
                        page.get_by_role('button',name='Start analysis').click()
                        wait_js(page,f'tasks.length === {amount} && !submitting')
                        assert page.evaluate("tasks.every(t => t.job.status === 'queued' && t.job.mail_status === 'disabled' && !t.failure)")
                        assert 'No email will be sent.' in page.locator('#tasks-description').inner_text()
                        assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
                        job_ids=page.evaluate('tasks.map(t => t.id)');return_link=page.url
                        button=page.get_by_role('button',name='Save private link' if amount==1 else 'Save private links',exact=True)
                        with page.expect_download() as pending:button.click()
                        download=pending.value;assert download.failure() is None
                        saved=Path(download.path()).read_text()
                        assert saved.count('/jobs/')==amount
                        for job_id in job_ids:
                            row=store.get(job_id)
                            assert row['email']=='' and row['mail_status']=='disabled' and row['mail_attempts']==0
                            assert job_id in saved
                            publish_report(ROOT/'example_report',settings.state/'reports'/job_id)
                            store.update(job_id,status='succeeded',progress=100)
                        notify(mail_config);sender.assert_not_called()
                        assert not bad,bad
                        assert not context.cookies()
                        context.close()
                        # Reopen in a clean browser context using the saved private URL.
                        context=browser.new_context(viewport={'width':width,'height':960},accept_downloads=True)
                        context.on('response',lambda response:bad.append(response.status) if response.status>=400 else None)
                        page=context.new_page();page.goto(return_link)
                        wait_js(page,f"tasks.length === {amount} && tasks.every(t => t.job?.status === 'succeeded')")
                        assert 'No email will be sent.' in page.locator('#tasks-description').inner_text()
                        images={}
                        for sample in ['HL120']+(['HL240'] if amount==2 else []):
                            with page.expect_popup() as opened:page.get_by_role('link',name='Open '+sample+' report').click()
                            report=opened.value;report.wait_for_load_state('load')
                            wait_js(report,'Array.from(document.images).every(img => img.complete && img.naturalWidth > 0)')
                            images[sample]=report.locator('img').count();report.close()
                        with page.expect_download() as pending:page.get_by_role('link',name='Download report ZIP').first.click()
                        assert pending.value.failure() is None
                        assert not context.cookies() and not bad,bad
                        page.screenshot(path=str(OUTPUT/f'{engine}-{width}-tasks.png'))
                        results.append({'browser':engine,'files':amount,'width':width,'blank_email_submission':True,
                            'invalid_email_blocked':True,'saved_links':True,'reopened_without_storage':True,
                            'report_images':images,'zip_download':True,'emails_sent':sender.call_count,'http_errors':bad})
                        context.close()
                finally:browser.close()
    finally:
        server.should_exit=True;thread.join(timeout=10);sock.close()
(OUTPUT/'validation.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
