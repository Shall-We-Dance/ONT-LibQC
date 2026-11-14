"""Report navigation remains usable without granting access to task operations."""
import io
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import zipfile
from urllib.parse import parse_qs, urlsplit
from fastapi.testclient import TestClient
from pydantic import ValidationError
from ont_portal.app import create_app
from ont_portal.report_entry import direct_report
from ont_portal.report_download import archive_filename
from ont_portal.schema import Options
from ont_portal.settings import Settings
from ont_portal.store import report_grant, token_for, valid_report_grant
from ont_portal.short_links import key_hash
from ont_portal.worker import publish_report
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'workflow/scripts'))
from custom_regions_ui import render_section, upgrade_html

JID='c'*32
SAMPLE='example'

class ReportAccessTests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        self.root=Path(tmp.name);self.now=time.time()
        self.settings=Settings(state=self.root/'state',public_url='https://testserver',allowed_hosts=['testserver'],
                               secret='test-only-'*8,submission_key='',public_catalog={},min_free_bytes=0)
        self.app=create_app(self.settings);self.store=self.app.state.store
        with self.store.connect() as db:
            db.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires) VALUES (?,?,?,?,?,?,?,?,?,?)',
                       (JID,self.now,self.now,'succeeded','','example.fastq',10,json.dumps({'sample':SAMPLE}),'test',self.now+7*86400))
        source=self.root/'source';(source/'samples').mkdir(parents=True);(source/'assets').mkdir()
        (source/'index.html').write_text('<a href="samples/example.html">Report</a>')
        (source/'samples/example.html').write_text('<html><body><a href="../assets/figure.svg" target="_blank" rel="noopener">Open figure ↗</a><footer class="report-footer"><p><a href="../index.html">All samples</a> · <a href="../data.tsv">Adapter table</a> · <a href="../data.tsv">Index assignments</a></p></footer></body></html>')
        (source/'assets/figure.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg"><text>Figure</text></svg>')
        (source/'data.tsv').write_text('a\tb\n1\t2\n')
        (source/'downloads.html').write_text('<html><body><a href="index.html">All samples</a><main class="report-section"><a href="data.tsv">data.tsv</a></main></body></html>')
        direct_report(source,SAMPLE)
        publish_report(source,self.settings.state/'reports'/JID)
        self.client=TestClient(self.app,base_url='https://testserver',headers={'X-Portal-Request':'1'})
        self.addCleanup(self.client.close)
        self.page=f'/reports/{JID}/samples/example.html'
        self.owner=token_for(self.settings.secret,JID)

    def task(self):
        response=self.client.get(f'/api/jobs/{JID}',headers={'Authorization':'Bearer '+self.owner})
        self.assertEqual(response.status_code,200)
        self.assertNotIn('set-cookie',response.headers)
        return response.json()

    def link_record(self,url):
        self.assertRegex(url,r'^/r/[A-Za-z0-9_-]{22}$')
        with self.store.connect() as db:
            row=db.execute('SELECT * FROM report_links WHERE key_hash=?',(key_hash(url.rsplit('/',1)[1]),)).fetchone()
        return dict(row) if row else None

    def test_private_task_issues_cookie_free_short_report_and_download_links(self):
        job=self.task()
        self.assertEqual(job['report_access_expires'],int(self.now+7*86400))
        for field in ['report_url','download_url']:
            url=job[field]
            self.assertNotIn(self.owner,url);self.assertNotIn(JID,url)
            self.assertEqual(self.link_record(url)['expires'],job['report_access_expires'])
            response=self.client.get(url)
            self.assertEqual(response.status_code,200)
            self.assertNotIn('set-cookie',response.headers)
        self.assertFalse(self.client.cookies)
        self.assertEqual(self.client.get(self.page).status_code,404)
        with patch('time.time',return_value=self.now+3700):
            self.assertEqual(self.task()['report_url'],job['report_url'])
            self.assertEqual(self.client.get(job['report_url']).status_code,200)

    def test_archive_download_uses_sample_and_analysis_date(self):
        root=self.settings.state/'reports'/JID
        path=root/'data/summary/provenance.json';path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps({'created_utc':'2026-01-01T01:00:00Z','analyses':{SAMPLE:{'analysis_completed_utc':'2026-02-03T04:00:00Z'}}}))
        response=self.client.get(self.task()['download_url'])
        self.assertEqual(response.status_code,200)
        self.assertIn('report_example_2026-02-03.zip',response.headers['content-disposition'])
        self.assertEqual(response.content,(root/'report.zip').read_bytes())
        # Other downloadable files keep their own names.
        grant=report_grant(self.settings.secret,JID,int(self.now)+100)
        response=self.client.get(f'/reports/{JID}/data.tsv?view={grant}')
        self.assertIn('data.tsv',response.headers['content-disposition'])

    def test_archive_name_has_safe_sample_and_missing_metadata_fallback(self):
        root=self.settings.state/'reports'/JID
        name=archive_filename(root,'../../bad\r\nInjected: true',0)
        self.assertEqual(name,'report_sample_1970-01-01.zip')
        self.assertEqual(archive_filename(root,'Valid_sample',0),'report_Valid_sample_1970-01-01.zip')

    def test_read_only_capabilities_cannot_access_or_change_tasks(self):
        key=self.task()['report_url'].rsplit('/',1)[1]
        grant=report_grant(self.settings.secret,JID,int(self.now)+100)
        for credential in [key,grant]:
            for method,path in [('GET',f'/api/jobs/{JID}'),('DELETE',f'/api/jobs/{JID}'),('POST',f'/api/jobs/{JID}/cancel')]:
                response=self.client.request(method,path,params={'view':credential},headers={'Authorization':'Bearer '+credential})
                self.assertEqual(response.status_code,404)
        self.assertEqual(self.client.post('/api/links/resolve',json={'key':key}).status_code,404)
        self.assertGreater(self.store.get(JID)['expires'],self.now)

    def test_legacy_navigation_redirects_without_extending_grant_lifetime(self):
        old=report_grant(self.settings.secret,JID,int(self.now)+100)
        for headers in ({},{'Authorization':'Bearer invalid'}):
            response=self.client.get(self.page,params={'view':old},headers=headers,follow_redirects=False)
            self.assertEqual(response.status_code,303)
            link=response.headers['location'];record=self.link_record(link)
            self.assertEqual(record['file'],'samples/example.html')
            self.assertEqual(record['expires'],int(self.now)+100)
            page=self.client.get(link)
            self.assertEqual(page.status_code,200)
            self.assertNotIn('?view=',page.text)
            children=re.findall(r'href="(/r/[A-Za-z0-9_-]{22})"',page.text)
            self.assertTrue(children)
            for child in children:
                self.assertEqual(self.link_record(child)['expires'],record['expires'])
                self.assertEqual(self.client.get(child).status_code,200)
            with patch('time.time',return_value=self.now+110):
                self.assertEqual(self.client.get(link).status_code,404)
                self.assertEqual(self.client.get(self.page,params={'view':old}).status_code,404)
                self.assertEqual(self.client.get(self.task()['report_url']).status_code,200)

    def test_browser_missing_access_gets_short_private_link_guidance(self):
        response=self.client.get(self.page,headers={'Accept':'text/html','Sec-Fetch-Dest':'document'},follow_redirects=False)
        self.assertEqual(response.status_code,303);self.assertEqual(response.headers['location'],'/r/')
        self.assertIn('private analysis link',self.client.get('/r/').text)
        self.assertNotIn('set-cookie',response.headers)
        self.assertEqual(self.client.get(self.page).status_code,404)
        self.assertEqual(self.client.get(f'/reports/{JID}/assets/figure.svg',headers={'Accept':'image/svg+xml','Sec-Fetch-Dest':'image'}).status_code,404)

    def test_tampered_and_other_report_credentials_are_rejected(self):
        grant=report_grant(self.settings.secret,JID,int(self.now)+100)
        other=report_grant(self.settings.secret,'d'*32,int(self.now)+300)
        tampered=grant[:-1]+('1' if grant[-1]!='1' else '2')
        for invalid in [other,tampered,self.owner,'bad',report_grant(self.settings.secret,JID,int(self.now)+7200)]:
            self.assertEqual(self.client.get(self.page,params={'view':invalid}).status_code,404)
        link=self.task()['report_url'];tampered=link[:-1]+('1' if link[-1]!='1' else '2')
        self.assertEqual(self.client.get(tampered).status_code,404)
        self.assertEqual(self.client.get(link+'/other-file').status_code,404)
        for file in ['../secret.json','.hidden.json','assets/../data.tsv']:
            response=self.client.get(f'/reports/{JID}/'+file.replace('/','%2F'),params={'view':grant})
            self.assertEqual(response.status_code,404)

    def test_deletion_and_retention_revoke_links_and_erase_alias_metadata(self):
        url=self.task()['report_url']
        for expires in [0,self.now-1]:
            self.store.update(JID,expires=expires)
            self.assertEqual(self.client.get(url).status_code,404)
            self.assertIsNone(self.link_record(url))
            job=self.task()
            for field in ['report_url','download_url','report_access_expires']:self.assertIsNone(job[field])

    def test_alias_expiry_never_exceeds_report_retention(self):
        self.store.update(JID,expires=int(self.now)+90)
        job=self.task()
        self.assertEqual(job['report_access_expires'],int(self.now)+90)
        page=self.client.get(job['report_url'])
        for child in re.findall(r'href="(/r/[A-Za-z0-9_-]{22})"',page.text):
            self.assertEqual(self.link_record(child)['expires'],int(self.now)+90)
        self.store.update(JID,status='expired',options='{}',filename='')
        self.assertIsNone(self.link_record(job['report_url']))

    def test_old_cookies_are_ignored_and_session_endpoint_is_removed(self):
        self.client.cookies.set('ont_report_'+JID,self.owner)
        self.assertEqual(self.client.get(self.page).status_code,404)
        response=self.client.post(f'/api/jobs/{JID}/session',json={},headers={'Authorization':'Bearer '+self.owner})
        self.assertEqual(response.status_code,404)
        self.assertNotIn('set-cookie',response.headers)

    def test_report_sandbox_and_privacy_headers_are_preserved(self):
        response=self.client.get(self.task()['report_url'])
        csp=response.headers['content-security-policy']
        self.assertIn('allow-popups',csp);self.assertNotIn('allow-same-origin',csp)
        self.assertIn('allow-modals',csp)
        self.assertIn('sandbox ',csp);self.assertIn("connect-src 'none'",csp)
        self.assertEqual(response.headers['referrer-policy'],'no-referrer')
        self.assertEqual(response.headers['cache-control'],'no-store, no-transform')
        self.assertIn('max-age=',response.headers['strict-transport-security'])
        self.assertNotIn('set-cookie',response.headers)

    def test_many_authenticated_report_images_are_not_rate_limited(self):
        self.settings.resource_request_burst=3
        self.settings.resource_requests_per_minute=1
        self.settings.global_resource_request_burst=4
        self.settings.global_resource_requests_per_minute=1
        page=self.client.get(self.task()['report_url'])
        image=re.search(r'href="(/r/[A-Za-z0-9_-]{22})"',page.text).group(1)
        legacy=f'/reports/{JID}/assets/figure.svg?view='+report_grant(self.settings.secret,JID,int(self.now)+100)
        for path in [image,legacy]*128:
            response=self.client.get(path,headers={'Sec-Fetch-Dest':'image'})
            self.assertEqual(response.status_code,200)
            self.assertEqual(response.headers['content-type'],'image/svg+xml')
        # A forged or expired capability still reaches ordinary access checks.
        self.assertEqual(self.client.get('/r/'+'x'*22,headers={'Sec-Fetch-Dest':'image'}).status_code,404)
        self.store.update(JID,expires=0)
        self.assertEqual(self.client.get(image,headers={'Sec-Fetch-Dest':'image'}).status_code,404)

    def test_image_destination_cannot_exempt_documents_downloads_or_invalid_links(self):
        self.settings.resource_request_burst=3
        self.settings.resource_requests_per_minute=1
        job=self.task();hint={'Sec-Fetch-Dest':'image'}
        self.assertEqual(self.client.get(job['report_url'],headers=hint).status_code,200)
        self.assertEqual(self.client.get(job['download_url'],headers=hint).status_code,200)
        self.assertEqual(self.client.get('/r/'+'x'*22,headers=hint).status_code,404)
        self.assertEqual(self.client.get(job['report_url'],headers=hint).status_code,429)
        self.assertEqual(self.client.get(job['download_url'],headers=hint).status_code,429)
        self.assertEqual(self.client.get('/r/'+'x'*22,headers=hint).status_code,429)
        self.assertEqual(self.client.get('/api/catalog').status_code,200)

    def test_bundle_buttons_are_short_and_preserve_offline_paths(self):
        job=self.task();grant=report_grant(self.settings.secret,JID,int(self.now)+100)
        for url in [self.page,f'/reports/{JID}/downloads.html',f'/reports/{JID}/index.html']:
            response=self.client.get(url,params={'view':grant});self.assertEqual(response.status_code,200)
            self.assertRegex(urlsplit(str(response.url)).path,r'^/r/[A-Za-z0-9_-]{22}$')
            self.assertIn('Download complete report (ZIP)',response.text)
            self.assertNotIn('?view=',response.text)
            self.assertNotIn('set-cookie',response.headers)
        with zipfile.ZipFile(io.BytesIO(self.client.get(job['download_url']).content)) as archive:
            for name in ['index.html','samples/example.html','assets/figure.svg','data.tsv','downloads.html']:
                self.assertIn('report/'+name,archive.namelist())
            sample=archive.read('report/samples/example.html').decode()
            self.assertIn('Browse report files',sample)
            self.assertNotIn('Adapter table',sample);self.assertNotIn('report.zip',sample);self.assertNotIn('?view=',sample)

class HttpsSettingsTests(unittest.TestCase):
    def test_https_enforcement_survives_legacy_config_migration(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'secrets.json').write_text(json.dumps({'signing_key':'CHANGE_ME'*8}))
            (root/'catalog.json').write_text('{"references":[]}')
            base={'state':str(root/'state'),'public_url':'http://testserver','allowed_hosts':['testserver'],
                  'secrets_file':str(root/'secrets.json'),'public_catalog_file':str(root/'catalog.json')}
            config=root/'web.yaml'
            with patch.dict(os.environ,PORTAL_CONFIG=str(config)):
                for extra in [{},{'secure_cookies':True},{'secure_cookies':False,'require_https':True}]:
                    with self.subTest(extra=extra):
                        config.write_text(json.dumps(base|extra))
                        with self.assertRaisesRegex(ValueError,'HTTPS is required'):Settings.load()
                for extra in [{'secure_cookies':False},{'require_https':False}]:
                    config.write_text(json.dumps(base|extra))
                    self.assertFalse(Settings.load().require_https)

class BedPresentationTests(unittest.TestCase):
    def test_filename_is_display_metadata_and_rejects_paths(self):
        base={'sample':'example','assay':'chromatin','reference':'dm6','custom_regions':True}
        self.assertEqual(Options(**base,regions_filename='HL280 peaks.bed').regions_filename,'HL280 peaks.bed')
        for name in ['../peaks.bed','folder/peaks.bed','folder\\peaks.bed','peaks\n.bed']:
            with self.subTest(name=name),self.assertRaises(ValidationError):Options(**base,regions_filename=name)

    def test_submitted_filename_survives_api_storage_and_pipeline_config(self):
        from ont_portal.catalog import job_config
        from test_custom_region_options import BASE, PUBLIC
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            settings=Settings(state=root/'state',public_url='https://testserver',allowed_hosts=['testserver'],
                              secret='test-only-'*8,submission_key='',public_catalog=PUBLIC,min_free_bytes=0)
            app=create_app(settings);app.state.store.heartbeat()
            name='HL280 peaks.bed'
            with TestClient(app,base_url='https://testserver',headers={'X-Portal-Request':'1'}) as client:
                response=client.post('/api/jobs',json={'email':'test@example.org','filename':'test.fastq','size':10,'consent':True,'data_policy_acknowledged':True,
                    'options':BASE | {'regions_filename':name,'regions_mode':'region','regions_assembly':'dm6'}})
                self.assertEqual(response.status_code,201,response.text)
                saved=Options.model_validate_json(app.state.store.get(response.json()['id'])['options'])
                self.assertEqual(saved.regions_filename,name)
                ref={key:str(root/key) for key in ('fasta','annotation','chrom_sizes','bowtie2_index','bwa_mem2_index')}
                ref.update(assembly='dm6',primary_chromosomes=['chr2L'])
                config,_=job_config(saved,{'references':{'dm6':ref}},root,2)
                self.assertEqual(config['samples'][saved.sample]['custom_regions']['bed_name'],name)

    def test_submitted_and_new_bed_are_distinct_and_names_are_escaped(self):
        payload={'reference':'dm6','sample':'example','scopes':{'total':{'name':'Total library','records':[]}}}
        initial="<p class='muted'>Uploaded regions</p><img class='custom-static' src='../assets/plot.svg'>"
        name='peaks <script>.bed'
        report=render_section(payload,initial,'region',name)
        self.assertIn('peaks &lt;script&gt;.bed',report);self.assertIn('Submitted BED',report)
        self.assertIn('Analyze another BED',report);self.assertIn('New BED result',report)
        self.assertIn("value='region' selected",report)
        self.assertEqual(upgrade_html(report,name),report)
        legacy=("<section class='report-section' id='custom-regions'><h2>Custom region metaplot</h2></div>"+initial+
                "<p class='muted'>Choose a BED file for dm6.</p><script id='custom-payload' type='application/json'>"+json.dumps(payload)+"</script><script>old();</script></section>")
        upgraded=upgrade_html(legacy)
        self.assertIn('original filename was not retained',upgraded);self.assertIn(initial,upgraded)
