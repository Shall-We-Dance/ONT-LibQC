"""Report navigation remains usable without granting access to task operations."""
import io
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
import zipfile
from fastapi.testclient import TestClient
from pydantic import ValidationError
from ont_portal.app import create_app
from ont_portal.report_entry import direct_report
from ont_portal.schema import Options
from ont_portal.settings import Settings
from ont_portal.store import report_grant, report_session_token, token_for, valid_report_grant
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

    def session(self):
        return self.client.post(f'/api/jobs/{JID}/session',json={},headers={'Authorization':'Bearer '+self.owner})

    def test_report_cookie_lasts_to_retention_and_cannot_manage_tasks(self):
        response=self.session();self.assertEqual(response.status_code,200)
        cookie=response.headers['set-cookie']
        for flag in ['HttpOnly','Secure','SameSite=lax',f'Path=/reports/{JID}/']:self.assertIn(flag,cookie)
        self.assertGreater(int(re.search(r'Max-Age=(\d+)',cookie)[1]),6*86400)
        credential=self.client.cookies.get('ont_report_'+JID)
        self.assertEqual(credential,report_session_token(self.settings.secret,JID));self.assertNotEqual(credential,self.owner)
        self.assertEqual(self.client.get(f'/api/jobs/{JID}').status_code,404)
        self.assertEqual(self.client.get(f'/api/jobs/{JID}',headers={'Authorization':'Bearer '+credential}).status_code,404)
        self.assertEqual(self.client.post(f'/api/jobs/{JID}/cancel',json={},headers={'Authorization':'Bearer '+credential}).status_code,404)

    def test_expired_resource_link_works_with_owner_report_cookie(self):
        self.session();old=report_grant(self.settings.secret,JID,int(self.now)+300)
        with patch('time.time',return_value=self.now+3700):
            page=self.client.get(self.page,params={'view':old});self.assertEqual(page.status_code,200)
            fresh=re.search(r'view=(\d+\.[a-f0-9]{64})',page.text)[1]
            self.assertTrue(valid_report_grant(self.settings.secret,JID,fresh));self.assertNotEqual(fresh,old)
            self.assertEqual(self.client.get(f'/reports/{JID}/assets/figure.svg',params={'view':old}).status_code,200)
            self.assertEqual(self.client.get(f'/reports/{JID}/report.zip',params={'view':old}).status_code,200)
            self.client.cookies.clear()
            self.assertEqual(self.client.get(self.page,params={'view':old}).status_code,404)

    def test_valid_view_grant_sets_only_report_cookie_and_refreshes_links(self):
        old=report_grant(self.settings.secret,JID,int(self.now)+100)
        response=self.client.get(self.page,params={'view':old})
        self.assertEqual(response.status_code,200)
        self.assertEqual(self.client.cookies.get('ont_report_'+JID),report_session_token(self.settings.secret,JID))
        fresh=re.search(r'view=(\d+\.[a-f0-9]{64})',response.text)[1]
        self.assertGreater(int(fresh.split('.')[0]),int(old.split('.')[0]))
        self.assertEqual(self.client.get(f'/api/jobs/{JID}').status_code,404)

    def test_browser_missing_access_gets_private_link_guidance(self):
        response=self.client.get(self.page,headers={'Accept':'text/html','Sec-Fetch-Dest':'document'},follow_redirects=False)
        self.assertEqual(response.status_code,303);self.assertEqual(response.headers['location'],f'/jobs/{JID}?report=1')
        self.assertEqual(self.client.get(self.page).status_code,404)
        self.assertEqual(self.client.get(f'/reports/{JID}/assets/figure.svg',headers={'Accept':'image/svg+xml','Sec-Fetch-Dest':'image'}).status_code,404)

    def test_expiry_and_other_report_credentials_still_block_access(self):
        self.client.cookies.set('ont_report_'+JID,report_session_token(self.settings.secret,'d'*32))
        self.assertEqual(self.client.get(self.page).status_code,404)
        self.client.cookies.clear();self.session();self.store.update(JID,expires=self.now-1)
        self.assertEqual(self.client.get(self.page).status_code,410)

    def test_existing_cookie_is_migrated_and_figures_can_open(self):
        self.client.cookies.set('ont_report_'+JID,self.owner)
        response=self.client.get(self.page);self.assertEqual(response.status_code,200)
        csp=response.headers['content-security-policy']
        self.assertIn('allow-popups',csp);self.assertNotIn('allow-same-origin',csp)
        self.assertIn('sandbox ',csp);self.assertIn('connect-src \'none\'',csp)
        self.assertIn(report_session_token(self.settings.secret,JID),response.headers['set-cookie'])

    def test_bundle_buttons_include_report_and_preserve_offline_paths(self):
        self.session()
        for url in [self.page,f'/reports/{JID}/downloads.html']:
            response=self.client.get(url);self.assertEqual(response.status_code,200)
            self.assertIn('Download complete report (ZIP)',response.text)
            self.assertIn('report.zip?view=',response.text)
        with zipfile.ZipFile(io.BytesIO(self.client.get(f'/reports/{JID}/report.zip').content)) as archive:
            for name in ['index.html','samples/example.html','assets/figure.svg','data.tsv','downloads.html']:
                self.assertIn('report/'+name,archive.namelist())
            sample=archive.read('report/samples/example.html').decode()
            self.assertIn('Browse report files',sample)
            self.assertNotIn('Adapter table',sample);self.assertNotIn('report.zip',sample);self.assertNotIn('?view=',sample)

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
