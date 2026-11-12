"""Private links remain compatible; public service history contains only counts."""
import json
import sqlite3
import tempfile
import time
import unittest
from pathlib import Path
from fastapi.testclient import TestClient
from ont_portal.app import create_app
from ont_portal.settings import Settings
from ont_portal.store import Store, token_for
from ont_portal.short_links import job_link, batch_link, key_hash
from ont_portal.mail_template import build_message

CATALOG={'references':[{'id':'test','ready':True,'mappers':['bowtie2'],'mycoplasma':False,'spike_in':False}]}
PAYLOAD={'filename':'private_sample.fastq','size':20,'consent':True,'data_policy_acknowledged':True,
         'email':'person@example.org','options':{'sample':'private_sample','assay':'chromatin','reference':'test'}}

class StatusLinksTests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        self.settings=Settings(state=Path(tmp.name)/'state',public_url='https://testserver',allowed_hosts=['testserver'],
            secret='test-only-'*8,submission_key='',public_catalog=CATALOG,min_free_bytes=0)
        self.app=create_app(self.settings);self.store=self.app.state.store;self.store.heartbeat()
        self.client=TestClient(self.app,base_url='https://testserver',headers={'X-Portal-Request':'1'})
        self.addCleanup(self.client.close)

    def create(self,n=1):
        response=self.client.post('/api/batches',json={'files':[PAYLOAD for _ in range(n)]})
        self.assertEqual(response.status_code,201,response.text)
        return response.json()['jobs']

    def resolve(self,url):
        return self.client.post('/api/links/resolve',json={'key':url.split('#')[1]})

    def test_migration_preserves_old_job_details_and_leaves_old_keys_unassigned(self):
        old=self.store.state/'legacy';(old/'control').mkdir(parents=True)
        db=sqlite3.connect(old/'control/jobs.sqlite3')
        db.execute("CREATE TABLE jobs (id TEXT PRIMARY KEY,created REAL,updated REAL,status TEXT,email TEXT,filename TEXT,expected_bytes INTEGER,options TEXT,client TEXT,expires REAL)")
        db.execute("INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?,?,?)",('a'*32,1,2,'succeeded','legacy@example.org','legacy.fastq',20,'{}','old-client',99))
        db.commit();db.close()
        migrated=Store(old);row=migrated.get('a'*32)
        self.assertEqual((row['email'],row['filename'],row['expires']),('legacy@example.org','legacy.fastq',99))
        self.assertEqual(row['outcome'],'succeeded');self.assertIsNone(row['short_link_hash'])
        self.assertEqual(migrated.statistics()['succeeded'],1)
        self.assertEqual(job_link(self.settings.secret,row),'/jobs/'+'a'*32+'#'+token_for(self.settings.secret,'a'*32))

    def test_new_single_link_is_short_private_and_isolated(self):
        one,two=self.create(2)
        self.assertRegex(one['private_link'],r'^/j/#[A-Za-z0-9_-]{22}$')
        self.assertNotEqual(one['private_link'],two['private_link'])
        result=self.resolve(one['private_link']);self.assertEqual(result.status_code,200)
        self.assertEqual([row['id'] for row in result.json()['jobs']],[one['id']])
        key=result.json()['jobs'][0]['token'];self.assertEqual(key,one['token'])
        url='/api/jobs/'+one['id']
        self.assertEqual(self.client.get(url).status_code,404)
        status=self.client.get(url,headers={'Authorization':'Bearer '+key})
        self.assertEqual(status.status_code,200)
        self.assertEqual(status.json()['private_link'],one['private_link'])
        stored=self.store.get(one['id']);self.assertEqual(stored['short_link_hash'],key_hash(one['private_link'].split('#')[1]))
        self.assertNotIn(one['private_link'].split('#')[1],str(stored))
        self.assertEqual(self.resolve('/j/#'+'z'*22).status_code,404)
        self.assertEqual(self.resolve('/j/#'+one['id']).status_code,404)

    def test_batch_capability_returns_only_its_members(self):
        batch=self.create(2);other=self.create()[0]
        self.assertRegex(batch[0]['batch_link'],r'^/b/#[A-Za-z0-9_-]{22}$')
        result=self.resolve(batch[0]['batch_link']).json()['jobs']
        self.assertEqual({row['id'] for row in result},{row['id'] for row in batch})
        self.assertNotIn(other['id'],str(result))
        with self.store.connect() as db:
            self.assertNotIn(batch[0]['batch_link'].split('#')[1],str([tuple(r) for r in db.execute('SELECT * FROM batch_links')]))

    def test_old_long_links_and_credentials_remain_unchanged(self):
        job=self.create()[0]
        with self.store.connect() as db:db.execute('UPDATE jobs SET short_link_hash=NULL WHERE id=?',(job['id'],))
        row=self.store.get(job['id']);old='/jobs/'+job['id']+'#'+token_for(self.settings.secret,job['id'])
        self.assertEqual(job_link(self.settings.secret,row),old)
        self.assertEqual(self.client.get(old).status_code,200)
        status=self.client.get('/api/jobs/'+job['id'],headers={'Authorization':'Bearer '+job['token']})
        self.assertEqual(status.status_code,200)
        self.assertEqual(status.json()['private_link'],old)
        self.assertEqual(self.resolve(job['private_link']).status_code,404)
        # Reopening/migrating the database never assigns new keys to old rows.
        Store(self.store.state)
        self.assertIsNone(self.store.get(job['id'])['short_link_hash'])

    def test_public_statistics_are_counts_only_and_survive_cleanup(self):
        jobs=self.create(5)
        for job,status in zip(jobs,('running','queued','succeeded','failed','cancelled')):
            self.store.update(job['id'],status=status)
        stats=self.client.get('/api/statistics')
        self.assertEqual(stats.json(),dict(total=5,running=1,queued=1,succeeded=1,failed=1,cancelled=1,archived=0))
        self.assertTrue(all(type(value) is int for value in stats.json().values()))
        for marker in (jobs[0]['id'],jobs[0]['token'],'private_sample','person@example.org','198.51'):
            self.assertNotIn(marker,stats.text)
        self.store.update(jobs[2]['id'],status='expired',options='{}',email='',filename='',client='',expires=0)
        self.assertEqual(self.client.get('/api/statistics').json(),stats.json())
        self.store.update(jobs[1]['id'],status='expired',options='{}',email='',filename='',client='')
        updated=self.client.get('/api/statistics').json()
        self.assertEqual((updated['archived'],updated['queued'],updated['succeeded']),(1,0,1))
        Store(self.store.state)
        self.assertEqual(self.store.statistics(),updated)

    def test_new_and_old_task_shells_keep_private_headers(self):
        for url in ('/j/','/b/','/jobs/'+'a'*32):
            response=self.client.get(url)
            self.assertEqual(response.status_code,200)
            self.assertIn('no-transform',response.headers['cache-control'])
            self.assertNotIn('cloudflare',response.headers['content-security-policy'])
            self.assertEqual(response.headers['referrer-policy'],'no-referrer')
            self.assertNotIn('set-cookie',response.headers)
        self.assertEqual(self.client.get('/api/links/resolve').status_code,405)
        self.assertEqual(self.client.post('/api/links/resolve',content='{',headers={'Content-Type':'application/json'}).status_code,404)
        self.assertEqual(self.client.post('/api/links/resolve',json={'key':None}).status_code,404)
        self.assertEqual(self.client.post('/api/links/resolve',json={'key':'x'*22},headers={'Origin':'https://other.example'}).status_code,403)

    def test_retention_defaults_to_thirty_days(self):
        before=time.time();job=self.create()[0];row=self.store.get(job['id'])
        self.assertEqual(self.client.get('/api/catalog').json()['retention_days'],30)
        self.assertGreaterEqual(row['expires'],before+30*86400)
        self.assertLess(row['expires'],time.time()+30*86400+1)

    def test_mail_uses_short_links_without_changing_legacy_rows(self):
        jobs=self.create(2);rows=[self.store.get(job['id']) for job in jobs]
        config={'public_url':'https://testserver'}
        smtp={'from_address':'reports@example.org'}
        text=build_message(config,smtp,rows[0],self.settings.secret,rows).get_body(preferencelist=('plain',)).get_content()
        self.assertIn('https://testserver'+jobs[0]['batch_link'],text)
        for job in jobs:self.assertIn('https://testserver'+job['private_link'],text)
        self.assertNotIn('/#jobs=',text)
        for row in rows:row['short_link_hash']=None
        legacy=build_message(config,smtp,rows[0],self.settings.secret,rows).get_body(preferencelist=('plain',)).get_content()
        self.assertIn('/#jobs=',legacy)
        for row in rows:self.assertIn('/jobs/'+row['id']+'#'+token_for(self.settings.secret,row['id']),legacy)

if __name__=='__main__':unittest.main()
