"""Optional email must preserve private access, admission limits and mail opt-in."""
import tempfile
import time
import unittest
from pathlib import Path

from fastapi.testclient import TestClient
from pydantic import ValidationError
from ont_portal.app import create_app
from ont_portal.schema import Submission
from ont_portal.settings import Settings
from ont_portal.store import client_key

BASE={'filename':'sample.fastq','size':20,'consent':True,'data_policy_acknowledged':True,
      'options':{'sample':'sample','assay':'chromatin','reference':'dm6'}}
CATALOG={'references':[{'id':'dm6','ready':True,'mappers':['bowtie2'],'mycoplasma':False,'spike_in':False}]}


class OptionalEmailTests(unittest.TestCase):
    def setUp(self):
        temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup)
        self.settings=Settings(state=Path(temp.name)/'state',public_url='https://testserver',
            allowed_hosts=['testserver'],secret='test-only-'*8,submission_key='',public_catalog=CATALOG,
            min_free_bytes=0,mail_enabled=True,max_files=2,max_pending=8,max_pending_per_client=2,max_jobs_per_day=4)
        self.app=create_app(self.settings);self.store=self.app.state.store;self.store.heartbeat()

    def client(self,number=1):
        client=TestClient(self.app,base_url='https://testserver',client=(f'198.51.100.{number}',50000),headers={'X-Portal-Request':'1'})
        self.addCleanup(client.close);return client

    def batch(self,client,email=None):
        submission=BASE if email is None else BASE|{'email':email}
        return client.post('/api/batches',json={'files':[submission,submission]})

    def expire(self,jobs):
        for job in jobs:
            self.store.update(job['id'],status='expired',email='',client='',options='{}',expires=0)

    def test_omitted_null_empty_and_whitespace_are_equivalent(self):
        for extra in [{},{'email':None},{'email':''},{'email':' \t\r\n '}]:
            with self.subTest(extra=extra):
                self.assertEqual(Submission(**(BASE|extra)).email,'')

    def test_nonempty_email_stays_validated_and_normalized(self):
        self.assertEqual(Submission(**(BASE|{'email':' user@EXAMPLE.org '})).email,'user@example.org')
        for address in ['invalid','user@','a@example.org,b@example.org','a@example.org\r\nBcc:b@example.org',123,{},[],False]:
            with self.subTest(address=address),self.assertRaises(ValidationError):Submission(**(BASE|{'email':address}))

    def test_single_submit_without_email_has_private_access_and_no_notification(self):
        for i,extra in enumerate([{},{'email':None},{'email':''},{'email':'   '}],1):
            client=self.client(i);response=client.post('/api/jobs',json=BASE|extra)
            self.assertEqual(response.status_code,201,response.text)
            job=response.json();saved=self.store.get(job['id'])
            self.assertEqual((saved['email'],saved['mail_status'],saved['mail_attempts']),('','disabled',0))
            self.assertEqual(job['mail_status'],'disabled')
            url='/api/jobs/'+job['id']
            self.assertEqual(client.get(url).status_code,404)
            status=client.get(url,headers={'Authorization':'Bearer '+job['token']})
            self.assertEqual(status.status_code,200)
            self.assertEqual(status.json()['mail_status'],'disabled')
            self.assertNotIn('set-cookie',status.headers)

    def test_empty_batch_uses_only_its_client_quota(self):
        client=self.client()
        response=client.post('/api/batches',json={'files':[BASE,BASE|{'email':None}]})
        self.assertEqual(response.status_code,201,response.text)
        jobs=response.json()['jobs']
        self.assertEqual(len({self.store.get(job['id'])['batch_id'] for job in jobs}),1)
        with self.store.connect() as db:
            events=[tuple(row) for row in db.execute('SELECT subject,amount FROM admission_events')]
        self.assertEqual(events,[(client_key(self.settings.secret,'198.51.100.1'),2)])

    def test_unrelated_clients_without_email_do_not_share_a_quota(self):
        # A legacy empty-email admission record must not affect new anonymous users.
        with self.store.connect() as db:
            db.execute('INSERT INTO admission_events VALUES (?,?,?)',(client_key(self.settings.secret,'email:'),time.time(),100))
        for i in range(1,4):
            response=self.batch(self.client(i))
            self.assertEqual(response.status_code,201,response.text)
        with self.store.connect() as db:
            self.assertEqual(db.execute('SELECT COUNT(*) FROM jobs').fetchone()[0],6)

    def test_client_active_limit_applies_with_or_without_email(self):
        client=self.client();self.assertEqual(self.batch(client).status_code,201)
        for extra in [{},{'email':'user@example.org'}]:
            response=client.post('/api/jobs',json=BASE|extra)
            self.assertEqual(response.status_code,429,response.text)
            self.assertIn('active task limit',response.json()['detail'])

    def test_daily_limit_survives_private_metadata_cleanup_without_email(self):
        client=self.client()
        for _ in range(2):
            response=self.batch(client);self.assertEqual(response.status_code,201,response.text)
            self.expire(response.json()['jobs'])
        response=client.post('/api/jobs',json=BASE)
        self.assertEqual(response.status_code,429,response.text)
        self.assertIn('Daily submission limit',response.json()['detail'])

    def test_provided_email_keeps_active_limit_across_clients_and_case(self):
        self.assertEqual(self.batch(self.client(1),'user@example.org').status_code,201)
        response=self.client(2).post('/api/jobs',json=BASE|{'email':'User@example.org'})
        self.assertEqual(response.status_code,429,response.text)
        self.assertIn('active task limit',response.json()['detail'])

    def test_provided_email_keeps_daily_limit_across_clients_after_cleanup(self):
        for i in [1,2]:
            response=self.batch(self.client(i),'user@example.org')
            self.assertEqual(response.status_code,201,response.text);self.expire(response.json()['jobs'])
        response=self.client(3).post('/api/jobs',json=BASE|{'email':'User@example.org'})
        self.assertEqual(response.status_code,429,response.text)
        self.assertIn('Daily submission limit',response.json()['detail'])

    def test_global_queue_limit_still_applies_without_email(self):
        for i in range(1,5):self.assertEqual(self.batch(self.client(i)).status_code,201)
        response=self.client(5).post('/api/jobs',json=BASE)
        self.assertEqual(response.status_code,503,response.text)
        self.assertIn('queue is full',response.json()['detail'])

    def test_mixed_recipients_are_rejected_atomically(self):
        client=self.client()
        for email in ['','other@example.org']:
            response=client.post('/api/batches',json={'files':[BASE|{'email':'user@example.org'},BASE|{'email':email}]})
            self.assertEqual(response.status_code,422,response.text)
        with self.store.connect() as db:
            self.assertEqual(db.execute('SELECT COUNT(*) FROM jobs').fetchone()[0],0)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM admission_events').fetchone()[0],0)
        self.assertEqual(list((self.settings.state/'uploads').iterdir()),[])

    def test_invalid_email_is_rejected_before_any_job_or_quota_is_reserved(self):
        response=self.client().post('/api/jobs',json=BASE|{'email':'not-an-address'})
        self.assertEqual(response.status_code,422)
        with self.store.connect() as db:
            self.assertEqual(db.execute('SELECT COUNT(*) FROM jobs').fetchone()[0],0)
            self.assertEqual(db.execute('SELECT COUNT(*) FROM admission_events').fetchone()[0],0)

    def test_notification_requires_both_email_and_enabled_delivery(self):
        client=self.client()
        first=client.post('/api/jobs',json=BASE|{'email':'user@example.org'})
        self.assertEqual(first.status_code,201);self.assertEqual(first.json()['mail_status'],'pending')
        self.settings.mail_enabled=False
        second=client.post('/api/jobs',json=BASE|{'email':'user@example.org'})
        self.assertEqual(second.status_code,201);self.assertEqual(second.json()['mail_status'],'disabled')


if __name__=='__main__':unittest.main()
