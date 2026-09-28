"""Batch delivery regression tests. SMTP is mocked; no mail is sent."""
import json
import smtplib
import sqlite3
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch
from fastapi.testclient import TestClient
from ont_portal.app import create_app
from ont_portal.notifications import notify
from ont_portal.settings import Settings
from ont_portal.store import Store, token_for

class BatchNotificationsTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.store=Store(self.root/'state')
        (self.root/'smtp.yaml').write_text('enabled: true\nfrom_address: reports@example.org\nprovider: resend\n')
        (self.root/'secrets.json').write_text(json.dumps({'signing_key':'CHANGE_ME'}))
        self.config=dict(state=str(self.store.state),mail_enabled=True,smtp_file=str(self.root/'smtp.yaml'),secrets_file=str(self.root/'secrets.json'),public_url='https://example.org')
        self.sender=patch('ont_portal.notifications.send');self.send=self.sender.start();self.addCleanup(self.sender.stop)

    def job(self,jid,status='succeeded',batch='batch-one',email='user@example.org',zip_ready=True):
        now=time.time()
        with self.store.connect() as c:
            c.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires,batch_id) VALUES (?,?,?,?,?,?,?,?,?,?,?)',(jid,now,now,status,email,jid+'.fastq',10,json.dumps(dict(sample=jid,assay='chromatin',reference='dm6')),'test',now+3600,batch))
        if status=='succeeded' and zip_ready:
            folder=self.store.state/'reports'/jid;folder.mkdir();(folder/'report.zip').touch()

    def test_waits_for_every_file_and_published_report_then_sends_once(self):
        self.job('a');self.job('b','running');self.job('c','created')
        notify(self.config);self.send.assert_not_called()
        self.store.update('b',status='failed');self.store.update('c',status='succeeded')
        notify(self.config);self.send.assert_not_called()
        folder=self.store.state/'reports/c';folder.mkdir();(folder/'report.zip').touch()
        notify(self.config);self.send.assert_called_once()
        message=self.send.call_args.args[1];text=message.get_body(preferencelist=('plain',)).get_content()
        self.assertIn('/#jobs=',text)
        for jid in 'abc':
            self.assertIn(jid+':'+token_for('CHANGE_ME',jid),text)
            self.assertEqual(self.store.get(jid)['mail_status'],'sent')
        self.assertIn('Failed',text)
        notify(self.config);self.send.assert_called_once()

    def test_cancelled_and_expired_members_do_not_block(self):
        self.job('a');self.job('b','cancelled');self.job('c','expired',email='')
        self.store.update('c',mail_status='expired')
        notify(self.config);self.send.assert_called_once()
        text=self.send.call_args.args[1].get_body(preferencelist=('plain',)).get_content()
        self.assertIn('Cancelled',text);self.assertIn('Expired',text)

    def test_same_email_separate_batches_and_legacy_remain_separate(self):
        self.job('a',batch='one');self.job('b',batch='two');self.job('c',batch=None)
        notify(self.config);self.assertEqual(self.send.call_count,3)

    def test_retry_is_one_identical_message_for_batch(self):
        self.job('a');self.job('b');self.send.side_effect=smtplib.SMTPServerDisconnected('test')
        notify(self.config)
        before=self.send.call_args.args[1].as_bytes()
        for jid in 'ab':self.store.update(jid,mail_updated=0)
        self.send.side_effect=None;notify(self.config)
        self.assertEqual(before,self.send.call_args.args[1].as_bytes())
        for jid in 'ab':self.assertEqual(self.store.get(jid)['mail_attempts'],2)
        notify(self.config);self.assertEqual(self.send.call_count,2)

    def test_permanent_failure_does_not_retry(self):
        self.job('a');self.job('b');self.send.side_effect=smtplib.SMTPRecipientsRefused({})
        notify(self.config)
        for jid in 'ab':self.assertEqual(self.store.get(jid)['mail_status'],'failed')
        notify(self.config);self.send.assert_called_once()

    def test_no_starvation_by_unfinished_batches(self):
        for i in range(7):self.job(f'a{i}',batch=f'wait{i}');self.job(f'b{i}','running',batch=f'wait{i}')
        self.job('ready',batch='ready');notify(self.config,limit=1);self.send.assert_called_once()

    def test_retry_remains_identical_when_first_member_expires(self):
        self.job('a');self.job('b');self.send.side_effect=smtplib.SMTPServerDisconnected('test')
        notify(self.config);before=self.send.call_args.args[1].as_bytes()
        self.store.update('a',status='expired',mail_status='expired',email='',options='{}')
        self.store.update('b',mail_updated=0);self.send.side_effect=None
        notify(self.config)
        self.assertEqual(before,self.send.call_args.args[1].as_bytes())
        self.assertEqual(self.store.get('b')['mail_status'],'sent')

    def test_old_ambiguous_attempt_stops_without_sending(self):
        self.job('a');self.job('b');self.send.side_effect=smtplib.SMTPServerDisconnected('test')
        notify(self.config)
        with self.store.connect() as c:
            c.execute('UPDATE notification_attempts SET first_attempt=?',(time.time()-24*3600,))
            c.execute('UPDATE jobs SET mail_updated=0')
        notify(self.config);self.send.assert_called_once()
        for jid in 'ab':self.assertEqual(self.store.get(jid)['mail_status'],'failed')

    def test_legacy_database_migration_preserves_existing_jobs(self):
        self.job('old',batch=None)
        with self.store.connect() as c:
            c.execute('DROP INDEX jobs_batch');c.execute('ALTER TABLE jobs DROP COLUMN batch_id')
        migrated=Store(self.store.state)
        self.assertIsNone(migrated.get('old')['batch_id'])
        self.assertEqual(migrated.get('old')['filename'],'old.fastq')

    def test_batch_api_persists_distinct_groups_and_requires_choices(self):
        catalog={'references':[dict(id='dm6',ready=True,mappers=['bowtie2'],mycoplasma=False,spike_in=False)]}
        settings=Settings(state=self.store.state,public_url='https://testserver',allowed_hosts=['testserver'],secret='test-only-'*8,submission_key='',public_catalog=catalog,min_free_bytes=0,mail_enabled=True)
        app=create_app(settings);app.state.store.heartbeat()
        submission=dict(email='test@example.org',filename='test.fastq',size=10,options=dict(sample='test',reference='dm6',assay='chromatin'),consent=True,data_policy_acknowledged=True)
        with TestClient(app,base_url='https://testserver',headers={'X-Portal-Request':'1'}) as client:
            result=client.post('/api/batches',json={'files':[submission,submission]})
            self.assertEqual(result.status_code,201,result.text)
            ids=[j['id'] for j in result.json()['jobs']]
            batch=self.store.get(ids[0])['batch_id'];self.assertTrue(batch)
            self.assertEqual(batch,self.store.get(ids[1])['batch_id'])
            one=client.post('/api/jobs',json=submission);self.assertEqual(one.status_code,201,one.text)
            self.assertNotEqual(batch,self.store.get(one.json()['id'])['batch_id'])
            for field in ('reference','assay'):
                invalid=json.loads(json.dumps(submission));invalid['options'][field]=''
                self.assertEqual(client.post('/api/batches',json={'files':[invalid]}).status_code,422)

if __name__=='__main__':unittest.main()
