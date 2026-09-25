"""Distinguish manual deletion from retention after private details are erased."""
import json
import tempfile
import time
import unittest
from pathlib import Path

from fastapi.testclient import TestClient
from ont_portal.app import create_app
from ont_portal.frontend_agent import Frontend
from ont_portal.settings import Settings
from ont_portal.store import token_for

JID = 'e' * 32

class DeletionStatusTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(); self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.settings = Settings(state=self.root / 'state', public_url='https://testserver',
                                 allowed_hosts=['testserver'], secret='test-only-' * 8,
                                 submission_key='', public_catalog={}, min_free_bytes=0)
        self.app = create_app(self.settings); self.store = self.app.state.store
        now = time.time()
        with self.store.connect() as c:
            c.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires) VALUES (?,?,?,?,?,?,?,?,?,?)',
                      (JID, now, now, 'succeeded', 'test@example.org', 'sample.fastq', 10,
                       json.dumps({'sample': 'sample', 'assay': 'chromatin', 'reference': 'dm6'}), 'test', now + 7 * 86400))
        report = self.settings.state / 'reports' / JID
        report.mkdir(); (report / 'report.zip').write_bytes(b'test')
        self.frontend = Frontend({'state': str(self.settings.state), 'retention_days': 7})
        self.client = TestClient(self.app, base_url='https://testserver', headers={'X-Portal-Request': '1'})
        self.addCleanup(self.client.close)
        self.auth = {'Authorization': 'Bearer ' + token_for(self.settings.secret, JID)}
        self.url = '/api/jobs/' + JID

    def get(self):
        response = self.client.get(self.url, headers=self.auth)
        self.assertEqual(response.status_code, 200, response.text)
        return response.json()

    def test_manual_delete_marker_survives_cleanup_and_privacy_erasure(self):
        self.assertIsNone(self.get()['deletion_reason'])
        self.assertFalse(self.get()['deletion_requested'])
        response = self.client.delete(self.url, headers=self.auth)
        self.assertEqual(response.status_code, 200)
        self.assertTrue(self.get()['deletion_requested'])
        self.assertIsNone(self.get()['deletion_reason'])
        self.frontend.maintain()
        result = self.get()
        self.assertEqual(result['status'], 'expired')
        self.assertEqual(result['deletion_reason'], 'manual')
        self.assertFalse(result['deletion_requested'])
        self.assertEqual(result['options'], {})
        self.assertEqual(result['filename'], '')
        self.assertFalse((self.settings.state / 'reports' / JID).exists())

    def test_natural_expiry_has_retention_reason(self):
        self.store.update(JID, expires=time.time() - 1)
        self.frontend.maintain()
        result = self.get()
        self.assertEqual(result['deletion_reason'], 'retention')
        self.assertEqual(result['options'], {})
        self.assertIsNone(result['report_url'])
        self.assertIsNone(result['download_url'])
        self.assertEqual(self.client.get('/api/catalog').json()['retention_days'], 7)

    def test_repeated_delete_cannot_relabel_natural_expiry(self):
        expired_at = time.time() - 1
        self.store.update(JID, expires=expired_at)
        self.frontend.maintain()
        response = self.client.delete(self.url, headers=self.auth)
        self.assertEqual(response.json()['status'], 'expired')
        self.assertEqual(self.get()['deletion_reason'], 'retention')
        self.assertEqual(self.store.get(JID)['expires'], expired_at)

    def test_delete_after_deadline_preserves_retention_before_cleanup(self):
        expired_at = time.time() - 1
        self.store.update(JID, expires=expired_at)
        self.assertEqual(self.client.delete(self.url, headers=self.auth).status_code, 200)
        self.assertEqual(self.store.get(JID)['expires'], expired_at)
        self.frontend.maintain()
        self.assertEqual(self.get()['deletion_reason'], 'retention')

    def test_existing_cleaned_records_retain_their_reason(self):
        for expires, expected in [(0, 'manual'), (time.time() - 1, 'retention')]:
            with self.subTest(reason=expected):
                self.store.update(JID, status='expired', expires=expires, options='{}', filename='', email='')
                self.assertEqual(self.get()['deletion_reason'], expected)

    def test_private_link_and_active_task_protections_remain(self):
        self.assertEqual(self.client.get(self.url).status_code, 404)
        self.assertEqual(self.client.delete(self.url).status_code, 404)
        self.store.update(JID, status='running')
        self.assertEqual(self.client.delete(self.url, headers=self.auth).status_code, 409)
        self.assertGreater(self.store.get(JID)['expires'], time.time())

if __name__ == '__main__': unittest.main()
