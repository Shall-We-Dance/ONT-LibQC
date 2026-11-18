import io
import json
import sys
import tempfile
import time
import unittest
import zipfile
from pathlib import Path

from fastapi.testclient import TestClient
from ont_portal.app import create_app
from ont_portal.settings import Settings
from ont_portal.store import token_for
from ont_portal.worker import publish_report

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow/scripts'))
from report_render import identity_panel_html


def write_report(root):
    (root/'samples').mkdir(parents=True)
    (root/'data/metrics').mkdir(parents=True)
    (root/'index.html').write_text('<a href="samples/sample.html">Report</a>')
    (root/'downloads.html').write_text('<main>Files</main>')
    scopes = [dict(scope_id='lib01', scope_name='lib01', i5='', i7='ACGT'),
              dict(scope_id='lib02', scope_name='B', i5=None, i7='TGCA')]
    pairs = [dict(scope, name=scope['scope_name']) for scope in scopes]
    metric = dict(sub_libraries=scopes, assay='rnaseq')
    panel = identity_panel_html(metric, dict(mode='multiplexed_library'), pairs)
    body = ('<!doctype html><html><head><title>Sample</title></head><body><main id="library-shape">' + panel +
            '<div class="table-scroll"><table><tr><th>lib01</th><th>B</th></tr></table></div></main></body></html>')
    for sample in ('sample', 'other'):
        (root/'data/metrics'/(sample+'.json')).write_text(json.dumps(metric))
        (root/'samples'/(sample+'.html')).write_text(body if sample == 'sample' else
            body.replace(' data-sublibrary="lib01"', '').replace(' data-sublibrary="lib02"', ''))


class NameTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        self.cfg = Settings(state=root/'state', public_url='https://testserver', allowed_hosts=['testserver'],
                            secret='test-only-' * 8, submission_key='', public_catalog={})
        app = create_app(self.cfg)
        self.store = app.state.store
        self.client = TestClient(app, base_url=self.cfg.public_url, headers={'X-Portal-Request': '1'})
        self.jid = '1' * 32
        now = time.time()
        src = root/'input'
        write_report(src)
        for jid in (self.jid, '2' * 32):
            with self.store.connect() as db:
                db.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires) VALUES(?,?,?,?,?,?,?,?,?,?)',
                    (jid, now, now, 'succeeded', '', 'sample.fastq', 100, '{"sample":"sample"}', 'local', now + 86400))
            publish_report(src, self.cfg.state/'reports'/jid)
        self.job = self.task(self.jid)
        self.api = self.job['report_url'].replace('/r/', '/api/sublibrary-names/')
        self.auth = {'Authorization': 'Bearer ' + self.job['name_edit_token']}

    def task(self, jid):
        return self.client.get('/api/jobs/' + jid, headers={'Authorization': 'Bearer ' + token_for(self.cfg.secret, jid)}).json()

    def save(self, names, revision=0, headers=None):
        return self.client.put(self.api, headers=self.auth if headers is None else headers,
                               json={'names': names, 'revision': revision})

    def test_save_reload_archive_and_clear(self):
        self.assertEqual(self.client.get(self.api, headers=self.auth).json()['revision'], 0)
        name = '对照组 1'
        self.assertEqual(self.save({'lib01': ' ' + name + ' '}).status_code, 200)
        response = self.client.get(self.job['report_url'])
        self.assertIn(name, response.text)
        self.assertIn('Edit names', response.text)
        self.assertIn('sandbox', response.headers['Content-Security-Policy'])
        self.assertNotIn(self.job['name_edit_token'], response.text)
        with zipfile.ZipFile(io.BytesIO(self.client.get(self.job['download_url']).content)) as z:
            page = z.read('report/samples/sample.html').decode()
            self.assertIn(name, page)
            self.assertNotIn('Edit names', page)
            self.assertNotIn(name, z.read('report/samples/other.html').decode())
            self.assertNotIn(name, z.read('report/data/metrics/sample.json').decode())
            self.assertIn(f'sample\tlib01\t{name}\t\tACGT\n', z.read('report/data/sublibrary_names.tsv').decode())
            self.assertIsNone(z.testzip())
        self.assertEqual(self.save({'lib02': 'Changed'}).status_code, 409)
        self.assertEqual(self.client.get(self.api, headers=self.auth).json()['names'], {'lib01': name})
        self.assertEqual(self.save({'lib01': ''}, revision=1).status_code, 200)
        self.assertNotIn(name, self.client.get(self.job['report_url']).text)

    def test_same_scope_in_two_samples_has_independent_names(self):
        owner = {'Authorization': 'Bearer ' + token_for(self.cfg.secret, self.jid)}
        second = self.client.get('/reports/' + self.jid + '/samples/other.html', headers=owner,
                                 follow_redirects=False).headers['location']
        self.assertEqual(self.save({'lib01': 'Control'}).status_code, 200)
        response = self.client.put(second.replace('/r/', '/api/sublibrary-names/'), headers=self.auth,
                                   json={'revision': 0, 'names': {'lib01': 'Treatment'}})
        self.assertEqual(response.status_code, 200)
        with zipfile.ZipFile(io.BytesIO(self.client.get(self.job['download_url']).content)) as z:
            mapping = z.read('report/data/sublibrary_names.tsv').decode()
            self.assertIn('other\tlib01\tTreatment\t', mapping)
            self.assertIn('sample\tlib01\tControl\t', mapping)
            self.assertNotIn('Treatment', z.read('report/samples/sample.html').decode())
            self.assertNotIn('Control', z.read('report/samples/other.html').decode())

    def test_permissions_validation_and_retention(self):
        self.assertEqual(self.save({'lib01':'A'}, headers={}).status_code, 403)
        other = self.task('2' * 32)
        self.assertEqual(self.save({'lib01':'A'}, headers={'Authorization':'Bearer ' + other['name_edit_token']}).status_code, 403)
        self.assertEqual(self.client.get('/api/jobs/' + self.jid, headers=self.auth).status_code, 404)
        self.assertEqual(self.save({'lib01':'A'}, headers=self.auth | {'Origin':'null'}).status_code, 403)
        for values in ({'unknown':'A'}, {'lib01':'A\tB'}, {'lib02':'=1'}, {'lib01':'a' * 49},
                       {'lib01':'A<bad>'}, {'lib01':'A\u202eZ'}, {'lib01':34}):
            self.assertEqual(self.save(values).status_code, 422)
        self.assertEqual(self.save({'lib01':'Rep-1'}).status_code, 200)
        self.store.update(self.jid, expires=0)
        self.assertEqual(self.client.get(self.api, headers=self.auth).status_code, 404)
        with self.store.connect() as db:
            self.assertIsNone(db.execute('SELECT 1 FROM sublibrary_names WHERE job_id=?', (self.jid,)).fetchone())

    def test_older_page_and_single_library(self):
        from ont_portal.sublibrary_names import annotate_page, scope_catalog
        root = self.cfg.state/'reports'/self.jid
        old = (root/'samples/other.html').read_text()
        self.assertIn('Alias A', annotate_page(old, scope_catalog(root, 'other'), {'lib01': 'Alias A'}))
        (root/'data/metrics/sample.json').write_text('{"sub_libraries":[]}')
        self.assertNotIn('Edit names', self.client.get(self.job['report_url']).text)
        self.assertEqual(self.save({'lib01':'A'}).status_code, 404)
