"""Custom BED reports must survive publication, transfer and private download.

Run: PYTHONPATH=portal python -m unittest discover -s portal/tests -v
Requires portal/requirements.txt and httpx for the HTTP tests.
"""
import json
from pathlib import Path
import shutil
import tempfile
import time
import unittest
import zipfile

from fastapi.testclient import TestClient
from ont_portal.app import create_app
from ont_portal.frontend_agent import Frontend
from ont_portal.settings import Settings
from ont_portal.store import Store, report_grant, token_for
from ont_portal.worker import publish_report

SAMPLE = 'DB6TSF_1_CR_HL280'
JOB_ID = 'a' * 32
PDF_NAME = f'assets/{SAMPLE}.custom_regions.pdf'
PDF = b'%PDF-1.4\n% custom region metaplot fixture\n%%EOF\n'


def create_job(store, status='running'):
    now = time.time()
    with store.connect() as connection:
        connection.execute(
            'INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,'
            'options,client,expires,bed_ready) VALUES (?,?,?,?,?,?,?,?,?,?,?)',
            (JOB_ID, now, now, status, '', SAMPLE + '.fastq', 100,
             json.dumps({'sample': SAMPLE, 'assay': 'chromatin', 'reference': 'dm6',
                         'custom_regions': True}), 'test', now + 3600, 1))


def archive_report(source, archive):
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as bundle:
        for path in source.rglob('*'):
            if path.is_file():
                bundle.write(path, 'report/' + path.relative_to(source).as_posix())


class ReportPublicationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / 'source'
        assets = f'{SAMPLE}.custom_regions'
        files = {
            'index.html': f'<a href="samples/{SAMPLE}.html">Sample</a>'.encode(),
            f'samples/{SAMPLE}.html': (
                f'<img src="../assets/{assets}.svg">'
                f'<a href="../assets/{assets}.pdf">PDF</a>'
                f'<a href="../assets/{assets}.bed">BED</a>'
                f'<a href="../data/custom_regions/{SAMPLE}.profile.tsv">Profile data</a>'
            ).encode(),
            PDF_NAME: PDF,
            f'assets/{assets}.svg': b'<svg xmlns="http://www.w3.org/2000/svg"/>',
            f'assets/{assets}.bed': b'chr2L\t100\t101\tsummit\n',
            f'data/custom_regions/{SAMPLE}.profile.tsv': b'scope\tbin\tmean_cpm\nTotal library\t0\t0\n',
            f'data/custom_regions/{SAMPLE}.json': b'{"retained_regions":1}',
        }
        for name, content in files.items():
            path = self.source / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)

    def test_worker_publishes_custom_bed_pdf_and_complete_zip(self):
        destination = self.root / 'published'
        (self.source / '.snakemake_timestamp').touch()
        publish_report(self.source, destination)
        self.assertEqual((destination / PDF_NAME).read_bytes(), PDF)
        with zipfile.ZipFile(destination / 'report.zip') as bundle:
            self.assertEqual(bundle.read('report/' + PDF_NAME), PDF)
            self.assertNotIn('report/.snakemake_timestamp', bundle.namelist())
            for path in self.source.rglob('*'):
                if path.is_file() and not path.name.startswith('.'):
                    self.assertEqual(bundle.read('report/' + path.relative_to(self.source).as_posix()),
                                     path.read_bytes())

    def test_frontend_installs_custom_bed_pdf_and_marks_success(self):
        frontend = Frontend({'state': str(self.root / 'frontend')})
        create_job(frontend.store)
        archive = self.root / 'incoming.zip'
        archive_report(self.source, archive)
        self.assertEqual(frontend.install_report(JOB_ID, archive, archive.stat().st_size), {'ok': True})
        job = frontend.store.get(JOB_ID)
        self.assertEqual((job['status'], job['progress'], job['error']), ('succeeded', 100, None))
        destination = frontend.state / 'reports' / JOB_ID
        self.assertEqual((destination / PDF_NAME).read_bytes(), PDF)
        self.assertEqual((destination / 'report.zip').read_bytes(), archive.read_bytes())

    def test_pdf_is_downloadable_only_with_private_report_access(self):
        settings = Settings(state=self.root / 'web', public_url='https://testserver',
                            allowed_hosts=['testserver'], secret='test-only-' * 8,
                            submission_key='', public_catalog={})
        store = Store(settings.state)
        create_job(store, 'succeeded')
        shutil.copytree(self.source, settings.state / 'reports' / JOB_ID)
        with TestClient(create_app(settings), base_url='https://testserver') as client:
            url = f'/reports/{JOB_ID}/{PDF_NAME}'
            self.assertEqual(client.get(url).status_code, 404)
            token = token_for(settings.secret, JOB_ID)
            response = client.get(url, headers={'Authorization': 'Bearer ' + token})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.content, PDF)
            self.assertEqual(response.headers['content-type'], 'application/pdf')
            self.assertEqual(response.headers['content-disposition'],
                             f'attachment; filename="{Path(PDF_NAME).name}"')
            self.assertEqual(response.headers['x-content-type-options'], 'nosniff')
            grant = report_grant(settings.secret, JOB_ID, int(time.time()) + 300)
            page = client.get(f'/reports/{JOB_ID}/samples/{SAMPLE}.html?view={grant}')
            self.assertIn(f'{SAMPLE}.custom_regions.pdf?view={grant}', page.text)
            self.assertEqual(client.get(url, params={'view': grant}).content, PDF)
            wrong_grant = report_grant(settings.secret, 'b' * 32, int(time.time()) + 300)
            self.assertEqual(client.get(url, params={'view': wrong_grant}).status_code, 404)

    def test_pdf_publication_still_enforces_size_limit(self):
        with self.assertRaisesRegex(ValueError, 'publication limits'):
            publish_report(self.source, self.root / 'published', max_bytes=len(PDF) - 1)
        self.assertFalse((self.root / 'published').exists())

    def test_worker_rejects_disallowed_and_disguised_intermediate_files(self):
        for name in ('unexpected.exe', 'reads.bam.pdf', '.hidden.pdf'):
            with self.subTest(name=name):
                invalid = self.source / 'assets' / name
                invalid.write_bytes(PDF)
                try:
                    with self.assertRaises(ValueError):
                        publish_report(self.source, self.root / 'published')
                    self.assertFalse((self.root / 'published').exists())
                finally:
                    invalid.unlink()

    def test_frontend_rejects_disallowed_and_disguised_intermediate_files(self):
        frontend = Frontend({'state': str(self.root / 'frontend')})
        create_job(frontend.store)
        for name in ('unexpected.exe', 'reads.bam.pdf', '.hidden.pdf'):
            with self.subTest(name=name):
                archive = self.root / 'incoming.zip'
                archive_report(self.source, archive)
                with zipfile.ZipFile(archive, 'a') as bundle:
                    bundle.writestr('report/assets/' + name, PDF)
                with self.assertRaises(ValueError):
                    frontend.install_report(JOB_ID, archive, archive.stat().st_size)
                self.assertEqual(frontend.store.get(JOB_ID)['status'], 'running')
                self.assertFalse((frontend.state / 'reports' / JOB_ID).exists())

    def test_worker_rejects_pdf_symlink(self):
        linked = self.source / 'assets/linked.pdf'
        linked.symlink_to(self.source / PDF_NAME)
        with self.assertRaisesRegex(ValueError, 'Unsupported file'):
            publish_report(self.source, self.root / 'published')
        self.assertFalse((self.root / 'published').exists())


if __name__ == '__main__':
    unittest.main()
