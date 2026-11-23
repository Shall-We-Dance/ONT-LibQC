"""Exercise both portal execution modes and the private single-cell report bundle."""
import gzip
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import zipfile

from fastapi.testclient import TestClient
from ont_portal.app import create_app
from ont_portal.run_job import normalize_fastq, run_phase
from ont_portal.settings import Settings
from ont_portal.store import token_for
from ont_portal.catalog import ReferenceMounts
from ont_portal.worker import check_reference_mounts, publish_report, Worker

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow/scripts'))
from scrna_report import render


class ScrnaRunnerTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.options = {'sample': 'control', 'reference': 'toy', 'assay': 'scrna', 'params': {'mapq': 32}}
        self.request = {'options': self.options, 'cores': 2, 'compressed': False, 'custom_regions': False,
                        'max_reads': 50, 'max_expanded_bytes': 10000}
        self.resources = {'assembly': 'toy', 'barcodes': '/refs/barcodes', 'star_index': '/refs/star',
                          'gtf': '/refs/genes.gtf', 'species_panel': '/job/scrna-species.json'}
        self.cfg = {'scrna': self.resources}
        self.upload = self.root / 'input.upload'
        self.upload.write_bytes(b'@synthetic\nACTG\n+\nIIII\n')
        self.counts = b'scope\tbarcode\tgene_id\tumis\texonic_umis\n'

    def normalize(self, path, *args):
        self.assertEqual(path, '/input/input.upload')
        return normalize_fastq(self.upload, *args)

    def analysis(self, args, *, cwd, check):
        self.assertEqual(args[2], 'scrna')
        self.assertTrue(check)
        flags = dict(zip(args[3::2], args[4::2]))
        self.assertEqual(flags['--mapq'], '32')
        self.assertEqual(flags['--threads'], '2')
        self.assertNotIn('--keep-work', args)
        self.assertNotIn('--single-capture', args)
        for key, value in self.resources.items():
            self.assertEqual(flags['--' + key.replace('_', '-')], value)
        with gzip.open(flags['--fastq'], 'rb') as handle:
            self.assertEqual(handle.read(), self.upload.read_bytes())
        out = Path(flags['--out'])
        out.mkdir()
        pair = {'scope_id': 'lib01', 'name': 'lib01', 'i5': '', 'i7': 'ACGT'}
        summary = {'sample': flags['--sample'], 'assembly': flags['--assembly'], 'raw_reads': 1,
                   'non_dimer_units': 1, 'metrics': {},
                   'analysis_provenance': {'pipeline_version': 'v1.2.0', 'pipeline_build': 'test',
                                           'analysis_completed_utc': '2026-01-12T10:00:00Z'},
                   'physical_library': {'library_units': 1, 'sub_libraries': [{'scope_id': 'lib01', 'library_units': 1}],
                                        'index_pairing': {'pairs': [pair]}},
                   'scopes': {'lib01': {'kit': None, 'units': 1, 'evidence': {}}}}
        (out / 'summary.json').write_text(json.dumps(summary))
        (out / 'report.html').write_text(render(summary))
        (out / 'barcode_gene_counts.tsv.gz').write_bytes(gzip.compress(self.counts))
        (out / 'transcript_regions.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
        # Only explicitly identified CLI deliverables may reach the portal ZIP.
        (out / 'work').mkdir()
        (out / 'work/private.fastq').write_bytes(self.upload.read_bytes())

    def execute(self, phases):
        with patch('ont_portal.run_job.normalize_fastq', side_effect=self.normalize) as norm, \
             patch('ont_portal.run_job.subprocess.run', side_effect=self.analysis) as cli:
            for phase in phases:
                result = run_phase(self.root, self.request, self.cfg, phase)
                (self.root / 'result.json').write_text(json.dumps(result))
                self.assertEqual((result['phase'], result['raw_reads']), (phase, 1))
                if phase == 'prepare':
                    self.assertFalse((self.root / 'scrna').exists())
                if phase in ('prepare', 'align'):
                    self.assertFalse((self.root / 'report').exists())
            self.assertEqual(result['validation']['sample_reports'], 1)
            norm.assert_called_once()
            cli.assert_called_once()
        self.assertFalse((self.root / 'workflow').exists())
        report = self.root / 'report'
        self.assertEqual(json.loads((report / 'manifest.json').read_text())['entrypoint'], 'samples/control.html')
        self.assertEqual(json.loads((report / 'data/submission.json').read_text())['options'], self.options)
        return report

    def test_full_runs_one_cli_without_species_panel(self):
        del self.resources['species_panel']
        self.execute(['full'])

    def test_missing_or_failed_phase_stops_report(self):
        for previous in ({'ok': False, 'phase': 'align'}, {'ok': True, 'phase': 'prepare'}):
            (self.root / 'result.json').write_text(json.dumps(previous))
            with self.assertRaisesRegex(ValueError, 'previous phase'):
                run_phase(self.root, self.request, self.cfg, 'report')
        (self.root / 'result.json').write_text('{"ok":true,"phase":"prepare","raw_reads":1}')
        with patch('ont_portal.run_job.subprocess.run', side_effect=subprocess.CalledProcessError(1, 'scrna')):
            with self.assertRaises(subprocess.CalledProcessError):
                run_phase(self.root, self.request, self.cfg, 'align')
        self.assertFalse((self.root / 'report').exists())

    def test_report_artifact_cannot_be_a_link_to_a_private_intermediate(self):
        out = self.root / 'scrna'
        out.mkdir()
        (out / 'summary.json').symlink_to(self.upload)
        (self.root / 'result.json').write_text('{"ok":true,"phase":"align","raw_reads":1}')
        with self.assertRaisesRegex(ValueError, 'file type'):
            run_phase(self.root, self.request, self.cfg, 'report')
        self.assertFalse((self.root / 'report').exists())

    def test_phased_report_sharing_editing_and_portable_archive(self):
        source = self.execute(['prepare', 'align', 'report'])
        settings = Settings(state=self.root/'state', public_url='https://testserver', allowed_hosts=['testserver'],
                            secret='test-only-'*8, submission_key='', public_catalog={}, min_free_bytes=0)
        app = create_app(settings)
        job_id = 'c'*32
        now = time.time()
        with app.state.store.connect() as db:
            db.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires) '
                       'VALUES(?,?,?,?,?,?,?,?,?,?)', (job_id, now, now, 'succeeded', '', 'control.fq', 20,
                                                     json.dumps(self.options), 'local', now+3600))
        publish_report(source, settings.state/'reports'/job_id)
        with TestClient(app, base_url=settings.public_url, headers={'X-Portal-Request': '1'}) as client:
            task_url = '/api/jobs/' + job_id
            self.assertEqual(client.get(task_url).status_code, 404)
            owner = {'Authorization': 'Bearer ' + token_for(settings.secret, job_id)}
            task = client.get(task_url, headers=owner).json()
            page_url = task['report_url']
            self.assertRegex(page_url, r'^/r/[\w-]{22}$')
            page = client.get(page_url)
            self.assertEqual(page.status_code, 200)
            self.assertIn("connect-src 'none'", page.headers['Content-Security-Policy'])
            self.assertIn('data:image/svg+xml;base64,', page.text)
            self.assertIn('Single-cell RNA QC', page.text)
            self.assertIn('Download complete report (ZIP)', page.text)
            self.assertIn('Edit names', page.text)
            file_list = re.search(r'<a href="([^"]+)">Browse report files</a>', page.text)[1]
            files = client.get(file_list)
            matrix = re.search(r'<a href="([^"]+)">Barcode and gene UMI counts', files.text)[1]
            self.assertEqual(gzip.decompress(client.get(matrix).content), self.counts)
            names = page_url.replace('/r/', '/api/sublibrary-names/')
            saved = client.put(names, json={'revision': 0, 'names': {'lib01': 'Rep A'}})
            self.assertEqual(saved.status_code, 200)
            self.assertIn('Rep A', client.get(page_url).text)
            archive = client.get(task['download_url'])
            self.assertIn('report_control_2026-01-12.zip', archive.headers['Content-Disposition'])
            with zipfile.ZipFile(io.BytesIO(archive.content)) as z:
                self.assertIn('Rep A', z.read('report/samples/control.html').decode())
                self.assertIn('control\tlib01\tRep A', z.read('report/data/sublibrary_names.tsv').decode())
                self.assertIn('report/assets/transcript_regions.svg', z.namelist())
                self.assertNotIn('report/assets/rna_families.svg', z.namelist())
                self.assertFalse(any('fastq' in name or '/work/' in name for name in z.namelist()))
                self.assertEqual(gzip.decompress(z.read('report/data/scrna/barcode_gene_counts.tsv.gz')), self.counts)

    def test_auto_star_peak_and_limit(self):
        worker = Worker.__new__(Worker)
        worker.memory = '2g'
        worker.cfg = {'max_memory': '8g'}
        reference = {'star_index': '/refs/star'}
        worker.catalog = {'references': {'toy': reference}}
        options = SimpleNamespace(mapper='auto', reference='toy', assay='scrna')
        with patch('ont_portal.worker.star_memory_bytes', return_value=6*1024**3) as required:
            self.assertEqual(worker.memory_for(options), 6*1024**3)
            required.assert_called_once_with(reference)
            worker.cfg['max_memory'] = '4g'
            with self.assertRaisesRegex(ValueError, 'memory budget'):
                worker.memory_for(options)

    def test_exact_reference_mounts_cannot_expose_private_inputs(self):
        state = self.root / 'state'
        state.mkdir()
        secret = self.root / 'secrets.json'
        secret.touch()
        private = state / 'another-job'
        private.touch()
        mounts = ReferenceMounts()
        mounts.bind(self.upload, '/refs/valid-file')
        check_reference_mounts(mounts, state, secret)
        for path in (private, secret):
            mounts.bind(path, '/refs/unavailable')
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, 'private service'):
                check_reference_mounts(mounts, state, secret)
        with self.assertRaisesRegex(ValueError, 'regular reference file'):
            mounts.bind(state, '/refs/directory')


if __name__ == '__main__':
    unittest.main()
