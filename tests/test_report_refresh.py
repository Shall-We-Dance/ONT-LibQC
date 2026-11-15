import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from report_refresh import refresh, rollback, build_archive, hashes
from report_profiles import render_profiles_card
from report_metrics import REPORT_CSS


class RefreshTests(unittest.TestCase):
    def fixture(self, root):
        for name in ('assets','samples','data/metrics','data/gene_profile'): (root/name).mkdir(parents=True,exist_ok=True)
        metric={'sample':'demo','assay':'chromatin','sub_libraries':[{'scope_id':'lib01','library_units':3},{'scope_id':'lib02','library_units':5}]}
        (root/'data/metrics/demo.json').write_text(json.dumps(metric))
        (root/'manifest.json').write_text(json.dumps({'schema_version':1,'samples':['demo']}))
        (root/'index.html').write_text('<html><body><a href="samples/demo.html">Demo</a><section id="rna-cohort"><h2>Retired</h2></section></body></html>')
        for sid in ('lib01','lib02'):
            (root/f'assets/demo.{sid}.full_library_profiles.tsv').write_text('size_bp\tlength_weighted_relative_signal\n0\t0\n100\t100\n2000\t0\n')
        scopes=[dict(metric,scope_id='total')]+metric['sub_libraries']
        sections=[f"<section class='profile-scope'><h3>{s['scope_id']}</h3><p>Keep this evidence</p></section>" for s in scopes]
        card=render_profiles_card('demo',scopes,sections,root/'assets',root/'data/gene_profile','../assets/')
        (root/'samples/demo.html').write_text('<html><head><style>'+REPORT_CSS+'</style></head><body>'+card+'</body></html>')
        build_archive(root)

    def test_preview_publish_archive_and_exact_rollback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'report';self.fixture(root);before=hashes(root)
            result=refresh(root)
            self.assertFalse(result['applied']);self.assertEqual(hashes(root),before)
            result=refresh(root,True)
            self.assertEqual(hashes(Path(result['backup'])),before)
            self.assertEqual(hashes(root)['data/metrics/demo.json'],before['data/metrics/demo.json'])
            self.assertNotIn('rna-cohort',(root/'index.html').read_text())
            with zipfile.ZipFile(root/'report.zip') as z:
                self.assertEqual(z.read('report/samples/demo.html'),(root/'samples/demo.html').read_bytes())
                self.assertNotIn('report/report.zip',z.namelist())
            rollback(root,result['backup']);self.assertEqual(hashes(root),before)

    def test_missing_numeric_input_leaves_original_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'report';self.fixture(root)
            (root/'assets/demo.lib02.full_library_profiles.tsv').unlink();before=hashes(root)
            with self.assertRaises(FileNotFoundError):refresh(root,True)
            self.assertEqual(hashes(root),before)

    def test_failed_publication_restores_original(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'report';self.fixture(root);before=hashes(root)
            rename=Path.rename
            def fail_stage(source,target):
                if '.refresh-' in str(source):raise OSError('injected publish failure')
                return rename(source,target)
            with patch.object(Path,'rename',fail_stage):
                with self.assertRaisesRegex(OSError,'publish failure'):refresh(root,True)
            self.assertEqual(hashes(root),before)

    def test_numeric_assets_are_immutable_during_refresh(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'report';self.fixture(root);before=hashes(root)
            from report_refresh import refresh_page
            def accidental_write(stage,sample):
                refresh_page(stage,sample)
                (stage/'assets/demo.lib01.full_library_profiles.tsv').write_text('changed')
            with patch('report_refresh.refresh_page',accidental_write):
                with self.assertRaisesRegex(ValueError,'Scientific record changed'):refresh(root,True)
            self.assertEqual(hashes(root),before)

    def test_symlink_and_wrong_backup_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'report';self.fixture(root)
            (root/'escape').symlink_to('/tmp')
            with self.assertRaisesRegex(ValueError,'Symlink'):refresh(root,True)
            with self.assertRaisesRegex(ValueError,'does not belong'):rollback(root,Path(tmp)/'unrelated')
