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
from report_metrics import COLLISION_ROWS, REPORT_CSS, scoped_table_html
from report_export import add_pdf_export
from report_refresh import refresh_collision_rows, refresh_presentation, Elements


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

    def test_collision_counts_and_low_depth_in_retained_report(self):
        scopes = [dict(scope_name='Total library', collision_observed_pairs=5, collision_oe=10,
                       collision_null_expected_pairs=0.5),
                  dict(scope_id='lib01', collision_observed_pairs=4, collision_oe=80,
                       collision_null_expected_pairs=0.05),
                  dict(scope_id='lib02', collision_observed_pairs=0, collision_null_expected_pairs=0),
                  dict(scope_id='lib03')]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'data/metrics').mkdir(parents=True)
            (root/'data/metrics/demo.json').write_text(json.dumps(dict(scopes[0], sub_libraries=scopes[1:])))
            old = scoped_table_html(scopes, ['collision_oe'], {})
            updated = refresh_collision_rows(root, 'demo', old)
            self.assertEqual(updated, scoped_table_html(scopes, COLLISION_ROWS, {}))
            self.assertEqual(updated, refresh_collision_rows(root, 'demo', updated))
            self.assertIn(">80</span></div><small class='metric-note'>Low depth</small>", updated)
            self.assertIn(">Not estimable</span></div><small class='metric-note'>Low depth</small>", updated)
            self.assertNotIn(">10</span></div><small class='metric-note'>", updated)
            self.assertEqual(updated.count("<small class='metric-note'>Low depth</small>"), 2)
            self.assertIn('>0.05</span>', updated)

    def test_archived_page_outside_the_manifest_keeps_its_original_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/'report'
            self.fixture(root)
            archival_page = root/'samples/archived.html'
            archival_page.write_text('<html><body>Prior report</body></html>')
            refresh(root, True)
            self.assertEqual(archival_page.read_text(), '<html><body>Prior report</body></html>')
            manifest = json.loads((root/'manifest.json').read_text())
            self.assertEqual(manifest['samples'], ['demo'])
            self.assertEqual(manifest['validation']['sample_reports'], 2)

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

    def test_export_styles_survive_refresh_and_repair_prior_corruption(self):
        for corrupted in (False, "named", "unnamed"):
            with self.subTest(corrupted=corrupted), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)/'report';self.fixture(root)
                page=root/'samples/demo.html'
                document=page.read_text().replace('<body>', '<body><header class="report-hero"><h1>Demo</h1></header>')
                document=add_pdf_export(document)
                # Another component can reference the same selectors; its ID
                # must keep it outside the base-style replacement boundary.
                component='<style id="custom-component">.report-hero .metric-tip{color:red}</style>'
                document=document.replace('</head>',component+'</head>')
                parsed=Elements(document)
                export=next(n for n in parsed.elements if n['attrs'].get('id')=='report-export-style')
                expected=parsed.fragment(export)
                if corrupted:
                    tag='<style id="report-export-style">' if corrupted=='named' else '<style>'
                    document=document.replace(expected,tag+REPORT_CSS+'</style>')
                page.write_text(document)
                refresh(root,True)
                updated=page.read_text()
                self.assertIn(expected,updated)
                self.assertIn(component,updated)
                self.assertEqual(updated.count('id="report-base-style"'),1)
                self.assertEqual(updated.count('id="report-export-style"'),1)
                self.assertEqual(updated.count('id="report-export-pdf"'),1)
                self.assertEqual(updated.count('id="report-export-script"'),1)
                self.assertEqual(refresh_presentation(updated),updated)
                with zipfile.ZipFile(root/'report.zip') as z:
                    self.assertEqual(z.read('report/samples/demo.html'),page.read_bytes())

    def test_single_library_preserves_embedded_card_layout(self):
        for wrapper in ("div class='library-profiles'", "section class='report-section'"):
            with self.subTest(wrapper=wrapper), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)/'report';self.fixture(root)
                metric=root/'data/metrics/demo.json'
                metric.write_text(json.dumps({'sample':'demo','assay':'chromatin','sub_libraries':[]}))
                scope="<section class='profile-scope'><h3>Total library</h3><p>Retained figure</p></section>"
                page=root/'samples/demo.html'
                page.write_text("<html><head><style>"+REPORT_CSS+"</style></head><body>"
                                "<section id='library-shape'><"+wrapper+" id='profiles'><h3>Size and genomic profiles</h3>"
                                +scope+"</"+wrapper.split()[0]+"></section></body></html>")
                refresh(root,True)
                nodes=Elements(page.read_text()).elements
                profile=next(n for n in nodes if n['attrs'].get('id')=='profiles')
                self.assertEqual(profile['tag'],'div')
                self.assertEqual(profile['attrs']['class'],'library-profiles')
                self.assertIn(scope,page.read_text())

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
