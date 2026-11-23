import fcntl
import gzip
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

    def test_manifest_precondition_uses_the_version_under_lock(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/'report'
            self.fixture(root)
            expected = hashes(root)['manifest.json']
            manifest = root/'manifest.json'
            newer = json.loads(manifest.read_bytes())
            newer['report_version'] = 'newer'
            lock = fcntl.flock

            def acquire_after_other_publisher(fd, flags):
                manifest.write_text(json.dumps(newer))
                lock(fd, flags)

            with patch('report_refresh.fcntl.flock', acquire_after_other_publisher):
                with self.assertRaisesRegex(ValueError, 'manifest has changed'):
                    refresh(root, apply=True, expected_manifest_sha256=expected)
            self.assertEqual(json.loads(manifest.read_bytes()), newer)

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

    def test_chromatin_profiles_are_sibling_cards_after_refresh(self):
        for count, wrapper in [(0, 'div'), (0, 'section'), (2, 'section')]:
            with self.subTest(scopes=count, wrapper=wrapper), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp)/'report'; self.fixture(root)
                page=root/'samples/demo.html'
                fixture=Elements(page.read_text())
                card=next(n for n in fixture.elements if n['attrs'].get('id')=='profiles')
                evidence="<section class='profile-scope'><h3>Total library</h3><p>Original figures</p></section>"
                if count:
                    inner=fixture.fragment(card)
                else:
                    (root/'data/metrics/demo.json').write_text(json.dumps(dict(sample='demo', assay='chromatin', sub_libraries=[])))
                    inner=f"<{wrapper} class='library-profiles' id='profiles'><h3>Size and genomic profiles</h3>{evidence}</{wrapper}>"
                nav="<nav><a href='#library-shape'>Library quality</a></nav>"
                quality="<section id='library-shape'><div>"+inner+"</div><p>Other evidence</p></section>"
                page.write_text("<html><head><style>"+REPORT_CSS+"</style></head><body>"+nav+quality+"</body></html>")
                refresh(root,True)
                updated=page.read_text()
                nodes=Elements(updated)
                profile=next(n for n in nodes.elements if n['attrs'].get('id')=='profiles')
                library=next(n for n in nodes.elements if n['attrs'].get('id')=='library-shape')
                self.assertEqual(profile['tag'],'section')
                self.assertIn('report-section',profile['attrs']['class'].split())
                self.assertIs(profile['parent'],library['parent'])
                self.assertGreaterEqual(profile['start'],library['end'])
                self.assertIn('<p>Other evidence</p>',nodes.fragment(library))
                self.assertIn('<h2>Size and genomic profiles</h2>',nodes.fragment(profile))
                self.assertEqual(sum(n['attrs'].get('href')=='#profiles' for n in nodes.elements),1)
                if count:
                    self.assertEqual(nodes.fragment(profile).count("class='profile-detail"),3)
                else:
                    self.assertIn(evidence,nodes.fragment(profile))
                refresh(root,True)
                self.assertEqual(page.read_text(),updated)

    def test_chromatin_region_identity_and_single_total_color(self):
        from report_regions import svg_region_bars
        from report_palette import POOLED, scope_color
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/'report'
            self.fixture(root)
            metric = dict(promoter_proximal_rate=.35, promoter_proximal_background_rate=.1)
            original = svg_region_bars(metric, 'Keep the title', '#2a6cbd')
            self.assertIn(POOLED, svg_region_bars(metric, 'Single library'))
            for sid in ('lib01', 'lib02', 'total'):
                (root/'assets'/f'demo.{sid}.regions.svg').write_text(original)
            refresh(root, True)
            for sid in ('lib01', 'lib02', 'total'):
                rendered = (root/'assets'/f'demo.{sid}.regions.svg').read_text()
                self.assertEqual(rendered, original.replace('#2a6cbd', scope_color(sid, 0)))
            self.assertNotEqual(POOLED, scope_color('lib01', 0))

    def test_rna_scope_palette_preserves_source_svg_and_evidence_colors(self):
        from report_sublibrary import scope_color
        from rna_report import rna_scope_card
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/'report'
            self.fixture(root)
            metrics = root/'data/metrics/demo.json'
            record = json.loads(metrics.read_text())
            record['assay'] = 'rnaseq'
            metrics.write_text(json.dumps(record))
            folder = root/'data/rna/demo'
            folder.mkdir(parents=True)
            old = '#23877f'
            trace = ('<svg xmlns="http://www.w3.org/2000/svg"><text>Observed coverage · equal gene weighting</text>'
                     '<g id="axes_1"><path d="M1 2 L3 4" style="fill: none; stroke: '+old+'; stroke-width: 1.7;"/>'
                     '<path d="M1 3 L3 5" style="fill: '+old+'; fill-opacity: 0.18"/></g></svg>')
            (folder/'lib02.profiles.svg').write_text(trace)
            categorical = trace.replace('</svg>', '<g id="axes_2"/></svg>')
            (folder/'lib01.profiles.svg').write_text(categorical)
            rna = {sid: dict(units=10, family_counts={}, plots={'profiles': sid+'.profiles.svg'})
                   for sid in ('lib01', 'lib02')}
            (folder/'summary.json').write_text(json.dumps(dict(scopes=rna)))
            page = root/'samples/demo.html'
            source_url = '../data/rna/demo/lib02.profiles.svg'
            evidence = '<p style="background:'+old+'">Reference evidence</p>'
            detail = rna_scope_card('lib02', 'Second', 10, old, f'{evidence}<a href="{source_url}">Plot</a>')
            page.write_text(page.read_text().replace('</body>', '<section id="rna-qc">'+detail+'</section></body>'))
            refresh(root, True)
            rendered = page.read_text()
            color = scope_color('lib02', 0)
            target = root/'assets/demo.lib02.rna_transcript.svg'
            self.assertEqual(target.read_text(), trace.replace(old, color))
            self.assertEqual((folder/'lib02.profiles.svg').read_text(), trace)
            self.assertIn('--scope-tint:'+color+'0f', rendered)
            self.assertIn(evidence, rendered)
            self.assertNotIn(source_url, rendered)
            self.assertGreaterEqual(rendered.count('../assets/'+target.name), 3)
            self.assertIn('../data/rna/demo/lib01.profiles.svg', rendered)
            self.assertFalse((root/'assets/demo.lib01.rna_transcript.svg').exists())
            refresh(root, True)
            self.assertEqual(page.read_text(), rendered)
            with zipfile.ZipFile(root/'report.zip') as archive:
                self.assertEqual(archive.read('report/assets/'+target.name), target.read_bytes())
                self.assertEqual(archive.read('report/data/rna/demo/lib02.profiles.svg'), trace.encode())

    def ligation_fixture(self, root, *, combined):
        """The original standalone section and the later two-plots-in-one SVG layout."""
        for folder in ('assets', 'samples', 'data/metrics', 'data/ligation'):
            (root/folder).mkdir(parents=True, exist_ok=True)
        (root/'manifest.json').write_text(json.dumps(dict(schema_version=1, samples=['demo'])))
        (root/'index.html').write_text('<a href="samples/demo.html">Report</a>')
        (root/'data/metrics/demo.json').write_text(json.dumps(dict(sample='demo', assay='microc')))
        from microc_ligation import summarize
        stats = summarize([dict(category='2', read_id='read1')])
        (root/'data/ligation/demo.json').write_text(json.dumps(dict(stats, scopes={})))
        with gzip.open(root/'data/ligation/demo.units.tsv.gz', 'wt') as stream:
            stream.write('scope_id\tcategory\tinsert_bp\nlib01\t2\t240\n')
        (root/'assets/demo.total.full_library_profiles.tsv').write_text(
            'size_bp\tlength_weighted_relative_signal\n0\t0\n100\t100\n2000\t0\n')
        legacy_svg = "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 500 160'><title>Retained ligation evidence</title></svg>"
        (root/'assets/demo.total.ligation.svg').write_text(legacy_svg)
        plot = "<img src='../assets/demo.total.ligation.svg'>"
        scope = "<section class='profile-scope'><h3>Total library</h3></section>"
        if combined:
            (root/'assets/demo.microc.ligation.svg').write_text(legacy_svg)
            inner = "<div class='microc-comparison' id='ligation'><img src='../assets/demo.microc.ligation.svg'></div>"
            profiles = "<section id='profiles'>" + inner + scope + "</section>"
        else:
            profiles = "<section id='profiles'>" + scope + "</section><section id='ligation'>" + plot + "</section>"
        page = root/'samples/demo.html'
        page.write_text("<html><head><style>"+REPORT_CSS+"</style></head><body><nav><a href='#library-shape'>Library quality</a></nav>"
                        "<section id='library-shape'><h2>Library quality</h2></section>" + profiles + "</body></html>")
        build_archive(root)
        return page

    def test_combined_ligation_upgrades_both_downloads_without_losing_pruned_panel(self):
        for source_to_remove in ('', 'demo.json', 'demo.units.tsv.gz'):
            with self.subTest(source_to_remove=source_to_remove), tempfile.TemporaryDirectory() as t:
                root = Path(t)/'report'
                page = self.ligation_fixture(root, combined=True)
                if source_to_remove:
                    (root/'data/ligation'/source_to_remove).unlink()
                old_svg = (root/'assets/demo.microc.ligation.svg').read_bytes()
                refresh(root, True)
                text = page.read_text()
                if source_to_remove:
                    self.assertNotIn('demo.microc.ligation_lengths.svg', text)
                    self.assertEqual((root/'assets/demo.microc.ligation.svg').read_bytes(), old_svg)
                else:
                    self.assertEqual(text.count("id='ligation'"), 1)
                    for name in ('demo.microc.ligation.svg', 'demo.microc.ligation_lengths.svg'):
                        self.assertIn("src='../assets/"+name+"'", text)
                        with zipfile.ZipFile(root/'report.zip') as archive:
                            self.assertEqual(archive.read('report/assets/'+name), (root/'assets'/name).read_bytes())
                    self.assertIn("href='#ligation'", text)
                assets = hashes(root/'assets')
                refresh(root, True)
                self.assertEqual(text, page.read_text())
                self.assertEqual(assets, hashes(root/'assets'))

    def test_standalone_ligation_with_one_source_keeps_scope_image_and_one_anchor(self):
        cases = [('demo.json', 'ligation_lengths', 'ligation'), ('demo.units.tsv.gz', 'ligation', 'ligation_lengths')]
        for missing, available, absent in cases:
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as t:
                root = Path(t)/'report'
                page = self.ligation_fixture(root, combined=False)
                (root/'data/ligation'/missing).unlink()
                refresh(root, True)
                first = page.read_text()
                self.assertEqual(first.count("id='ligation'"), 1)
                self.assertEqual(first.count("href='#ligation'"), 1)
                self.assertIn("href='#profiles'", first)
                self.assertIn("src='../assets/demo.total.ligation.svg'", first)
                self.assertIn(f"src='../assets/demo.microc.{available}.svg'", first)
                self.assertNotIn(f"src='../assets/demo.microc.{absent}.svg'", first)
                # A previously split one-panel report can refresh its present
                # chart from its own source without demanding the missing other input.
                chart = root/'assets'/f'demo.microc.{available}.svg'
                expected = chart.read_bytes()
                chart.write_text('<svg><title>Earlier rendering</title></svg>')
                refresh(root, True)
                self.assertEqual(first, page.read_text())
                self.assertEqual(chart.read_bytes(), expected)

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
