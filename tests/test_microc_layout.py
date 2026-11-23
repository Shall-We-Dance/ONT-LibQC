"""Render and reopen a portable Micro-C overview with and without sublibraries."""
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from microc_ligation import summarize
from report_microc_layout import collect_ligation_assets, overview_html, prepare_microc_profiles, write_overview_figure
from report_profiles import render_profiles_card
from report_refresh import Elements, refresh_profile_layout


class MicrocLayoutTests(unittest.TestCase):
    def test_ligation_links_and_export_follow_available_panels(self):
        for counts, lengths in ((True, True), (True, False), (False, True)):
            with self.subTest(counts=counts, lengths=lengths), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                evidence, assets = root/'data/ligation', root/'assets'
                evidence.mkdir(parents=True)
                assets.mkdir()
                ids = [f'lib{i:02d}' for i in range(1, 10)]
                scopes = [dict(scope_id=sid, unique_valid_contacts=0) for sid in ['total', *ids]]
                if counts:
                    stats = summarize([{'category': '1', 'read_id': 'test'}])
                    (evidence/'s.json').write_text(json.dumps(dict(stats, scopes={sid: stats for sid in ids})))
                if lengths:
                    with gzip.open(evidence/'s.units.tsv.gz', 'wt') as units:
                        units.write('scope_id\tcategory\tinsert_bp\n')
                        for sid in ids:
                            units.write(f'{sid}\t1\t120\n')
                (assets/'s.total.ligation.svg').write_text('<svg/>')
                (assets/'s.ligation_examples.svg').write_text('<svg/>')
                profiles = prepare_microc_profiles('s', scopes, root/'data', assets, '../assets/', '../data')
                self.assertEqual(profiles.section_ids, ('ligation',))
                self.assertIn('s.total.ligation.svg', profiles.details['total'])
                self.assertIn('s.ligation_examples.svg', profiles.details['total'])
                html = overview_html('../assets/', profiles)
                parsed = Elements(html)
                self.assertEqual(html.count("id='ligation'"), 1)
                self.assertEqual(html.count('<summary>'), 1)
                urls = [n['attrs']['href'] for n in parsed.elements if n['tag'] == 'a'
                        and n['attrs']['href'].endswith('.svg')]
                expected = (['../assets/s.microc.ligation.svg'] if counts else []) + (
                    ['../assets/s.microc.ligation_lengths.svg'] if lengths else [])
                self.assertEqual(urls, expected)
                self.assertEqual([n['attrs']['src'] for n in parsed.elements if n['tag'] == 'img'], expected)

                for name, width in [('size.svg', 780), ('gel.svg', 190)]:
                    (assets/name).write_text(f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} 308'></svg>")
                name = write_overview_figure('s', assets, 'size.svg', 'gel.svg', profiles, scopes[1:])
                panels = ET.parse(assets/name).getroot().findall('{*}svg')
                ligation_row = panels[2:-1]
                self.assertEqual(len(ligation_row), len(expected))
                self.assertEqual(len({p.get('y') for p in ligation_row}), 1)
                for panel in ligation_row:
                    end = float(panel.get('y')) + float(panel.get('height'))
                    self.assertLess(end, float(panels[-1].get('y')))
                if counts and lengths:
                    left, right = ligation_row
                    self.assertLess(float(left.get('x')) + float(left.get('width')), float(right.get('x')))

    def test_reused_work_directory_does_not_publish_old_ligation_scope_images(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            plots, assets, data = out/'plots/ligation', out/'reports/assets', out/'qc/ligation'
            for folder in (plots, assets, data):
                folder.mkdir(parents=True)
            scopes = [{'scope_id': sid} for sid in ('total', 'lib01', 'lib02')]
            summary = data/'s.json'
            summary.write_text(json.dumps({'scopes': {'lib01': {}, 'lib02': {}}}))
            for scope in scopes:
                for suffix in ('svg', 'mnase.svg', 'single_block.svg'):
                    (plots/f"s.{scope['scope_id']}.{suffix}").write_text('first run')
            (plots/'s.examples.svg').write_text('old examples')
            collect_ligation_assets(out, 's', scopes)

            # The plotting tool enumerates observed scopes on rerun. It does not
            # remove source files for lib02 when that scope now has zero units.
            summary.write_text(json.dumps({'scopes': {'lib01': {}}}))
            (plots/'s.lib01.svg').write_text('second run')
            (plots/'s.lib01.mnase.svg').unlink()
            (plots/'s.examples.svg').unlink()
            collect_ligation_assets(out, 's', scopes)
            self.assertEqual((assets/'s.lib01.ligation.svg').read_text(), 'second run')
            self.assertFalse((assets/'s.lib01.mnase.svg').exists())
            self.assertFalse((assets/'s.ligation_examples.svg').exists())
            self.assertFalse(list(assets.glob('s.lib02.*')))
            self.assertEqual(len(list(plots.glob('s.lib02.*'))), 3)
            self.assertTrue((assets/'s.total.ligation.svg').is_file())

            # An older run without the JSON may only retain a report-side plot.
            summary.unlink()
            (plots/'s.lib01.svg').unlink()
            collect_ligation_assets(out, 's', scopes)
            self.assertEqual((assets/'s.lib01.ligation.svg').read_text(), 'second run')

    def test_historical_multiscope_report_with_missing_size_table_keeps_existing_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'assets').mkdir()
            (root/'data/metrics').mkdir(parents=True)
            subs = [dict(scope_id=sid, library_units=1) for sid in ('lib01', 'lib02')]
            metric = dict(assay='microc', sample='s', library_units=2, sub_libraries=subs)
            (root/'data/metrics/s.json').write_text(json.dumps(metric))
            (root/'assets/s.lib01.full_library_profiles.tsv').write_text(
                'size_bp\tlength_weighted_relative_signal\n0\t0\n100\t100\n')
            gallery = ("<nav><a href='#profiles'>Profiles</a><a href='#ligation'>Ligation</a></nav>"
                       "<section id='profiles'><h2>Size profiles</h2><img src='../assets/s.sizes.svg'></section>"
                       "<section id='ligation'><img src='../assets/s.ligation.svg'></section>")
            page = "<section id='sublibrary-composition'>Old mixture</section>" + gallery
            updated = refresh_profile_layout(root, 's', page)
            self.assertTrue(updated.endswith(gallery))
            self.assertIn('class="mixture-bar"', updated)
            self.assertEqual(refresh_profile_layout(root, 's', updated), updated)

    def test_scope_disclosures_download_and_refresh(self):
        for ids in [('total',), ('total', 'lib01', 'lib02')]:
            with self.subTest(scopes=ids), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                assets, data = root / 'assets', root / 'data'
                assets.mkdir()
                (data / 'ligation').mkdir(parents=True)
                (data / 'metrics').mkdir()
                scopes = [dict(scope_id=sid, library_units=3) for sid in ids]
                scopes[0]['scope_name'] = 'Total library'
                sample = 'test'
                with gzip.open(data / 'ligation/test.mnase_fragments.tsv.gz', 'wt') as table:
                    table.write('scope_id\treference_bp\nlib01\t128\n')
                metric = dict(assay='microc', sample=sample, library_units=3, sub_libraries=scopes[1:])
                (data / 'metrics/test.json').write_text(json.dumps(metric))
                sections = []
                for sid in ids:
                    (assets / f'test.{sid}.full_library_profiles.tsv').write_text(
                        'size_bp\tlength_weighted_relative_signal\n0\t0\n100\t100\n2000\t0\n')
                    (assets / f'test.{sid}.single_block.svg').write_text('<svg/>')
                    sections.append(f"<section class='profile-scope'><p>{sid}</p></section>")
                prepared = prepare_microc_profiles(sample, scopes, data, assets, '../assets/', '../data')
                page = render_profiles_card(sample, scopes, sections, assets, data/'gene_profile',
                                            '../assets/', microc=prepared)
                parsed = Elements(page)
                disclosures = [n for n in parsed.elements if 'profile-detail' in n['attrs'].get('class', '').split()]
                self.assertEqual(len(disclosures), len(ids))
                self.assertEqual([n['attrs'].get('data-scope', 'total') for n in disclosures],
                                 list(ids[1:]) + ['total'])
                for n in disclosures:
                    sid = n['attrs'].get('data-scope', 'total')
                    self.assertIn(f'test.{sid}.single_block.svg', parsed.fragment(n))
                    self.assertEqual('open' in n['attrs'], len(ids) == 1)
                self.assertIn("href='../assets/test.microc_qc.svg' download", page)
                self.assertNotIn('test.mnase.json', page)
                combined = ET.parse(assets / 'test.microc_qc.svg')
                self.assertEqual(len(combined.findall('{*}svg')), 3)
                self.assertEqual(refresh_profile_layout(root, sample, page), page)

                # A compact archive can retain the chart without its source table.
                # Refresh must preserve the image, the scope details, and its navigation target.
                (data/'ligation/test.mnase_fragments.tsv.gz').unlink()
                retained = page.replace("<p><a href='../data/ligation/test.mnase_fragments.tsv.gz'>"
                                        "Fragment lengths and parent molecules</a></p>", '')
                self.assertNotIn('mnase_fragments.tsv.gz', retained)
                self.assertEqual(refresh_profile_layout(root, sample, retained), retained)
