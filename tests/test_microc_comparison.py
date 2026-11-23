"""Population accounting and missing measurements in the Micro-C overview."""
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow/scripts'))
from microc_ligation import CATEGORIES as BLOCK_CATEGORIES
from microc_ligation import summarize as summarize_ligation
from microc_pairtools import CATEGORIES as PAIR_CATEGORIES
from report_microc import native_contact_summary, write_microc_comparisons
from report_microc_layout import color_retained_ligation
from report_sublibrary import bar_ink, scope_color


class MicroCComparisonTests(unittest.TestCase):
    def test_ligation_figures_have_separate_titles_and_independent_inputs(self):
        names = [f'long_library_label_{i}' for i in range(9)]
        scopes = [{'scope_id': name} for name in names]
        rows = [dict(scope_id=name, category='2', read_id=name) for name in names]
        ligation = dict(summarize_ligation(rows),
                        scopes={name: summarize_ligation([row]) for name, row in zip(names, rows)})
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            data, assets = root / 'ligation', root / 'assets'
            data.mkdir()
            record, units = data / 'x.json', data / 'x.units.tsv.gz'
            record.write_text(json.dumps(ligation))
            self.assertEqual(write_microc_comparisons('x', scopes, root, assets),
                             {'ligation': 'x.microc.ligation.svg'})
            with gzip.open(units, 'wt') as stream:
                stream.write('scope_id\tcategory\tinsert_bp\n')
                for i, name in enumerate(names):
                    stream.write(f'{name}\t2\t{150+i}\n')
            for group in ([{'scope_id': 'total'}], scopes[:2], scopes[:4], scopes):
                result = write_microc_comparisons('x', group, root, assets)
                self.assertEqual(result, {'ligation': 'x.microc.ligation.svg',
                                          'ligation_lengths': 'x.microc.ligation_lengths.svg'})
                counts = ET.parse(assets / result['ligation'])
                lengths = ET.parse(assets / result['ligation_lengths'])
                self.assertEqual(counts.find('{*}title').text, 'Genomic blocks per library unit')
                self.assertEqual(lengths.find('{*}title').text, 'Insert length by block count')
                self.assertIsNone(counts.find(".//{*}g[@data-boxplot]"))
                self.assertIsNone(lengths.find(".//{*}rect[@data-category]"))
                keys = lengths.findall(".//{*}rect[@data-library-key]")
                self.assertEqual([k.get('data-library-key') for k in keys], [sc['scope_id'] for sc in group])
                labels = [t.text for t in lengths.findall('.//{*}text')]
                self.assertIn('Adapter-bounded insert length (bp)', labels)
                self.assertNotIn('Left to right', labels)
                header = next(t for t in counts.findall('{*}text') if t.text == 'Library')
                left, previous = float(header.get('x')), float(header.get('y'))
                for sc in group:
                    row = counts.find(f"{{*}}g[@data-scope='{sc['scope_id']}']")
                    bar = row.find('{*}rect')
                    self.assertEqual(float(bar.get('x')), left)
                    self.assertGreaterEqual(float(row.find('{*}text').get('y')) - previous, 26)
                    previous = float(bar.get('y')) + float(bar.get('height'))
                legend_top = min(float(r.get('y')) for r in counts.findall('{*}rect') if r.get('y'))
                self.assertGreaterEqual(legend_top - previous, 18)
            record.unlink()
            self.assertEqual(write_microc_comparisons('x', scopes, root, assets),
                             {'ligation_lengths': 'x.microc.ligation_lengths.svg'})

    def test_retained_matplotlib_detail_changes_only_box_facecolors(self):
        original = (ROOT / 'example_report/assets/HL120.total.ligation.svg').read_bytes()
        expected = ET.fromstring(original)
        axes = expected.find(".//{*}g[@id='axes_2']")
        boxes = [g.find('{*}path') for g in axes.findall('{*}g') if g.get('id', '').startswith('patch_')]
        color = scope_color('lib02', 1)
        changed = 0
        for box in boxes:
            style = box.get('style', '')
            if box.get('clip-path') and style.startswith('fill: #') and 'stroke: #000000' in style:
                box.set('style', 'fill: ' + color + ';' + style.split(';', 1)[1])
                changed += 1
        self.assertGreater(changed, 0)
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 's.lib02.ligation.svg'
            path.write_bytes(original)
            color_retained_ligation(path, color)
            updated = path.read_bytes()
            self.assertEqual(len(updated), len(original))
            self.assertEqual(ET.tostring(ET.fromstring(updated)), ET.tostring(expected))
            color_retained_ligation(path, color)
            self.assertEqual(path.read_bytes(), updated)
            # Some portable reports were recolored once before boxes gained stable ids.
            path.write_bytes(original)
            color_retained_ligation(path, '#bd805c')
            self.assertNotEqual(path.read_bytes(), original)
            color_retained_ligation(path, color)
            self.assertEqual(path.read_bytes(), updated)

    def test_retained_box_colors_accept_xml_quoting_namespaces_and_reordered_attributes(self):
        def chart(old, marked):
            return f"""<s:svg xmlns:s='http://www.w3.org/2000/svg'>
              <s:g id='axes_1'><s:g id='patch_4'><s:path clip-path='url(#a)' style='fill: #8293a4; stroke: #000000' d='M0,0 z'/></s:g></s:g>
              <s:g class="plot" id = 'axes_2'>
                <s:g data-note='quoted > and "id"' id = 'patch_14' class="legacy">
                  <s:path style = 'stroke: #000000; fill : {old}' d='M0,0 z' clip-path='url(#a)'/>
                </s:g>
                <s:g id="ligation-box-1"><s:path d='M0,0 z' style='fill: {marked}; stroke: #000000'/></s:g>
                <s:g id="patch_16"><s:path d='M 0 0 L 1 1 L 0 2 z' clip-path='url(#a)' style='fill: #bd805c; stroke: #000000'/></s:g>
                <s:g id='line2d_3'><s:path d="M0,0 L5,0" style='fill: none; stroke: #ff7f0e'/></s:g>
              </s:g></s:svg>""".encode()
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'detail.svg'
            path.write_bytes(chart('#8293a4', '#654321'))
            color_retained_ligation(path, '#123456')
            self.assertEqual(path.read_bytes(), chart('#123456', '#123456'))
            color_retained_ligation(path, '#123456')
            self.assertEqual(path.read_bytes(), chart('#123456', '#123456'))

    def test_ligation_insert_boxes_use_assigned_units_and_pool_unassigned_only_for_total(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            data = root / 'ligation'
            data.mkdir()
            units = [('lib01', '1', x) for x in (10, 20, 30, 40, 1000)]
            units += [('lib01', '2', 101), ('lib02', '1', 100), ('lib02', '1', 108),
                      ('unassigned', 'unresolved', 222)]
            stats = {}
            for sid in ('total', 'lib01', 'lib02', 'empty'):
                selected = [r for r in units if sid == 'total' or r[0] == sid]
                counts = {c: sum(r[1] == c for r in selected) for c in BLOCK_CATEGORIES}
                stats[sid] = dict(category_counts=counts, ligation_screened_units=len(selected))
            (data / 'x.json').write_text(json.dumps(dict(stats.pop('total'), scopes=stats)))
            table = data / 'x.units.tsv.gz'
            with gzip.open(table, 'wt') as stream:
                stream.write('scope_id\tcategory\tinsert_bp\n')
                for sid, category, value in units:
                    stream.write(f'{sid}\t{category}\t{value}\n')
            assets = root / 'assets'
            def render(ids):
                files = write_microc_comparisons('x', [{'scope_id': sid} for sid in ids], root, assets)
                chart = ET.parse(assets / files['ligation_lengths'])
                return chart.find("{*}g[@data-panel='insert-length-by-blocks']")
            def statistics_by_box(panel):
                return {(b.get('data-scope'), b.get('data-category')):
                        tuple(float(b.get('data-' + k)) for k in ('count', 'whislo', 'q1', 'med', 'q3', 'whishi'))
                        for b in panel.findall("{*}g[@data-boxplot='insert-length']")}
            panel = render(['lib01', 'lib02', 'empty'])
            split = statistics_by_box(panel)
            centers = {}
            for box in panel.findall("{*}g[@data-boxplot='insert-length']"):
                sid = box.get('data-scope')
                rect = box.find('{*}rect')
                whiskers, median = box.findall('{*}path')
                key = panel.find(f"{{*}}rect[@data-library-key='{sid}']")
                self.assertEqual(rect.get('fill'), scope_color(sid, 0))
                self.assertEqual(rect.get('fill'), key.get('fill'))
                self.assertEqual(rect.get('stroke'), '#000000')
                self.assertEqual(whiskers.get('stroke'), '#000000')
                self.assertEqual(median.get('stroke'), '#ff7f0e')
                centers[sid, box.get('data-category')] = float(rect.get('x')) + float(rect.get('width')) / 2
            self.assertLess(centers['lib01', '1'], centers['lib02', '1'])
            second_tick = next(t for t in panel.findall('{*}text') if t.text == '2')
            self.assertLess(centers['lib01', '2'], float(second_tick.get('x')))
            self.assertEqual(split, {('lib01', '1'): (5, 10, 20, 30, 40, 40),
                                     ('lib01', '2'): (1, 101, 101, 101, 101, 101),
                                     ('lib02', '1'): (2, 100, 102, 104, 106, 108)})
            pooled = statistics_by_box(render(['total']))
            self.assertEqual(set(pooled), {('total', '1'), ('total', '2'), ('total', 'unresolved')})
            self.assertEqual(pooled['total', '1'], (7, 10, 25, 40, 104, 108))
            self.assertEqual(pooled['total', 'unresolved'], (1, 222, 222, 222, 222, 222))
            files = write_microc_comparisons('x', [{'scope_id': 'empty'}], root, assets)
            empty_lengths = ET.parse(assets / files['ligation_lengths']).getroot()
            self.assertIn('No units', ''.join(empty_lengths.itertext()))
            self.assertIsNone(empty_lengths.find(".//{*}g[@data-boxplot]"))

    def test_native_contact_denominators_zero_and_distance_threshold(self):
        self.assertIsNone(native_contact_summary([{'scope_id': 'total'}]))
        library = dict(scope_id='lib01', unique_valid_contacts=20, mappable_library_units=40,
                       valid_contact_rate=.5, cis_contacts=15, cis_ratio=.75,
                       long_range_cis_contacts=10, long_range_cis_fraction=.5,
                       long_range_cis_min_bp=12000)
        empty = dict(scope_id='lib02', unique_valid_contacts=0, mappable_library_units=14,
                     valid_contact_rate=0, cis_contacts=0, long_range_cis_contacts=0,
                     long_range_cis_min_bp=12000)
        chart = ET.fromstring(native_contact_summary([library, empty]))
        self.assertIn('Long-range cis ≥12,000 bp.', ''.join(chart.itertext()))
        rows = [[node.text for node in row.findall('{*}text')] for row in chart.findall('{*}g')]
        self.assertEqual(rows[0], ['lib01', '20', '50.0%', '20 / 40', '75.0%', '15 / 20',
                                   '50.0%', '10 / 20'])
        self.assertEqual(rows[1], ['lib02', '0', '0.0%', '0 / 14', '—', '0 / 0', '—', '0 / 0'])
        empty['long_range_cis_min_bp'] = None
        chart = ET.fromstring(native_contact_summary([library, empty]))
        rows = [list(row.itertext()) for row in chart.findall('{*}g')]
        self.assertIn('≥12,000 bp', rows[0])
        self.assertNotIn('≥12,000 bp', rows[1])

    def test_pooled_fragments_include_unassigned_and_off_scale_observations(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'ligation').mkdir()
            with gzip.open(root / 'ligation/a.v2.mnase_fragments.tsv.gz', 'wt') as stream:
                stream.write('scope_id\treference_bp\nlib01\t100\nlib01\t199\n'
                             'lib02\t2000\nunassigned\t200\n')
            scopes = [{'scope_id': sid} for sid in ['total', 'lib01', 'lib02', 'lib03']]
            assets = root / 'assets'
            files = write_microc_comparisons('a.v2', scopes, root, assets)
            self.assertEqual(set(files), {'mnase'})
            chart = ET.parse(assets / files['mnase'])
            expected = [('Total library', '4', '199.5', '50.0%', '1'),
                        ('lib01', '2', '149.5', '100.0%', '0'),
                        ('lib02', '1', '2,000', '0.0%', '1'),
                        ('lib03', '0', '—', '—', '0')]
            for scope, values in zip(scopes, expected):
                row = chart.find(f"{{*}}g[@data-scope='{scope['scope_id']}']")
                self.assertEqual(tuple(t.text for t in row.findall('{*}text')), values)
                self.assertEqual(row.find('{*}polyline') is not None,
                                 scope['scope_id'] in {'total', 'lib01'})

    def test_categories_use_their_own_denominators_and_retain_missing_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            scopes = [dict(scope_id='lib01', scope_name='A & B'), dict(scope_id='lib02', ligation_screened_units=0),
                      dict(scope_id='custom_sublibrary_id_03', scope_name='Long name'),
                      dict(scope_id='lib04', ligation_screened_units=10)]
            (root / 'pairtools').mkdir()
            pair = dict(zip(PAIR_CATEGORIES, [1, 1, 1, 0, 1]))
            zero = dict.fromkeys(PAIR_CATEGORIES, 0)
            results = [dict(scope_id='lib01', parsed_pairs=4, after={'counts': pair}),
                       dict(scope_id='lib02', parsed_pairs=0, after={'counts': zero})]
            (root / 'pairtools/x.json').write_text(json.dumps({'scopes': results}))
            (root / 'ligation').mkdir()
            counts = dict(zip(BLOCK_CATEGORIES, [1, 2, 0, 0, 0, 3]))
            ligation = dict(ligation_screened_units=6, category_counts=counts)
            (root / 'ligation/x.json').write_text(json.dumps(dict(ligation, scopes={'lib01': ligation})))
            assets = root / 'assets'
            files = write_microc_comparisons('x', scopes, root, assets)
            for key, denominator, stat in [('pairtools', '4', '1.00'), ('ligation', '6', '33.3%')]:
                chart = ET.parse(assets / files[key])
                row = chart.find("{*}g[@data-scope='lib01']")
                self.assertEqual(row.find('{*}title').text, 'A & B')
                self.assertEqual([t.text for t in row.findall('{*}text')][1:3], [denominator, stat])
                bars = row.findall('{*}rect')
                percentage_labels = [t for t in row.findall('{*}text')[3:] if t.text.endswith('%')]
                self.assertEqual([t.get('fill') for t in percentage_labels],
                                 [bar_ink(r.get('fill')) for r in bars if float(r.get('width')) >= 49])
                self.assertEqual(sum(int(r.get('data-count')) for r in bars), int(denominator))
                width = sum(float(r.get('width')) for r in bars)
                for bar in bars:
                    self.assertAlmostEqual(float(bar.get('width')) / width,
                                           int(bar.get('data-count')) / int(denominator), places=5)
                missing = chart.find("{*}g[@data-scope='custom_sublibrary_id_03']")
                self.assertEqual(missing.find('{*}text').text, 'custom…id_03')
                self.assertIn('custom_sublibrary_id_03', missing.find('{*}title').text)
                self.assertIn('Unavailable', [t.text for t in missing.findall('{*}text')])
                empty = chart.find("{*}g[@data-scope='lib02']")
                self.assertEqual([t.text for t in empty.findall('{*}text')],
                                 ['lib02', '0', '—', 'No observations'])
                incomplete = chart.find("{*}g[@data-scope='lib04']")
                self.assertIn('Unavailable', [t.text for t in incomplete.findall('{*}text')])


if __name__ == '__main__':
    unittest.main()
