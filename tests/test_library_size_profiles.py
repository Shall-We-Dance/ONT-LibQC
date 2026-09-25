"""Validate physical length models, weighting and the default report disclosure."""
import csv
import gzip
from html.parser import HTMLParser
import json
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).parents[1] / 'workflow/scripts'))
import ont_qc as qc
from qc_plot_metrics import size_trace_data, svg_size_trace


def unit(name, insert, scope='lib01', family='truseq', assigned=True):
    return dict(unit_id=name, read_id=name.split('|')[0], insert_length=str(insert),
                scope_id=scope, assignment_status='assigned' if assigned else 'unassigned',
                family=family, p5_detected='0', p7_detected='0', raw_length='900', mean_q='25')


def summary(pairs=None, length=8):
    return dict(i5={'length': length}, i7={'length': length},
                index_pairing={'pairs': pairs or []})


class LengthWeightingTests(unittest.TestCase):
    def test_equal_molecule_counts_give_twice_the_weight_to_twice_the_length(self):
        values = [300] * 100 + [600] * 100
        for weighting, expected in [('count', 1), ('length', 2)]:
            data = size_trace_data(values, weighting=weighting)
            short = sum(y for x, y in zip(data['lengths'], data['signal']) if 225 < x < 399)
            long = sum(y for x, y in zip(data['lengths'], data['signal']) if 450 < x < 850)
            self.assertAlmostEqual(long / short, expected, places=10)
            self.assertEqual(data['n'], 200)
            self.assertEqual(data['shown'], 200)
        self.assertEqual(data['shown_weight'], 90000)

    def test_length_weighted_shape_is_invariant_to_uniform_replication(self):
        values = [136, 250, 300, 300, 600, 2100]
        a = size_trace_data(values, weighting='length')
        b = size_trace_data(values * 7, weighting='length')
        for x, y in zip(a['signal'], b['signal']):
            self.assertAlmostEqual(x, y, places=11)
        self.assertEqual(b['above_range'], 7)

    def test_long_outliers_do_not_change_the_visible_weighted_curve(self):
        a = size_trace_data([300] * 10, weighting='length')
        b = size_trace_data([300] * 10 + [1000000], weighting='length')
        self.assertEqual(a['signal'], b['signal'])
        self.assertEqual(b['shown'], 10)
        self.assertEqual(b['above_range_weight'], 1000000)

    def test_zero_length_has_no_mass_but_remains_an_observation(self):
        data = size_trace_data([0], weighting='length')
        self.assertEqual(data['n'], 1)
        self.assertEqual(data['weight_sum'], 0)
        self.assertFalse(any(data['signal']))
        self.assertEqual(size_trace_data([0])['signal'][0], 100)

    def test_svg_lane_identifies_the_same_weighting_as_the_curve(self):
        for weighting, label in [('count', 'count-based'), ('length', 'length-weighted')]:
            root = ET.fromstring(svg_size_trace([136, 300, 600], 'Test <library>', weighting=weighting))
            self.assertEqual(root.attrib['data-weighting'], weighting)
            text = ''.join(root.itertext())
            self.assertIn('darkness follows the ' + label, text)
            self.assertIn('Not measured TapeStation fluorescence', text)


class FullLibraryModelTests(unittest.TestCase):
    def test_model_matches_synthetic_complete_library_with_nonstandard_indexes(self):
        for family, (r1, r2) in qc.FAMILIES.items():
            i5, i7, insert = 'A' * 10, 'C' * 12, 'ACGT' * 50
            molecule = qc.P5 + i5 + r1 + insert + qc.rc(r2) + i7 + qc.rc(qc.P7)
            s = summary([dict(scope_id='lib01', i5=i5, i7=i7)])
            p = qc.full_library_profile([unit('r|u1', len(insert), family=family)], s)
            self.assertEqual(p['lengths'], [len(molecule)])
            self.assertEqual(p['adapter_offsets_bp'], [142])
            self.assertEqual(p['outer_incomplete_modeled_units'], 1)

    def test_dimer_and_multiplexed_indexes_keep_their_actual_lengths(self):
        rows = [unit('r|u1', 0, 'lib01'), unit('r|u2', 50, 'lib02')]
        s = summary([dict(scope_id='lib01', i5='A' * 8, i7='C' * 8),
                     dict(scope_id='lib02', i5='G' * 10, i7='T' * 10)])
        p = qc.full_library_profile(rows, s)
        self.assertEqual(p['lengths'], [136, 190])
        self.assertEqual(p['modeled_units'], 2)
        self.assertEqual(p['adapter_offsets_bp'], [136, 140])

    def test_unassigned_units_use_sample_length_consensus_without_a_fake_assignment(self):
        row = unit('r|u1', 100, assigned=False)
        s = summary([dict(scope_id='lib01', i5='A' * 12, i7='C' * 12)])
        p = qc.full_library_profile([row], s)
        self.assertEqual(p['lengths'], [236])
        self.assertEqual(p['rows'][0]['model_basis'], 'sample_index_length_consensus')
        self.assertEqual(row['assignment_status'], 'unassigned')

    def test_unknown_indexes_are_excluded_instead_of_assuming_136bp(self):
        p = qc.full_library_profile([unit('r|u1', 100)], summary(length=None))
        self.assertEqual(p['lengths'], [])
        self.assertEqual(p['source_units'], 1)
        self.assertEqual(p['exclusion_reasons'], {'unknown_index_length': 1})

    def test_total_includes_unassigned_units_while_scopes_keep_assigned_units(self):
        rows = [unit('r1|u1', 100), unit('r1|u2', 200), unit('r2|u1', 300, assigned=False)]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'reads.tsv.gz'
            with gzip.open(path, 'wt') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
                writer.writeheader()
                writer.writerows(rows)
            old = qc.read_scope_profiles(path)
            p = qc.read_scope_profiles(path, include_total=True)
        self.assertNotIn('total', old)
        self.assertEqual(p['total']['insert_lengths'], [100, 200, 300])
        self.assertEqual(len(p['total']['raw_lengths']), 2)
        self.assertEqual(p['lib01']['insert_lengths'], [100, 200])
        self.assertEqual(len(p['lib01']['raw_lengths']), 1)


class ProfileParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.details = []
        self.images = []
        self.count_panels = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'details':
            self.details.append(attrs)
            if attrs.get('class') == 'count-based-details':
                self.count_panels.append(attrs)
        if tag == 'img':
            self.images.append((attrs, list(self.details)))

    def handle_endtag(self, tag):
        if tag == 'details':
            self.details.pop()


class LibraryReportTests(unittest.TestCase):
    def test_main_plot_and_collapsed_counts_share_exported_unit_lengths(self):
        with tempfile.TemporaryDirectory() as tmp:
            assets = qc.write_library_size_assets(tmp, 's.total', 'Total library',
                [unit('r1|u1', 164), unit('r2|u1', 464)], summary())
            model = json.loads((Path(tmp) / assets['model']).read_text())
            self.assertEqual(model['modeled_units'], 2)
            self.assertEqual(model['shown_model_bases'], 900)
            with (Path(tmp) / assets['profiles']).open() as handle:
                data = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(len(data), 650)
            self.assertNotEqual([d['length_weighted_relative_signal'] for d in data],
                                [d['count_based_relative_signal'] for d in data])
            page = qc.profile_scope_html('Total library', 'raw.svg', 'insert.svg',
                'region.svg', 'karyotype.svg', full_library_assets=assets)
        parser = ProfileParser()
        parser.feed(page)
        self.assertEqual(len(parser.count_panels), 1)
        self.assertNotIn('open', parser.count_panels[0])
        weighted = [(img, parents) for img, parents in parser.images if 'length-weighted' in img['alt']]
        count = [(img, parents) for img, parents in parser.images if 'count-based' in img['alt']]
        self.assertEqual(len(weighted), 1)
        self.assertFalse(weighted[0][1])
        self.assertEqual(len(count), 3)
        self.assertTrue(all(any(p.get('class') == 'count-based-details' for p in parents) for _, parents in count))

    def test_unavailable_model_still_emits_valid_empty_profile_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            assets = qc.write_library_size_assets(tmp, 's.total', 'Total library',
                [unit('r|u1', 164)], summary(length=None))
            for key in ('weighted', 'count'):
                root = ET.fromstring((Path(tmp) / assets[key]).read_text())
                self.assertIn('no observations', ''.join(root.itertext()))
            self.assertEqual(assets['model_summary']['excluded_units'], 1)


if __name__ == '__main__':
    unittest.main()
