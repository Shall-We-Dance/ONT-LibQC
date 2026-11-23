import sys
import unittest
from html.parser import HTMLParser
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from report_palette import UNASSIGNED_COLOR, scope_color
from scrna_composition import sublibrary_composition


class Chart(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.segments = []
        self.labels = []
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'mixture-segment' in attrs.get('class', '').split():
            self.segments.append(attrs)
        if attrs.get('role') == 'img':
            self.labels.append(attrs['aria-label'])


class SublibraryComposition(unittest.TestCase):
    def test_index_share_uses_all_units_including_dimers(self):
        summary = {'physical_library': {
            'library_units': 110, 'sub_libraries': [{'scope_id': 'lib01', 'library_units': 100}],
            'index_pairing': {'unassigned_units': 10, 'ambiguous_units': 0}},
            'non_dimer_units': 70, 'scopes': {'lib01': {'units': 60}}}
        markup = sublibrary_composition(summary)
        chart = Chart(markup)
        self.assertEqual(chart.labels, ['lib01: 90.9%, 100 units; Unassigned: 9.1%, 10 units'])
        self.assertIn('mixture-deferred', chart.segments[1]['class'])
        self.assertIn('Unassigned <strong>9.1%</strong>', markup.split('<details class="mixture-other">')[1])
        self.assertIn(f'background:{scope_color("lib01")};color:#ffffff', chart.segments[0]['style'])
        self.assertIn(f'background:{UNASSIGNED_COLOR};color:#000000', chart.segments[1]['style'])

    def test_native_and_ambiguous_are_distinct(self):
        physical = {'library_units': 9,
                    'sub_libraries': [{'scope_name': 'Native unindexed', 'scope_id': 'native-unindexed', 'library_units': 7}],
                    'index_pairing': {'ambiguous_units': 1, 'unassigned_units': 1}}
        rows = Chart(sublibrary_composition({'physical_library': physical})).segments
        self.assertEqual([r['title'] for r in rows],
                         ['Native unindexed: 77.8% (7 units)', 'Ambiguous: 11.1% (1 unit)', 'Unassigned: 11.1% (1 unit)'])
        self.assertEqual(sublibrary_composition({'scopes': {'lib01': {'units': 65}}}), '')
        self.assertEqual(sublibrary_composition({'physical_library': {'library_units': 0}}), '')
