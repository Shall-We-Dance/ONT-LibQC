import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from scrna_library_report import library_identity, naming_assets


class LibraryIdentityTests(unittest.TestCase):
    def test_editor_lists_only_stable_physical_library_ids(self):
        summary = {'sample': 'embryos & cells', 'physical_indexes': {
            'lib01': {'i5': 'AA', 'i7': 'TT'}, 'unassigned': {}, 'native-unindexed': {}}}
        page = naming_assets(summary)
        self.assertIn('data-scope-id="lib01"', page)
        self.assertNotIn('data-scope-id="unassigned"', page)
        self.assertNotIn('data-scope-id="native-unindexed"', page)
        self.assertIn('embryos &amp; cells', page)
        self.assertIn('data-sublibrary="lib01"', library_identity(summary))
        self.assertEqual(naming_assets({'physical_indexes': {'unassigned': {}}}), '')

    def test_index_support_is_among_informative_units(self):
        summary = {'assembly': 'test-assembly', 'physical_library': {
            'library_units': 2000, 'sub_libraries': [{'scope_id': 'lib01', 'library_units': 1600}],
            'index_pairing': {'mode': 'single_library', 'top_pair_fraction': 1080 / 1125,
                              'single_library_fraction': .9, 'minimum_informative_units': 20,
                              'informative_units': 1125, 'pairs': [{
                                  'scope_id': 'lib01', 'name': '<library>', 'i5': 'AACCAACC', 'i7': 'TTGGTTGG',
                                  'adapter_family_counts': {'truseq': 1600}, 'informative_support': 1080,
                                  'support_fraction': 1080 / 1125}]}}}
        block = library_identity(summary)
        for text in ('test-assembly', 'scRNA-seq', 'Single Library', '96.0%', '1,080 / 1,125',
                     'AACCAACC', 'TTGGTTGG', 'TruSeq', '&lt;library&gt;'):
            self.assertIn(text, block)
        for text in ('<library>', 'Edit names', '80.0%'):
            self.assertNotIn(text, block)

    def test_older_summaries_keep_index_sequences_and_do_not_impute_shares(self):
        block = library_identity({'physical_indexes': {'lib01': {'i5': '', 'i7': 'TTGGTTGG'},
                                                         'lib02': {'i5': None, 'i7': 'CCTTCCTT'}}})
        for text in ('Not Recorded', 'No index', 'Not detected', 'TTGGTTGG', 'CCTTCCTT'):
            self.assertIn(text, block)
        self.assertNotIn('0.0%', block)
        self.assertIn('Physical indexes were not recorded', library_identity({}))

    def test_native_unindexed_units_do_not_create_an_inferred_pair(self):
        block = library_identity({'physical_library': {
            'library_units': 57, 'sub_libraries': [{'scope_id': 'native-unindexed', 'library_units': 57}],
            'index_pairing': {'mode': 'inconclusive_low_index_information', 'informative_units': 0,
                              'top_pair_fraction': None, 'pairs': []}}})
        self.assertIn('No confident index combination resolved', block)
        self.assertNotIn('0.0%', block)
