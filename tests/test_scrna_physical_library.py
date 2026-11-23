import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from report_sublibrary import mixture_counts
from scrna_qc import physical_library_summary


class PhysicalLibrary(unittest.TestCase):
    def test_distinct_assignment_and_informative_denominators(self):
        pair = dict(scope_id='lib01', name='lib01', i5='', i7='AACG', support_units=3,
                    informative_support=2, support_fraction=2/3, unit_ids=['internal-read-id'])
        source = dict(mode='single_library', top_pair_fraction=2/3, single_library_fraction=.9,
                      minimum_informative_units=20, informative_units=3, ambiguous_units=1,
                      unassigned_units=1, pairs=[pair])
        adapter = dict(index_pairing=source, library_units=5, adapter_dimer_units=1)
        rows = [('lib01', 'assigned', 'truseq', '1'), ('lib01', 'assigned', 'truseq', '0'),
                ('lib01', 'assigned', 'nextera', '0'), ('', 'ambiguous', 'truseq', '0'),
                ('', 'unassigned', 'nextera', '0')]
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            with gzip.open(work / 'reads.tsv.gz', 'wt') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow(['scope_id', 'assignment_status', 'family', 'is_adapter_dimer'])
                writer.writerows(rows)
            result = physical_library_summary(adapter, work, native_mgi_units=2)
        self.assertEqual(result['adapter_dimer_units'], 1)
        recorded = result['index_pairing']['pairs'][0]
        self.assertEqual(recorded['adapter_family_counts'], {'truseq': 2, 'nextera': 1})
        self.assertEqual(recorded['support_fraction'], 2/3)
        self.assertNotIn('internal-read-id', json.dumps(result))
        self.assertEqual(recorded['i5'], '')
        total, composition = mixture_counts(result)
        self.assertEqual(total, 7)
        self.assertEqual([(label, n) for label, n, _ in composition],
                         [('lib01', 3), ('Native unindexed', 2), ('Ambiguous', 1), ('Unassigned', 1)])
