import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
import rna_report as report

class LayoutTests(unittest.TestCase):
    def test_source_precedence_and_conserved_counts(self):
        base = dict(family='rRNA', mapping_status='unique', region_bases='{"CDS": 12}')
        rows = [dict(base, source_status='assigned', source_taxon=taxon) for taxon in ('Homo sapiens', 'Escherichia coli', 'Homo sapiens,Mus musculus')]
        rows.append(dict(base, mapping_status='unmapped', source_status='unassigned', region_bases='{}'))
        families, regions = report.source_breakdown(rows, 'Homo sapiens')
        self.assertEqual(sum(families['rRNA'].values()), 4)
        self.assertEqual(sum(regions['CDS'].values()), 36)
        self.assertEqual(families['rRNA']['E. coli candidate'], 1)
        self.assertEqual(families['rRNA']['Uncertain / unassigned'], 2)
        self.assertEqual(report.source_category(rows[1], None), 'Uncertain / unassigned')

    def test_zero_categories_and_missing_denominator(self):
        text = report.evidence_bars({k:0 for k in report.FAMILIES}, 10)
        for family in report.FAMILIES:
            self.assertIn(family, text)
        self.assertEqual(text.count('>0.0%'), len(report.FAMILIES))
        empty = report.evidence_bars({'Spliced':0}, 0)
        self.assertIn('Not evaluable', empty)
        self.assertNotIn('0.0%', empty)

if __name__ == '__main__':
    unittest.main()
