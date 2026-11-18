import sys
import unittest
import tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
import rna_report as report
from report_profiles import render_profiles_card

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

    def test_scope_overview_uses_local_denominators_and_links(self):
        def metrics(units, bases):
            return dict(units=units, region_denominator_bases=bases, family_counts={'rRNA': 2},
                        region_bases={'CDS': 2}, plots={'profiles': 'profile.svg'})
        scopes = {'total': metrics(20, 50), 'lib01': metrics(4, 10), 'lib02': metrics(4, 0)}
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            for scope in scopes:
                (folder/f's.{scope}.full_library_profiles.tsv').write_text(
                    'size_bp\tlength_weighted_relative_signal\n0\t0\n100\t100\n2000\t0\n')
            for ids in [('total', 'lib01', 'lib02'), ('total',)]:
                scope_rows = [dict(scope_id=s, scope_name='A & B' if s=='lib01' else s) for s in ids]
                overview = render_profiles_card('s', scope_rows, ['<p>Details</p>']*len(ids), folder, folder,
                                                rna_summary={'scopes': scopes}, rna_prefix='../data/rna/s')
                regions = (folder/'s.rna_region_distribution.svg').read_text()
                self.assertLess(overview.index('virtual_gel.svg'), overview.index('rna_region_distribution.svg'))
                self.assertLess(overview.index('rna_region_distribution.svg'), overview.index('<h3>Transcript metaplots'))
                self.assertIn('../data/rna/s/profile.svg', overview)
                self.assertEqual(overview.count(' transcript coverage'), max(1, len(ids)-1))
                if len(ids)>1:
                    self.assertIn('20.0%', regions)
                    self.assertIn('—', regions)
                    self.assertIn('50.0%', overview)
                    self.assertIn('A &amp; B', overview)
                else:
                    self.assertIn('4.0%', regions)

if __name__ == '__main__':
    unittest.main()
