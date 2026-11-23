import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from scrna_species import analyze, make_queries
from species_composition import render


class SingleCellSpeciesTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.args = SimpleNamespace(star_index=self.root / 'star', species_panel=None, assembly='danRer11',
                                    bowtie2='bowtie2', threads=2)

    def panel(self):
        for name in ('fish', 'human'):
            for part in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
                (self.root / f'{name}.{part}.bt2').touch()
        cfg = {'indexes': {'GRCz12ab': 'fish', 'hg38': 'human'}, 'assembly_aliases': {'danRer11': 'GRCz12ab'}}
        (self.root / 'species-panel.json').write_text(json.dumps(cfg))

    def inputs(self):
        names = ['u3', 'u1', 'u4', 'u2']
        with gzip.open(self.root / 'cdna.fq.gz', 'wt') as out:
            for name in names:
                out.write(f"@{name}\n{'A' * 60}\n+\n{'I' * 60}\n")
        # Deliberately different ordering between files.
        with gzip.open(self.root / 'rna_families.tsv.gz', 'wt') as out:
            out.write('unit_id\tfamily\nu1\trRNA\nu2\tprotein_coding\nu3\tmt-rRNA\nu4\tunknown\n')
        return {n: {'scope_id': 'lib01' if n != 'u4' else 'unassigned'} for n in sorted(names)}

    def test_trimmed_cdna_counts_partition_and_non_rrna_control(self):
        self.panel()
        tagged = self.inputs()

        def mock_map(fasta, size, indexes, *args):
            self.assertEqual(fasta.read_text().splitlines()[::2], ['>0', '>1', '>2', '>3'])
            self.assertEqual(size, 4)
            self.assertEqual(indexes['hg38'], str(self.root / 'human'))
            return [[('GRCz12ab', 70, True, ('chr1', 1))],
                    [('GRCz12ab', 70, True, None), ('hg38', 69, True, None)],
                    [], [('hg38', 60, True, ('chr2', 3))]]

        with patch('scrna_species.align_hits', side_effect=mock_map), patch(
                'scrna_species.subprocess.check_output', return_value='bowtie2-align-s version 2.5.5'):
            result = analyze(self.args, self.root, tagged)
        self.assertEqual(result['mapper_version'], '2.5.5')
        rows = {r['scope']: r for r in result['rows']}
        self.assertEqual(set(rows), {'Total inserts', 'lib01', 'Unassigned index'})
        self.assertEqual(rows['Total inserts']['n'], len(tagged))
        self.assertEqual(sum(rows['lib01']['counts'].values()), 3)
        self.assertEqual(rows['Unassigned index']['counts'], {'No confident match': 1})
        self.assertEqual(result['non_rrna_row']['counts'], {'No confident match': 1, 'hg38': 1})
        self.assertFalse((self.root / 'species_queries.fa').exists())
        html = render(result)
        self.assertEqual(html.count("class='species-bar'"), 2)
        self.assertIn('Excluding rRNA and mt-rRNA', html)
        self.assertIn('Zebrafish (expected)', html)
        self.assertNotIn('Reads without a reconstructed insert', html)

    def test_absence_of_configuration_or_queries_is_explicit(self):
        missing = analyze(self.args, self.root, {})
        self.assertEqual(missing['status'], 'not_configured')
        self.assertIn('was not configured', render(missing))
        self.args.species_panel = self.root / 'requested.json'
        with self.assertRaises(FileNotFoundError):
            analyze(self.args, self.root, {})
        self.args.species_panel = None
        self.panel()
        empty = analyze(self.args, self.root, {})
        self.assertEqual(empty['status'], 'no_qualified_cdna')
        self.assertNotIn('class=\'species-bar\'', render(empty))
        self.args.assembly = 'dm6'
        with self.assertRaisesRegex(ValueError, 'expected species'):
            analyze(self.args, self.root, {})
        self.args.assembly = 'danRer11'
        (self.root / 'fish.3.bt2').unlink()
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            analyze(self.args, self.root, {})

    def test_rejects_mismatched_units(self):
        tagged = self.inputs()
        tagged['extra'] = {'scope_id': 'lib01'}
        with self.assertRaisesRegex(ValueError, 'missing qualified'):
            make_queries(self.root / 'cdna.fq.gz', tagged, self.root / 'species.fa')


if __name__ == '__main__':
    unittest.main()
