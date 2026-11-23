import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from scrna_families import analyze, family_counts, repeat_loci
from scrna_report import family_plot, family_section


class TestSingleCellRNAFamilies(unittest.TestCase):
    def test_family_gtf_with_bam_header_list_and_missing_mapper(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            gtf = work / 'families.gtf'
            gtf.write_text('chr1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "a"; gene_biotype "rRNA";\n')
            with gzip.open(work / 'cdna.fq.gz', 'wt') as f:
                f.write('@unit\n' + 'A' * 30 + '\n+\n' + 'I' * 30 + '\n')
            args = SimpleNamespace(rna_family_gtf=gtf, rna_source_panel=None, minimap2=str(work / 'missing'),
                                   repeatmasker=None, threads=1, mapq=20)
            result = analyze(args, work, [], {}, ['chr1'], {'unit': {'scope_id': 'lib01'}}, {})
            self.assertEqual(result['panel']['status'], 'not_available')
            self.assertEqual(result['cohorts']['all'], {'units': 1, 'counts': {'unknown': 1}})
            with self.assertRaisesRegex(ValueError, 'reference names'):
                analyze(args, work, [], {}, ['chr2'], {}, {})

    def test_intersection_uses_aligned_blocks_and_union_and_retains_missing_hits(self):
        def row(chrom, a, b, cls):
            return '\t'.join(['0'] * 5 + [chrom, str(a), str(b), '0', '+', 'name', cls, 'family'])
        def rec(chrom, blocks, nh):
            return dict(chrom=chrom, blocks=blocks, nh=nh, supplementary=False)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'rmsk.txt'
            path.write_text('\n'.join([row('chr1', 0, 3, 'rRNA'), row('chr1', 1, 4, 'rRNA'),
                row('chr1', 8, 100, 'LINE'), row('chr1', 6, 8, 'Simple_repeat'), row('chr2', 0, 20, 'rRNA')]))
            alignments = {'a': [rec('chr1', [(0, 8), (100, 104)], 2), rec('chr2', [(0, 9)], 2)],
                          'b': [rec('chr1', [(5, 13)], 2)]}
            self.assertEqual(repeat_loci(path, alignments, {'chr1', 'chr2'}), {'a': (2, 1, 0), 'b': (2, 0, 1)})
            with self.assertRaisesRegex(ValueError, 'No shared chromosomes'):
                repeat_loci(path, alignments, {'other'})

    def test_family_priority_one_vote_per_fragment_and_separate_denominators(self):
        families = ['ambiguous_family', 'unknown', 'protein_coding', 'ambiguous_family', 'unknown', 'mt-rRNA']
        mapping = ['multimapped'] * 4 + ['unmapped', 'unique']
        rows = [dict(unit_id=str(i), scope_id='lib01', mapping_status=mapping[i], family=family)
                for i, family in enumerate(families)]
        hits = [dict(unit_id=str(i), reference_name='ref' if i in (0, 5) else '', matched_fraction='', identity='',
                     status='assigned' if i in (0, 5) else 'unassigned', source_class='rRNA' if i in (0, 5) else '')
                for i in range(6)]
        repeats = {'1': (2, 2, 0), '2': (2, 1, 0), '3': (3, 0, 3)}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'units.tsv.gz'
            cohorts, evidence = family_counts(rows, hits, {'ref': {}}, repeats, path)
            self.assertEqual(cohorts['all']['units'], 6)
            self.assertEqual(cohorts['multimapped'], {'units': 4, 'counts': {'rRNA': 2, 'protein_coding': 1, 'ambiguous_family': 1}})
            self.assertEqual(cohorts['all']['counts']['mt-rRNA'], 1)
            self.assertEqual(evidence, {'sequence_rRNA': 1, 'repeat_rRNA_all_loci': 1, 'rRNA_some_loci': 1, 'other_repeat_all_loci': 1})
            with gzip.open(path, 'rt') as f:
                self.assertEqual(len(list(csv.DictReader(f, delimiter='\t'))), 6)
            self.assertEqual(rows[0]['family'], 'ambiguous_family')
            with self.assertRaisesRegex(ValueError, 'different input fragments'):
                family_counts(rows, hits[:-1], {}, {}, path)

    def test_no_multimappers_and_unavailable_independent_screen(self):
        cohorts = {'all': {'units': 1, 'counts': {'unknown': 1}}, 'multimapped': {'units': 0, 'counts': {}}}
        data = dict(cohorts=cohorts, panel={'status': 'not_available'}, repeatmasker_sha256=None,
                    multimapping_evidence={})
        svg = family_plot(cohorts)
        self.assertIn(b'Unclassified', svg)
        html = family_section(data, svg)
        self.assertIn('independent rRNA sequence screen was unavailable', html)
        self.assertNotIn('Evidence for multimapping', html)
        self.assertIsNone(family_plot(dict.fromkeys(cohorts, {'units': 0, 'counts': {}})))


if __name__ == '__main__':
    unittest.main()
