"""RNA source evidence must reject short hits and expose competing taxa."""
import argparse
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('rna_sources', ROOT / 'workflow/scripts/rna_sources.py')
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)


def paf(unit, target, qlen=100, span=100, matches=98, block=100, score=190):
    return '\t'.join(map(str, [unit, qlen, 0, span, '+', target, 1000, 0, span,
                              matches, block, 0, 'AS:i:' + str(score)])) + '\n'


class RNASourceTests(unittest.TestCase):
    def setUp(self):
        self.refs = {'human': {'source_class': 'rRNA', 'taxon': 'Homo sapiens'},
                     'mouse': {'source_class': 'rRNA', 'taxon': 'Mus musculus'},
                     'ecoli': {'source_class': 'rRNA', 'taxon': 'Escherichia coli'}}

    def test_short_exact_matches_do_not_classify_long_inserts(self):
        hits = source.parse_paf([paf('short', 'human', qlen=500, span=30, matches=30, block=30),
                                 paf('low_identity', 'human', matches=84),
                                 paf('tiny', 'human', qlen=25, span=25, matches=25, block=25)], self.refs)
        self.assertEqual(hits, {})

    def test_near_equal_competing_taxa_are_ambiguous(self):
        hits = source.parse_paf([paf('conserved', 'human', score=190),
                                 paf('conserved', 'mouse', matches=97, score=185),
                                 paf('conserved', 'ecoli', matches=86, score=160)], self.refs)
        row = hits['conserved']
        self.assertEqual(row['source_class'], 'rRNA')
        self.assertEqual(row['taxon_ambiguous'], 1)
        self.assertEqual(row['taxon'], 'Homo sapiens,Mus musculus')
        self.assertEqual(row['reference_name'], 'human')

    def test_threshold_boundary_and_indel_identity(self):
        rows = source.parse_paf([paf('boundary', 'human', qlen=100, span=80, matches=68, block=80),
                                 paf('indels', 'human', span=100, matches=90, block=110)], self.refs)
        self.assertIn('boundary', rows)
        self.assertNotIn('indels', rows)

    def test_bundled_panel_provenance_and_tamper_check(self):
        panel = source.load_panel(source.RESOURCE / 'panel.fa', source.RESOURCE / 'panel.json')
        self.assertEqual(panel['reference_count'], 66)
        self.assertEqual(len({v['taxon'] for v in panel['references'].values()}), 6)
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / 'bad.fa'
            bad.write_text('>invented\nACGT\n')
            with self.assertRaisesRegex(ValueError, 'provenance'):
                source.load_panel(bad, source.RESOURCE / 'panel.json')

    def test_unavailable_panel_preserves_denominator_and_all_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            fq = p / 'reads.fastq'
            fq.write_text('@one\nACGT\n+\nIIII\n@two\nTGCA\n+\nIIII\n')
            args = argparse.Namespace(input=str(fq), reference_fasta=str(p/'missing.fa'),
                                      reference_manifest=str(p/'missing.json'), minimap2='missing',
                                      output=str(p/'units.tsv.gz'), summary=str(p/'summary.json'), threads=1)
            summary = source.run(args)
            self.assertEqual(summary['input_units'], 2)
            self.assertEqual(summary['not_available_units'], 2)
            self.assertEqual(summary['assigned_units'], 0)
            with gzip.open(args.output, 'rt') as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual([r['unit_id'] for r in rows], ['one', 'two'])
            self.assertTrue(all(r['status'] == 'not_available' for r in rows))

    def test_all_input_units_remain_in_complete_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            assigned = source.parse_paf([paf('one', 'ecoli')], self.refs)
            summary = source.write_results(['one', 'two'], assigned, p/'units.tsv.gz', p/'summary.json')
            self.assertEqual(summary['assigned_units'], 1)
            self.assertEqual(summary['unassigned_units'], 1)
            self.assertEqual(summary['source_counts'], {'rRNA': 1})
            with gzip.open(p/'units.tsv.gz', 'rt') as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(rows[1]['status'], 'unassigned')
            self.assertEqual(rows[1]['source_class'], '')


if __name__ == '__main__':
    unittest.main()
