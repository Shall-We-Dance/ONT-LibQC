"""GRIP preparation regressions: independent scopes, strand, and observed tag boundaries."""
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
import rna_prepare as prep


def tag(number, length):
    return ''.join('ACGT'[(number >> (2 * i)) & 3] for i in range(length))


def molecule(i, insert='GATTACAGACTGACTGACTGACTGAC', barcode='CCTATAT', prefix=None):
    return (prefix if prefix is not None else tag(i, 5)) + barcode + insert + tag(12345 + i, 10)


class RNAPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def run_prepare(self, entries, **kwargs):
        fastq, table = self.root / 'in.fastq', self.root / 'table.tsv'
        with fastq.open('w') as reads, table.open('w') as rows:
            writer = csv.DictWriter(rows, delimiter='\t', fieldnames=['unit_id', 'scope_id', 'assignment_status', 'orientation', 'insert_length', 'umi_status', 'umi_sequence'])
            writer.writeheader()
            for name, normalized, scope, reverse, umi in entries:
                seq = prep.reverse_complement(normalized) if reverse else normalized
                quality = ''.join(chr(40 + (i % 30)) for i in range(len(seq)))
                reads.write(f'@{name}\n{seq}\n+\n{quality}\n')
                writer.writerow(dict(unit_id=name, scope_id=scope, assignment_status='assigned' if scope else 'unassigned',
                    orientation='p7_to_p5' if reverse else 'p5_to_p7', insert_length=len(seq),
                    umi_status='complete' if umi else 'not_configured', umi_sequence=umi))
        counts = prep.prepare(fastq, table, self.root / 'out.fastq.gz', self.root / 'summary.json', units_tsv=self.root / 'units.tsv.gz', **kwargs)
        with gzip.open(self.root / 'units.tsv.gz', 'rt') as handle:
            units = {r['unit_id']: r for r in csv.DictReader(handle, delimiter='\t')}
        out = {name: (seq, quality) for name, seq, quality in prep.fastq_records(self.root / 'out.fastq.gz')}
        return counts, units, out

    def test_supported_scope_orientation_and_observed_umi(self):
        entries = [(f'u{i}', molecule(i), 'positive', bool(i % 2), '') for i in range(24)]
        counts, units, out = self.run_prepare(entries)
        self.assertTrue(counts['auto_grip']['scope_evidence']['positive']['supported'])
        self.assertEqual(counts['umi_complete_units'], 24)
        for i in range(24):
            row = units[f'u{i}']
            self.assertEqual(row['p5_random_tag'], tag(i, 5))
            self.assertEqual(row['p7_umi'], tag(12345 + i, 10))
            self.assertEqual(row['umi_sequence'], row['p7_umi'])
            self.assertEqual(row['p5_trim_bp'], '12')
            self.assertEqual(row['p7_trim_bp'], '10')
            biological = molecule(i)[12:-10]
            self.assertEqual(out[f'u{i}'][0], prep.reverse_complement(biological) if i % 2 else biological)
            raw_quality = ''.join(chr(40 + (j % 30)) for j in range(len(molecule(i))))
            self.assertEqual(out[f'u{i}'][1], raw_quality[10:-12] if i % 2 else raw_quality[12:-10])

    def test_mixed_scopes_do_not_propagate_and_uncertain_units_are_preserved(self):
        entries = [(f'p{i}', molecule(i), 'positive', False, '') for i in range(20)]
        entries += [('uncertain', molecule(44, barcode='GGGGGGG'), 'positive', False, '')]
        entries += [(f'n{i}', molecule(i, barcode='GGGGGGG'), 'negative', False, '') for i in range(25)]
        entries += [(f'x{i}', molecule(i), '', False, '') for i in range(10)]
        counts, units, out = self.run_prepare(entries)
        self.assertTrue(counts['auto_grip']['scope_evidence']['positive']['supported'])
        self.assertFalse(counts['auto_grip']['scope_evidence']['negative']['supported'])
        self.assertFalse(counts['auto_grip']['scope_evidence']['unassigned']['supported'])
        self.assertEqual(units['uncertain']['status'], 'uncertain_barcode')
        self.assertEqual(units['uncertain']['umi_status'], 'structure_uncertain')
        self.assertEqual(out['uncertain'][0], entries[20][1])
        self.assertEqual(out['n0'][0], molecule(0, barcode='GGGGGGG'))
        self.assertEqual(out['x0'][0], molecule(0))
        self.assertEqual(counts['umi_complete_units'], 20)

    def test_shifted_and_error_tag_cuts_at_observed_barcode_end(self):
        entries = [(f'p{i}', molecule(i), 'positive', False, '') for i in range(20)]
        entries += [('shift', molecule(30, prefix='ACGTAC'), 'positive', True, ''),
                    ('mismatch', molecule(31, barcode='CCTACAT'), 'positive', False, ''),
                    ('deletion', molecule(32, barcode='CCTTAT'), 'positive', False, '')]
        _, units, out = self.run_prepare(entries)
        self.assertEqual(units['shift']['p5_trim_bp'], '13')
        self.assertEqual(out['shift'][0], prep.reverse_complement(molecule(0)[12:-10]))
        self.assertEqual(units['mismatch']['barcode_edits'], '1')
        self.assertEqual(units['deletion']['p5_trim_bp'], '11')
        self.assertEqual(out['deletion'][0], molecule(0)[12:-10])

    def test_clonal_fixed_tags_do_not_establish_structure(self):
        entries = [(f'u{i}', molecule(0), 'one', False, '') for i in range(25)]
        counts, _, out = self.run_prepare(entries)
        ev = counts['auto_grip']['scope_evidence']['one']
        self.assertEqual(ev['decision'], 'insufficient_random_tag_diversity')
        self.assertEqual(len(out['u0'][0]), len(molecule(0)))
        self.assertEqual(counts['umi_complete_units'], 0)

    def test_custom_umi_or_manual_trim_bypasses_auto_and_keeps_custom_tag(self):
        entries = [(f'u{i}', molecule(i), 'one', False, 'ACTG') for i in range(20)]
        counts, units, out = self.run_prepare(entries)
        self.assertEqual(counts['auto_grip']['bypass_reason'], 'configured_umi_extraction')
        self.assertEqual(units['u0']['umi_sequence'], 'ACTG')
        self.assertEqual(out['u0'][0], molecule(0))
        counts, units, out = self.run_prepare([(n, s, c, r, '') for n, s, c, r, _ in entries], p5=2, p7=3)
        self.assertEqual(counts['auto_grip']['bypass_reason'], 'manual_end_trimming')
        self.assertEqual(out['u0'][0], molecule(0)[2:-3])
        self.assertEqual(counts['umi_complete_units'], 0)
        counts, _, _ = self.run_prepare([], structure_json='{"umis":[{"position":"before_insert","length":5}]}')
        self.assertEqual(counts['auto_grip']['bypass_reason'], 'configured_umi_extraction')

    def test_empty_analysis_and_empty_input_are_audited(self):
        entries = [(f'p{i}', molecule(i), 'positive', False, '') for i in range(20)]
        entries += [('empty', molecule(30, insert=''), 'positive', False, '')]
        counts, units, out = self.run_prepare(entries)
        self.assertEqual(counts['empty_after_trim'], 1)
        self.assertEqual(counts['umi_complete_units'], 21)
        self.assertNotIn('empty', out)
        self.assertEqual(units['empty']['analysis_length'], '0')
        self.assertEqual(units['empty']['umi_status'], 'complete')
        counts, units, out = self.run_prepare([])
        self.assertEqual(counts['input_units'], 0)
        self.assertEqual(units, {})
        self.assertEqual(out, {})

    def test_invalid_umi_in_supported_scope_is_not_labeled_complete(self):
        entries = [(f'p{i}', molecule(i), 'positive', False, '') for i in range(20)]
        entries += [('invalid', molecule(30)[:-1] + 'N', 'positive', False, '')]
        _, units, out = self.run_prepare(entries)
        self.assertEqual(units['invalid']['status'], 'uncertain_random_tag')
        self.assertEqual(units['invalid']['umi_sequence'], '')
        self.assertEqual(out['invalid'][0], entries[-1][1])

if __name__ == '__main__':
    unittest.main()
