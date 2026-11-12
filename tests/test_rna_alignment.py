"""Boundary checks for ONT RNA preparation and STAR input handling."""
import gzip
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from rna_prepare import trim_insert
from run_star_rna import prepare_fastq, star_command
from ont_qc import assay_reference_ranges


class RNAAlignmentTests(unittest.TestCase):
    def test_physical_end_trimming_is_invariant_to_ont_direction(self):
        seq = 'AACGTGATTACAATTT'
        qual = ''.join(chr(65 + i) for i in range(len(seq)))
        forward = trim_insert(seq, qual, 'p5_to_p7', 3, 2)
        reverse = trim_insert(seq[::-1], qual[::-1], 'p7_to_p5', 3, 2)
        self.assertEqual(forward, tuple(x[::-1] for x in reverse))
        self.assertEqual(forward, (seq[3:-2], qual[3:-2]))
        self.assertEqual(trim_insert('ACG', 'ABC', 'p5_to_p7', 2, 2), ('', ''))
        with self.assertRaises(ValueError):
            trim_insert(seq, qual, 'unknown', 3, 2)

    def test_variable_length_fastq_preserves_unit_names_and_quality(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / 'reads.fastq.gz', Path(tmp) / 'star.fastq'
            records = [('unit/origin.1', 'A' * 30, 'I' * 30), ('unit/origin.2', 'AC' * 4500, 'HI' * 4500)]
            with gzip.open(source, 'wt') as handle:
                for name, seq, quality in records:
                    handle.write(f'@{name} extra description\n{seq}\n+\n{quality}\n')
            self.assertEqual(prepare_fastq(source, target), (2, 9000))
            lines = target.read_text().splitlines()
            self.assertEqual(lines[0], '@unit/origin.1')
            self.assertEqual(lines[4], '@unit/origin.2')
            self.assertEqual(lines[7], records[1][2])
            command = star_command('STARlong', 'index', target, 'out.', 4, 9000)
            self.assertEqual(command[command.index('--readNameSeparator') + 1], ' ')
            self.assertEqual(command[command.index('--outSAMtype') + 2], 'Unsorted')

    def test_malformed_and_empty_fastq(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, target = Path(tmp) / 'reads.fastq', Path(tmp) / 'out.fastq'
            source.write_text('')
            self.assertEqual(prepare_fastq(source, target), (0, 0))
            source.write_text('@read\nACGT\n+\nIII\n')
            with self.assertRaises(ValueError):
                prepare_fastq(source, target)

    def test_rna_has_no_implicit_dna_mapping_threshold(self):
        self.assertIsNone(assay_reference_ranges({}, 'rnaseq')['primary_mapping_rate'])
        self.assertIsNotNone(assay_reference_ranges({}, 'chromatin')['primary_mapping_rate'])
        explicit = {'qc_reference_ranges': {'primary_mapping_rate': {'min': .4}}}
        self.assertEqual(assay_reference_ranges(explicit, 'clip')['primary_mapping_rate']['min'], .4)


if __name__ == '__main__':
    unittest.main()
