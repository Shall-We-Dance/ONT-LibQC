"""End-to-end fixed-length UMI extraction on both physical orientations."""
import contextlib
import csv
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
import ont_qc as qc
from library_structure import normalize_structure, POSITIONS, handles

I5, I7 = 'ACGTACGT', 'TGCATGCA'
INSERT = 'GATCACCGTTGACGATTCGG' * 5
UMI = 'GCTAACGATGACTCGAGTCAAGGTACCTAGCA'


def molecule(layout, family='truseq', insert=INSERT, umi_values=None):
    lengths = {u['position']: u['length'] for u in layout['umis']}
    values = {p: (umi_values or {}).get(p, UMI[:n]) for p, n in lengths.items()}
    r1, r2 = qc.FAMILIES[family]; p5, p7 = handles(layout)
    blocks = [p5, I5, r1, insert, qc.rc(r2), I7, p7]
    sequence = ''.join(block + values.get(POSITIONS[i], '') if i < len(POSITIONS) else block for i, block in enumerate(blocks))
    return sequence, values


def prepare(root, sequences, layout=None, assay='chromatin'):
    path = root / 'input.fastq'
    qualities = [''.join(chr(35 + i % 40) for i in range(len(seq))) for seq in sequences]
    path.write_text(''.join(f'@r{i}\n{s}\n+\n{q}\n' for i, (s, q) in enumerate(zip(sequences, qualities))))
    command = ['prepare', '--sample', 'test', '--assay', assay, '--fastq', str(path), '--minimum-index-cluster-units', '1', '--minimum-index-informative-units', '1']
    for key in ['inserts', 'legs', 'leg1', 'leg2', 'read-table', 'index-assignments', 'summary', 'raw-lengths', 'insert-lengths']:
        command += ['--' + key, str(root / (key + ('.json' if key == 'summary' else '.gz')))]
    if layout is not None: command += ['--library-structure', json.dumps(layout)]
    with contextlib.redirect_stdout(io.StringIO()): qc.cmd_prepare(qc.parser().parse_args(command))
    summary = json.loads((root / 'summary.json').read_text())
    with gzip.open(root / 'read-table.gz', 'rt') as stream: rows = list(csv.DictReader(stream, delimiter='\t'))
    inserts = list(qc.fastq_records([root / 'inserts.gz']))
    return summary, rows, inserts, qualities


class UMIExtractionTests(unittest.TestCase):
    def test_every_position_both_strands_and_adapter_families(self):
        for family in qc.FAMILIES:
            for position in POSITIONS:
                with self.subTest(family=family, position=position), tempfile.TemporaryDirectory() as tmp:
                    layout = normalize_structure({'mode': 'custom', 'umis': [{'position': position, 'length': 12}]})
                    seq, umis = molecule(layout, family)
                    summary, rows, inserts, quals = prepare(Path(tmp), [seq, qc.rc(seq)], layout)
                    self.assertEqual(summary['library_units'], 2)
                    self.assertEqual(summary['umi_summary']['complete_units'], 2)
                    self.assertEqual(summary['index_pairing']['mode'], 'single_library')
                    self.assertEqual(summary['index_pairing']['pairs'][0]['i5'], I5)
                    self.assertEqual(summary['index_pairing']['pairs'][0]['i7'], I7)
                    self.assertEqual([r[1] for r in inserts], [INSERT, qc.rc(INSERT)])
                    for i, row in enumerate(rows):
                        self.assertEqual(row['umi_sequence'], umis[position])
                        self.assertEqual(row['umi_' + position], umis[position])
                        self.assertEqual(row['umi_status'], 'complete')
                        self.assertEqual(inserts[i][2], quals[i][int(row['insert_start']):int(row['insert_end'])])
                    model = qc.full_library_profile(rows, summary)
                    self.assertEqual(model['lengths'], [len(seq), len(seq)])
                    self.assertEqual([row['model_umi_bases'] for row in model['rows']], [12, 12])

    def test_all_six_umis_and_concatemers_keep_unit_sequences_separate(self):
        layout = normalize_structure({'mode': 'custom', 'umis': [{'position': p, 'length': 4 + i * 2} for i, p in enumerate(POSITIONS)]})
        seq, values = molecule(layout)
        with tempfile.TemporaryDirectory() as tmp:
            summary, rows, inserts, _ = prepare(Path(tmp), [seq + qc.rc(seq)], layout, 'microc')
            self.assertEqual(summary['library_units'], 2)
            self.assertEqual(summary['concatemer_reads'], 1)
            self.assertEqual(summary['umi_summary']['complete_units'], 2)
            self.assertEqual([r['umi_sequence'] for r in rows], ['-'.join(values[p] for p in POSITIONS)] * 2)
            self.assertEqual([r[1] for r in inserts], [INSERT, qc.rc(INSERT)])
            self.assertEqual(len(list(qc.fastq_records([Path(tmp) / 'leg1.gz']))), 2)
            self.assertEqual(qc.full_library_profile(rows, summary)['lengths'], [len(seq)] * 2)

    def test_random_umis_do_not_create_false_index_pairs(self):
        layout = normalize_structure({'mode': 'custom', 'umis': [{'position': 'after_i5', 'length': 8}, {'position': 'before_i7', 'length': 12}]})
        sequences = [molecule(layout, umi_values={'after_i5': value[:8], 'before_i7': value[:12]})[0] for value in [UMI, qc.rc(UMI), 'AGTC' * 4, 'TCGA' * 4]]
        with tempfile.TemporaryDirectory() as tmp:
            summary, rows, _, _ = prepare(Path(tmp), sequences, layout)
            self.assertEqual(len(summary['index_pairing']['pairs']), 1)
            self.assertEqual(sum(row['assignment_status'] == 'assigned' for row in rows), 4)
            self.assertEqual(len({row['umi_sequence'] for row in rows}), 4)

    def test_missing_outer_anchors_keep_umi_unknown_and_still_trim_inner_umi(self):
        layout = normalize_structure({'mode': 'custom', 'umis': [{'position': 'before_i5', 'length': 8}, {'position': 'before_insert', 'length': 12}]})
        complete, _ = molecule(layout)
        partial = qc.TRU_R1 + UMI[:12] + INSERT + qc.rc(qc.TRU_R2)
        with tempfile.TemporaryDirectory() as tmp:
            summary, rows, inserts, _ = prepare(Path(tmp), [complete, partial, qc.rc(partial)], layout)
            self.assertEqual(summary['umi_summary']['complete_units'], 1)
            self.assertEqual(summary['umi_summary']['partial_units'], 2)
            self.assertEqual([r['umi_before_i5'] for r in rows[1:]], ['', ''])
            self.assertEqual([r['umi_sequence'] for r in rows[1:]], ['', ''])
            self.assertEqual([r[1] for r in inserts], [INSERT, INSERT, qc.rc(INSERT)])

    def test_layout_too_long_is_counted_and_not_mapped_as_a_dimer(self):
        layout = normalize_structure({'mode': 'custom', 'umis': [{'position': 'before_insert', 'length': 12}]})
        short = qc.P5 + I5 + qc.TRU_R1 + 'GCTA' + qc.rc(qc.TRU_R2) + I7 + qc.rc(qc.P7)
        with tempfile.TemporaryDirectory() as tmp:
            summary, rows, inserts, _ = prepare(Path(tmp), [short], layout)
            self.assertEqual(summary['umi_summary']['rejected_short_units'], 1)
            self.assertEqual(summary['library_units'], 0)
            self.assertEqual(summary['adapter_dimer_units'], 0)
            self.assertEqual(inserts, [])

    def test_zero_insert_umi_library_is_a_dimer_after_umi_removal(self):
        layout = normalize_structure({'mode': 'custom', 'umis': [{'position': 'before_insert', 'length': 12}, {'position': 'after_insert', 'length': 8}]})
        seq, _ = molecule(layout, insert='')
        with tempfile.TemporaryDirectory() as tmp:
            summary, rows, inserts, _ = prepare(Path(tmp), [seq, qc.rc(seq)], layout)
            self.assertEqual(summary['adapter_dimer_units'], 2)
            self.assertEqual(summary['unreconstructed_dimer_screen'], 'not_assessed_for_umi_layout')
            self.assertIsNone(summary['unreconstructed_probable_dimer_reads'])
            self.assertEqual(inserts, [])
            self.assertEqual(qc.full_library_profile(rows, summary)['lengths'], [len(seq)] * 2)

    def test_custom_handles_and_expected_index_lengths_reach_length_model(self):
        layout = normalize_structure({'mode': 'custom', 'adapter_family': 'nextera', 'p5_sequence': 'AGGTCCTGACGATCGGTAACTG', 'p7_sequence': 'TGGACTGAGCCTAGCATTCGGATC', 'i5_length': 8, 'i7_length': 8, 'umis': [{'position': 'after_i7', 'length': 8}]})
        seq, _ = molecule(layout, 'nextera')
        with tempfile.TemporaryDirectory() as tmp:
            summary, rows, inserts, _ = prepare(Path(tmp), [seq, qc.rc(seq)], layout)
            self.assertEqual(summary['outer_complete_units'], 2)
            self.assertEqual([r[1] for r in inserts], [INSERT, qc.rc(INSERT)])
            self.assertEqual(qc.full_library_profile(rows, summary)['lengths'], [len(seq)] * 2)

    def test_long_configured_indexes_and_single_index_libraries(self):
        for left, right in [('ACGTACGTACGTACGTAC', 'TGCATGCATGCATGCATG'), ('', I7)]:
            layout = normalize_structure({'mode': 'custom', 'i5_length': len(left), 'i7_length': len(right), 'umis': [{'position': 'after_i5', 'length': 8}]})
            seq = qc.P5 + left + UMI[:8] + qc.TRU_R1 + INSERT + qc.rc(qc.TRU_R2) + right + qc.rc(qc.P7)
            with self.subTest(i5_length=len(left)), tempfile.TemporaryDirectory() as tmp:
                summary, rows, _, _ = prepare(Path(tmp), [seq, qc.rc(seq)], layout)
                self.assertEqual(summary['i5']['sequence'], left)
                self.assertEqual(summary['i7']['sequence'], right)
                self.assertEqual(rows[0]['i5_observed'], left)
                self.assertEqual(rows[1]['i5_observed'], qc.rc(left))
                self.assertEqual(qc.full_library_profile(rows, summary)['lengths'], [len(seq)] * 2)

    def test_auto_and_empty_custom_do_not_change_sequences_or_assignments(self):
        seq, _ = molecule(normalize_structure())
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b, tempfile.TemporaryDirectory() as c:
            default = prepare(Path(a), [seq, qc.rc(seq)])
            auto = prepare(Path(b), [seq, qc.rc(seq)], {'mode': 'auto'})
            custom = prepare(Path(c), [seq, qc.rc(seq)], {'mode': 'custom'})
            self.assertEqual(default, auto)
            self.assertEqual(default[1:], custom[1:])
            self.assertEqual(default[0]['index_pairing'], custom[0]['index_pairing'])


class StructureValidationTests(unittest.TestCase):
    def test_invalid_settings_are_rejected_by_the_pipeline(self):
        for bad in [{'mode': 'auto', 'umis': [{'position': 'before_insert', 'length': 8}]}, {'mode': 'custom', 'p5_sequence': 'x;echo unsafe'}, {'mode': 'custom', 'i5_length': True}, {'mode': 'custom', 'i7_length': 3.5}, {'mode': 'custom', 'i5_length': 25}, {'mode': 'custom', 'umis': [{'position': 'before_insert', 'length': 0}]}, {'mode': 'custom', 'umis': [{'position': 'before_insert', 'length': 65}]}, {'mode': 'custom', 'umis': [{'position': 'anywhere', 'length': 8}]}, {'mode': 'custom', 'umis': [{'position': 'before_insert', 'length': 8}] * 2}]:
            with self.subTest(bad=bad), self.assertRaises(ValueError): normalize_structure(bad)

if __name__ == '__main__': unittest.main()
