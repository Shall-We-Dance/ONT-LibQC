import gzip
import random
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from library_sequences import rc
from scrna_chemistry import (C4, C4_COMPACT, MGI_R1_END, MGI_R2_START, PARSE, SMART_HANDLE, TENX,
                             TENX_TSO, Unit, detect, extract, native_mgi_units, round_match, trim_smart_handle)
from scrna_qc import load_units, make_queries, molecule_summary, read_well_map
from scrna_report import render

rng = random.Random(17)
def dna(n):
    return ''.join(rng.choice('ACGT') for _ in range(n))
P = [[dna(8) for _ in range(5)] for _ in range(3)]
MG = [[dna(10) for _ in range(5)] for _ in range(2)]
X = {book: [dna(16) for _ in range(5)] for book, _, _ in TENX.values()}
REF = {'parse': [set(x) for x in P], 'mgi': [set(x) for x in MG], 'tenx': {k: set(v) for k, v in X.items()}}
LINKS = ('CGACCGTATCAG', 'TGGACCACTTCA')
UMI = 'ACGTATCCAG'
CDNA = 'CGCACCTCAAGGACTCCGATGCGACTACTAGCACGAGGACTGC'


def fixture(kit, i=0, parent=None):
    j = i % 5
    if kit == PARSE:
        seq = SMART_HANDLE + CDNA + rc(UMI + P[0][j] + LINKS[0] + P[1][j] + LINKS[1] + P[2][j] + 'T' * 8)
    elif kit in (C4, C4_COMPACT):
        if kit == C4:
            prefix = 'TCTGCG' + MG[0][j] + 'CCTTCC' + MG[1][j] + 'CGATG'
        else:
            prefix = MG[0][j] + MG[1][j]
        seq = prefix + UMI + 'T' * 9 + CDNA
    else:
        book, n, anchor = TENX[kit]
        seq = X[book][j] + (UMI + 'GA')[:n] + ('T' * 9 if anchor == 'polyT' else TENX_TSO) + CDNA
    return Unit(f'u{i}', parent or f'parent{i}', 'lib01', seq, 'I' * len(seq))


class TestChemistry(unittest.TestCase):
    def test_all_supported_layouts_and_exact_cdna_boundary(self):
        for kit in [PARSE, C4, C4_COMPACT, *TENX]:
            with self.subTest(kit=kit):
                group = [fixture(kit, i) for i in range(15)]
                call, linkers = detect(group, REF)
                self.assertEqual(call['kit'], kit)
                self.assertEqual(call['evidence'][kit]['supported_units'], len(group))
                tag, cdna, quality = extract(group[0], kit, REF, linkers)
                self.assertEqual(cdna, CDNA)
                self.assertEqual(len(quality), len(cdna))
                self.assertTrue(tag.umi.startswith(UMI))

    def test_reject_spurious_tags_and_mixed_libraries(self):
        group = [fixture(PARSE, i, 'same-read') for i in range(25)]
        self.assertIsNone(detect(group, REF)[0]['kit'])
        group = [fixture('10x-3p-v2', i) for i in range(15)]
        for unit in group:
            unit.seq = unit.seq.replace('TTTTTTTTT', 'ACTGCACAC')
        call, _ = detect(group, REF)
        self.assertIsNone(call['kit'])
        self.assertEqual(call['evidence']['10x-3p-v2']['exact_barcode_units'], 15)
        group = [fixture('10x-3p-v3', i) for i in range(15)] + [fixture(C4, i + 15) for i in range(15)]
        self.assertEqual(detect(group, REF)[0]['status'], 'mixed_or_ambiguous')
        self.assertEqual(detect([], REF)[0]['status'], 'insufficient_evidence')

    def test_one_low_quality_substitution_and_linker_indel(self):
        u = fixture(PARSE)
        seq = rc(u.seq)
        pos = 11
        other = next(b for b in 'ACGT' if b != seq[pos])
        changed = seq[:pos] + other + seq[pos + 1:25] + 'A' + seq[25:]
        u.seq = rc(changed)
        u.qual = ('I' * pos + '!' + 'I' * (len(changed) - pos - 1))[::-1]
        tag, cdna, _ = extract(u, PARSE, REF, LINKS)
        self.assertEqual((tag.barcode, tag.edits, cdna), ('-'.join(p[0] for p in P), 1, CDNA))
        u.qual = 'I' * len(changed)
        self.assertIsNone(extract(u, PARSE, REF, LINKS))
        self.assertIsNone(round_match('CCAA', '!!!!', {'CAAA', 'ACAA'}, True))

    def test_parse_random_priming_and_anchored_handle(self):
        seq = UMI + P[0][0] + LINKS[0] + P[1][0] + LINKS[1] + P[2][0]
        prefix = (SMART_HANDLE + 'GAATGGG')[8:]
        u = Unit('u', 'r', 'lib01', prefix + CDNA + rc(seq), 'I' * (len(prefix) + len(CDNA) + len(seq)))
        self.assertEqual(extract(u, PARSE, REF, LINKS)[1], CDNA)
        self.assertEqual(trim_smart_handle(CDNA + SMART_HANDLE), 0)
        self.assertEqual(trim_smart_handle(SMART_HANDLE[-6:] + CDNA), 0)

    def test_native_mgi_orientation_boundaries(self):
        inner = fixture(C4).seq
        seq = MGI_R1_END + inner + MGI_R2_START
        for read in (seq, rc(seq)):
            unit, = native_mgi_units('raw', read, 'I' * len(read))
            self.assertEqual(unit.seq, inner)
            self.assertEqual(unit.scope, 'native-unindexed')
        self.assertEqual(len(native_mgi_units('raw', seq + seq, 'I' * (len(seq) * 2))), 2)
        self.assertFalse(native_mgi_units('raw', MGI_R1_END + inner, 'I' * (len(MGI_R1_END) + len(inner))))
        both = seq + MGI_R2_START
        self.assertFalse(native_mgi_units('raw', both, 'I' * len(both)))

    def test_restore_reverse_fastq_qualities_and_keep_native_units_separate(self):
        seq = fixture('10x-3p-v3').seq
        qual = '!' * 3 + 'I' * (len(seq) - 3)
        native = MGI_R1_END + fixture(C4).seq + MGI_R2_START
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            with gzip.open(d / 'reads.tsv.gz', 'wt') as f:
                f.write('unit_id\tread_id\torientation\tis_adapter_dimer\nu\tr\tp7_to_p5\t0\n')
            with gzip.open(d / 'assign.tsv.gz', 'wt') as f:
                f.write('unit_id\tassignment_status\tscope_id\nu\tassigned\tlib01\n')
            with gzip.open(d / 'inserts.fq.gz', 'wt') as f:
                f.write(f'@u\n{rc(seq)}\n+\n{qual[::-1]}\n')
            (d / 'raw.fq').write_text(f'@r\n{native}\n+\n{"I" * len(native)}\n@n\n{rc(native)}\n+\n{"I" * len(native)}\n')
            standard, mgi = load_units(d / 'raw.fq', d)
            self.assertEqual((standard.seq, standard.qual, standard.scope), (seq, qual, 'lib01'))
            self.assertEqual((mgi.seq, mgi.parent), (fixture(C4).seq, 'n'))

    def test_per_scope_calls_and_parse_well_map(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            path = d / 'wells.tsv'
            path.write_text('round1_barcode\twell\n' + ''.join(f'{b}\tW{i}\n' for i, b in enumerate(P[2])))
            wells = read_well_map(path, set(P[2]))
            group = [fixture(PARSE, i) for i in range(15)] + [fixture(C4, i) for i in range(15, 30)]
            for u in group[15:]:
                u.scope = 'lib02'
            scopes, units = make_queries(group, REF, wells, d / 'q.fq.gz')
            self.assertEqual([v['kit'] for v in scopes.values()], [PARSE, C4])
            self.assertTrue(units['u0']['barcode'].endswith('-W0'))
            with self.assertRaises(ValueError):
                make_queries(group, REF, wells, d / 'q.fq.gz', single_capture=True)
            with self.assertRaises(ValueError):
                read_well_map(path, {*P[2], 'TTTTTTTT'})

    def test_dedup_strand_and_matrix_denominators(self):
        row = {'mapping_status': 'unique', 'eligible_gene': True, 'gene_id': 'G',
               'strand_evidence': 'sense', 'region_bases': {'intron': 35},
               'family': 'unknown', 'junctions': [], 'known_junctions': 0}
        scopes = {name: {'funnel': Counter()} for name in ('lib01', 'unassigned')}
        tagged = {str(i): {'scope_id': 'lib01', 'molecule_scope': 'lib01', 'kit': PARSE, 'countable': True,
                           'umi': UMI, 'barcode': 'BC'} for i in range(6)}
        tagged['4'].update(scope_id='unassigned', countable=False)
        tagged['5'].update(umi='G' + UMI[1:])
        rows = [dict(row, unit_id=str(i)) for i in range(6)]
        rows[1]['region_bases'] = {'CDS': 24, 'intron': 0}
        rows[2]['strand_evidence'] = 'antisense'
        rows[3]['mapping_status'] = 'multimapped'
        m, entries = molecule_summary(rows, tagged, scopes)
        self.assertEqual((m['assigned_units'], m['molecules'], m['genes'], m['barcodes']), (3, 2, 1, 1))
        self.assertEqual((m['exonic_molecules'], m['intronic_only_molecules']), (1, 1))
        self.assertEqual(m['rarefaction'][-1]['molecules'], 2)
        self.assertEqual(sum(v['umis'] for v in entries), 2)

    def test_confirmed_single_capture_and_unassessed_native_counts(self):
        with tempfile.TemporaryDirectory() as t:
            q = Path(t) / 'q.fq.gz'
            group = [fixture(C4, i) for i in range(30)]
            for u in group[15:]:
                u.scope = 'unassigned'
            _, tagged = make_queries(group, REF, {}, q, single_capture=True)
            self.assertEqual({u['molecule_scope'] for u in tagged.values()}, {'lib01'})
            self.assertTrue(all(u['countable'] for u in tagged.values()))
            for u in group:
                u.scope = 'native-unindexed'
            scopes, tagged = make_queries(group, REF, {}, q)
            metrics, _ = molecule_summary([], tagged, scopes)
            self.assertEqual(metrics['molecule_assessment'], 'no_assigned_capture')
            data = {'sample': 'RNA <A&B>', 'assembly': 'fish', 'raw_reads': 30,
                    'reconstructed_parents': 30, 'non_dimer_units': 30, 'scopes': scopes,
                    'analysis_provenance': {'pipeline_version': 'dev', 'pipeline_build': 'test'},
                    'metrics': metrics, 'mapq': 20, 'annotation': None, 'single_capture': False,
                    'parse_well_map_sha256': None}
            page = render(data)
            self.assertIn('RNA &lt;A&amp;B&gt;', page)
            self.assertIn('<span class="rna-beta">Beta</span>', page)
            self.assertIn('Molecule counting requires a confirmed sample or capture', page)
            self.assertNotIn('Barcode combinations with a UMI', page)


if __name__ == '__main__':
    unittest.main()
