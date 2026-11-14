"""Physical single-index junctions, mixed pools and incomplete-read ambiguity."""
import tempfile
import unittest
from pathlib import Path
from test_library_structure import qc, prepare, INSERT, I5, I7, normalize_structure


def molecule(i7=I7, i5='', overlap=4, family='truseq'):
    r1, r2 = qc.FAMILIES[family]
    p5 = qc.P5[:-overlap] if overlap else qc.P5
    return p5 + i5 + r1 + INSERT + qc.rc(r2) + i7 + qc.rc(qc.P7)


class SingleIndexTests(unittest.TestCase):
    def test_truseq_single_junction_both_strands_and_full_size(self):
        seq = molecule()
        with tempfile.TemporaryDirectory() as tmp:
            summary, rows, inserts, _ = prepare(Path(tmp), [seq, qc.rc(seq)])
            self.assertEqual(summary['index_architecture'], 'classic_single_index')
            self.assertEqual(summary['index_pairing']['informative_units'], 2)
            self.assertEqual(summary['index_pairing']['pairs'][0]['i5'], '')
            self.assertEqual(summary['index_pairing']['pairs'][0]['i7'], I7)
            self.assertEqual([r['p5_adapter_overlap'] for r in rows], ['4', '4'])
            self.assertEqual([x[1] for x in inserts], [INSERT, qc.rc(INSERT)])
            self.assertEqual(qc.full_library_profile(rows, summary)['lengths'], [len(seq)]*2)
            self.assertIn('No index', qc.identity_panel_html({},summary['index_pairing'],summary['index_pairing']['pairs'], rows))

    def test_multiplexed_single_indexes_and_mixed_dual_pool(self):
        single = [molecule(i7=v) for v in ('GATCAG','TAGCTT','CTTGTA','GGCTAC')]
        dual = molecule(i7='GATCAG',i5=I5,overlap=0)
        with tempfile.TemporaryDirectory() as tmp:
            summary, rows, _, _ = prepare(Path(tmp), single + [dual])
            pairing = summary['index_pairing']
            self.assertEqual(pairing['mode'],'multiplexed_library')
            self.assertEqual(len(pairing['pairs']),5)
            self.assertEqual(summary['index_architecture'],'mixed_index_architecture')
            self.assertEqual(sum(p['support_units'] for p in pairing['pairs']),5)
            self.assertEqual(qc.full_library_profile(rows,summary)['lengths'],list(map(len,single+[dual])))

    def test_i5_only_and_dual_nextera_still_work(self):
        for family in ('truseq','nextera'):
            for i7 in ('',I7):
                with self.subTest(family=family,i7=i7),tempfile.TemporaryDirectory() as tmp:
                    seq=molecule(i7=i7,i5=I5,overlap=0,family=family)
                    summary, rows, _, _=prepare(Path(tmp),[seq,qc.rc(seq)])
                    self.assertEqual(summary['index_architecture'],'physical_dual_index' if i7 else 'i5_single_index')
                    self.assertEqual(qc.full_library_profile(rows,summary)['lengths'],[len(seq)]*2)

    def test_junction_errors_are_tolerated_without_changing_insert(self):
        seq=molecule()
        variants=[seq[:7]+'T'+seq[8:],seq[:7]+'G'+seq[7:],seq[:7]+seq[8:]]
        for variant in variants:
            for oriented in (variant,qc.rc(variant)):
                units,_=qc.find_units(oriented,5000);qc.add_indexes(oriented,units)
                self.assertEqual(len(units),1)
                self.assertEqual(units[0]['i5'],'')
                self.assertEqual(units[0]['p5_adapter_overlap'],4)
                self.assertEqual(units[0]['i7'],I7)

    def test_truncated_outer_end_is_unknown_not_absent(self):
        # No outer P5 bases: the inner adapter and i7 cannot prove a single-index architecture.
        seq=qc.TRU_R1+INSERT+qc.rc(qc.TRU_R2)+I7+qc.rc(qc.P7)
        with tempfile.TemporaryDirectory() as tmp:
            summary,rows,_,_=prepare(Path(tmp),[seq,qc.rc(seq)])
            self.assertEqual(summary['index_pairing']['informative_units'],0)
            self.assertEqual(summary['index_pairing']['pairs'],[])
            self.assertEqual(summary['index_architecture'],'partial_or_custom')
            self.assertTrue(all(r['p5_detected']=='0' for r in rows))

    def test_partial_index_reads_require_discriminating_evidence(self):
        def unit(name,i5,i7):
            return dict(unit_id=name,read_id=name,i5=i5,i7=i7,p5_detected=i5 is not None,
                        p7_detected=i7 is not None,is_adapter_dimer=False)
        # Two dual-index libraries share i7: a missing i5 must remain ambiguous.
        units=[unit('a',I5,I7),unit('b','GATCTAGC',I7),unit('partial',None,I7),unit('none',None,None)]
        result,assigned=qc.cluster_index_pairs(units,minimum_cluster_units=1)
        self.assertEqual(assigned[2]['assignment_status'],'ambiguous')
        self.assertEqual(assigned[3]['assignment_status'],'unassigned')
        units=[unit('single','',I7),unit('missing_index','',None),unit('partial',None,I7)]
        result,assigned=qc.cluster_index_pairs(units,minimum_cluster_units=1)
        self.assertEqual(assigned[1]['assignment_status'],'unassigned')
        self.assertEqual(assigned[2]['assignment_status'],'assigned')
        # Unknown after configured UMI extraction is not an empty index or a TypeError.
        units[2]['p5_detected']=True
        qc.cluster_index_pairs(units,minimum_cluster_units=1)

    def test_empty_and_short_indexes_do_not_collapse_as_indel_variants(self):
        self.assertFalse(qc.index_close('', 'T', 2))
        self.assertGreater(qc.index_pair_distance({'i5':'','i7':I7},{'i5':'T','i7':I7},2),2)

    def test_custom_outer_umis_do_not_trigger_overlapping_p5_rescue(self):
        layout=normalize_structure({'mode':'custom','umis':[{'position':'after_i5','length':8}]})
        seq=molecule();units,_=qc.find_units(seq,5000)
        self.assertIsNone(qc.single_index_p5_hit(seq,units[0],layout))

    def test_single_index_dimer_template_uses_real_junction(self):
        seq=molecule().replace(INSERT,'')
        self.assertEqual(qc.zero_insert_template('truseq','',I7,p5_adapter_overlap=4),seq)
        with tempfile.TemporaryDirectory() as tmp:
            summary,rows,inserts,_=prepare(Path(tmp),[seq,qc.rc(seq)])
            self.assertEqual(summary['adapter_dimer_units'],2)
            self.assertEqual(inserts,[])
            self.assertEqual(qc.full_library_profile(rows,summary)['lengths'],[len(seq)]*2)

if __name__=='__main__':unittest.main()
