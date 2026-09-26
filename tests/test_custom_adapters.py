"""Custom adapter sequences must drive reconstruction and every derived display."""
import tempfile
import unittest
from pathlib import Path
from test_library_structure import prepare, INSERT, I5, I7, UMI, qc, normalize_structure
from library_structure import adapter_families, handles
LEFT='GCTACAGTCGATGACCTGATCGTAGCGA'
RIGHT='TCAGGATCGACTAGGCTACGTAGACTGACG'

class CustomAdapterTests(unittest.TestCase):
    def test_one_or_both_sides_both_strands_and_presets(self):
        for preset in ('auto','truseq','nextera'):
            for side in ('left','right','both'):
                with self.subTest(preset=preset,side=side),tempfile.TemporaryDirectory() as tmp:
                    layout=normalize_structure(dict(mode='custom',adapter_family=preset,
                        adapter_i5_sequence=LEFT if side!='right' else None,
                        adapter_i7_sequence=RIGHT if side!='left' else None,
                        umis=[dict(position='before_insert',length=8),dict(position='after_insert',length=12)]))
                    families=adapter_families(layout,qc.FAMILIES)
                    family='custom' if side=='both' else 'custom_'+('nextera' if preset=='nextera' else 'truseq')
                    r1,r2=families[family];p5,p7=handles(layout)
                    seq=p5+I5+r1+UMI[:8]+INSERT+UMI[:12]+qc.rc(r2)+I7+p7
                    summary,rows,inserts,_=prepare(Path(tmp),[seq,qc.rc(seq)],layout)
                    self.assertEqual(summary['library_units'],2)
                    self.assertEqual([r['family'] for r in rows],[family]*2)
                    self.assertEqual([x[1] for x in inserts],[INSERT,qc.rc(INSERT)])
                    self.assertEqual(summary['umi_summary']['complete_units'],2)
                    self.assertEqual(summary['i5']['sequence'],I5);self.assertEqual(summary['i7']['sequence'],I7)
                    self.assertEqual(qc.full_library_profile(rows,summary)['lengths'],[len(seq)]*2)
                    metric={'library_units':2,'library_structure':layout,'index_pairing':summary['index_pairing']}
                    page=qc.library_structure_html(metric,summary,rows)
                    self.assertIn(r1,page);self.assertIn(qc.rc(r2),page)
                    self.assertIn('Custom',page)
    def test_custom_zero_insert_templates_and_dimers(self):
        layout=normalize_structure(dict(mode='custom',adapter_i5_sequence=LEFT,adapter_i7_sequence=RIGHT))
        seq=qc.P5+I5+LEFT+RIGHT+I7+qc.rc(qc.P7)
        self.assertEqual(qc.zero_insert_template('custom',I5,I7,layout),seq)
        for orientation in (seq,qc.rc(seq)):
            self.assertEqual(qc.full_length_template_identity(orientation,seq),1)
        with tempfile.TemporaryDirectory() as tmp:
            summary,rows,inserts,_=prepare(Path(tmp),[seq,qc.rc(seq)],layout)
            self.assertEqual(summary['adapter_dimer_units'],2)
            self.assertEqual(inserts,[])
            self.assertEqual(qc.full_library_profile(rows,summary)['lengths'],[len(seq)]*2)
    def test_custom_sequences_without_evidence_are_labelled_configured(self):
        layout=normalize_structure(dict(mode='custom',adapter_i5_sequence=LEFT,adapter_i7_sequence=RIGHT))
        text=qc.library_structure_html({'library_structure':layout},unit_rows=[])
        self.assertEqual(text.count('Configured · not detected'),2)
        self.assertIn(LEFT,text);self.assertIn(RIGHT,text)
    def test_preset_not_mutated_and_empty_custom_keeps_auto_results(self):
        before=dict(qc.FAMILIES)
        self.assertEqual(adapter_families(normalize_structure(),qc.FAMILIES),before)
        adapter_families(normalize_structure(dict(mode='custom',adapter_i5_sequence=LEFT)),qc.FAMILIES)
        self.assertEqual(qc.FAMILIES,before)
    def test_sequences_are_bounded_and_auto_cannot_override(self):
        for key in ('adapter_i5_sequence','adapter_i7_sequence'):
            for value in ('','A'*15,'A'*81,'ACGT'*4+'N',123,True,'ACGT;unsafeACGTACGT'):
                with self.subTest(key=key,value=value),self.assertRaises(ValueError):normalize_structure({'mode':'custom',key:value})
            with self.assertRaises(ValueError):normalize_structure({'mode':'auto',key:LEFT})
            for value in ('ACGT'*4,'ACGT'*20):self.assertEqual(normalize_structure({'mode':'custom',key:value})[key],value)

if __name__=='__main__':unittest.main()
