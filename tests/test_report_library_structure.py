"""Check that report diagrams reflect saved evidence, including missing data."""
import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
import ont_qc as qc

class LibraryDiagramTests(unittest.TestCase):
    def setUp(self):
        self.metric={'library_units':3,'adapter_family':'truseq','i5':'ACGTACGT','i7':'TGCAAAAA'}
        self.rows=[dict(family='truseq',p5_detected='1',p7_detected='0',scope_id='a',assignment_status='assigned'),
                   dict(family='truseq',p5_detected='1',p7_detected='1',scope_id='b',assignment_status='assigned'),
                   dict(family='nextera',p5_detected='0',p7_detected='1',scope_id='b',assignment_status='ambiguous')]
    def render(self,metric=None,rows=None,summary=None):
        return qc.library_structure_html(metric or self.metric, summary, self.rows if rows is None else rows)
    def test_auto_includes_sequences_and_evidence(self):
        text=self.render()
        for value in ('Auto detection','ACGTACGT','TGCAAAAA','Matched · 2 units',qc.P5,qc.rc(qc.P7),qc.TRU_R1,qc.rc(qc.TRU_R2),'matched reference motifs'):
            self.assertIn(value,text)
        self.assertIn('TruSeq',text);self.assertNotIn('illumina',text.lower())
    def test_overview_only_shows_index_sequences(self):
        text=self.render()
        diagram=text.split("<div class='library-diagram'>",1)[1].split("<p class='library-diagram-note'>",1)[0]
        self.assertIn('ACGTACGT',diagram);self.assertIn('TGCAAAAA',diagram)
        for motif in (qc.P5,qc.rc(qc.P7),qc.TRU_R1,qc.rc(qc.TRU_R2),qc.NEX_R1):self.assertNotIn(motif,diagram)
        self.assertIn('Standard',diagram);self.assertIn('TruSeq / Nextera',diagram)
        self.assertIn("<details class='library-sequence-details'>",text)
        self.assertIn(qc.TRU_R1,text.split("class='library-sequence-details'",1)[1])

    def test_single_custom_side_does_not_relabel_its_preset_partner(self):
        m=dict(self.metric,library_structure={'mode':'custom','adapter_i5_sequence':'ACGT'*6})
        rows=[dict(self.rows[0],family='custom_truseq')]
        text=self.render(m,rows)
        diagram=text.split("<div class='library-diagram'>",1)[1].split("<p class='library-diagram-note'>",1)[0]
        self.assertIn("class='library-type'>Custom<",diagram)
        self.assertIn("class='library-type'>TruSeq<",diagram)
        self.assertNotIn('ACGT'*6,diagram)
        self.assertIn('ACGT'*6,text)

    def test_undetected_motifs_are_not_claimed(self):
        text=self.render(rows=[])
        self.assertIn('Not detected',text)
        for value in (qc.P5,qc.TRU_R1,qc.rc(qc.TRU_R2)):
            self.assertNotIn(value,text)
    def test_custom_not_detected_is_explicit(self):
        m=copy.deepcopy(self.metric);m['library_structure']={'mode':'custom','p5_sequence':'ACGT'*6}
        text=self.render(m,[])
        self.assertIn('Configured · not detected',text);self.assertIn('ACGT'*6,text)
    def test_umi_rows_and_counts(self):
        m=copy.deepcopy(self.metric);m['library_structure']={'mode':'custom','umis':[{'position':'before_insert','length':12},{'position':'after_insert','length':8}]};m['umi_summary']={'complete_units':2,'partial_units':1,'rejected_short_units':0}
        text=self.render(m)
        for value in ('UMI 1','UMI 2','12 bp','8 bp','grid-column:7/10','<b>2</b>Complete units','excluded from mapped inserts'):
            self.assertIn(value,text)
    def test_subdiagrams_use_their_own_evidence(self):
        m=copy.deepcopy(self.metric);m['index_pairing']={'mode':'multiplexed','pairs':[{'name':'A','scope_id':'a','i5':'AAAA','i7':'CCCC','support_units':1},{'name':'B','scope_id':'b','i5':'TTTT','i7':'GGGG','support_units':1}]}
        text=self.render(m);subs=[s.split('</section>',1)[0] for s in text.split("<section class='library-sublibrary'>")[1:]]
        self.assertEqual(len(subs),2)
        for actual,expected,absent in ((subs[0],'AAAA','TTTT'),(subs[1],'TTTT','AAAA')):
            self.assertIn(expected,actual);self.assertNotIn('>'+absent+'<',actual)
            self.assertIn('33.3% of all units',actual)
        detail_a=text.split("data-scope='a'",1)[1].split('</section>',1)[0]
        detail_b=text.split("data-scope='b'",1)[1].split('</section>',1)[0]
        self.assertNotIn(qc.rc(qc.P7),detail_a);self.assertIn(qc.rc(qc.P7),detail_b)
        self.assertNotIn(qc.NEX_R1,detail_b) # ambiguous units are not sublibrary evidence
    def test_single_library_shows_dominant_index_pair(self):
        m=copy.deepcopy(self.metric);m['index_pairing']={'mode':'single_library','pairs':[{'i5':'AAAA','i7':'CCCC','informative_support':20},{'i5':'TTTT','i7':'GGGG','informative_support':1}]}
        text=self.render(m);self.assertIn('>AAAA<',text);self.assertNotIn('>TTTT<',text);self.assertNotIn('Resolved sublibraries',text)
    def test_legacy_sublibrary_metrics_render_without_pairs(self):
        m=copy.deepcopy(self.metric);m['sub_libraries']=[{'scope_id':'a','scope_name':'A','i5':'AAAA','i7':'CCCC','library_units':1}]
        self.assertEqual(self.render(m).count("class='library-sublibrary'"),1)
    def test_zero_index_differs_from_missing_index(self):
        m=dict(self.metric,i5='',i7=None);text=self.render(m)
        self.assertIn('No index · 0 bp',text);self.assertIn('Not detected',text)
    def test_names_and_sequences_are_escaped(self):
        m=dict(self.metric,i5='<script>bad()</script>');text=self.render(m)
        self.assertNotIn('<script>',text);self.assertIn('&lt;script&gt;',text)
    def test_summary_fallback_and_generic_family_label(self):
        text=qc.library_structure_html(self.metric,{'i5':{'sequence':'CCCC','all_flanked_observations':2},'i7':{'sequence':'AAAA','all_flanked_observations':0}})
        self.assertIn(qc.P5,text);self.assertNotIn(qc.rc(qc.P7),text);self.assertIn('>CCCC<',text)
        self.assertEqual(qc.display_value('adapter_family','truseq',{}),'TruSeq')

if __name__=='__main__':unittest.main()
