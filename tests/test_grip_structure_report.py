import copy
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from report_library_structure import library_structure_html

class GripStructureTests(unittest.TestCase):
    def setUp(self):
        self.metric = {'assay':'clip', 'library_units':200, 'sub_libraries':[
            {'scope_id':'lib01', 'scope_name':'A', 'i5':'AAAA', 'i7':'CCCC', 'library_units':100},
            {'scope_id':'lib02', 'scope_name':'B', 'i5':'AAAA', 'i7':'GGGG', 'library_units':100}]}
        self.rna = {'preparation':{'mode':'auto_grip', 'auto_grip':{'scope_evidence':{
            'lib01':{'supported':True}, 'lib02':{'supported':False}}}},
            'scopes':{s:{'units':n, 'umi':{'auto_grip_p7_complete_units':count}} for s,n,count in
                [('total',200,80), ('lib01',100,80), ('lib02',100,0)]}}
    def render(self, metric=None, rna=None):
        return library_structure_html(metric or self.metric, adapter_templates={'truseq':('ACGT','TGCA')}, rna_summary=rna if rna is not None else self.rna)
    def test_observed_and_expected_are_scope_specific(self):
        page=self.render()
        self.assertEqual(page.count("grip-umi detected"),2)
        self.assertEqual(page.count("grip-umi expected"),1)
        self.assertIn('80 / 100 recovered units (80.0%)',page)
        self.assertIn('expected position; not confirmed in this sublibrary',page)
        self.assertEqual(page.count("data-position='after_insert'"),3)
        self.assertIn('RNA-derived DNA',page)
        self.assertNotIn('Sample DNA',page)
        self.assertNotIn('UMI 2',page)
    def test_protocol_only_has_no_detected_claim(self):
        rna=copy.deepcopy(self.rna)
        for scope in rna['preparation']['auto_grip']['scope_evidence'].values():scope['supported']=False
        for scope in rna['scopes'].values():scope['umi']['auto_grip_p7_complete_units']=0
        self.assertEqual(self.render(rna=rna).count('grip-umi expected'),3)
        self.assertNotIn('grip-umi detected',self.render(rna=rna))
        self.assertNotIn('grip-umi',self.render(metric=dict(self.metric,assay='rnaseq'),rna=rna))
    def test_configured_umi_takes_precedence(self):
        metric=dict(self.metric,library_structure={'mode':'custom','umis':[{'position':'after_insert','length':12}]})
        page=self.render(metric=metric)
        self.assertNotIn('grip-umi',page)
        self.assertIn('UMI 1',page)
        self.assertIn('12 bp',page)
    def test_disabled_and_missing_inference_preserve_generic_diagram(self):
        for bypass in ['configured_umi_extraction','manual_end_trimming','auto_disabled']:
            rna=copy.deepcopy(self.rna);rna['preparation']['auto_grip']['bypass_reason']=bypass
            self.assertNotIn('grip-umi',self.render(rna=rna))
        self.assertNotIn('grip-umi',self.render(rna={}))

if __name__=='__main__':unittest.main()
