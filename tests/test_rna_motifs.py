"""Validate real motif recovery, incomplete tools, and visible RNA screening."""
import os
from pathlib import Path
import random
import shutil
import sys
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from rna_qc import dinucleotide_shuffle, motif_aligned_sequence, reverse_complement
from rna_motifs import discover_homer, parse_motifs, logo_svg, adapter_match
from rna_report import render_motifs
from ont_qc import rna_screening_html

class MotifTests(unittest.TestCase):
    def sequences(self):
        rng=random.Random(42)
        return [(str(i), ''.join(rng.choices('ACGT', k=25))+'TGCATG'+''.join(rng.choices('ACGT', k=25))) for i in range(180)]

    @unittest.skipUnless(shutil.which(os.environ.get('HOMER2', 'homer2')), 'HOMER integration executable is not installed')
    def test_homer_recovers_planted_rna_motif(self):
        result=discover_homer(self.sequences(),dinucleotide_shuffle,shuffles=4)
        self.assertEqual(result['status'],'complete')
        motif=next(m for m in result['motifs'] if m['motif']=='UGCAUG')
        self.assertGreater(motif['enrichment'],10)
        self.assertEqual(motif['observed_units'],180)
        ET.fromstring(logo_svg(motif['pwm']))
        self.assertIn('UGCAUG',render_motifs(result))

    def test_missing_tool_is_not_a_negative_motif_result(self):
        result=discover_homer(self.sequences(),dinucleotide_shuffle,executable='/nonexistent/homer2')
        self.assertEqual(result['status'],'engine_unavailable')
        self.assertIn('unavailable',render_motifs(result))
        self.assertNotIn('No recurrent',render_motifs(result))

    def test_motif_input_omits_softclips_and_respects_rna_strand(self):
        aligned='GACTGGCA';query='TTCCGATCT'+aligned+'AAAAA'
        record={'cigar':'9S8M5S','reverse':False}
        self.assertEqual(motif_aligned_sequence(query,record,'+'),aligned)
        self.assertEqual(motif_aligned_sequence(query,record,'-'),reverse_complement(aligned))
        record['reverse']=True
        self.assertEqual(motif_aligned_sequence(reverse_complement(query),record,'+'),aligned)
        self.assertEqual(motif_aligned_sequence(reverse_complement(query),record,'-'),reverse_complement(aligned))
        self.assertEqual(motif_aligned_sequence(query,{'cigar':'1H22M','reverse':False},'+'),'')

    def test_adapter_matches_are_flagged_in_both_orientations(self):
        self.assertTrue(adapter_match('GAUCGGAA'))
        self.assertTrue(adapter_match('UUCCGAUC'))
        self.assertFalse(adapter_match('UGCAUG'))
        self.assertFalse(adapter_match('UUCCG'))

    def test_incomplete_pwm_is_rejected(self):
        with self.assertRaises(ValueError):
            parse_motifs('>AC\t1-AC\t2\t-3\t0\tT:2(20.0%),B:1(10.0%)\n.25 .25 .25 .25\n',10,10)

    def test_screen_preserves_not_run_and_zero_support(self):
        scope={'scope_id':'total','scope_name':'Total library','mycoplasma_status':'NOT_CONFIGURED','spike_in_status':'NOT_CONFIGURED'}
        rendered=rna_screening_html([scope])
        self.assertIn('id=\'contamination\'',rendered)
        self.assertIn('NOT_CONFIGURED',rendered)
        scope.update(mycoplasma_status='NOT_DETECTED_AT_THIS_DEPTH',mycoplasma_screened_reads=1000,mycoplasma_qualifying_reads=0,mycoplasma_qualifying_fraction=0)
        self.assertIn('NOT_DETECTED_AT_THIS_DEPTH',rna_screening_html([scope]))

if __name__=='__main__':unittest.main()
