import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from report_spike_matrix import spike_matrix_html
class MatrixTests(unittest.TestCase):
 def test_independent_denominators_and_zero(self):
  screen={'candidates':{'e':{'label':'E. coli'},'x':{'label':'Other'}}}
  libs=[{'scope_id':'lib01','spike_in_screened_reads':100,'spike_in_species_counts':{'e':10}}, {'scope_id':'lib02','spike_in_screened_reads':400,'spike_in_species_counts':{'e':10}}]
  s=spike_matrix_html(screen,libs,'sample');self.assertIn('10.00%',s);self.assertIn('2.50%',s);self.assertEqual(s.count('>0%</td>'),2);self.assertEqual(s.count('light amber'),2);self.assertNotIn('Total library',s)
 def test_missing_and_host_are_not_zero(self):
  s=spike_matrix_html({'candidates':{'e':{'status':'HOST_NOT_DISTINGUISHABLE'},'x':{}}},[{'scope_id':'<lib>','spike_in_screened_reads':0}], 'x')
  self.assertEqual(s.count('>N/A</td>'),2);self.assertIn('&lt;lib&gt;',s);self.assertNotIn('>0%</td>',s)
 def test_small_positive_not_rounded_to_zero(self):
  s=spike_matrix_html({'candidates':{'e':{}}},[{'scope_id':'lib01','spike_in_screened_reads':100000,'spike_in_species_counts':{'e':1}}],'s')
  self.assertIn('&lt;0.01%',s);self.assertIn('light amber',s)
