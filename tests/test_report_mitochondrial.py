from pathlib import Path
import sys, unittest
import xml.etree.ElementTree as ET
import json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from report_mitochondrial import coverage_html

class MitochondrialCoverageTests(unittest.TestCase):
    def test_mitochondrial_midpoints_and_denominator_exclude_low_quality_and_other_contigs(self):
        rows=[('chrM',0,4,20),('chrM',96,100,60),('chrM',2,8,19),('chrM',2,8,255),('chr1',10,20,60),('other',0,10,60)]
        records={str(i):dict(chrom=c,start=s,end=e,mapq=q) for i,(c,s,e,q) in enumerate(rows)}
        rendered=coverage_html(records,{'chr1':1000,'chrM':100,'other':500},{'chrM'},20,primary=['chr1'])
        root=ET.fromstring(rendered)
        m=json.loads(root.find('.//{*}metadata').text)
        self.assertEqual(m['mapped_alignments'],2);self.assertEqual(m['denominator'],3)
        self.assertEqual(m['midpoint_counts'][1],1);self.assertEqual(m['midpoint_counts'][49],1)
        self.assertIn('66.7%',rendered)

    def test_missing_positions_differ_from_zero_coverage_and_use_reference_name(self):
        saved=coverage_html(None,{'MT':100},{'MT'},20,{'mitochondrial_fraction':.1})
        self.assertIn('MT',saved);self.assertIn('10.0%',saved);self.assertNotIn('<svg',saved)
        zero=coverage_html({}, {'MT':100},{'MT'},20)
        self.assertIn('0 alignments',zero);self.assertIn('<svg',zero)
        self.assertEqual(coverage_html({}, {'chr1':100},{'chrM'},20),'')

class IntegratedCoverageTests(unittest.TestCase):
    def test_integrated_row_preserves_counts_nuclear_geometry_and_footer(self):
        from report_mitochondrial import integrate_karyotype
        svg = '<svg xmlns="http://www.w3.org/2000/svg" width="1008" height="720" viewBox="0 0 1008 720"><path id="nuclear" d="M10 10L20 30"/><g id="caption"><use y="711"/></g></svg>'
        coverage = coverage_html({'a':dict(chrom='chrM',start=2,end=8,mapq=60)}, {'chrM':100}, {'chrM'},20)
        rendered = integrate_karyotype(svg, coverage)
        root = ET.fromstring(rendered)
        self.assertEqual(root.get('viewBox'), '0 0 1008 774')
        self.assertEqual(root.find('.//{*}path[@id="nuclear"]').get('d'), 'M10 10L20 30')
        self.assertEqual(root.find('.//{*}g[@id="caption"]').get('transform'), 'translate(0 54)')
        self.assertEqual(json.loads(root.find('.//{*}metadata').text)['mapped_alignments'],1)
        self.assertEqual(integrate_karyotype(rendered, coverage), rendered)
        zero = coverage_html({}, {'chrM':100}, {'chrM'},20)
        updated = integrate_karyotype(rendered, zero)
        self.assertEqual(updated.count('id="mitochondrial-tracks"'),1)
        self.assertEqual(json.loads(ET.fromstring(updated).find('.//{*}metadata').text)['mapped_alignments'],0)

    def test_missing_positions_remain_unavailable_inside_figure(self):
        from report_mitochondrial import integrate_karyotype
        svg = '<svg xmlns="http://www.w3.org/2000/svg" width="1008" height="720" viewBox="0 0 1008 720"/>'
        rendered = integrate_karyotype(svg,coverage_html(None,{'MT':100},{'MT'},20,{'mitochondrial_fraction':.1}))
        self.assertIn('10.0%', rendered)
        self.assertIn('Position data unavailable',rendered)
        self.assertIsNone(ET.fromstring(rendered).find('.//{*}metadata'))
        self.assertEqual(integrate_karyotype(svg,''),svg)
