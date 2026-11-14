import math
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from report_profiles import read_gene_profile, render_profiles_card, gene_comparison, size_comparison, virtual_lanes


class ProfileComparisonTests(unittest.TestCase):
    def rows(self):
        return [dict(id='lib01',label='First <library>',color='#2a6cbd',sizes=([0,100,2000],[0,100,0]),gene=(['TSS','','TES'],[1,2,3],[1,2,3])),
                dict(id='lib02',label='Second',color='#23877f',sizes=([0,100,2000],[0,50,100]),gene=(['TSS','','TES'],[1,2,3],[10,20,30]))]

    def test_shared_cpm_scale_keeps_tenfold_signal_difference(self):
        root=ET.fromstring(gene_comparison(self.rows()))
        curves=root.findall('{*}polyline')
        self.assertEqual([c.attrib['stroke'] for c in curves],['#2a6cbd','#23877f'])
        low=float(curves[0].attrib['points'].split()[-1].split(',')[1])
        high=float(curves[1].attrib['points'].split()[-1].split(',')[1])
        self.assertAlmostEqual((350-high)/(350-low),10)
        incompatible=self.rows();incompatible[1]['gene']=(['TSS','','TES'],[1,2,4],[10,20,30])
        with self.assertRaises(ValueError):gene_comparison(incompatible)

    def test_lanes_share_size_coordinates_and_match_trace_peaks(self):
        root=ET.fromstring(virtual_lanes(self.rows()))
        lanes=[node for node in root.findall('{*}g') if node.get('data-scope')]
        self.assertEqual(len(lanes),2)
        peak_y=[]
        for lane in lanes:
            dark=[r for r in lane.findall('{*}rect') if r.get('fill')=='rgb(43,43,43)']
            self.assertEqual(len(dark),1);peak_y.append(float(dark[0].get('y')))
        self.assertLess(peak_y[1],peak_y[0])  # Larger fragments run closer to the top.
        traces=ET.fromstring(size_comparison(self.rows())).findall('{*}polyline')
        self.assertEqual(len(traces),2)

    def test_profile_data_missing_bins_and_empty_libraries_are_not_fabricated(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'profile.tsv'
            p.write_text('No MAPQ-pass primary mapped reads\n')
            self.assertIsNone(read_gene_profile(p))
            p.write_text('bin labels\t\tTSS\t\tTES\nbins\t\t1\t2\t3\nlib\tgenes\t1\tnan\t3\n')
            data=read_gene_profile(p);self.assertTrue(math.isnan(data[2][1]))
            rows=self.rows();rows[0]['gene']=data;rows[1]['gene']=None
            curves=ET.fromstring(gene_comparison(rows)).findall('{*}polyline')
            self.assertEqual(len(curves),2)  # No connecting line across the missing bin.

    def test_comparison_replaces_pooled_overview_and_details_start_collapsed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            for sid in ('lib01','lib02'):
                (p/f's.{sid}.full_library_profiles.tsv').write_text('size_bp\tlength_weighted_relative_signal\n0\t0\n100\t100\n2000\t0\n')
            scopes=[{'scope_id':'total'},{'scope_id':'lib01','scope_name':'A & B','library_units':12},{'scope_id':'lib02','library_units':2}]
            sections=["<section class='profile-scope'><h3>"+s['scope_id']+'</h3><p>Detail</p></section>' for s in scopes]
            doc=render_profiles_card('s',scopes,sections,p,p)
            self.assertEqual(doc.count("class='profile-detail'"),2)
            self.assertIn('A &amp; B',doc)
            self.assertNotIn(' open',doc)
            self.assertIn("<section class='report-section' id='profiles'>",doc)
            self.assertLess(doc.index('profile-comparison-grid\'>'),doc.index('<details'))
            self.assertIn('pooled-profile',doc)
            self.assertTrue((p/'s.sublibraries.virtual_gel.svg').is_file())
            single=render_profiles_card('s',scopes[:1],sections[:1],p,p)
            self.assertNotIn('<details',single)
            self.assertNotIn('sublibraries.size_profiles.svg',single)
            self.assertIn(sections[0],single)

if __name__=='__main__':unittest.main()
