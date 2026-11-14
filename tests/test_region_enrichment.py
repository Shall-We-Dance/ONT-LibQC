import itertools
from pathlib import Path
import sys,tempfile,unittest
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from report_regions import OverlapBackground,region_counts,svg_region_bars,load_te_bed

class RegionEnrichmentTests(unittest.TestCase):
    def test_exact_background_matches_exhaustive_valid_placements(self):
        for intervals in [[],[(0,2)],[(8,10)],[(2,4),(3,6),(8,9)],[(0,10)]]:
            bg=OverlapBackground(intervals,10)
            for length in range(1,11):
                hits=sum(any(s<b and s+length>a for a,b in intervals) for s in range(11-length))
                self.assertAlmostEqual(bg.probability(length),hits/(11-length))

    def test_exclusive_gene_regions_and_independent_te(self):
        records={str(i):dict(chrom='chr1',start=s,end=e,mapq=q) for i,(s,e,q) in enumerate([(1,4,60),(5,6,60),(8,9,60),(1,4,0)])}
        m=region_counts(records,{'chr1':[(2,4)]},{'chr1':[(3,7)]},['chr1'],20,2000,{'chr1':10},{'chr1':[(3,6)]})
        self.assertEqual([m[k+'_reads'] for k in ['promoter_proximal','gene_body','intergenic','te']],[1,1,1,2])
        self.assertAlmostEqual(sum(m[k+'_background_rate'] for k in ['promoter_proximal','gene_body','intergenic']),1)
        self.assertAlmostEqual(m['promoter_proximal_background_rate'],(4/8+2/10+2/10)/3)
        self.assertAlmostEqual(m['te_background_rate'],(5/8+3/10+3/10)/3)
        self.assertAlmostEqual(m['te_enrichment'],(2/3)/m['te_background_rate'])
        self.assertEqual(m['gene_body_genome_rate'],.3)

    def test_missing_annotation_zero_background_and_empty_library(self):
        r={'r':dict(chrom='c',start=0,end=1,mapq=60)}
        missing=region_counts(r,{}, {},['c'],20,2000,{'c':10})
        self.assertIsNone(missing['te_rate']);self.assertIsNone(missing['gene_body_background_rate'])
        empty=region_counts(r,{'c':[(1,3)]},{},['c'],20,2000,{'c':10},{})
        self.assertEqual(empty['te_rate'],0);self.assertEqual(empty['te_background_rate'],0);self.assertIsNone(empty['te_enrichment'])
        no_reads=region_counts({}, {'c':[(1,3)]},{},['c'],20,2000,{'c':10},{'c':[(1,3)]})
        self.assertIsNone(no_reads['te_rate']);self.assertIsNone(no_reads['te_enrichment']);self.assertEqual(no_reads['te_genome_rate'],.2)
        svg=svg_region_bars(missing,'Test <sample>');ET.fromstring(svg)
        self.assertNotIn('TE overlap', [node.text for node in ET.fromstring(svg).iter('{http://www.w3.org/2000/svg}text')]);self.assertNotIn('nan',svg)

    def test_background_matches_chromosome_and_length_distribution(self):
        r={'r':dict(chrom='a',start=0,end=5,mapq=60)}
        m=region_counts(r,{'a':[(5,6)]},{},['a','b'],20,2000,{'a':10,'b':1000},{})
        self.assertAlmostEqual(m['promoter_proximal_background_rate'],5/6)
        self.assertAlmostEqual(m['promoter_proximal_genome_rate'],1/1010)

    def test_te_bed_merges_overlap_and_validates_bounds(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'te.bed';p.write_text('1\t1\t4\nchr1\t3\t6\n')
            self.assertEqual(load_te_bed(p,{'chr1':10}),{'chr1':[(1,6)]})
            p.write_text('chr1\t1\t11\n')
            with self.assertRaises(ValueError):load_te_bed(p,{'chr1':10})

if __name__=='__main__':unittest.main()

class TEPreparationTests(unittest.TestCase):
    def test_repeatmasker_classes_and_coordinate_conversion(self):
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'portal/scripts'))
        from prepare_te_annotation import prepare
        import gzip
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);src=root/'repeats.out';sizes=root/'sizes';source_sizes=root/'source.sizes';aliases=root/'aliases';out=root/'te.bed.gz'
            sizes.write_text('chr1\t100\n');source_sizes.write_text('NC_1\t100\n');aliases.write_text('NC_1\t1\tchr1\n')
            src.write_text('  100 1.0 0.0 0.0 NC_1 3 8 (92) + L1 LINE/L1 1 6 (0) 1\n'
                           '  100 1.0 0.0 0.0 NC_1 7 12 (88) + L2 SINE/Alu 1 6 (0) 2\n'
                           '  100 1.0 0.0 0.0 NC_1 30 40 (60) + simple Simple_repeat 1 11 (0) 3\n')
            m=prepare(src,'out',sizes,source_sizes,out,'test','https://example.org/repeats',aliases,['chr1'])
            with gzip.open(out,'rt') as f:self.assertEqual(f.read(),'chr1\t2\t12\n')
            self.assertEqual(m['included_classes'],{'LINE':1,'SINE':1});self.assertEqual(m['excluded_classes'],{'Simple_repeat':1})
            source_sizes.write_text('NC_1\t101\n')
            with self.assertRaisesRegex(ValueError,'length mismatch'):prepare(src,'out',sizes,source_sizes,out,'test','https://example.org/repeats',aliases,['chr1'])
