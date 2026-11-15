from pathlib import Path
import sys,tempfile,unittest,gzip,json,xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from report_similarity import calculate,ranks,pearson,cluster_order,matrix_for,write_assets,report_html

class SimilarityTests(unittest.TestCase):
    def test_ties_and_correlations(self):
        self.assertEqual(ranks([3,0,0,2]),[4,1.5,1.5,3])
        self.assertAlmostEqual(pearson([0,1,2],[0,10,20]),1)
        self.assertAlmostEqual(pearson([0,1,2],[2,1,0]),-1)
        self.assertIsNone(pearson([1,1],[2,3]))
        self.assertIsNone(pearson([],[]))
        self.assertIsNone(pearson([1],[2]))

    def test_bins_filters_splices_and_one_count_per_alignment(self):
        scopes=[dict(scope_id='a'),dict(scope_id='b')]
        def record(start,cigar,**kw):return dict(chrom='chr1',start=start,end=start+20,mapq=20,cigar=cigar,**kw)
        records={'a':{'r1':record(5,'10M'),'r2':record(0,'2M2I2M20N2M'),
                      'low':dict(record(1,'1M'),mapq=19),'mito':dict(record(0,'10M'),chrom='chrM'),
                      'supp':record(1,'1M',flag=2048),'bad':record(0,'10S')},
                 'b':{'r':record(90,'30M')}}
        data,tables=calculate(scopes,records,{'chr1':100,'chrM':20},['chr1'],20,[10])
        self.assertEqual(data['eligible_alignments'],[2,1])
        self.assertEqual(data['scales']['10']['covered_bins'],[3,1])
        self.assertEqual(data['scales']['10']['bins_used'],10)
        self.assertEqual(data['scales']['10']['all_zero_bins'],6)
        self.assertEqual(tables[10].splitlines(),['chrom\tstart\tend\ta\tb','chr1\t0\t10\t2\t0','chr1\t10\t20\t1\t0','chr1\t20\t30\t1\t0','chr1\t90\t100\t0\t1'])

    def test_joint_zero_filter_keeps_disjoint_signal_and_empty_scope(self):
        scopes=[dict(scope_id=k) for k in 'abc']
        r=lambda x:dict(chrom='c',start=x,end=x+1,mapq=60)
        data,_=calculate(scopes,{'a':{'a':r(0)},'b':{'b':r(10)}},{'c':100},['c'],bin_sizes=[10])
        for method in ('spearman','pearson'):
            m=data['scales']['10']['nonzero_correlations'][method]['matrix']
            self.assertEqual(m,[[1,-1,None],[-1,1,None],[None,None,None]])

    def test_implicit_shared_zeros_equal_full_matrix(self):
        columns=[[0,2,1,0,4],[1,0,0,2,3],[0,0,0,0,0]]
        for method in ('spearman','pearson'):
            compact=matrix_for(columns,50,method)['matrix']
            expanded=matrix_for([c+[0]*50 for c in columns],0,method)['matrix']
            for a,b in zip(compact,expanded):
                for x,y in zip(a,b):
                    if x is None:self.assertIsNone(y)
                    else:self.assertAlmostEqual(x,y,places=12)
        # Disjoint sparse reads are nearly uncorrelated with all genomic bins.
        m=matrix_for([[1,0],[0,1]],98,'spearman')['matrix']
        self.assertAlmostEqual(m[0][1],-1/99)

    def test_average_linkage_and_deterministic_order(self):
        m=[[1,.1,.9,None],[.1,1,.2,None],[.9,.2,1,None],[None]*4]
        order,merges=cluster_order(m)
        self.assertEqual(merges[0]['left'],[0]);self.assertEqual(merges[0]['right'],[2])
        self.assertEqual(order,[0,2,1,3])
        self.assertAlmostEqual(merges[1]['distance'],.85)

    def test_downloads_and_html_for_sparse_and_empty_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);scopes=[dict(scope_id='lib01',scope_name='<unsafe>'),dict(scope_id='lib02')]
            result=write_assets('sample',scopes,{}, {'chr1':5000000},['chr1'],20,p)
            self.assertEqual(result['scales']['100000']['bins_used'],50)
            for svg in p.glob('*.svg'):
                ET.fromstring(svg.read_text());self.assertNotIn('nan',svg.read_text())
            for counts in p.glob('*.gz'):
                with gzip.open(counts,'rt') as f:self.assertEqual(len(f.readlines()),1)
            content=report_html('sample',p,'../assets/')
            self.assertIn('data-bin=\'100000\' data-method=\'spearman\' data-zeros=\'all\'><',content)
            self.assertIn('No downsampling',content);self.assertIn('&lt;unsafe&gt;',content)
            self.assertNotIn('<unsafe>',content);self.assertEqual(content.count("class='similarity-panel'"),12)
            self.assertEqual(report_html('absent',p),'')
            self.assertIsNone(write_assets('single',scopes[:1],{}, {},[],20,p))

if __name__=='__main__':unittest.main()
