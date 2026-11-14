"""Validate terminal correction and one shared downstream eligibility filter."""
import io
from pathlib import Path
import sys,tempfile,unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from chromatin_alignment import foldback_evidence, rc, collect_candidates, verified_remap, read_sam, eligible
from ont_qc import sam_primary_records, write_sparse_cpm_bedgraph
from report_provenance import analysis_label

class Process:
    def __init__(self,text):self.stdout=io.StringIO(text)
    def wait(self):return 0

class ChromatinAlignmentTests(unittest.TestCase):
    def test_requires_opposite_end_not_arbitrary_internal_match(self):
        end='ACGTGATCGTACGA';clip=rc(end)
        self.assertIsNotNone(foldback_evidence(clip,'GCTAGCTAGCTAGCTA'+end,'left'))
        self.assertIsNotNone(foldback_evidence(clip,end+'GCTAGCTAGCTAGCTA','right'))
        self.assertIsNone(foldback_evidence(clip,'G'*15+end+'T'*15,'left'))
        self.assertIsNone(foldback_evidence('A'*14,'T'*40,'left'))
        self.assertIsNone(foldback_evidence(clip[:13],end,'right'))

    def test_remap_must_keep_locus_strand_and_confidence(self):
        prior=dict(chrom='chr1',strand='+',start=10,end=60)
        corrected=dict(prior,mapq=20,aligned=50,coverage=1,identity=.98)
        self.assertTrue(verified_remap(prior,corrected,20))
        for change in [dict(mapq=19),dict(mapq=255),dict(start=16),dict(strand='-'),dict(chrom='chr2'),dict(identity=.89),dict(aligned=29),dict(coverage=.79)]:
            self.assertFalse(verified_remap(prior,dict(corrected,**change),20))

    def test_candidate_clip_fraction_and_other_terminal(self):
        mid='GCTAGCTAGCTAGCTA'+'ACGTGATCGTACGA';clip=rc(mid[-14:]);seq=clip+mid
        record=dict(chrom='chr1',mapq=20,aligned=len(mid),left=14,right=0)
        self.assertEqual(len(collect_candidates({'unit':record},{'unit':(seq,'I'*len(seq))},{'chr1'},20)),1)
        self.assertEqual(collect_candidates({'unit':dict(record,mapq=19)},{'unit':(seq,'I'*len(seq))},{'chr1'},20),{})
        self.assertEqual(collect_candidates({'unit':dict(record,right=4)},{'unit':(seq,'I'*len(seq))},{'chr1'},20),{})

    def test_reverse_strand_terminal_and_query_coverage(self):
        line='read\t16\tchr1\t11\t30\t30M14S\t*\t0\t0\t'+'A'*44+'\t'+'I'*44+'\tNM:i:0\n'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.sam';p.write_text(line);_,records=read_sam(p)
        r=records['read'];self.assertEqual((r['left'],r['right']),(14,0));self.assertAlmostEqual(r['coverage'],30/44);self.assertFalse(eligible(r))

    def test_metrics_and_coverage_exclude_same_rejected_alignment(self):
        lines=''.join(f'{name}\t0\tchr1\t1\t30\t40M\t*\t0\t0\t'+ 'A'*40+'\t'+'I'*40+f'\tZQ:i:{valid}\n' for name,valid in [('keep',1),('reject',0)])
        with patch('ont_qc.subprocess.Popen',side_effect=lambda *a,**k:Process(lines)):
            self.assertEqual(set(sam_primary_records('unused','samtools')),{'keep'})
            self.assertEqual(len(sam_primary_records('unused','samtools',include_ineligible=True)),2)
            with tempfile.TemporaryDirectory() as d:
                path=Path(d)/'a.bedgraph';count,covered=write_sparse_cpm_bedgraph('samtools','unused',20,{'chr1':100},{'chr1'},path)
                self.assertEqual(count,1);self.assertTrue(covered);self.assertIn('1000000',path.read_text())

    def test_report_revision_does_not_invent_a_release(self):
        self.assertEqual(analysis_label(dict(pipeline_version='v1.1.0',analysis_revision='r2',analysis_completed_utc='2026-01-02T00:00:00Z')),'v1.1.0 · r2 · 2026-01-02 UTC')

if __name__=='__main__':unittest.main()
