"""The full-insert CLI must record the BWA variant that produced its SAM."""
from pathlib import Path
import json,subprocess,sys,tempfile,unittest
SCRIPT=Path(__file__).resolve().parents[1]/'workflow/scripts/microc_ligation.py'
class BwaLigationProvenanceTests(unittest.TestCase):
 def test_cli_records_actual_mapper_and_preserves_default(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'in.fastq').write_text('');(p/'in.sam').write_text('@HD\tVN:1.6\n')
   for aligner in [None,'bwa-mem','bwa-mem2']:
    cmd=[sys.executable,str(SCRIPT),'--sample','empty','--fastq',str(p/'in.fastq'),'--sam',str(p/'in.sam'),'--primary','chr1','--output',str(p/'result.json')]
    if aligner:cmd+=['--aligner',aligner]
    subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL)
    result=json.loads((p/'result.json').read_text());expected=aligner or 'bwa-mem2'
    self.assertEqual(result['aligner'],expected)
    self.assertEqual(result['method'],'full_adapter_bounded_insert_'+expected.upper().replace('-','_')+'_split_alignment')
    self.assertEqual(result['ligation_screened_units'],0)
if __name__=='__main__':unittest.main()
