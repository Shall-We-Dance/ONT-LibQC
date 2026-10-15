"""Low-memory BWA fallback stays consistent across validation, mounts and budgeting."""
import json,struct,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from ont_portal.catalog import public_catalog,job_config,bwa_mem_ready
from ont_portal.schema import Options
from ont_portal.worker import Worker

class BwaFallbackTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
  for name in ['genome.fa','genes.gtf','sizes']:(self.root/name).write_text('chr1\t100\n')
  for part in ['1','2','3','4','rev.1','rev.2']:(self.root/('bt.'+part+'.bt2')).touch()
  for part in ['0123','amb','ann','bwt.2bit.64','pac']:(self.root/('mem2.'+part)).touch()
  prefix=self.root/'classic';header=struct.pack('<5Q',1,50,100,150,200)
  Path(str(prefix)+'.ann').write_text('100 1 11\n');Path(str(prefix)+'.amb').touch()
  Path(str(prefix)+'.bwt').write_bytes(header);Path(str(prefix)+'.sa').write_bytes(header+struct.pack('<2Q',32,200)+bytes(48))
  Path(str(prefix)+'.pac').write_bytes(bytes(27))
  self.ref=dict(fasta=str(self.root/'genome.fa'),annotation=str(self.root/'genes.gtf'),chrom_sizes=str(self.root/'sizes'),
   bowtie2_index=str(self.root/'bt'),bwa_mem2_index=str(self.root/'mem2'),bwa_index=str(prefix),primary_chromosomes=['chr1'],
   microc_mapper='bwa-mem',disabled_mappers=['bwa-mem2'],availability_note='Memory limitation')
  self.catalog={'references':{'hg38':self.ref}}
 def options(self,**kwargs):
  return Options(sample='test',reference='hg38',assay='microc',mycoplasma=False,spike_in=False,**kwargs)
 def test_mem2_remains_default_and_classic_is_opt_in(self):
  self.ref.pop('microc_mapper');self.ref['disabled_mappers']=[]
  public=public_catalog(self.catalog)
  self.assertEqual(public['references'][0]['microc_mapper'],'bwa-mem2')
  self.options().check_catalog(public)
  default,_=job_config(self.options(),self.catalog,self.root,2)
  self.assertEqual(default['samples']['test']['mapper'],'bwa-mem2')
  self.options(mapper='bwa-mem').check_catalog(public)
  classic,_=job_config(self.options(mapper='bwa-mem'),self.catalog,self.root,2)
  self.assertEqual(classic['samples']['test']['mapper'],'bwa-mem')
  self.assertEqual(classic['references']['hg38']['microc_mapper'],'bwa-mem')
 def test_missing_mem2_requires_explicit_classic_selection(self):
  self.ref.pop('microc_mapper')
  public=public_catalog(self.catalog)
  self.assertEqual(public['references'][0]['microc_mapper'],'bwa-mem2')
  with self.assertRaisesRegex(ValueError,'Mapper is unavailable'):self.options().check_catalog(public)
  self.options(mapper='bwa-mem').check_catalog(public)
 def test_disabled_mem2_cannot_be_selected_and_auto_uses_classic(self):
  public=public_catalog(self.catalog);row=public['references'][0]
  self.assertTrue(row['ready']);self.assertEqual(row['mappers'],['bowtie2','bwa-mem']);self.assertEqual(row['microc_mapper'],'bwa-mem')
  self.options().check_catalog(public)
  with self.assertRaisesRegex(ValueError,'Mapper is unavailable'):self.options(mapper='bwa-mem2').check_catalog(public)
  cfg,mounts=job_config(self.options(),self.catalog,self.root,2)
  self.assertEqual(cfg['samples']['test']['mapper'],'bwa-mem');self.assertEqual(cfg['references']['hg38']['microc_mapper'],'bwa-mem')
  self.assertTrue(cfg['references']['hg38']['bwa_index'].startswith('/refs/'))
  self.assertTrue(all('readonly' in a for a in mounts.arguments() if a!='--mount'))
 def test_corrupt_packed_index_is_not_advertised(self):
  Path(self.ref['bwa_index']+'.ann').write_text('200 1 11\n')
  self.assertFalse(bwa_mem_ready(self.ref['bwa_index']))
  self.assertFalse(public_catalog(self.catalog)['references'][0]['ready'])
 def test_truncated_suffix_array_is_not_advertised(self):
  Path(self.ref['bwa_index']+'.sa').write_bytes(Path(self.ref['bwa_index']+'.sa').read_bytes()[:56])
  self.assertFalse(bwa_mem_ready(self.ref['bwa_index']))
 def test_classic_memory_and_minimap_auxiliary_use_classic_sizes(self):
  worker=Worker.__new__(Worker);worker.memory='5g';worker.cfg={'max_memory':'6g'};worker.catalog=self.catalog
  self.ref['minimap2_index']='short.mmi'
  sizes={'.bwt':3209286200,'.sa':1604643104,'.pac':802321528,'.mmi':1024}
  def fake_path(name):return SimpleNamespace(stat=lambda:SimpleNamespace(st_size=next(v for k,v in sizes.items() if name.endswith(k))))
  with patch('ont_portal.worker.Path',side_effect=fake_path):
   self.assertEqual(worker.memory_for(self.options()),6*1024**3)
   self.assertEqual(worker.memory_for(self.options(mapper='minimap2')),6*1024**3)
   sizes.update({'.bwt':2728222532,'.sa':1364111280,'.pac':682055614})
   self.assertEqual(worker.memory_for(self.options()),6*1024**3)
   worker.cfg['max_memory']='5g'
   with self.assertRaisesRegex(ValueError,'memory budget'):worker.memory_for(self.options())

if __name__=='__main__':unittest.main()
