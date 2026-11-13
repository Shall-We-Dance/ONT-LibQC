"""TE annotation is an operator reference, mounted read-only for the report stage."""
from pathlib import Path
import tempfile,unittest
from ont_portal.catalog import job_config
from ont_portal.schema import Options

class TEAnnotationTests(unittest.TestCase):
    def test_te_reference_reaches_pipeline_as_readonly_mount(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);annotations=root/'annotations';annotations.mkdir()
            te=annotations/'te.bed.gz';te.write_bytes(b'test')
            ref={k:str(root/k) for k in ('fasta','annotation','chrom_sizes','bowtie2_index')}
            ref.update(assembly='test',primary_chromosomes=['chr1'],te_bed=str(te),te_source='https://example.org/repeats')
            options=Options(sample='synthetic',assay='chromatin',reference='test',mapper='bowtie2',mycoplasma=False,spike_in=False)
            cfg,mounts=job_config(options,{'references':{'test':ref}},root,2)
            self.assertTrue(cfg['references']['test']['te_bed'].startswith('/refs/'))
            self.assertEqual(cfg['references']['test']['te_source'],'https://example.org/repeats')
            self.assertTrue(any(str(annotations) in arg and arg.endswith(',readonly') for arg in mounts.arguments()))
            ref.pop('te_bed');ref.pop('te_source')
            cfg,_=job_config(options,{'references':{'test':ref}},root,2)
            self.assertNotIn('te_bed',cfg['references']['test'])

if __name__=='__main__':unittest.main()
