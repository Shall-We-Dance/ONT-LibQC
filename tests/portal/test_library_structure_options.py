"""Verify structure settings are validated, persisted and sent to the workflow."""
import copy
import json
import re
from pathlib import Path
import sys
import tempfile
import unittest
from fastapi.testclient import TestClient
from pydantic import ValidationError
from ont_portal.app import create_app
from ont_portal.schema import Options, schema_payload
from ont_portal.structure import LibraryStructure
from ont_portal.catalog import job_config
from ont_portal.settings import Settings

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'workflow/scripts'))
from library_structure import normalize_structure, handles
PUBLIC = {'references': [{'id': 'dm6', 'assembly': 'dm6', 'ready': True, 'mappers': ['bowtie2', 'bwa-mem2'], 'mycoplasma': False, 'spike_in': False, 'chromosomes': ['chr2L']}]}
BASE = {'sample': 'umi_test', 'assay': 'chromatin', 'reference': 'dm6'}
CUSTOM = {'mode': 'custom', 'adapter_i5_sequence': 'GCTACAGTCGATGACCTGATCGTAGCGA', 'adapter_i7_sequence': 'TCAGGATCGACTAGGCTACGTAGACTGACG', 'umis': [{'position': 'before_insert', 'length': 8}, {'position': 'after_i7', 'length': 12}]}

class StructureOptionsTests(unittest.TestCase):
    def test_old_jobs_default_to_auto(self):
        options = Options(**BASE)
        self.assertEqual(options.library_structure.mode, 'auto')
        self.assertEqual(options.library_structure.umis, [])
        self.assertEqual(normalize_structure(options.library_structure.model_dump()), normalize_structure())

    def test_api_and_pipeline_validation_agree(self):
        for value in [None, {'mode': 'auto'}, CUSTOM, {'mode': 'custom', 'i5_length': 0, 'i7_length': 24, 'adapter_family': 'nextera', 'p5_sequence': 'ACGT' * 8}]:
            with self.subTest(value=value):
                api = LibraryStructure(**(value or {})).model_dump()
                self.assertEqual(api, normalize_structure(value))
        for value in [{'mode': 'auto', 'umis': CUSTOM['umis']}, {'mode': 'custom', 'i5_length': True}, {'mode': 'custom', 'i7_length': 1.5}, {'mode': 'custom', 'p5_sequence': 'ACGT;unsafe'}, {'mode': 'custom', 'umis': [{'position': 'before_insert', 'length': 65}]}, {'mode': 'custom', 'umis': [{'position': 'before_insert', 'length': 8}] * 2}, {'mode': 'custom', 'umis': [{'position': 'inside_adapter', 'length': 8}]}]:
            with self.subTest(value=value):
                with self.assertRaises(ValidationError): LibraryStructure(**value)
                with self.assertRaises(ValueError): normalize_structure(value)

    def test_custom_adapter_api_and_pipeline_validation(self):
        for key in ('adapter_i5_sequence', 'adapter_i7_sequence'):
            for size in (16, 80):
                value={'mode':'custom', key:'A'*size}
                self.assertEqual(LibraryStructure(**value).model_dump(), normalize_structure(value))
            for value in ('', 'A'*15, 'A'*81, 'ACGT'*4+'N', 123, True):
                settings={'mode':'custom',key:value}
                with self.assertRaises(ValidationError):LibraryStructure(**settings)
                with self.assertRaises(ValueError):normalize_structure(settings)
            with self.assertRaises(ValidationError):LibraryStructure(mode='auto', **{key:'ACGT'*6})

    def test_schema_defaults_match_pipeline_handles(self):
        values=schema_payload()['library_structure']
        self.assertEqual((values['p5_sequence'], values['p7_sequence']), handles(normalize_structure()))

    def test_expected_index_pair_lengths_cannot_conflict(self):
        options=Options(**BASE, library_structure={'mode':'custom','i5_length':10}, index_pairs=[{'name':'a','i5':'ACGTACGT','i7':'TGCATGCA'}])
        with self.assertRaisesRegex(ValueError, 'configured length'):options.check_catalog(PUBLIC)

    def test_custom_layout_reaches_worker_sample_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ref={key:str(root/key) for key in ['fasta','annotation','chrom_sizes','bowtie2_index','bwa_mem2_index']};ref.update(assembly='dm6',primary_chromosomes=['chr2L'])
            for layout in [{'mode':'auto'},CUSTOM]:
                options=Options(**BASE,library_structure=layout)
                saved=Options.model_validate_json(options.model_dump_json())
                config,_=job_config(saved,{'references':{'dm6':ref}},root,2)
                self.assertEqual(config['samples']['umi_test']['library_structure'],normalize_structure(layout))

    def test_homepage_static_assets_are_served_by_the_real_app(self):
        with tempfile.TemporaryDirectory() as tmp:
            settings=Settings(state=Path(tmp),public_url='https://testserver',allowed_hosts=['testserver'],secret='test-only-'*8,submission_key='',public_catalog=PUBLIC,min_free_bytes=0)
            app=create_app(settings)
            with TestClient(app,base_url='https://testserver') as client:
                page=client.get('/');self.assertEqual(page.status_code,200)
                assets=set(re.findall(r'(?:src|href)="(/static/[^" ]+)"',page.text))
                self.assertTrue(any('library-structure.js' in path for path in assets))
                for path in assets:
                    with self.subTest(asset=path):
                        response=client.get(path);self.assertEqual(response.status_code,200)
                        if '.js' in path:self.assertIn('javascript',response.headers['content-type'])
                self.assertEqual(client.get('/static/config.yaml').status_code,404)

    def test_single_and_batch_api_store_per_file_layouts(self):
        with tempfile.TemporaryDirectory() as tmp:
            settings=Settings(state=Path(tmp),public_url='https://testserver',allowed_hosts=['testserver'],secret='test-only-'*8,submission_key='',public_catalog=PUBLIC,min_free_bytes=0)
            app=create_app(settings);app.state.store.heartbeat()
            with TestClient(app,base_url='https://testserver',headers={'X-Portal-Request':'1'}) as client:
                base={'email':'test@example.org','filename':'test.fastq','size':20,'consent':True,'data_policy_acknowledged':True}
                invalid=client.post('/api/jobs',json=base|{'options':BASE|{'library_structure':{'mode':'auto','umis':CUSTOM['umis']}}})
                self.assertEqual(invalid.status_code,422)
                batch=client.post('/api/batches',json={'files':[base|{'options':BASE|{'library_structure':CUSTOM}},base|{'options':BASE|{'sample':'auto_test'}}]})
                self.assertEqual(batch.status_code,201,batch.text)
                saved=[Options.model_validate_json(app.state.store.get(job['id'])['options']) for job in batch.json()['jobs']]
                self.assertEqual(saved[0].library_structure.model_dump(),normalize_structure(CUSTOM))
                self.assertEqual(saved[1].library_structure.mode,'auto')
                result=client.post('/api/jobs',json=base|{'options':BASE|{'library_structure':CUSTOM}})
                self.assertEqual(result.status_code,201,result.text)
                job=app.state.store.get(result.json()['id']);self.assertEqual(Options.model_validate_json(job['options']).library_structure.umis[1].length,12)

if __name__=='__main__':unittest.main()
