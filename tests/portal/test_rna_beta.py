"""RNA/CLIP admission must select a validated spliced mapper and reserve its RAM."""
import copy
import json
import os
from pathlib import Path
import tempfile
import unittest

from fastapi.testclient import TestClient
from ont_portal.app import create_app
from ont_portal.catalog import job_config, public_catalog
from ont_portal.schema import Options, schema_payload
from ont_portal.mappers import resolve_mapper
from ont_portal.settings import Settings
from ont_portal.star_resources import STAR_FILES, GIB, star_ready
from ont_portal.worker import Worker


def make_catalog(root):
    for name, text in [('genome.fa', '>chr1\nACGT\n'), ('genes.gtf', 'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g";\n'), ('sizes', 'chr1\t4\n')]:
        (root / name).write_text(text)
    for part in ('1', '2', '3', '4', 'rev.1', 'rev.2'):
        (root / ('bt.' + part + '.bt2')).touch()
    for part in ('0123', 'amb', 'ann', 'bwt.2bit.64', 'pac'):
        (root / ('bwa.' + part)).touch()
    index = root / 'star'; index.mkdir()
    for name in STAR_FILES:
        (index / name).write_text('synthetic-index\n')
    (index / 'genomeParameters.txt').write_text('versionGenome 2.7.4a\ngenomeSAsparseD 4\ngenomeSAindexNbases 2\n')
    marker = dict(star_version='2.7.11b', parameters=dict(genomeSAsparseD=4, genomeSAindexNbases=2),
                  index_files={name: (index / name).stat().st_size for name in STAR_FILES}, star_alignment_memory_bytes=4*GIB)
    (index / 'ontlibqc-star-index.json').write_text(json.dumps(marker))
    ref = dict(fasta=str(root/'genome.fa'), annotation=str(root/'genes.gtf'), chrom_sizes=str(root/'sizes'),
               bowtie2_index=str(root/'bt'), bwa_mem2_index=str(root/'bwa'), star_index=str(index),
               primary_chromosomes=['chr1'], assembly='synthetic')
    return {'references': {'synthetic': ref}}


class RnaBetaTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(); self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name); self.catalog = make_catalog(self.root)
        self.ref = self.catalog['references']['synthetic']

    def options(self, assay='rnaseq', **values):
        return Options(sample='Synthetic', reference='synthetic', assay=assay, **values)

    def test_screens_default_on_with_panels_and_respect_opt_out(self):
        public = public_catalog(self.catalog)
        public['references'][0].update(mycoplasma=True, spike_in=True)
        for assay in ('chromatin', 'microc', 'rnaseq', 'clip'):
            options = self.options(assay).check_catalog(public)
            self.assertTrue(options.mycoplasma)
            self.assertTrue(options.spike_in)
            self.assertFalse(options.declared_spike_in)
            disabled = self.options(assay, mycoplasma=False, spike_in=False).check_catalog(public)
            self.assertFalse(disabled.mycoplasma)
            self.assertFalse(disabled.spike_in)
        configured = copy.deepcopy(self.catalog)
        configured['references']['synthetic'].update(mycoplasma_index=str(self.root/'bt'),
            spike_in={'combined_bowtie2_index':str(self.root/'bt'), 'mode':'configured', 'target_contig_regex':'^spike_'})
        enabled = self.options().check_catalog(public_catalog(configured))
        cfg, _ = job_config(enabled, configured, self.root, 2)
        self.assertTrue(cfg['mycoplasma_screen']['enabled'])
        self.assertTrue(cfg['spike_in']['enabled'])
        unavailable = self.options().check_catalog(public_catalog(self.catalog))
        self.assertFalse(unavailable.mycoplasma)
        self.assertFalse(unavailable.spike_in)
        with self.assertRaisesRegex(ValueError, 'Mycoplasma reference is unavailable'):
            self.options(mycoplasma=True).check_catalog(public_catalog(self.catalog))
        with self.assertRaisesRegex(ValueError, 'Spike-in reference is unavailable'):
            self.options(spike_in=True).check_catalog(public_catalog(self.catalog))

    def test_production_microc_fallback_is_preserved(self):
        self.assertEqual(resolve_mapper('auto', 'microc', {'mappers': ['bowtie2', 'bwa-mem'], 'microc_mapper': 'bwa-mem'}), 'bwa-mem')
        self.assertEqual(resolve_mapper('auto', 'microc', {'mappers': ['bwa-mem', 'bwa-mem2']}), 'bwa-mem2')
        self.assertEqual(resolve_mapper('auto', 'microc', {'mappers': ['bwa-mem2'], 'microc_mapper': 'bwa-mem'}), 'bwa-mem')

    def test_ready_catalog_accepts_both_beta_assays_without_metadata(self):
        public = public_catalog(self.catalog); row = public['references'][0]
        self.assertTrue(row['ready']); self.assertTrue(row['rna_ready'])
        self.assertEqual(row['assays'], ['chromatin', 'microc', 'rnaseq', 'clip'])
        self.assertNotIn(str(self.root), json.dumps(public))
        for assay in ('rnaseq', 'clip', 'chromatin', 'microc'):
            self.options(assay).check_catalog(public)
        for assay in ('rnaseq', 'clip'):
            with self.assertRaisesRegex(ValueError, 'Mapper is unavailable'):
                self.options(assay, mapper='bowtie2').check_catalog(public)
        with self.assertRaisesRegex(ValueError, 'Mapper is unavailable'):
            self.options('chromatin', mapper='starlong').check_catalog(public)

    def test_missing_or_stale_index_blocks_only_rna(self):
        (Path(self.ref['star_index'])/'SA').write_text('truncated')
        public = public_catalog(self.catalog)
        self.assertFalse(public['references'][0]['rna_ready'])
        self.options('chromatin').check_catalog(public)
        with self.assertRaisesRegex(ValueError, 'RNA/CLIP beta is unavailable'):
            self.options().check_catalog(public)

    def test_wrong_version_sparse_layout_and_missing_marker_are_rejected(self):
        directory = Path(self.ref['star_index']); path = directory/'ontlibqc-star-index.json'
        original = json.loads(path.read_text())
        for mutation in ({'star_version': '2.7.10a'}, {'parameters': {'genomeSAsparseD': 1, 'genomeSAindexNbases': 2}}, {'index_files': {}}):
            path.write_text(json.dumps(original | mutation)); self.assertFalse(star_ready(directory))
        path.unlink(); self.assertFalse(star_ready(directory))

    def test_rna_reference_does_not_require_dna_mapper_indexes(self):
        self.ref.pop('bowtie2_index'); self.ref.pop('bwa_mem2_index')
        public = public_catalog(self.catalog)
        self.assertFalse(public['references'][0]['ready']); self.assertTrue(public['references'][0]['rna_ready'])
        self.options().check_catalog(public)
        with self.assertRaisesRegex(ValueError, 'Reference is unavailable'):
            self.options('chromatin').check_catalog(public)

    def test_rna_job_mounts_index_and_preserves_explicit_overrides(self):
        for assay in ('rnaseq', 'clip'):
            cfg, mounts = job_config(self.options(assay).check_catalog(public_catalog(self.catalog)), self.catalog, self.root, 2)
            self.assertEqual(cfg['samples']['Synthetic']['mapper'], 'starlong')
            self.assertEqual(cfg['samples']['Synthetic']['type'], assay)
            self.assertTrue(cfg['references']['synthetic']['star_index'].startswith('/refs/'))
            self.assertFalse(cfg['spike_in']['enabled']); self.assertFalse(cfg['mycoplasma_screen']['enabled'])
            for metric in ('primary_mapping_rate', 'usable_fragment_rate', 'mitochondrial_fraction'):
                self.assertIsNone(cfg['qc_reference_ranges'][metric])
            self.assertTrue(all('readonly' in arg for arg in mounts.arguments() if arg != '--mount'))
        cfg, _ = job_config(self.options(qc_reference_ranges={'primary_mapping_rate': {'min': .4}}).check_catalog(public_catalog(self.catalog)), self.catalog, self.root, 2)
        self.assertEqual(cfg['qc_reference_ranges']['primary_mapping_rate']['min'], .4)
        cfg, _ = job_config(self.options('chromatin').check_catalog(public_catalog(self.catalog)), self.catalog, self.root, 2)
        self.assertEqual(cfg['samples']['Synthetic']['mapper'], 'bowtie2')
        self.assertEqual(cfg['qc_reference_ranges']['mitochondrial_fraction']['max'], .2)

    def test_resource_budget_cannot_be_reduced_below_index_requirement(self):
        worker = Worker.__new__(Worker); worker.memory='2g'; worker.cfg={'max_memory':'8g'}; worker.catalog=self.catalog
        self.assertEqual(worker.memory_for(self.options()), 4*GIB)
        self.ref['star_alignment_memory_bytes'] = 6*GIB
        self.assertEqual(worker.memory_for(self.options('clip')), 6*GIB)
        worker.cfg['max_memory']='5g'
        with self.assertRaisesRegex(ValueError, 'memory budget'):
            worker.memory_for(self.options())
        self.ref['star_alignment_memory_bytes'] = -1
        with self.assertRaisesRegex(ValueError, 'Invalid STAR'):
            worker.memory_for(self.options())

    def test_api_acceptance_keeps_private_access_and_data_policy(self):
        settings = Settings(state=self.root/'state', public_url='https://testserver', allowed_hosts=['testserver'],
                            secret='test-only-'*8, submission_key='', public_catalog=public_catalog(self.catalog), min_free_bytes=0)
        app=create_app(settings); app.state.store.heartbeat()
        with TestClient(app, base_url='https://testserver', headers={'X-Portal-Request':'1'}) as client:
            for assay in ('rnaseq', 'clip'):
                payload={'filename':'synthetic.fastq','size':32,'consent':True,'data_policy_acknowledged':True,
                         'options':{'sample':'Synthetic','reference':'synthetic','assay':assay}}
                response=client.post('/api/jobs', json=payload)
                self.assertEqual(response.status_code, 201, response.text)
                job=response.json(); saved=json.loads(app.state.store.get(job['id'])['options'])
                self.assertEqual(saved['assay'], assay); self.assertEqual(saved['mapper'], 'auto')
                self.assertEqual(client.get('/api/jobs/'+job['id']).status_code, 404)
                self.assertEqual(client.get('/api/jobs/'+job['id'],headers={'Authorization':'Bearer '+job['token']}).status_code, 200)
                payload['data_policy_acknowledged']=False
                self.assertEqual(client.post('/api/jobs',json=payload).status_code, 422)


class RnaBetaBrowserTests(unittest.TestCase):
    def test_form_beta_defaults_availability_and_dna_restore(self):
        try:
            from playwright.sync_api import sync_playwright
        except ImportError:
            self.skipTest('Install Playwright to run the RNA portal browser acceptance test')
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); source=Path(__file__).resolve().parents[2]/'portal/static'
            catalog=public_catalog(make_catalog(root)) | schema_payload() | dict(max_files=8,max_upload_bytes=1000000,
                max_expanded_bytes=2000000,max_reads=100000,retention_days=7,mail_enabled=False,worker_online=True,
                submission_key_required=False,pipeline_url='https://github.com/Shall-We-Dance/ONT-LibQC')
            unavailable=copy.deepcopy(catalog['references'][0]); unavailable.update(id='dna_only',label='DNA only',rna_ready=False,
                assays=['chromatin','microc'],mappers=['bowtie2','bwa-mem2'])
            catalog['references'].append(unavailable)
            catalog['references'][0].update(mycoplasma=True,spike_in=True)
            errors=[]
            with sync_playwright() as playwright:
                # Use a downloaded browser through Playwright's normal installation path.
                engine = os.environ.get('ONT_TEST_BROWSER', 'chromium')
                executable = os.environ.get('PLAYWRIGHT_CHROMIUM_EXECUTABLE' if engine == 'chromium' else 'ONT_WEBKIT_EXECUTABLE')
                browser = getattr(playwright, engine).launch(headless=True, **({'executable_path': executable} if executable else {}))
                try:
                    page=browser.new_page(); page.on('pageerror', lambda error: errors.append(str(error)))
                    summaries={}; authorizations=[]
                    def serve(route):
                        from urllib.parse import urlsplit
                        path=urlsplit(route.request.url).path
                        if path=='/api/catalog': route.fulfill(json=catalog)
                        elif path in summaries:
                            authorizations.append((path, route.request.headers.get('authorization')))
                            route.fulfill(json=summaries[path])
                        else:
                            target=source/('index.html' if path=='/' else path.removeprefix('/static/'))
                            if target.is_file(): route.fulfill(path=target)
                            else: route.fulfill(status=404,body='')
                    page.route('**/*', serve); page.goto('http://portal.test/')
                    page.locator('#fastq-input').set_input_files({'name':'synthetic.fastq','mimeType':'text/plain','buffer':b'@synthetic\nACGT\n+\nIIII\n'})
                    page.get_by_label('Reference genomeRequired', exact=True).select_option('synthetic')
                    page.get_by_label('Library assayRequired', exact=True).select_option('rnaseq')
                    result=page.evaluate('optionsFor(files[0])')
                    self.assertEqual(result['assay'], 'rnaseq'); self.assertEqual(result['mapper'], 'auto')
                    self.assertTrue(result['mycoplasma']); self.assertTrue(result['spike_in'])
                    self.assertIsNone(result['qc_reference_ranges']['mitochondrial_fraction'])
                    self.assertIn('RNA-seq (Beta)', page.get_by_label('Library assayRequired', exact=True).inner_text())
                    self.assertEqual(page.get_by_label('Alignment method').locator('option').all_text_contents(), ['Default · STARlong','STARlong'])
                    page.evaluate('files[0].myco.input.click(); files[0].spike.input.click()')
                    page.get_by_label('Library assayRequired', exact=True).select_option('clip')
                    self.assertFalse(page.evaluate('optionsFor(files[0]).mycoplasma'))
                    self.assertFalse(page.evaluate('optionsFor(files[0]).spike_in'))
                    self.assertEqual(page.evaluate('optionsFor(files[0]).assay'), 'clip')
                    page.get_by_label('Library assayRequired', exact=True).select_option('chromatin')
                    restored=page.evaluate('optionsFor(files[0])')
                    self.assertTrue(restored['mycoplasma']); self.assertTrue(restored['spike_in'])
                    self.assertEqual(restored['qc_reference_ranges']['mitochondrial_fraction']['max'], .2)
                    page.get_by_label('Reference genomeRequired', exact=True).select_option('dna_only')
                    self.assertTrue(page.get_by_label('Library assayRequired', exact=True).locator('option[value="rnaseq"]').is_disabled())
                    self.assertEqual(errors, [])
                    genes={f'g{i}':i for i in range(1,31)}
                    jobs=[]
                    for i, sample in enumerate(('SyntheticA','SyntheticB','OtherAnnotation'),1):
                        jid=f'{i:032x}'; token=f'{i:064x}'
                        counts={gene:count*i for gene,count in genes.items()}
                        summaries[f'/reports/{jid}/data/rna/{sample}/summary.json']={
                            'schema_version':1,'sample':sample,'assay':'rnaseq','assembly':'synthetic',
                            'annotation':{'sha256':('a' if i<3 else 'b')*64},
                            'scopes':{'total':{'gene_assigned_units':sum(counts.values()),'gene_counts':counts}}}
                        jobs.append({'id':jid,'token':token,'job':{'id':jid,'sample':sample,'assay':'rnaseq','reference':'synthetic',
                            'status':'succeeded','expires':4102444800,'mail_status':'disabled','filename':sample+'.fastq'}})
                    page.evaluate('(jobs) => { tasks=jobs; renderTasks(); }', jobs)
                    self.assertEqual(page.locator('#rna-consistency').count(), 0)
                    self.assertEqual(page.evaluate('typeof window.RnaConsistency'), 'undefined')
                    self.assertEqual(authorizations, [])
                    self.assertEqual(errors, [])
                finally:
                    browser.close()


if __name__=='__main__':
    unittest.main(verbosity=2)
