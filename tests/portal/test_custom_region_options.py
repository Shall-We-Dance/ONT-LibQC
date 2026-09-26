"""Keep upload controls, persisted options and report plot modes consistent."""
import copy
from pathlib import Path
import tempfile
import unittest

from fastapi.testclient import TestClient
from pydantic import ValidationError
from ont_portal.app import create_app
from ont_portal.catalog import job_config
from ont_portal.schema import Options
from ont_portal.settings import Settings

PUBLIC = {'references': [{'id': 'dm6', 'assembly': 'dm6', 'ready': True,
                          'mappers': ['bowtie2', 'bwa-mem2'], 'mycoplasma': False,
                          'spike_in': False, 'chromosomes': ['chr2L']}]}
BASE = dict(sample='metaplot_test', assay='chromatin', reference='dm6', custom_regions=True)


class CustomRegionOptionsTests(unittest.TestCase):
    def test_new_submissions_require_matching_coordinates(self):
        for assembly in (None, 'hg38'):
            with self.subTest(assembly=assembly), self.assertRaisesRegex(ValueError, 'dm6 coordinates'):
                Options(**BASE, regions_assembly=assembly).check_catalog(PUBLIC, require_regions_confirmation=True)
        options = Options(**BASE, regions_assembly='dm6')
        self.assertIs(options.check_catalog(PUBLIC, require_regions_confirmation=True), options)

    def test_existing_jobs_keep_auto_without_new_confirmation(self):
        legacy = Options.model_validate(BASE)
        self.assertEqual(legacy.regions_mode, 'auto')
        self.assertIs(legacy.check_catalog(PUBLIC), legacy)

    def test_explicit_confirmation_is_checked_again_by_worker(self):
        with self.assertRaisesRegex(ValueError, 'dm6 coordinates'):
            Options(**BASE, regions_assembly='hg38').check_catalog(PUBLIC)

    def test_no_custom_bed_needs_no_confirmation(self):
        options = Options(**(BASE | {'custom_regions': False}))
        self.assertIs(options.check_catalog(PUBLIC, require_regions_confirmation=True), options)

    def test_reference_alias_uses_assembly_for_confirmation(self):
        catalog = copy.deepcopy(PUBLIC)
        catalog['references'][0]['id'] = 'fly'
        options = Options(**(BASE | {'reference': 'fly'}), regions_assembly='dm6')
        self.assertIs(options.check_catalog(catalog, require_regions_confirmation=True), options)

    def test_unknown_modes_are_rejected(self):
        with self.assertRaises(ValidationError):
            Options(**BASE, regions_mode='unsupported')

    def test_each_plot_mode_reaches_pipeline_configuration(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = {key: str(root / key) for key in ('fasta', 'annotation', 'chrom_sizes', 'bowtie2_index', 'bwa_mem2_index')}
            ref.update(assembly='dm6', primary_chromosomes=['chr2L'])
            for mode in ('auto', 'reference_point', 'region'):
                with self.subTest(mode=mode):
                    options = Options(**BASE, regions_mode=mode, regions_assembly='dm6')
                    options = Options.model_validate_json(options.model_dump_json())
                    config, _ = job_config(options, {'references': {'dm6': ref}}, root, 2)
                    self.assertEqual(config['samples']['metaplot_test']['custom_regions'],
                                     {'bed': '/job/regions.bed', 'assembly': 'dm6',
                                      'label': 'Uploaded regions', 'mode': mode})

    def test_single_and_batch_api_validate_and_persist_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            settings = Settings(state=Path(tmp), public_url='https://testserver',
                                allowed_hosts=['testserver'], secret='test-only-' * 8,
                                submission_key='', public_catalog=PUBLIC, min_free_bytes=0)
            app = create_app(settings)
            app.state.store.heartbeat()
            with TestClient(app, base_url='https://testserver', headers={'X-Portal-Request': '1'}) as client:
                for endpoint in ('/api/jobs', '/api/batches'):
                    for assembly in (None, 'hg38', 'dm6'):
                        with self.subTest(endpoint=endpoint, assembly=assembly):
                            options = BASE | {'regions_mode': 'reference_point', 'regions_assembly': assembly}
                            submission = {'email': 'test@example.org', 'filename': 'test.fastq', 'size': 10,
                                          'options': options, 'consent': True}
                            body = {'files': [submission]} if endpoint.endswith('batches') else submission
                            response = client.post(endpoint, json=body)
                            if assembly != 'dm6':
                                self.assertEqual(response.status_code, 422)
                                self.assertIn('dm6 coordinates', response.json()['detail'])
                            else:
                                self.assertEqual(response.status_code, 201, response.text)
                                result = response.json()
                                job = result['jobs'][0] if 'jobs' in result else result
                                saved = Options.model_validate_json(app.state.store.get(job['id'])['options'])
                                self.assertEqual((saved.regions_mode, saved.regions_assembly), ('reference_point', 'dm6'))


if __name__ == '__main__':
    unittest.main()
