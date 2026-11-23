"""The portal advertises and mounts only supported single-cell RNA resources."""
import hashlib
import json
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch

from ont_portal.catalog import BT2_PARTS, job_config, public_catalog
from ont_portal.schema import Options
from ont_portal.scrna_catalog import barcodes_ready, file_digest
from test_rna_beta import make_catalog


class ScrnaCatalogTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.catalog = make_catalog(self.root)
        self.ref = self.catalog['references']['synthetic']
        self.index = Path(self.ref['star_index'])
        self.manifest = self.index / 'ontlibqc-star-index.json'
        self.meta = json.loads(self.manifest.read_text())
        self.meta.update(reference_id='synthetic', source_sha256={'annotation': file_digest(self.ref['annotation'])})
        self.manifest.write_text(json.dumps(self.meta))
        self.barcodes = self.root / 'barcodes'
        self.barcodes.mkdir()
        public_data = b'synthetic barcode list'
        (self.barcodes / 'fixture.gz').write_bytes(public_data)
        self.catalog['scrna_barcodes'] = str(self.barcodes)
        sums = patch('ont_portal.scrna_catalog.barcode_hashes',
                     return_value={'fixture.gz': hashlib.sha256(public_data).hexdigest()})
        sums.start()
        self.addCleanup(sums.stop)

    def supports(self):
        public = public_catalog(self.catalog)
        self.assertNotIn(str(self.root), json.dumps(public))
        row = public['references'][0]
        self.assertEqual(row['scrna_ready'], 'scrna' in row['assays'])
        return row['scrna_ready']

    def build(self):
        opts = Options(sample='Synthetic', assay='scrna', reference='synthetic', spike_in=False, mycoplasma=False)
        return job_config(opts, self.catalog, self.root, 2)

    def test_readiness_requires_matching_barcode_reference_and_gtf(self):
        self.assertTrue(self.supports())
        target = self.barcodes / 'fixture.gz'
        target.write_bytes(b'corrupt synthetic list')
        self.assertFalse(barcodes_ready(self.barcodes))
        self.assertFalse(self.supports())
        with self.assertRaises(ValueError):
            self.build()
        target.write_bytes(b'synthetic barcode list')
        self.assertTrue(self.supports())
        self.manifest.write_text(json.dumps(dict(self.meta, reference_id='different')))
        self.assertFalse(self.supports())
        self.manifest.write_text(json.dumps(self.meta))
        Path(self.ref['annotation']).write_text('different source GTF')
        self.assertFalse(self.supports())
        self.assertIn('rnaseq', public_catalog(self.catalog)['references'][0]['assays'])

    def test_allowlisted_readonly_files_and_species_aliases(self):
        (self.index / 'do-not-mount').write_text('unrelated')
        (self.barcodes / 'do-not-mount').write_text('unrelated')
        host = self.ref['bowtie2_index']
        competitor = self.root / 'other'
        for part in BT2_PARTS:
            (self.root / f'other.{part}.bt2l').touch()
        self.catalog['species_composition'] = {
            'indexes': {'host': host, 'competitor': str(competitor)},
            'assembly_aliases': {'synthetic': 'host'},
        }
        cfg, mounts = self.build()
        settings = cfg['scrna']
        self.assertEqual(settings['assembly'], 'synthetic')
        self.assertEqual(settings['species_panel'], '/job/scrna-species.json')
        panel = json.loads((self.root / 'scrna-species.json').read_text())
        self.assertNotIn(str(self.root), json.dumps([cfg, panel]))
        self.assertEqual(panel['assembly_aliases'], {'synthetic': 'host'})
        self.assertEqual(set(panel['indexes']), {'host', 'competitor'})
        for field, suffix in [('host', 'bt2'), ('competitor', 'bt2l')]:
            for part in BT2_PARTS:
                self.assertIn(f"{panel['indexes'][field]}.{part}.{suffix}", mounts.exact)
        self.assertTrue(all(source.is_file() for source in mounts.exact.values()))
        self.assertIn(settings['star_index'] + '/Genome', mounts.exact)
        self.assertIn(settings['barcodes'] + '/fixture.gz', mounts.exact)
        self.assertFalse(mounts.parents)
        args = mounts.arguments()
        self.assertTrue(all(arg.endswith(',readonly') for arg in args if arg != '--mount'))
        self.assertNotIn('do-not-mount', str(args))
        (self.root / 'other.1.bt2l').unlink()
        self.assertNotIn('species_panel', self.build()[0]['scrna'])

    def test_missing_optional_resource_disables_only_scrna(self):
        self.ref['rna_family_gtf'] = str(self.root / 'missing.gtf')
        self.assertFalse(self.supports())
        self.ref.pop('rna_family_gtf')
        self.ref['rna_source_panel'] = str(self.root)
        (self.root / 'panel.fa').write_bytes(b'>r\nACGT\n')
        self.assertFalse(self.supports())
        (self.root / 'panel.json').write_text(json.dumps({'fasta_sha256': file_digest(self.root / 'panel.fa'),
                                                         'references': {'r': {'source_class': 'rRNA', 'taxon': 'Example'}}}))
        self.assertTrue(self.supports())
        cfg, mounts = self.build()
        self.assertIn(cfg['scrna']['rna_source_panel'] + '/panel.json', mounts.exact)

    def test_enable_command_writes_both_catalogs_and_preserves_source(self):
        preparation = Path(__file__).resolve().parents[2] / 'portal/scripts/prepare_references.py'
        enable = runpy.run_path(str(preparation))['enable_scrna']
        del self.catalog['scrna_barcodes']
        original = self.root / 'original.json'
        original.write_text(json.dumps(self.catalog))
        out = self.root / 'configured'
        enable(original, self.barcodes, out)
        admin = json.loads((out / 'catalog.admin.json').read_text())
        public = (out / 'catalog.public.json').read_text()
        self.assertEqual(admin['scrna_barcodes'], str(self.barcodes))
        self.assertTrue(json.loads(public)['references'][0]['scrna_ready'])
        self.assertNotIn(str(self.root), public)
        self.assertEqual(self.catalog, json.loads(original.read_text()))
        (self.barcodes / 'fixture.gz').unlink()
        with self.assertRaisesRegex(ValueError, 'barcode'):
            enable(original, self.barcodes, out)
        self.assertEqual(public, (out / 'catalog.public.json').read_text())
