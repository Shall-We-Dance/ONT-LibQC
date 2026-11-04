"""Canonical reference defaults retain sex chromosomes present in the manifest."""
import json
from pathlib import Path
import unittest

import yaml
import ont_portal


SOURCE = Path(ont_portal.__file__).resolve().parents[2]


class ReferenceChromosomeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.common = yaml.safe_load(
            (SOURCE / 'portal/deploy/references.common.yaml').read_text()
        )['references']

    def test_common_xy_references_include_both_sex_chromosomes(self):
        # These are assembly defaults, not requirements for custom user subsets.
        for assembly in ('hg38', 'mm39', 'rn8', 'dm6'):
            with self.subTest(assembly=assembly):
                chromosomes = self.common[assembly]['primary_chromosomes']
                self.assertIn('chrX', chromosomes)
                self.assertIn('chrY', chromosomes)

    def test_example_chromosomes_match_common_reference_manifest(self):
        for name in ('config.example.yaml', 'portal/deploy/references.example.yaml'):
            references = yaml.safe_load((SOURCE / name).read_text())['references']
            for assembly, reference in references.items():
                with self.subTest(config=name, assembly=assembly):
                    self.assertCountEqual(
                        reference['primary_chromosomes'],
                        self.common[assembly]['primary_chromosomes'],
                    )

    def test_metric_chromosome_metadata_matches_common_reference_manifest(self):
        definitions = json.loads(
            (SOURCE / 'workflow/scripts/metric_definitions.json').read_text()
        )
        for assembly, chromosomes in definitions['primary_chromosomes'].items():
            with self.subTest(assembly=assembly):
                self.assertCountEqual(
                    chromosomes, self.common[assembly]['primary_chromosomes'],
                )


if __name__ == '__main__':
    unittest.main()
