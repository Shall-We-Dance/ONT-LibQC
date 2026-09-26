"""Micro-C's BWA analyses need enough memory even when native mapping uses minimap2."""
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from ont_portal.worker import Worker

class MicrocWorkerBudgetTests(unittest.TestCase):
    def worker(self):
        worker=Worker.__new__(Worker)
        worker.memory='1g';worker.cfg={'max_memory':'8g'}
        worker.catalog={'references':{'reference':{'bwa_mem2_index':'bwa','minimap2_index':'minimap.mmi'}}}
        return worker

    @staticmethod
    def path(name):
        size=1024**2 if name.endswith('.mmi') else 1024**3
        return SimpleNamespace(stat=lambda:SimpleNamespace(st_size=size))

    @patch('ont_portal.worker.Path',side_effect=path)
    def test_microc_budget_covers_bwa_with_minimap_native_mapping(self,_path):
        worker=self.worker()
        options=SimpleNamespace(mapper='minimap2',assay='microc',reference='reference')
        self.assertEqual(worker.memory_for(options),6*1024**3)
        options.assay='chromatin'
        self.assertEqual(worker.memory_for(options),4*1024**3)

    @patch('ont_portal.worker.Path',side_effect=path)
    def test_bwa_requirement_cannot_bypass_configured_memory_cap(self,_path):
        worker=self.worker();worker.cfg['max_memory']='5g'
        options=SimpleNamespace(mapper='minimap2',assay='microc',reference='reference')
        with self.assertRaisesRegex(ValueError,'memory budget'):worker.memory_for(options)

if __name__=='__main__':unittest.main()
