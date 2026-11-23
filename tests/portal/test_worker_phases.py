import unittest
from multiprocessing import Pipe
from unittest.mock import Mock

from ont_portal.run_job import phase_targets
from ont_portal.worker_phases import PhasePool


class PhaseTests(unittest.TestCase):
    def test_three_preparations_can_wait_and_all_finish(self):
        pool = PhasePool(None)
        children = {}
        for jid in ('a', 'b', 'c'):
            parent, child = Pipe()
            self.addCleanup(parent.close)
            self.addCleanup(child.close)
            pool.start(jid, 8, parent)
            children[jid] = child
            child.send(12)
        fits = lambda reservations, budget: sum(reservations.values()) + budget + 4 <= 31
        pool.advance(fits)
        self.assertEqual(pool.reserved, {'a':12, 'b':12, 'c':0})
        self.assertIsNone(children['a'].recv())
        self.assertIsNone(children['b'].recv())
        self.assertFalse(children['c'].poll())
        children['a'].send(8)
        children['b'].send(8)
        pool.advance(fits)
        self.assertIsNone(children['c'].recv())
        self.assertEqual(sum(pool.reserved.values()), 20)
        pool.release('c')
        pool.advance(fits)
        self.assertEqual(pool.reserved, {'a':8, 'b':8})

    def test_waiting_child_releases_cpu_and_can_be_cancelled(self):
        cpu = Mock()
        cpu.owners = {'a':1}
        parent, child = Pipe()
        self.addCleanup(child.close)
        pool = PhasePool(cpu)
        pool.start('a', 8, parent)
        child.send(24)
        pool.advance(lambda reservations, budget: False)
        cpu.release.assert_called_once_with('a')
        self.assertFalse(child.poll())
        cpu.owners.clear()
        pool.release('a')
        self.assertEqual((pool.waiting, pool.reserved), ({}, {}))

    def test_index_consumers_are_in_the_alignment_phase(self):
        cfg = {'output_dir':'.work', 'species_composition':{'indexes':{'dm6':'index'}}}
        for assay in ('chromatin', 'microc', 'rnaseq', 'clip'):
            options = {'sample':'sample', 'assay':assay}
            targets = phase_targets('align',cfg,options)
            self.assertIn('.work/qc/species_composition/sample.json', targets)
            self.assertIn('.work/bam/sample.bam.bai', targets)
            self.assertEqual('.work/qc/pairtools/sample.json' in targets, assay=='microc')
            self.assertEqual('.work/qc/ligation/sample.json' in targets, assay=='microc')
            self.assertEqual('.work/qc/rna/sample' in targets, assay in ('rnaseq','clip'))
            self.assertEqual(phase_targets('report',cfg,options), [])
        self.assertNotIn('.work/qc/species_composition/sample.json',
                         phase_targets('align', {'output_dir':'.work'}, options))
