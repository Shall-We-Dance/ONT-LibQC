import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
import ont_qc as report


class LibraryQualityTests(unittest.TestCase):
    def setUp(self):
        self.scope = dict(scope_id='lib01', scope_name='lib01', raw_reads=2500,
                          library_units=2500, sequencing_yield_call='pass',
                          adapter_dimer_rate=.01, primary_mapping_rate=.975,
                          usable_fragment_rate=.28, mitochondrial_fraction=.60)

    def render(self, scope=None, ranges=None, assay='chromatin'):
        return report.traffic_html(scope or self.scope, ranges, assay)

    def test_chromatin_mito_failure_overrides_passing_technical_metrics(self):
        text = self.render()
        self.assertIn('<b>Library quality</b><em>Fail</em>', text)
        self.assertIn('Mitochondrial fraction: 60.0% (outside range)', text)
        self.assertIn('Usable fragment rate: 28.0% (borderline)', text)
        self.assertIn('<b>Sequencing</b><em>Pass</em>', text)

    def test_total_warning_and_good_sublibraries(self):
        for usable, mito, expected in [(.4763, .3367, 'Warning'), (.77, .002, 'Pass'),
                                       (.73, .012, 'Pass'), (.717, .020, 'Pass'),
                                       (.149, .01, 'Fail'), (.15, .01, 'Warning'),
                                       (.30, .20, 'Pass'), (.30, .40, 'Warning'),
                                       (.30, .4001, 'Fail')]:
            with self.subTest(usable=usable, mito=mito):
                scope = dict(self.scope, usable_fragment_rate=usable, mitochondrial_fraction=mito)
                self.assertIn(f'<b>Library quality</b><em>{expected}</em>', self.render(scope))

    def test_custom_cutoffs_and_disabled_metrics_are_respected(self):
        ranges = copy.deepcopy(report.DEFAULT_REFERENCE_RANGES)
        ranges.update(usable_fragment_rate=None, mitochondrial_fraction=None)
        self.assertIn('<b>Library quality</b><em>Pass</em>', self.render(ranges=ranges))
        ranges['mitochondrial_fraction'] = {'max': .8, 'warning_max': .9}
        self.assertIn('<b>Library quality</b><em>Pass</em>', self.render(ranges=ranges))
        self.assertNotIn('<b>Library quality</b><em>Pass</em>', self.render(ranges={}))

    def test_other_assays_do_not_inherit_chromatin_thresholds(self):
        for assay in ('microc', 'rnaseq', 'clip'):
            with self.subTest(assay=assay):
                self.assertIn('<b>Library quality</b><em>Pass</em>', self.render(assay=assay))

    def test_low_depth_stays_insufficient_and_missing_metrics_are_not_failures(self):
        low = dict(self.scope, raw_reads=10, library_units=10)
        self.assertIn('<em>Insufficient evidence</em>', self.render(low))
        missing = dict(self.scope, usable_fragment_rate=None, mitochondrial_fraction=None)
        self.assertIn('<b>Library quality</b><em>Pass</em>', self.render(missing))


class EcoliNoteTests(unittest.TestCase):
    def test_note_depends_on_screened_candidate_not_detection(self):
        self.assertEqual(report.ecoli_calibration_note_html({}), '')
        self.assertEqual(report.ecoli_calibration_note_html({'candidates': {'yeast': {'label': 'Yeast'}}}), '')
        for status in ('CANDIDATE_SIGNAL', 'TRACE', 'NOT_DETECTED'):
            note = report.ecoli_calibration_note_html({'candidates': {'e_coli': {'label': 'Escherichia coli', 'status': status}}})
            self.assertIn('pA/G-MNase', note)
            self.assertIn('pA/G-Tn5', note)
            self.assertIn('same enzyme batch and amount', note)
            self.assertIn('10.7554/eLife.46314', note)
            self.assertIn('10.1038/s41467-019-09982-5', note)

    def test_note_is_in_rendered_candidate_section(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            (out / 'qc/spike_in').mkdir(parents=True)
            (out / 'qc/spike_in/sample.json').write_text(json.dumps(dict(mode='auto', screened_reads=100,
                candidates={'e_coli': dict(label='Escherichia coli', role='spikein_candidate', reads=0,
                    fraction=0, loci=0, status='NOT_DETECTED')})))
            text = report.extension_report_html(out, 'sample', [], {'assay': 'chromatin'})
            self.assertEqual(text.count("class='note ecoli-calibration-note'"), 1)
            self.assertNotIn('No spike-in was added', text)


if __name__ == '__main__':
    unittest.main()
