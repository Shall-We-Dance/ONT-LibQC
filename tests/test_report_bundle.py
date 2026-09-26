import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from report_bundle import build, validate, initialize_work, cleanup_work, check_locations, cleaned


class ReportBundleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.work, self.report = self.root/'work', self.root/'report'
        self.config = self.root/'config.json'
        self.config.write_text(json.dumps({'samples': {'example': {'type': 'chromatin'}}}))
        initialize_work(self.work, self.report, {})
        for name, text in {
            'summary/index.html': '<html><a href="../reports/example.html">Example</a></html>',
            'reports/example.html': '<html><h1 id="sample">Example</h1><a href="../summary/index.html">All</a><a href="#sample">Sample</a><img src="assets/plot.svg"><a href="../qc/adapters/example.reads.tsv.gz">Data</a></html>',
            'reports/assets/plot.svg': '<svg><defs><path id="glyph"/></defs><use href="#glyph"/></svg>',
            'qc/adapters/example.reads.tsv.gz': 'fixture',
            'qc/adapters/example.summary.json': '{}',
            'qc/metrics/example.json': '{"collision_oe":2,"collision_sensitivity":{"bins":[]},"collision_top1_signal_fraction":0.5}',
            'qc/mycoplasma/example.json': '{}',
            'qc/spike_in/example.json': '{}',
        }.items():
            p = self.work/name; p.parent.mkdir(parents=True, exist_ok=True); p.write_text(text)

    def test_relocation_and_cleanup_preserve_all_links(self):
        (self.work/'bam').mkdir();(self.work/'bam/example.bam').write_text('transient')
        build(self.work, self.report, self.config)
        metrics=json.loads((self.report/'data/metrics/example.json').read_text())
        self.assertEqual(metrics, {'collision_oe':2})
        self.assertTrue(cleanup_work(self.work,self.report));self.assertFalse(self.work.exists())
        moved=self.root/'moved';shutil.move(self.report,moved)
        self.assertEqual(validate(moved)['sample_reports'],1)

    def test_invalid_dependency_preserves_previous_delivery(self):
        build(self.work,self.report,self.config)
        before=(self.report/'manifest.json').read_bytes()
        (self.work/'reports/assets/plot.svg').unlink()
        with self.assertRaises(ValueError):build(self.work,self.report,self.config)
        self.assertEqual((self.report/'manifest.json').read_bytes(),before)

    def test_bam_link_cannot_enter_delivery(self):
        (self.work/'bam').mkdir();(self.work/'bam/example.bam').write_text('alignment')
        (self.work/'reports/example.html').write_text('<a href="../bam/example.bam">BAM</a>')
        with self.assertRaises(ValueError):build(self.work,self.report,self.config)
        self.assertFalse(self.report.exists())

    def test_work_cleanup_refuses_unknown_files(self):
        build(self.work,self.report,self.config)
        (self.work/'user-notes.txt').write_text('keep')
        with self.assertRaises(ValueError):cleanup_work(self.work,self.report)
        self.assertTrue((self.work/'user-notes.txt').is_file())

    def test_partial_run_does_not_delete_work(self):
        build(self.work,self.report,self.config)
        initialize_work(self.work,self.report,{})
        self.assertFalse(cleanup_work(self.work,self.report))
        self.assertTrue(self.work.is_dir())

    def test_changed_report_blocks_cleanup(self):
        build(self.work,self.report,self.config)
        (self.report/'assets/plot.svg').write_text('<svg><image href="missing.svg"/></svg>')
        with self.assertRaises(ValueError):cleanup_work(self.work,self.report)
        self.assertTrue(self.work.is_dir())

    def test_output_paths_cannot_enclose_inputs_or_each_other(self):
        for work,report,config in [
            (self.work,self.work/'report',{}),
            (self.work,self.report,{'samples':{'s':{'fastq':str(self.work/'raw.fastq')}}}),
            (self.work,self.report,{'references':{'ref':{'fasta':str(self.report/'ref.fa')}}}),
        ]:
            with self.assertRaises(ValueError):check_locations(work,report,config)

    def test_unowned_work_is_never_adopted(self):
        unowned=self.root/'unowned';unowned.mkdir();(unowned/'data.txt').write_text('keep')
        with self.assertRaises(ValueError):initialize_work(unowned,self.report,{})

    def test_external_link_and_missing_anchor_are_distinct(self):
        page=self.work/'reports/example.html'
        page.write_text('<a href="https://example.org/methods">Methods</a><h1 id="valid">Example</h1>')
        build(self.work,self.report,self.config)
        page.write_text('<a href="#missing">Invalid</a>')
        with self.assertRaises(ValueError):build(self.work,self.root/'second-report',self.config)


if __name__=='__main__':unittest.main()
