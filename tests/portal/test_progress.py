import json
from pathlib import Path
import tempfile
import unittest
from ont_portal.progress import analysis_stage

class ProgressTests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        self.work=Path(tmp.name)
        self.metrics=self.work/'.work/qc/metrics';self.metrics.mkdir(parents=True)
        self.assets=self.work/'.work/reports/assets';self.assets.mkdir(parents=True)

    def test_labels_follow_latest_recognized_rule_without_exposing_logs(self):
        self.assertEqual(analysis_stage(self.work,''),'Validating input reads')
        log='private.fastq secret@example.org /private/path\nlocalrule prepare_library_units:\n    input: private.fastq\nrule align:\n'
        self.assertEqual(analysis_stage(self.work,log),'Mapping inserts to the reference genome')
        self.assertEqual(analysis_stage(self.work,log+'rule unknown_rule:\n'), 'Mapping inserts to the reference genome')
        self.assertEqual(analysis_stage(self.work,log+'localrule report_bundle:\n'), 'Packaging report files')

    def test_report_counts_completed_sublibraries_then_chromosome_plots(self):
        (self.metrics/'hidden.json').write_text(json.dumps({'sample':'hidden','assay':'chromatin','sub_libraries':[{'scope_id':'lib01'},{'scope_id':'lib02'}]}))
        def stage():return analysis_stage(self.work,'localrule reports:\n')
        self.assertEqual(stage(),'Plotting sublibrary gene coverage (0 of 2 complete)')
        (self.assets/'hidden.lib01.TSS_TES.svg').write_text('<svg/>')
        (self.assets/'hidden.lib02.TSS_TES.svg').touch()
        self.assertEqual(stage(),'Plotting sublibrary gene coverage (1 of 2 complete)')
        (self.assets/'hidden.lib02.TSS_TES.svg').write_text('<svg/>')
        self.assertEqual(stage(),'Plotting chromosome coverage (0 of 3 complete)')
        for scope in ('total','lib01','lib02'):(self.assets/f'hidden.{scope}.karyotype.svg').write_text('<svg/>')
        self.assertEqual(stage(),'Preparing plots and the report')
        self.assertNotIn('hidden',stage())

    def test_partial_metrics_and_unsafe_names_do_not_break_progress(self):
        path=self.metrics/'partial.json';path.write_text('{')
        self.assertEqual(analysis_stage(self.work,'rule reports:\n'),'Preparing plots and the report')
        path.write_text(json.dumps({'sample':'../../private','sub_libraries':[]}))
        self.assertEqual(analysis_stage(self.work,'rule reports:\n'),'Preparing plots and the report')

if __name__=='__main__':unittest.main()
