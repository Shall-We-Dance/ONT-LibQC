import datetime
from pathlib import Path
import sys,tempfile,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from report_provenance import capture_analysis, analysis_record, analysis_label

class AnalysisProvenanceTests(unittest.TestCase):
    def test_build_is_reproducible_across_locations_but_changes_with_code(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            for root in [Path(a),Path(b)]:
                (root/'workflow').mkdir();(root/'Snakefile').write_text('rule all: pass\n')
                (root/'workflow/a.py').write_text('x=1\n')
            stamp=datetime.datetime(2026,1,2,tzinfo=datetime.timezone.utc)
            first=capture_analysis(a,stamp);second=capture_analysis(b,stamp)
            self.assertEqual(first,second)
            (Path(b)/'workflow/a.py').write_text('x=2\n')
            self.assertEqual(first['pipeline_version'],'v1.0.0')
            self.assertEqual(capture_analysis(b,stamp)['pipeline_version'],'v1.0.0')
            self.assertNotEqual(first['pipeline_sha256'],capture_analysis(b,stamp)['pipeline_sha256'])
            self.assertNotIn(a,str(first))

    def test_saved_analysis_survives_later_report_generation(self):
        old={'pipeline_version':'v1.0.0','analysis_completed_utc':'2026-01-01T23:00:00-05:00'}
        record=analysis_record({'sample':'example','analysis_provenance':old},{'created_utc':'2026-03-01T00:00:00Z'})
        self.assertEqual(record,old)
        self.assertEqual(analysis_label(record),'v1.0.0 · 2026-01-02 UTC')

    def test_historical_date_does_not_invent_a_build(self):
        record=analysis_record({}, {'created_utc':'2026-09-10T17:52:51+00:00'})
        self.assertEqual(analysis_label(record),'v1.0.0 · 2026-09-10 UTC')
        self.assertEqual(analysis_label(analysis_record({})),'v1.0.0 · Date unrecorded')

    def test_legacy_build_is_retained_beneath_release_label(self):
        record=analysis_record({'analysis_provenance':{'pipeline_version':'build.123456789abc','pipeline_sha256':'a'*64,'analysis_completed_utc':'2026-11-13T12:00:00Z'}})
        self.assertEqual(record['pipeline_version'],'v1.0.0')
        self.assertEqual(record['pipeline_build'],'build.123456789abc')
        self.assertEqual(record['pipeline_sha256'],'a'*64)
        self.assertEqual(analysis_record({'analysis_provenance':{'pipeline_version':'v2.0.0'}})['pipeline_version'],'v2.0.0')
