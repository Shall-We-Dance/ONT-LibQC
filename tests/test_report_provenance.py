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
            self.assertNotEqual(first['pipeline_version'],capture_analysis(b,stamp)['pipeline_version'])
            self.assertNotIn(a,str(first))

    def test_saved_analysis_survives_later_report_generation(self):
        old={'pipeline_version':'build.123456789abc','analysis_completed_utc':'2026-01-01T23:00:00-05:00'}
        record=analysis_record({'sample':'example','analysis_provenance':old},{'created_utc':'2026-03-01T00:00:00Z'})
        self.assertEqual(record,old)
        self.assertEqual(analysis_label(record),'build.123456789abc · 2026-01-02 UTC')

    def test_historical_date_does_not_invent_a_build(self):
        record=analysis_record({}, {'created_utc':'2026-09-10T17:52:51+00:00'})
        self.assertEqual(analysis_label(record),'Version unrecorded · 2026-09-10 UTC')
        self.assertEqual(analysis_label(analysis_record({})),'Version unrecorded · Date unrecorded')
