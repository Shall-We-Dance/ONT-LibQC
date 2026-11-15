"""Guard module direction and report orchestration as the CLI remains compatible."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
SCRIPTS=Path(__file__).resolve().parents[1]/'workflow/scripts'
sys.path.insert(0,str(SCRIPTS))
import report_render


class MaintenanceBoundaries(unittest.TestCase):
    def test_sequence_consumers_do_not_import_report_stack(self):
        script="""import sys
sys.path.insert(0, sys.argv[1])
import chromatin_alignment, rna_motifs, library_structure, library_sequences
assert library_structure.P5 == library_sequences.P5
assert library_structure.P7_FORWARD == library_sequences.rc(library_sequences.P7)
assert 'ont_qc' not in sys.modules
assert 'report_metrics' not in sys.modules
assert 'report_render' not in sys.modules
"""
        subprocess.run([sys.executable,'-c',script,str(SCRIPTS)],check=True)

    def test_cli_orchestrates_samples_and_preserves_prior_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=root/'config.json';out=root/'work'
            config.write_text(json.dumps({'samples':{'first':{},'second':{}}}))
            (out/'summary').mkdir(parents=True)
            (out/'summary/provenance.json').write_text('{"created_utc":"2020-01-01"}')
            with patch.object(report_render,'render_sample_report',side_effect=[{'sample':'first'},{'sample':'second'}]) as sample,patch.object(report_render,'write_report_summary') as summary:
                report_render.cmd_report(SimpleNamespace(config=config,out=out))
            self.assertEqual([c.args[0] for c in sample.call_args_list],['first','second'])
            self.assertIs(sample.call_args_list[0].args[5],sample.call_args_list[1].args[5])
            self.assertEqual(summary.call_args.args[0],[{'sample':'first'},{'sample':'second'}])
            self.assertEqual(summary.call_args.args[3],{'created_utc':'2020-01-01'})

    def test_live_renderer_has_no_similarity_import_or_call(self):
        import ast
        for name in ('report_render.py','report_profiles.py','rna_report.py'):
            tree=ast.parse((SCRIPTS/name).read_text())
            for n in ast.walk(tree):
                if isinstance(n,ast.ImportFrom):self.assertNotIn(n.module,('report_similarity','rna_cohort_prototype'))
                if isinstance(n,ast.Call) and isinstance(n.func,ast.Name):self.assertNotIn(n.func.id,('cohort_analysis','render_cohort_qc'))
