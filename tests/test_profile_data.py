"""Numerical sidecars match the plotted matrix, never geometry or invented SEM."""
import gzip
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'workflow/scripts'))
from profile_data import save_profile, read_numeric_profile, read_gene_profile, write_data
from report_profiles import read_region_metrics


class ProfileDataTests(unittest.TestCase):
    def fixture(self, root, values):
        matrix, table, output = (root/name for name in ('matrix.gz','profile.tsv','profile.json'))
        ma = np.ma.masked_invalid(np.array(values, dtype=float))
        means = np.ma.mean(ma, axis=0).filled(np.nan)
        width = len(means)
        with gzip.open(matrix, 'wt') as f:
            f.write('@'+json.dumps({'sample_boundaries':[0,width], 'group_boundaries':[0,len(values)]})+'\n')
            for row in values: f.write('chr1\t0\t100\tgene\t0\t+\t'+'\t'.join(map(str,row))+'\n')
        table.write_text('bin labels\t\t'+'\t'.join(['TSS']+['']*(width-2)+['TES'])+'\nbins\t\t'+'\t'.join(map(str,range(1,width+1)))+'\nsample\tgenes\t'+'\t'.join(map(str,means))+'\n')
        return matrix, table, output, ma

    def test_exact_deeptools_sem_with_masked_bins(self):
        with tempfile.TemporaryDirectory() as tmp:
            matrix, table, output, ma = self.fixture(Path(tmp), [[0,2,float('nan')],[4,float('nan'),float('nan')],[2,8,float('nan')]])
            save_profile(matrix, table, output)
            gene, band = read_numeric_profile(output, read_gene_profile(table))
            expected = np.std(ma, axis=0)/np.sqrt(ma.shape[0])
            for i in (0,1):
                self.assertAlmostEqual(gene[2][i]-band[0][i], expected[i])
                self.assertAlmostEqual(band[1][i]-gene[2][i], expected[i])
            self.assertTrue(math.isnan(band[0][2]))
            self.assertNotIn('NaN',output.read_text())
            self.assertEqual(json.loads(output.read_text())['data']['valid_counts'],[3,2,0])

    def test_mismatched_table_or_multiple_samples_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            matrix, table, output, _ = self.fixture(Path(tmp), [[0,2],[4,8]])
            table.write_text(table.read_text().replace('2.0\t5.0','9.0\t5.0'))
            with self.assertRaisesRegex(ValueError,'means disagree'): save_profile(matrix,table,output)
            with gzip.open(matrix,'rt') as f: content=f.read()
            with gzip.open(matrix,'wt') as f:f.write(content.replace('[0, 2]', '[0, 1, 2]'))
            with self.assertRaisesRegex(ValueError,'one profile sample'):save_profile(matrix,table,output)

    def test_region_json_preferred_and_invalid_schema_does_not_fall_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            svg=Path(tmp)/'regions.svg';svg.write_text('not XML')
            write_data(svg.with_suffix('.json'),'genomic_regions',{'gene_body_rate':.4})
            self.assertEqual(read_region_metrics(svg),{'gene_body_rate':.4})
            svg.with_suffix('.json').write_text('{"schema_version":9}')
            with self.assertRaises(ValueError):read_region_metrics(svg)

    def test_missing_profile_removes_stale_sidecar(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);table=root/'p.tsv';out=root/'p.json'
            table.write_text('No mapped reads\n');out.write_text('{}')
            save_profile(root/'absent.gz',table,out)
            self.assertFalse(out.exists())
