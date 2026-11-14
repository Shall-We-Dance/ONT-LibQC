"""Chromosome-prefix compatibility across upload, pipeline and offline browser paths."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'workflow/scripts'), str(ROOT / 'portal')]
from custom_regions import parse_bed, compute_profile
from ont_portal.run_job import InputError, normalize_bed


class CustomChromosomeTests(unittest.TestCase):
    def browser(self, bed, sizes, primary):
        if not shutil.which('node'):
            self.skipTest('Node is required for browser-parser parity')
        code = '''const c=require(process.argv[1]),p=JSON.parse(require('fs').readFileSync(0,'utf8'));
try {console.log(JSON.stringify(c.parseBed(p.bed,p.sizes,p.primary)))}
catch(e){console.log(JSON.stringify({error:e.message}))}'''
        return json.loads(subprocess.check_output(['node', '-e', code, str(ROOT/'workflow/scripts/custom_regions.js')],
            input=json.dumps(dict(bed=bed, sizes=sizes, primary=primary)), text=True))

    def compare(self, bed, sizes, primary):
        rows, info = parse_bed(bed, sizes, primary)
        js = self.browser(bed, sizes, primary)
        self.assertEqual([{k:v for k,v in row.items() if k!='name'} for row in rows], js['rows'])
        for py, key in [('renamed_rows','renamed_rows'), ('renamed_chromosomes','renamed_chromosomes'),
                        ('duplicate_intervals_removed','duplicates'), ('skipped_nonprimary','skipped'),
                        ('excluded_unknown_contigs','unknown')]:
            self.assertEqual(info[py], js[key])
        return rows, info

    def test_narrowpeak_prefix_and_reverse_direction_preserve_intervals(self):
        bed='1\t10\t25\tpeak_a\t10\t.\t2.3\t4\t5\t7\nX\t30\t42\tpeak_b\t20\t-\n'
        for sizes, expected in [({'chr1':100,'chrX':100},['chr1','chrX']), ({'1':100,'X':100},['1','X'])]:
            text=bed if 'chr1' in sizes else bed.replace('1\t10','chr1\t10').replace('X\t30','chrX\t30')
            rows,info=self.compare(text,sizes,list(sizes))
            self.assertEqual([r['chrom'] for r in rows],expected)
            self.assertEqual([(r['start'],r['end'],r['strand']) for r in rows],[(10,25,'.'),(30,42,'-')])
            self.assertEqual(info['renamed_rows'],2)

    def test_normalization_precedes_deduplication_and_profile(self):
        sizes={'chr1':100};bed='1 10 20 a 0 +\nchr1 10 20 b 0 +\n1 30 40 c 0 -\n'
        rows,info=self.compare(bed,sizes,list(sizes))
        canonical,_=parse_bed('chr1 10 20 a 0 +\nchr1 30 40 c 0 -\n',sizes,list(sizes))
        self.assertEqual(info['duplicate_intervals_removed'],1)
        records=[{'chrom':'chr1','start':12,'end':18},{'chrom':'chr1','start':28,'end':45}]
        self.assertEqual(compute_profile(records,rows,sizes),compute_profile(records,canonical,sizes))

    def test_exact_matches_remain_distinct(self):
        rows,info=self.compare('1 1 2\nchr1 1 2\n',{'1':10,'chr1':20},['1','chr1'])
        self.assertEqual(len(rows),2)
        self.assertEqual(info['renamed_rows'],0)

    def test_unknown_and_nonprimary_are_still_excluded(self):
        rows,info=self.compare('1 1 2\nUn 2 3\nMT 1 2\nmissing 1 2\n',{'chr1':10,'chrUn':10,'chrM':20},['chr1'])
        self.assertEqual(len(rows),1)
        self.assertEqual(info['skipped_nonprimary'],1)
        self.assertEqual([r['chrom'] for r in info['excluded_unknown_contigs']],['MT','missing'])

    def test_coordinate_bounds_checked_after_alias_resolution(self):
        for bed in ('1 1 101','chr1 1 101','1 4 4','1 -1 4'):
            with self.subTest(bed=bed):
                with self.assertRaises(ValueError):parse_bed(bed,{'chr1':100},['chr1'])
                self.assertIn('error',self.browser(bed,{'chr1':100},['chr1']))

    def test_upload_normalizer_uses_reference_names_and_bounds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'input';dest=root/'output';sizes=root/'sizes'
            for reference,uploaded in [('chr1','1'),('1','chr1')]:
                sizes.write_text(reference+'\t100\n')
                source.write_text(uploaded+'\t4\t9\tpeak\t5\t-\t8\t9\t10\t3\n')
                normalize_bed(source,dest,sizes)
                self.assertEqual(dest.read_text(),reference+'\t4\t9\tregion_1\t0\t-\n')
                for bad in (uploaded+' 4 101','missing 4 9'):
                    source.write_text(bad)
                    with self.assertRaises(InputError):normalize_bed(source,dest,sizes)
            sizes.write_text('1\t10\nchr1\t100\n');source.write_text('1 4 50')
            with self.assertRaises(InputError):normalize_bed(source,dest,sizes)

if __name__=='__main__':unittest.main()
