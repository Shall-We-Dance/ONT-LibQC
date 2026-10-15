"""Population accounting, endpoint geometry, missing evidence and report compatibility."""
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'workflow/scripts'))
import microc_pairtools as pt
import ont_qc as qc

class PairtoolsTests(unittest.TestCase):
    def test_report_names_actual_bwa_variant(self):
        result=pt.summarize('classic',{}, {},set(),set(),{})
        result['settings']['aligner']='bwa-mem mem -SP'
        with tempfile.TemporaryDirectory() as d:
            page=pt.report_html(result,d,[{'scope_id':'total'}])
        self.assertIn('aligned with BWA-MEM and',page)
        self.assertIn('bwa-mem mem -SP',page)
        self.assertNotIn('BWA-MEM2',page)

    def test_full_denominator_and_inherited_duplicate_marks(self):
        before=dict(a='UU', b='UU', c='RU', d='WW', e='NN')
        after=dict(before,b='DD')
        eligible=set('abcdef')
        assignments={u:dict(unit_id=u,read_id='parent',scope_id='A' if u in 'abcf' else 'B',
                            assignment_status='assigned') for u in 'abcdefg'}
        result=pt.summarize('example.v2',before,after,eligible,{'g'},assignments,{'empty':'No insert'})
        total=result['scopes'][0]
        self.assertEqual(total['parsed_pairs'],5)
        self.assertEqual(total['fastp_excluded_units'],1)
        self.assertEqual(total['short_inserts_excluded'],1)
        self.assertEqual(total['after']['counts'],{'UU':1,'UR/RU':1,'DD':1,'WW':1,'Other':1})
        self.assertEqual(total['after']['fractions'],dict.fromkeys(pt.CATEGORIES,.2))
        self.assertEqual(total['after']['uu_to_ur_ru'],1)
        self.assertNotIn('before',total)
        subset=next(s for s in result['scopes'] if s['scope_id']=='A')
        self.assertEqual(subset['parsed_pairs'],3)
        self.assertEqual(subset['after']['counts']['DD'],1)
        empty=next(s for s in result['scopes'] if s['scope_id']=='empty')
        self.assertTrue(all(x is None for x in empty['after']['fractions'].values()))
        with tempfile.TemporaryDirectory() as d:
            prefix=Path(d)/'example.v2'
            pt.save_results(result,after,eligible,{'g'},assignments,prefix)
            exported=json.loads(Path(str(prefix)+'.json').read_text())
            self.assertTrue(all('before' not in scope for scope in exported['scopes']))
            with Path(str(prefix)+'.summary.tsv').open() as f:
                summaries=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual({row['stage'] for row in summaries},{'after'})
            with gzip.open(str(prefix)+'.units.tsv.gz','rt') as f:
                rows={r['unit_id']:r for r in csv.DictReader(f,delimiter='\t')}
            self.assertEqual(len(rows),7)
            self.assertTrue(all(not key.startswith('before') for row in rows.values() for key in row))
            self.assertEqual(rows['f']['status'],'fastp_excluded')
            self.assertEqual(rows['g']['status'],'short_insert')
            self.assertEqual(rows['b']['after_pair_type'],'DD')
            self.assertEqual(rows['b']['read_id'],'parent')

    def test_population_changes_fail_instead_of_changing_denominator(self):
        for before,after,eligible in [({'a':'UU'},{}, {'a'}),({'a':'UU'},{'a':'UU'},set())]:
            with self.assertRaises(ValueError):
                pt.summarize('x',before,after,eligible,set(),{})
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'duplicated.pairs'
            path.write_text('#columns: readID pair_type\na\tUU\na\tRU\n')
            with self.assertRaises(ValueError):pt.pair_types(path)

    def test_overlapping_ends_and_exclusion_boundary(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);source=d/'inserts.fastq'
            seq='ACGT'*30;qual=''.join(chr(33+i%40) for i in range(120))
            source.write_text('@overlap\n'+seq+'\n+\n'+qual+'\n@short\n'+'A'*49+'\n+\n'+'I'*49+
                              '\n@boundary\n'+'C'*50+'\n+\n'+'J'*50+'\n@long\n'+'A'*160+'C'*160+'\n+\n'+'I'*320+'\n')
            eligible,short=pt.emit_ends(source,d/'r1.gz',d/'r2.gz')
            self.assertEqual(eligible,{'overlap','boundary','long'})
            self.assertEqual(short,{'short'})
            left={u:(s,q) for u,s,q in pt.fastq(d/'r1.gz')}
            right={u:(s,q) for u,s,q in pt.fastq(d/'r2.gz')}
            self.assertEqual(left['overlap'],(seq,qual))
            self.assertEqual(right['overlap'],(seq,qual[::-1]))
            self.assertEqual(right['boundary'][0],'G'*50)
            self.assertEqual(left['long'][0],'A'*150)
            self.assertEqual(right['long'][0],'G'*150)

    def test_no_pairs_is_missing_evidence_and_pies_are_optional(self):
        result=pt.summarize('empty',{}, {},set(),set(),{})
        with tempfile.TemporaryDirectory() as d:
            page=pt.report_html(result,d,[{'scope_id':'total'}])
            self.assertIn('0 parsed pairs',page)
            self.assertNotIn('Before:',page)
            self.assertNotIn('before marking',page)
            self.assertNotIn('Before marking',(Path(d)/'empty.total.pairtools.svg').read_text())
            self.assertIn('<td>0</td><td>—</td>',page)
            self.assertIn('No parsed pairs',(Path(d)/'empty.total.pairtools.svg').read_text())
            with self.assertRaises(ValueError):pt.report_html(result,d,[{'scope_id':'unknown'}])
            assets=dict(model_summary=dict(adapter_offsets_bp=[],modeled_units=0,excluded_units=0),
                        weighted='weight.svg',count='count.svg',profiles='profiles.tsv',units='units.tsv',model='model.json')
            microc=qc.profile_scope_html('Micro-C','raw.svg','insert.svg',None,'k.svg',full_library_assets=assets)
            chromatin=qc.profile_scope_html('Chromatin','raw.svg','insert.svg','regions.svg','k.svg',full_library_assets=assets)
            self.assertNotIn('Genomic regions',microc)
            self.assertNotIn('analysis-plots',microc)
            self.assertIn('Genomic regions',chromatin)

if __name__=='__main__': unittest.main()
