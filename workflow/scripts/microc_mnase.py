"""Length evidence conditional on resolved multi-block Micro-C library units."""
import csv,gzip,json,statistics
from pathlib import Path
from qc_plot_metrics import svg_size_trace

BANDS=[('short_fraction',0,100),('mono_range_fraction',100,200),('intermediate_fraction',200,300),('long_fraction',300,float('inf'))]

def summarize(rows):
    sizes=[int(r['reference_bp']) for r in rows]; n=len(sizes)
    # These are size ranges, not experimentally calibrated digestion calls.
    out={'mnase_resolved_fragments':n,'mnase_contributing_units':len({r['unit_id'] for r in rows}),
         'mnase_contributing_molecules':len({r['read_id'] for r in rows}),
         'mnase_median_bp':statistics.median(sizes) if n else None}
    for key,lo,hi in BANDS:out['mnase_'+key]=sum(lo<=v<hi for v in sizes)/n if n else None
    return out

def extract(units,segments):
    eligible={r['unit_id']:r for r in units if str(r['resolved']).lower()=='true' and int(r['inferred_blocks'])>=2}
    rows=[]
    for s in segments:
        if s['unit_id'] not in eligible:continue
        u=eligible[s['unit_id']];reference_bp=int(s['end'])-int(s['start']);query_bp=int(s['qend'])-int(s['qstart'])
        if min(reference_bp,query_bp)<=0:raise ValueError('Nonpositive resolved segment length')
        rows.append(dict(unit_id=s['unit_id'],read_id=u['read_id'],scope_id=u['scope_id'],block=s['block'],
                         inferred_blocks=int(u['inferred_blocks']),reference_bp=reference_bp,query_bp=query_bp,
                         chrom=s['chrom'],start=s['start'],end=s['end'],strand=s['strand']))
    assert len(rows)==sum(int(u['inferred_blocks']) for u in eligible.values()),'Missing resolved segment'
    return rows

def write_mnase(ligation_path,plot_prefix):
    path=Path(ligation_path);lig=json.loads(path.read_text());unitpath=path.with_suffix('.units.tsv.gz');segpath=path.with_suffix('.segments.tsv.gz')
    with gzip.open(unitpath,'rt') as f:units=list(csv.DictReader(f,delimiter='\t'))
    with gzip.open(segpath,'rt') as f:segments=list(csv.DictReader(f,delimiter='\t'))
    rows=extract(units,segments);table=path.with_suffix('.mnase_fragments.tsv.gz')
    fields=['unit_id','read_id','scope_id','block','inferred_blocks','reference_bp','query_bp','chrom','start','end','strand']
    with gzip.open(table,'wt') as f:w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
    result=dict(sample=lig['sample'],method='reference_span_of_resolved_blocks_in_multi_block_units',
                denominator='inferred genomic blocks from resolved units with at least two blocks',
                ranges_bp={'short':'<100','mono_range':'100 <= length < 200','intermediate':'200 <= length < 300','long':'>=300'},
                fragment_table=str(table),**summarize(rows),scopes={},
                input_paths=[str(unitpath),str(segpath)],
                limitations='Conditional on successful ligation, library selection, sequencing and alignment resolution. Short/repetitive fragments may be lost. Adjacent collinear fragments can merge; reference span estimates protected genomic length, not exact physical pre-ligation length. No calibrated digestion efficiency or nucleosome identity is inferred.')
    for scope in ['total']+list(lig['scopes']):
        selected=rows if scope=='total' else [r for r in rows if r['scope_id']==scope]
        if scope!='total':result['scopes'][scope]=summarize(selected)
        title=f"{lig['sample']} · {scope}: inferred pre-ligation fragments"
        note='Resolved multi-block units only; reference spans. Conditional length evidence, not a direct digestion trace.'
        Path(str(plot_prefix)+f'.{scope}.mnase.svg').write_text(svg_size_trace([r['reference_bp'] for r in selected],title,1000,.045,footnote=note))
        single=[int(r['insert_bp']) for r in units if r['category']=='1' and (scope=='total' or r['scope_id']==scope)]
        Path(str(plot_prefix)+f'.{scope}.single_block.svg').write_text(svg_size_trace(single,f"{lig['sample']} · {scope}: single-block library units",1000,.045,footnote='Whole insert lengths for resolved one-block units; one block does not establish unligated DNA.'))
    path.with_suffix('.mnase.json').write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--ligation',required=True);p.add_argument('--plot-prefix',required=True);a=p.parse_args();r=write_mnase(a.ligation,a.plot_prefix);print(r['sample'],r['mnase_resolved_fragments'],r['mnase_median_bp'])
