#!/usr/bin/env python3
"""Conservative whole-insert split-alignment evidence for blunt-end Micro-C ligation.
Counts are resolved genomic blocks, not direct biochemical junction measurements.
"""
import argparse,collections,csv,gzip,json,math,re,subprocess
from pathlib import Path
CATEGORIES=['1','2','3','4','5+','unresolved']

def fastq_lengths(path):
    result={}
    with gzip.open(path,'rt') if str(path).endswith('.gz') else open(path) as f:
        while True:
            name=f.readline()
            if not name:break
            seq=f.readline().strip();plus=f.readline();qual=f.readline().strip()
            if not plus.startswith('+') or len(seq)!=len(qual):raise ValueError('Invalid FASTQ')
            key=name[1:].split()[0]
            if key in result:raise ValueError('Duplicate insert ID')
            result[key]=len(seq)
    return result

def parse_segment(fields,query_length,min_bp=40,mapq=20,min_identity=.85):
    flag=int(fields[1]);mq=int(fields[4])
    if flag & (4|256|512|1024) or mq<mapq or mq==255:return None
    ops=[(int(n),op) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',fields[5])]
    full=sum(n for n,o in ops if o in 'MISH=X')
    if full!=query_length:return None
    left=0
    for n,o in ops:
        if o in 'SH':left+=n
        else:break
    aligned_q=sum(n for n,o in ops if o in 'MI=X');aligned=sum(n for n,o in ops if o in 'M=X')
    if aligned<min_bp:return None
    tags={x[:2]:x[5:] for x in fields[11:] if len(x)>5}
    nm=int(tags.get('NM',999999));identity=max(0.,1-nm/max(1,aligned))
    if identity<min_identity:return None
    qstart,qend=left,left+aligned_q
    if flag&16:qstart,qend=query_length-qend,query_length-qstart
    start=int(fields[3])-1;end=start+sum(n for n,o in ops if o in 'MDN=X')
    return dict(qstart=qstart,qend=qend,chrom=fields[2],start=start,end=end,strand='-' if flag&16 else '+',mapq=mq,identity=identity,aligned_bp=aligned,long_gap=any(n>=50 and o in 'IDN' for n,o in ops),cigar=fields[5])

def chain_segments(segments,max_overlap=10):
    """Maximum-covered query chain; alternatives covering the same query are not extra inserts."""
    segments=sorted(segments,key=lambda x:(x['qend'],x['qstart'],-x['mapq']))
    paths=[];scores=[]
    for i,s in enumerate(segments):
        best=s['qend']-s['qstart'];path=[s]
        for j,t in enumerate(segments[:i]):
            if s['qstart']<t['qend']-max_overlap or s['qend']-t['qend']<30:continue
            score=scores[j]+s['qend']-max(s['qstart'],t['qend'])
            if score>best:best=score;path=paths[j]+[s]
        scores.append(best);paths.append(path)
    return paths[max(range(len(scores)),key=lambda i:(scores[i],sum(s['mapq'] for s in paths[i]),-len(paths[i])))] if scores else []

def coalesce_collinear(chain,tolerance=20):
    out=[]
    for source in chain:
        s=dict(source)
        if out:
            prev=out[-1];qgap=s['qstart']-prev['qend']
            rgap=s['start']-prev['end'] if s['strand']=='+' else prev['start']-s['end']
            if s['chrom']==prev['chrom'] and s['strand']==prev['strand'] and abs(rgap-qgap)<=tolerance and abs(rgap)<=30:
                prev.update(qend=s['qend'],start=min(prev['start'],s['start']),end=max(prev['end'],s['end']),mapq=min(prev['mapq'],s['mapq']),identity=min(prev['identity'],s['identity']),long_gap=prev['long_gap'] or s['long_gap'])
                continue
        out.append(s)
    return out

def classify(segments,length,min_coverage=.8,max_gap=30):
    chain=coalesce_collinear(chain_segments(segments))
    covered=0;end=0
    for s in chain:covered+=max(0,s['qend']-max(end,s['qstart']));end=max(end,s['qend'])
    coverage=covered/length if length else 0
    gaps=[chain[0]['qstart'],length-chain[-1]['qend']]+[max(0,b['qstart']-a['qend']) for a,b in zip(chain,chain[1:])] if chain else [length]
    reasons=[]
    if not chain:reasons.append('no_qualifying_nuclear_alignment')
    if coverage<min_coverage:reasons.append('partial_query_coverage')
    if max(gaps)>max_gap:reasons.append('unexplained_query_gap')
    if any(s['long_gap'] for s in chain):reasons.append('large_indel_requires_resolution')
    resolved=not reasons;n=len(chain)
    return dict(category=str(n) if n<5 else '5+',resolved=resolved,segments=n,query_coverage=coverage,reason=';'.join(reasons) or 'resolved_genomic_blocks',chain=chain) | ({} if resolved else {'category':'unresolved'})

def summarize(rows):
    counts=collections.Counter(x['category'] for x in rows);total=len(rows);resolved=total-counts['unresolved'];ligated=sum(counts[c] for c in ['2','3','4','5+'])
    return dict(ligation_screened_units=total,ligation_resolved_units=resolved,ligation_unresolved_units=counts['unresolved'],ligation_resolved_fraction=resolved/total if total else None,
                ligation_two_units=counts['2'],ligation_three_units=counts['3'],ligation_four_units=counts['4'],ligation_five_plus_units=counts['5+'],ligation_single_block_units=counts['1'],
                ligation_multifragment_fraction=ligated/total if total else None,ligation_multifragment_fraction_resolved=ligated/resolved if resolved else None,
                category_counts={c:counts[c] for c in CATEGORIES},raw_molecules_with_units=len({x['read_id'] for x in rows}))

def analyze(args):
    lengths=fastq_lengths(args.fastq);assignment={}
    if args.assignments:
        with gzip.open(args.assignments,'rt') as f:assignment={r['unit_id']:r for r in csv.DictReader(f,delimiter='\t')}
    primary=set(args.primary.split(','));by_id=collections.defaultdict(list)
    with open(args.sam) as f:
        for line in f:
            if line.startswith('@'):continue
            fields=line.rstrip().split('\t');name=fields[0]
            if name not in lengths:raise ValueError('Alignment query not in input FASTQ: '+name)
            if fields[2] not in primary:continue
            seg=parse_segment(fields,lengths[name],args.minimum_segment_bp,args.mapq,args.minimum_identity)
            if seg:by_id[name].append(seg)
    rows=[];segment_rows=[]
    for name,length in lengths.items():
        result=classify(by_id[name],length,args.minimum_query_coverage,args.maximum_query_gap)
        a=assignment.get(name,{});scope=a.get('scope_id','unassigned') if a.get('assignment_status')=='assigned' else 'unassigned'
        row=dict(unit_id=name,read_id=a.get('read_id',name.rsplit('|u',1)[0]),scope_id=scope,insert_bp=length,category=result['category'],resolved=result['resolved'],inferred_blocks=result['segments'],query_coverage=result['query_coverage'],reason=result['reason']);rows.append(row)
        for i,s in enumerate(result['chain'],1):segment_rows.append(dict(unit_id=name,block=i,**s))
    out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    for path,table,fields in [(out.with_suffix('.units.tsv.gz'),rows,['unit_id','read_id','scope_id','insert_bp','category','resolved','inferred_blocks','query_coverage','reason']),(out.with_suffix('.segments.tsv.gz'),segment_rows,['unit_id','block','qstart','qend','chrom','start','end','strand','mapq','identity','aligned_bp','long_gap','cigar'])]:
        with gzip.open(path,'wt') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(table)
    scopes={s:summarize([x for x in rows if x['scope_id']==s]) for s in sorted({x['scope_id'] for x in rows})}
    result=dict(sample=args.sample,method='full_adapter_bounded_insert_'+getattr(args,'aligner','bwa-mem2').upper().replace('-','_')+'_split_alignment',aligner=getattr(args,'aligner','bwa-mem2'),status='EXPLORATORY',unit_of_analysis='adapter_bounded_non_dimer_library_unit',**summarize(rows),scopes=scopes,parameters={k:getattr(args,k) for k in ['minimum_segment_bp','mapq','minimum_identity','minimum_query_coverage','maximum_query_gap']},units_tsv=str(out.with_suffix('.units.tsv.gz')),segments_tsv=str(out.with_suffix('.segments.tsv.gz')),limitations=['Resolved blocks are evidence-supported minimum fragment counts, not exact biochemical insert counts.','Collinear or immediately adjacent genomic fragments may be indistinguishable from one continuous block.','Short, repetitive, low-quality or unaligned internal fragments may be missed.','ONT adapter-linked library units are not counted as Micro-C junctions; units sharing a raw molecule are dependent.','No positive wet-lab ligation standards or no-ligation control have been supplied.'])
    out.write_text(json.dumps(result,indent=2)+'\n');return result

def plots(result,units_path,segments_path,output):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    with gzip.open(units_path,'rt') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    with gzip.open(segments_path,'rt') as f:segs=list(csv.DictReader(f,delimiter='\t'))
    base=Path(output);base.parent.mkdir(parents=True,exist_ok=True)
    for scope in ['total']+list(result['scopes']):
        selected=rows if scope=='total' else [r for r in rows if r['scope_id']==scope];stat=result if scope=='total' else result['scopes'][scope]
        fig,axs=plt.subplots(1,2,figsize=(11,4.8),gridspec_kw={'width_ratios':[1.35,1]});colors=['#8293a4','#3286a1','#479976','#c49a42','#926b99','#c1c7cc'];vals=[stat['category_counts'][c] for c in CATEGORIES]
        bars=axs[0].bar(range(6),vals,color=colors);axs[0].set_xticks(range(6),['1 block','2','3','4','≥5','Unresolved']);axs[0].set_xlabel('Inferred genomic fragments per library unit');axs[0].set_ylabel('Library units');axs[0].set_ylim(0,max(vals+[1])*1.25)
        for bar,v in zip(bars,vals):axs[0].text(bar.get_x()+bar.get_width()/2,v,f'{v:,}\n{v/max(1,len(selected)):.1%}',ha='center',va='bottom',fontsize=9)
        for i,c in enumerate(CATEGORIES):
            lengths=[int(x['insert_bp']) for x in selected if x['category']==c]
            if lengths:
                b=axs[1].boxplot([lengths],positions=[i],widths=.55,showfliers=False,patch_artist=True);b['boxes'][0].set_facecolor(colors[i])
        axs[1].set_xticks(range(6),['1','2','3','4','≥5','Unres.']);axs[1].set_ylabel('Adapter-bounded insert length (bp)');axs[1].set_xlabel('Inferred fragment count');axs[1].set_ylim(bottom=0)
        for ax in axs:ax.spines[['top','right']].set_visible(False)
        fig.suptitle(f"{result['sample']} · {scope}: Micro-C ligation evidence",fontsize=16)
        fig.text(.5,.015,'Counts use all non-dimer reconstructed units. One block does not prove absence of ligation; unresolved units remain separate.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.055,1,.94))
        for ext in ['svg','png','pdf']:fig.savefig(str(base)+f'.{scope}.{ext}',dpi=180,bbox_inches='tight')
        plt.close(fig)
    # One representative per category: linked query blocks in original molecule orientation.
    examples=[]
    for cat in ['2','3','4','5+','1','unresolved']:
        candidates=[r for r in rows if r['category']==cat]
        if candidates:examples.append(sorted(candidates,key=lambda r:(-float(r['query_coverage']),int(r['insert_bp']),r['unit_id']))[len(candidates)//2])
    fig,ax=plt.subplots(figsize=(12,max(3,len(examples)*.9+1.4)))
    for i,row in enumerate(examples):
        length=int(row['insert_bp']);ax.plot([0,length],[i,i],color='#c5cbd0',lw=9,zorder=1)
        for j,s in enumerate(x for x in segs if x['unit_id']==row['unit_id']):
            qs,qe=int(s['qstart']),int(s['qend']);ax.plot([qs,qe],[i,i],lw=10,color=plt.get_cmap('tab10')(j%10),zorder=2)
            ax.text((qs+qe)/2,i+.2,s['chrom']+s['strand'],fontsize=8,ha='center')
    ax.set_yticks(range(len(examples)),[r['category']+' · '+r['unit_id'] for r in examples],fontsize=8);ax.set_xlabel('Query position within one adapter-bounded unit (bp)');ax.set_title(result['sample']+': representative alignment block structures');ax.spines[['top','right']].set_visible(False);fig.tight_layout()
    for ext in ['svg','png','pdf']:fig.savefig(str(base)+f'.examples.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)

def main():
    p=argparse.ArgumentParser();p.add_argument('--aligner',choices=('bwa-mem','bwa-mem2'),default='bwa-mem2');p.add_argument('--sample',required=True);p.add_argument('--fastq',required=True);p.add_argument('--sam',required=True);p.add_argument('--assignments');p.add_argument('--primary',required=True);p.add_argument('--output',required=True);p.add_argument('--plot-prefix');p.add_argument('--minimum-segment-bp',type=int,default=40);p.add_argument('--mapq',type=int,default=20);p.add_argument('--minimum-identity',type=float,default=.85);p.add_argument('--minimum-query-coverage',type=float,default=.8);p.add_argument('--maximum-query-gap',type=int,default=30);a=p.parse_args();result=analyze(a)
    if a.plot_prefix:
        plots(result,result['units_tsv'],result['segments_tsv'],a.plot_prefix)
        from microc_mnase import write_mnase
        write_mnase(a.output,a.plot_prefix)
    print(json.dumps({k:v for k,v in result.items() if k.startswith('ligation_') or k=='sample'},indent=2))
if __name__=='__main__':main()
