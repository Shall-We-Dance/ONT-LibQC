"""Offline custom BED profiles: exact mean coverage around reference points or in scaled bins, with CPM."""
import collections,html,json,math,re
from pathlib import Path
import numpy as np
from custom_regions_ui import render_section

PARAMS=dict(mode='auto',upstream_bp=5000,downstream_bp=5000,flank_bins=50,body_bins=100,reference_bins=100)

def resolve_parameters(regions,params=None):
    resolved={**PARAMS,**(params or {})}
    requested=resolved['mode']
    if requested not in ('auto','reference_point','region'):raise ValueError('Custom region mode must be auto, reference_point or region')
    resolved['requested_mode']=resolved.get('requested_mode',requested) if requested!='auto' else 'auto'
    resolved['mode']=('reference_point' if all(r['end']-r['start']==1 for r in regions) else 'region') if requested=='auto' else requested
    return resolved

def mode_description(params):
    return 'Reference point (interval midpoint)' if params['mode']=='reference_point' else 'Region (scaled body)'


def parse_bed(text,sizes,primary):
    rows=[];seen=set();skipped=duplicates=0;unknown=[]
    for number,line in enumerate(text.splitlines(),1):
        if not line.strip() or line.startswith(('#','track','browser')):continue
        v=line.split()
        if len(v)<3:raise ValueError(f'BED line {number}: at least three columns required')
        try:start,end=int(v[1]),int(v[2])
        except ValueError:raise ValueError(f'BED line {number}: coordinates must be integers')
        if start<0 or end<=start:raise ValueError(f'BED line {number}: invalid interval')
        if v[0] not in sizes:
            unknown.append(dict(line=number,chrom=v[0],start=start,end=end,reason='chromosome_absent_from_reference'));continue
        if start<0 or end<=start or end>sizes[v[0]]:raise ValueError(f'BED line {number}: invalid or out-of-reference interval')
        strand=v[5] if len(v)>=6 else '.'
        if strand not in ['+','-','.']:raise ValueError(f'BED line {number}: strand must be +, - or .')
        if v[0] not in primary:skipped+=1;continue
        key=(v[0],start,end,strand)
        if key in seen:duplicates+=1;continue
        seen.add(key);rows.append(dict(chrom=v[0],start=start,end=end,strand=strand,name=v[3] if len(v)>3 else f'region_{number}'))
    if not rows:raise ValueError('No regions on configured primary chromosomes')
    return rows,dict(retained_regions=len(rows),skipped_nonprimary=skipped,duplicate_intervals_removed=duplicates,excluded_unknown_contigs=unknown)

def aligned_blocks(record):
    if not record.get('cigar'):return [(record['start'],record['end'])]
    result=[];pos=record['start']
    for length,op in re.findall(r'(\d+)([MIDNSHP=X])',record['cigar']):
        length=int(length)
        if op in 'M=X':result.append((pos,pos+length))
        if op in 'MDN=X':pos+=length
    return result

def track(records):
    events=collections.defaultdict(int)
    for r in records:
        for start,end in aligned_blocks(r):events[start]+=1;events[end]-=1
    x=np.array(sorted(events),float)
    if not len(x):return x,x,x
    y=np.cumsum([events[v] for v in x]);area=np.r_[0,np.cumsum(np.diff(x)*y[:-1])]
    return x,y,area

def integral(t,v):
    x,y,area=t
    if not len(x):return np.zeros_like(v)
    i=np.searchsorted(x,v,side='right')-1;j=np.maximum(i,0)
    return np.where(i<0,0,area[j]+(v-x[j])*y[j])

def region_edges(r,params=PARAMS):
    if params.get('mode')=='reference_point':
        center=(r['start']+r['end'])/2
        u,d=params['upstream_bp'],params['downstream_bp']
        left,right=(d,u) if r['strand']=='-' else (u,d)
        return np.linspace(center-left,center+right,params['reference_bins']+1)
    u,d=params['upstream_bp'],params['downstream_bp'];a,b=params['flank_bins'],params['body_bins'];s,e=r['start'],r['end']
    left,right=(d,u) if r['strand']=='-' else (u,d)
    return np.r_[np.linspace(s-left,s,a+1)[:-1],np.linspace(s,e,b+1)[:-1],np.linspace(e,e+right,a+1)]

def compute_profile(records,regions,sizes,params=PARAMS):
    params=resolve_parameters(regions,params)
    by=collections.defaultdict(list)
    for r in records:by[r['chrom']].append(r)
    tracks={c:track(rs) for c,rs in by.items()};matrix=[]
    for r in regions:
        edges=np.clip(region_edges(r,params),0,sizes[r['chrom']]);widths=np.diff(edges)
        vals=np.diff(integral(tracks.get(r['chrom'],(np.array([]),)*3),edges))
        vals=np.divide(vals,widths,out=np.full(len(widths),np.nan),where=widths>0)
        if r['strand']=='-':vals=vals[::-1]
        matrix.append(vals)
    matrix=np.asarray(matrix);counts=np.sum(np.isfinite(matrix),axis=0)
    if not len(records):return dict(mean=[None]*matrix.shape[1],sem=[None]*matrix.shape[1],regions_per_bin=counts.tolist(),eligible_alignments=0)
    matrix*=1e6/len(records);sums=np.nansum(matrix,axis=0);mean=np.divide(sums,counts,out=np.full(len(counts),np.nan),where=counts>0)
    ss=np.nansum((matrix-mean)**2,axis=0);sem=np.sqrt(np.divide(ss,counts*(counts-1),out=np.full(len(counts),np.nan),where=counts>1))
    clean=lambda v:[float(x) if math.isfinite(x) else None for x in v]
    return dict(mean=clean(mean),sem=clean(sem),regions_per_bin=counts.tolist(),eligible_alignments=len(records))

def static_plot(profiles,label,path,params=PARAMS):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(9,5));colors=['#176e87','#c47a35','#7970a8','#4f9477']
    for i,(scope,p) in enumerate(profiles.items()):
        y=np.array([np.nan if v is None else v for v in p['mean']]);se=np.array([np.nan if v is None else v for v in p['sem']]);x=np.arange(len(y))+.5
        if np.isfinite(y).any():ax.plot(x,y,label=scope,color=colors[i%len(colors)],lw=2);ax.fill_between(x,y-se,y+se,color=colors[i%len(colors)],alpha=.13)
    point=params.get('mode')=='reference_point'
    n=params['reference_bins'] if point else 2*params['flank_bins']+params['body_bins']
    u,d=params['upstream_bp'],params['downstream_bp']
    ax.set(xlim=(0,n),ylim=(0,None),ylabel='Mean coverage (CPM)',title=label)
    if point:
        center=n*u/(u+d)
        ax.set_xticks([0,center,n],[f'−{u/1000:g} kb','Reference point (0)',f'+{d/1000:g} kb']);ax.axvline(center,color='#7793a6',ls='--',lw=1)
    else:
        a,b=params['flank_bins'],params['body_bins']
        ax.set_xticks([0,a,a+b,n],[f'−{u/1000:g} kb','Region start','Region end',f'+{d/1000:g} kb']);ax.axvspan(a,a+b,color='#7793a6',alpha=.06)
    ax.spines[['top','right']].set_visible(False)
    if any(p['eligible_alignments'] for p in profiles.values()):ax.legend(frameon=False)
    else:ax.text(.5,.5,'No eligible nuclear alignments',ha='center',transform=ax.transAxes)
    fig.text(.5,.015,'Equal region weights; '+('distance from interval midpoint' if point else 'scaled region bodies')+'; bands = SEM across regions.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.04,1,1))
    for ext in ['svg','pdf','png']:fig.savefig(Path(path).with_suffix('.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)

def report_html(out,sample,reference,sizes,primary,scope_records,scopes,sample_cfg,config_root,mapq):
    eligible={s:[r for r in rs.values() if r['chrom'] in primary and mapq<=r['mapq']<255] for s,rs in scope_records.items()}
    payload=dict(sample=sample,reference=reference,sizes=sizes,primary=list(primary),parameters=PARAMS,
                 scopes={s:{'name':next(x.get('scope_name',s) for x in scopes if x['scope_id']==s),'records':[[r['chrom'],r['start'],r['end'],aligned_blocks(r)] for r in rs]} for s,rs in eligible.items()})
    initial='';configured=sample_cfg.get('custom_regions');requested_mode='auto';bed_name=None
    if configured:
        if isinstance(configured,str):configured=dict(bed=configured,assembly=reference)
        if configured.get('assembly')!=reference:raise ValueError(f'{sample}: custom BED assembly mismatch')
        source=Path(configured['bed']);source=source if source.is_absolute() else Path(config_root)/source
        bed_name=configured.get('bed_name') or (source.name if str(source)!='/job/regions.bed' else None)
        bedtext=source.read_text();regions,info=parse_bed(bedtext,sizes,primary)
        requested_mode=configured.get('mode','auto')
        params=resolve_parameters(regions,dict(mode=requested_mode))
        profiles={payload['scopes'][s]['name']:compute_profile(rs,regions,sizes,params) for s,rs in eligible.items()}
        folder=out/'qc/custom_regions';folder.mkdir(exist_ok=True)
        result=dict(sample=sample,assembly=reference,label=configured.get('label',source.stem),bed=str(source.resolve()),bed_name=bed_name,region_origin=configured.get('region_origin',{}),parameters=params,**info,profiles=profiles)
        (folder/f'{sample}.json').write_text(json.dumps(result,indent=2)+'\n')
        bedasset=out/'reports/assets'/f'{sample}.custom_regions.bed';bedasset.write_text(bedtext)
        plot=out/'reports/assets'/f'{sample}.custom_regions.svg';static_plot(profiles,sample+' · '+result['label'],plot,params)
        with (folder/f'{sample}.profile.tsv').open('w') as f:
            f.write('scope\tbin\tmean_cpm\tsem_across_regions\tregions_in_bin\n')
            for scope,p in profiles.items():
                for i,(mean,sem,n) in enumerate(zip(p['mean'],p['sem'],p['regions_per_bin'])):f.write(f'{scope}\t{i}\t{mean if mean is not None else "NA"}\t{sem if sem is not None else "NA"}\t{n}\n')
        exclusion=f" · {len(info['excluded_unknown_contigs'])} regions excluded on absent contigs: "+', '.join(sorted({r['chrom'] for r in info['excluded_unknown_contigs']})) if info['excluded_unknown_contigs'] else ''
        if info['skipped_nonprimary']:exclusion+=f" · {info['skipped_nonprimary']} non-primary regions excluded"
        note=result['region_origin'].get('interpretation','')
        tip=f" <span class='metric-tip info-icon' tabindex='0' data-tip='{html.escape(note,quote=True)}'>?</span>" if note else ''
        initial=f"<p class='muted'>{html.escape(result['label'])} · {info['retained_regions']:,} regions · {mode_description(params)} · {reference}{html.escape(exclusion)}{tip}</p><img class='custom-static' src='assets/{plot.name}' alt='Configured custom region metaplot'><p><a href='assets/{bedasset.name}'>BED</a> · <a href='../qc/custom_regions/{sample}.profile.tsv'>Profile data</a> · <a href='assets/{plot.with_suffix('.pdf').name}'>PDF</a></p>"
    return render_section(payload,initial,requested_mode,bed_name)
