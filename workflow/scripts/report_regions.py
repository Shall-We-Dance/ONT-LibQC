"""Region overlap fractions and exact chromosome/length-matched genomic backgrounds."""
import bisect
import collections
import gzip
import html
import json
from pathlib import Path


def merge(intervals):
    out=[]
    for start,end in sorted(intervals):
        if out and start<=out[-1][1]:out[-1]=(out[-1][0],max(out[-1][1],end))
        else:out.append((start,end))
    return out


def overlaps(intervals,start,end):
    i=bisect.bisect_left(intervals,(end,))
    return i>0 and intervals[i-1][1]>start


def load_te_bed(path,sizes):
    """Operator-prepared TE-only BED; missing annotation is distinct from zero TE."""
    if not path:return None
    intervals=collections.defaultdict(list)
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as handle:
        for number,line in enumerate(handle,1):
            if not line.strip() or line.startswith(('#','track','browser')):continue
            fields=line.split();chrom=fields[0]
            if chrom not in sizes:
                alternate=chrom[3:] if chrom.startswith('chr') else 'chr'+chrom
                if alternate in sizes:chrom=alternate
                else:raise ValueError(f'TE BED line {number}: chromosome absent from reference')
            start,end=map(int,fields[1:3])
            if not 0<=start<end<=sizes[chrom]:raise ValueError(f'TE BED line {number}: invalid coordinates')
            intervals[chrom].append((start,end))
    return {c:merge(rows) for c,rows in intervals.items()}


class OverlapBackground:
    """Exact chance of any overlap for a uniformly placed interval of length L.

    Non-overlapping placements fit wholly inside an annotation-free gap. Sorted
    gap lengths and prefix sums give O(log number_of_gaps) queries, without
    random sampling or repeatedly expanding millions of TE intervals.
    """
    def __init__(self,intervals,size):
        self.size=size;end=0;gaps=[]
        for start,stop in merge((max(0,a),min(size,b)) for a,b in intervals if a<size and b>0):
            if start>end:gaps.append(start-end)
            end=stop
        if end<size:gaps.append(size-end)
        self.gaps=sorted(gaps);self.sums=[0]
        for gap in self.gaps:self.sums.append(self.sums[-1]+gap)

    def probability(self,length):
        if not 1<=length<=self.size:raise ValueError('Alignment span exceeds chromosome bounds')
        i=bisect.bisect_left(self.gaps,length);n=len(self.gaps)-i
        empty=self.sums[-1]-self.sums[i]-(length-1)*n
        return 1-empty/(self.size-length+1)


def region_counts(records,promoters,genes,primary_names,mapq,promoter_window,sizes=None,te=None,te_source=None,background=None):
    primary=set(primary_names)
    eligible=[r for r in records.values() if r['mapq']>=mapq and r['chrom'] in primary]
    keys=('promoter_proximal','gene_body','intergenic','te')
    count=dict.fromkeys(keys,0);lengths=collections.Counter()
    for r in eligible:
        c,s,e=r['chrom'],r['start'],r['end'];lengths[c,e-s]+=1
        if overlaps(promoters.get(c,[]),s,e):count['promoter_proximal']+=1
        elif overlaps(genes.get(c,[]),s,e):count['gene_body']+=1
        else:count['intergenic']+=1
        if te is not None and overlaps(te.get(c,[]),s,e):count['te']+=1
    total=len(eligible);annotated=bool(promoters or genes)
    result=dict(annotation_mapped_reads=total,promoter_window_bp=promoter_window,
        annotation_region_definition='mutually exclusive: promoter proximal takes precedence over gene body; TE overlap is independent',
        region_background_definition='Uniform valid interval starts on the same chromosome, matching each observed alignment span; no mappability correction or IP/input inference.',
        te_annotation_available=te is not None,te_source=te_source)
    for key in keys:
        available=te is not None if key=='te' else annotated
        result[key+'_reads']=count[key] if available else None
        result[key+'_rate']=count[key]/total if available and total else None
        result[key+'_background_rate']=result[key+'_enrichment']=result[key+'_genome_rate']=None
    if not sizes:return result
    if background is None:background=prepare_background(promoters,genes,primary,sizes,te)
    expected=dict.fromkeys(keys,0.0)
    for (chrom,length),n in lengths.items():
        p=background[chrom]['promoter'].probability(length)
        union=background[chrom]['union'].probability(length)
        expected['promoter_proximal']+=n*p;expected['gene_body']+=n*max(0,union-p);expected['intergenic']+=n*(1-union)
        if te is not None:expected['te']+=n*background[chrom]['te'].probability(length)
    genome_total=sum(sizes[c] for c in primary)
    for key in keys:
        available=te is not None if key=='te' else annotated
        if not available:continue
        result[key+'_genome_rate']=sum(background[c]['bases'][key] for c in primary)/genome_total if genome_total else None
        bg=expected[key]/total if total else None
        result[key+'_background_rate']=bg
        rate=result[key+'_rate'];result[key+'_enrichment']=rate/bg if rate is not None and bg else None
    return result


def prepare_background(promoters,genes,primary,sizes,te=None):
    result={}
    for c in primary:
        p=merge(promoters.get(c,[]));union=merge(p+genes.get(c,[]));repeat=(te or {}).get(c,[])
        pb=sum(b-a for a,b in p);ub=sum(b-a for a,b in union)
        result[c]=dict(promoter=OverlapBackground(p,sizes[c]),union=OverlapBackground(union,sizes[c]),
            te=OverlapBackground(repeat,sizes[c]) if te is not None else None,
            bases=dict(promoter_proximal=pb,gene_body=ub-pb,intergenic=sizes[c]-ub,te=sum(b-a for a,b in repeat)))
    return result


def svg_region_bars(metric,title):
    rows=[('Promoter','promoter_proximal'),('Gene body','gene_body'),('Distal','intergenic'),('TE overlap','te')]
    w,h,bx,bw=820,340,145,265;ink='#34495e';muted='#65758c'
    svg=[f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {w} {h}' role='img' aria-label='{html.escape(title,quote=True)}' style='font-family:Arial,Helvetica,sans-serif'>",
        f'<title>{html.escape(title)}</title>',"<metadata id='genomic-region-metrics'>"+html.escape(json.dumps(metric,sort_keys=True))+'</metadata>',
        "<rect width='100%' height='100%' fill='white'/>",f"<text x='22' y='28' font-size='17' fill='#203040'>{html.escape(title)}</text>"]
    for x,label in [(540,'Sample'),(665,'Background'),(798,'Enrichment')]:svg.append(f"<text x='{x}' y='61' text-anchor='end' font-size='12' fill='{muted}'>{label}</text>")
    percent=lambda x:('—' if x is None else f'{x:.1%}' if x==0 or x>=.001 else f'{100*x:.3g}%')
    for i,(label,key) in enumerate(rows):
        y=80+i*47+(8 if key=='te' else 0);rate=metric.get(key+'_rate');bg=metric.get(key+'_background_rate');fold=metric.get(key+'_enrichment');count=metric.get(key+'_reads')
        if key=='te':svg.append("<line x1='22' x2='798' y1='215' y2='215' stroke='#e3e9ef'/>")
        tooltip=f'{label}: sample {percent(rate)}, background {percent(bg)}.'
        if key=='te':tooltip+=' TE overlaps are independent of the three genomic regions.'
        svg.append(f"<g><title>{html.escape(tooltip)}</title><text x='22' y='{y+15}' font-size='14' fill='{ink}'>{label}</text>")
        for j,(value,color) in enumerate([(rate,'#3276b5'),(bg,'#aebccc')]):
            yy=y+j*12;svg.append(f"<rect x='{bx}' y='{yy}' width='{bw}' height='9' rx='2' fill='#eef2f6'/>")
            if value is not None:svg.append(f"<rect x='{bx}' y='{yy}' width='{bw*value:.3f}' height='9' rx='2' fill='{color}'/>")
        label_value=percent(rate)+(f' ({count:,})' if count is not None and rate is not None else '')
        for x,value in [(540,label_value),(665,percent(bg)),(798,f'{fold:,.2f}×' if fold is not None else '—')]:svg.append(f"<text x='{x}' y='{y+15}' text-anchor='end' font-size='13' fill='{ink}'>{value}</text>")
        svg.append('</g>')
    total=metric.get('annotation_mapped_reads',0);window=metric.get('promoter_window_bp',2000)
    footer=f'n={total:,} MAPQ-pass alignments · Promoter: TSS ±{window/1000:g} kb'
    svg.append(f"<text x='22' y='288' font-size='11' fill='{muted}'>{html.escape(footer)}</text>")
    note='Background matches chromosome and fragment length. TE overlaps may also be promoter or gene body.'
    if not metric.get('te_annotation_available'):note='Background matches chromosome and fragment length. TE annotation unavailable for this reference.'
    svg.append(f"<text x='22' y='309' font-size='11' fill='{muted}'>{note}</text></svg>")
    return ''.join(svg)
