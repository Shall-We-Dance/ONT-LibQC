"""Region overlap fractions and fixed genomic base-pair backgrounds."""
import bisect
import collections
import gzip
import html
import json
import re
from report_sublibrary import TOTAL_COLOR
from report_svg import svg_start, text


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


BACKGROUND_DEFINITION = ('Fraction of primary nuclear genome bases in each annotation category; '
                         'shared across libraries using the same reference, annotation and promoter window.')


def with_fixed_background(metric):
    result = dict(metric, region_background_definition=BACKGROUND_DEFINITION,
                  region_background_method='genome-bases-v1')
    for key in ('promoter_proximal', 'gene_body', 'intergenic', 'te'):
        bg = result.get(key + '_genome_rate')
        rate = result.get(key + '_rate')
        result[key + '_background_rate'] = bg
        result[key + '_enrichment'] = rate / bg if rate is not None and bg else None
    return result


def region_counts(records,promoters,genes,primary_names,mapq,promoter_window,sizes=None,te=None,te_source=None,background=None):
    primary=set(primary_names)
    eligible=[r for r in records.values() if r['mapq']>=mapq and r['chrom'] in primary]
    keys=('promoter_proximal','gene_body','intergenic','te')
    count=dict.fromkeys(keys,0)
    for r in eligible:
        c,s,e=r['chrom'],r['start'],r['end']
        if overlaps(promoters.get(c,[]),s,e):count['promoter_proximal']+=1
        elif overlaps(genes.get(c,[]),s,e):count['gene_body']+=1
        else:count['intergenic']+=1
        if te is not None and overlaps(te.get(c,[]),s,e):count['te']+=1
    total=len(eligible);annotated=bool(promoters or genes)
    result=dict(annotation_mapped_reads=total,promoter_window_bp=promoter_window,
        annotation_region_definition='mutually exclusive: promoter proximal takes precedence over gene body; TE overlap is independent',
        region_background_definition=BACKGROUND_DEFINITION,
        te_annotation_available=te is not None,te_source=te_source)
    for key in keys:
        available=te is not None if key=='te' else annotated
        result[key+'_reads']=count[key] if available else None
        result[key+'_rate']=count[key]/total if available and total else None
        result[key+'_background_rate']=result[key+'_enrichment']=result[key+'_genome_rate']=None
    if not sizes:return with_fixed_background(result)
    if background is None:background=prepare_background(promoters,genes,primary,sizes,te)
    genome_total=sum(sizes[c] for c in primary)
    for key in keys:
        available=te is not None if key=='te' else annotated
        if not available:continue
        result[key+'_genome_rate']=sum(background[c]['bases'][key] for c in primary)/genome_total if genome_total else None
    return with_fixed_background(result)


def prepare_background(promoters,genes,primary,sizes,te=None):
    result={}
    for c in primary:
        clipped=lambda rows: merge((max(0,a),min(sizes[c],b)) for a,b in rows if a<sizes[c] and b>0)
        p=clipped(promoters.get(c,[]));union=clipped(p+genes.get(c,[]));repeat=clipped((te or {}).get(c,[]))
        pb=sum(b-a for a,b in p);ub=sum(b-a for a,b in union)
        result[c]=dict(bases=dict(promoter_proximal=pb,gene_body=ub-pb,intergenic=sizes[c]-ub,te=sum(b-a for a,b in repeat)))
    return result


REGION_ROWS = (('Promoter', 'promoter_proximal'), ('Gene body', 'gene_body'), ('Distal', 'intergenic'))
REGION_LEGEND = 'Bar: sample proportion · Tick: genome background'


def region_percent(value):
    if value is None:
        return '—'
    return f'{value:.1%}' if value == 0 or value >= .001 else f'{100 * value:.3g}%'


def region_track(rate, background, color, x, y, width):
    parts = [f"<rect x='{x}' y='{y-7}' width='{width}' height='11' rx='2' fill='#eef2f6'/>"]
    if rate is not None:
        parts.append(f"<rect data-value='{rate}' x='{x}' y='{y-7}' width='{width*rate:.3f}' height='11' rx='2' fill='{color}'/>")
    if background is not None:
        tick = x + width * background
        parts.append(f"<line x1='{tick:.3f}' x2='{tick:.3f}' y1='{y-10}' y2='{y+7}' stroke='#536475' stroke-width='1.3'/>")
    return ''.join(parts)


def refresh_region_color(path, color):
    if not path.is_file():
        return
    svg = path.read_text()
    # Only the observed proportions have a data-value. Older category plots
    # use colors for the region classes and have no library identity to update.
    updated = re.sub(r"""(<rect\b[^>]*\bdata-value=['"][^'"]+['"][^>]*\bfill=['"])#[0-9a-fA-F]{6}""",
                     lambda m: m[1] + color, svg)
    if updated != svg:
        path.write_text(updated)


def svg_region_bars(metric, title, color=TOTAL_COLOR):
    svg = svg_start(820, 235, title)
    svg.append("<metadata id='genomic-region-metrics'>" + html.escape(json.dumps(metric, sort_keys=True)) + '</metadata>')
    for x, label in [(535, 'Sample'), (658, 'Background'), (796, 'Enrichment')]:
        svg.append(text(x, 60, label, 'end', 13))
    for index, (label, key) in enumerate(REGION_ROWS):
        y = 87 + index * 36
        rate = metric.get(key + '_rate')
        background = metric.get(key + '_background_rate')
        fold = metric.get(key + '_enrichment')
        count = metric.get(key + '_reads')
        value = region_percent(rate)
        if count is not None and rate is not None:
            value += f' ({count:,})'
        svg.append(f"<g data-region='{key}'>")
        svg.append(text(24, y + 4, label, 'start', 15, '#203040'))
        svg.append(region_track(rate, background, color, 125, y, 272))
        for x, label_value in [(535, value), (658, region_percent(background)), (796, '—' if fold is None else f'{fold:,.2f}×')]:
            svg.append(text(x, y + 4, label_value, 'end', 14, '#34495e'))
        svg.append('</g>')
    svg.append(text(24, 198, REGION_LEGEND, 'start', 12))
    total = metric.get('annotation_mapped_reads', 0)
    window = metric.get('promoter_window_bp', 2000)
    svg.append(text(24, 220, f'n={total:,} MAPQ-pass alignments · Promoter: TSS ±{window/1000:g} kb', 'start', 12))
    return ''.join(svg) + '</svg>'
