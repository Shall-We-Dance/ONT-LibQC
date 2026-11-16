"""Within-sample genomic-bin correlations from assigned chromatin alignments.

No pooled-library column, random sampling or coordinate smoothing. Shared-zero
bins are retained by default. Alignment blocks count once per overlapping bin, regardless of depth.
"""
import collections
import csv
import gzip
import html
import io
import json
import math
from pathlib import Path
import re

BIN_SIZES = (10000, 100000, 1000000)
METHOD_VERSION = 'sublibrary-bin-correlation-v1'


def aligned_blocks(record):
    """Exclude skipped introns, deletions and soft clips from bin overlaps."""
    pos = int(record['start'])
    cigar = record.get('cigar')
    if not cigar:
        yield pos, int(record['end'])
        return
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(length)
        if op in 'M=X':
            yield pos, pos + length
        if op in 'MDN=X':
            pos += length


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + j - 1) / 2 + 1
        for index in order[i:j]:
            result[index] = rank
        i = j
    return result


def pearson(a, b, repeats=0, repeated_a=0, repeated_b=0):
    n = len(a) + repeats
    if n < 2:
        return None
    ma, mb = (sum(a)+repeats*repeated_a)/n, (sum(b)+repeats*repeated_b)/n
    va = sum((x-ma)**2 for x in a) + repeats*(repeated_a-ma)**2
    vb = sum((x-mb)**2 for x in b) + repeats*(repeated_b-mb)**2
    if not va or not vb:
        return None
    covariance = sum((x-ma)*(y-mb) for x,y in zip(a,b)) + repeats*(repeated_a-ma)*(repeated_b-mb)
    return max(-1.0, min(1.0, covariance / math.sqrt(va*vb)))


def matrix_for(columns, shared_zeros, method):
    if method == 'spearman':
        # Average ranks include the implicit shared-zero rows, without allocating them.
        zero_ranks = [(c.count(0)+shared_zeros+1)/2 for c in columns]
        values = [[rank+shared_zeros if x else zero_ranks[i] for x,rank in zip(c,ranks(c))]
                  for i,c in enumerate(columns)]
    else:
        values=columns;zero_ranks=[0]*len(columns)
    matrix=[[pearson(a,b,shared_zeros,zero_ranks[i],zero_ranks[j]) for j,b in enumerate(values)] for i,a in enumerate(values)]
    order,merges=cluster_order(matrix)
    return dict(matrix=matrix,order=order,merges=merges)


def cluster_order(matrix):
    """Deterministic average linkage on 1-r; undefined columns appear last."""
    valid = [i for i in range(len(matrix)) if matrix[i][i] is not None]
    missing = [i for i in range(len(matrix)) if matrix[i][i] is None]
    groups = [(i,) for i in valid]
    merges = []
    while len(groups) > 1:
        options = []
        for i, a in enumerate(groups):
            for j in range(i+1,len(groups)):
                b = groups[j]
                distance = sum(1-matrix[x][y] for x in a for y in b)/(len(a)*len(b))
                options.append((distance,a,b,i,j))
        distance,a,b,i,j = min(options)
        merges.append(dict(left=list(a),right=list(b),distance=distance))
        groups = [g for k,g in enumerate(groups) if k not in (i,j)] + [a+b]
        groups.sort()
    return (list(groups[0]) if groups else []) + missing, merges


def calculate(scopes, scope_records, sizes, primary, mapq=20, bin_sizes=BIN_SIZES):
    labels = [dict(id=s['scope_id'],label=s.get('scope_name',s['scope_id'])) for s in scopes]
    primary = [c for c in primary if c in sizes and sizes[c] > 0]
    eligible = []
    for scope in labels:
        records=[]
        for r in scope_records.get(scope['id'],{}).values():
            if r['chrom'] not in primary or r['mapq'] < mapq or r.get('flag',0) & (4|256|512|1024|2048):
                continue
            blocks=[(max(0,a),min(sizes[r['chrom']],b)) for a,b in aligned_blocks(r) if b>0 and a<sizes[r['chrom']] and a<b]
            if blocks:records.append((r['chrom'],blocks))
        eligible.append(records)
    result=dict(method_version=METHOD_VERSION,mapq=mapq,scopes=labels,
                eligible_alignments=[len(x) for x in eligible],default_bin_size=100000,
                default_method='spearman',primary_chromosomes=primary,
                chromosome_sizes={c:sizes[c] for c in primary},scales={})
    tables={}
    for size in bin_sizes:
        counts=[]
        for records in eligible:
            counter=collections.Counter()
            for chrom,blocks in records:
                bins={index for start,end in blocks for index in range(start//size,(end-1)//size+1)}
                counter.update((chrom,index) for index in bins)
            counts.append(counter)
        chrom_order={c:i for i,c in enumerate(primary)}
        union=sorted(set().union(*(set(x) for x in counts)), key=lambda k:(chrom_order[k[0]],k[1]))
        columns=[[c.get(key,0) for key in union] for c in counts]
        total_bins=sum((sizes[c]+size-1)//size for c in primary)
        shared_zeros=total_bins-len(union)
        matrices={method:matrix_for(columns,shared_zeros,method) for method in ('spearman','pearson')}
        nonzero={method:matrix_for(columns,0,method) for method in ('spearman','pearson')}
        result['scales'][str(size)]=dict(bin_size=size,bins_used=total_bins,signal_bins=len(union),
            all_zero_bins=shared_zeros,covered_bins=[len(c) for c in counts],correlations=matrices,
            nonzero_correlations=nonzero)
        stream=io.StringIO();writer=csv.writer(stream,delimiter='\t',lineterminator='\n')
        writer.writerow(['chrom','start','end']+[s['id'] for s in labels])
        for chrom,index in union:
            writer.writerow([chrom,index*size,min((index+1)*size,sizes[chrom])]+[c.get((chrom,index),0) for c in counts])
        tables[size]=stream.getvalue()
    return result,tables


def color(value):
    if value is None:return '#eef1f5'
    start=(247,249,252);end=(43,108,176) if value >= 0 else (182,66,57)
    t=abs(value)
    return '#'+''.join(f'{round(a+(b-a)*t):02x}' for a,b in zip(start,end))


def heatmap(result, size, method, skip_zeros=False):
    scale=result['scales'][str(size)];data=scale['nonzero_correlations' if skip_zeros else 'correlations'][method];order=data['order'];n=len(order)
    bins_used=scale['signal_bins'] if skip_zeros else scale['bins_used']
    cell=48;left=100;top=135;width=max(410,left+n*cell+32);height=top+n*cell+82
    esc=html.escape
    parts=[f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' style='font-family:Arial,Helvetica,sans-serif'>",
        f"<title>Sublibrary similarity · {method.title()} · {size//1000:,} kb</title>",
        "<metadata id='sublibrary-similarity'>"+esc(json.dumps(dict(method=method,bin_size=size,skip_all_zero_bins=skip_zeros,bins_used=bins_used,scopes=result['scopes'],**data)))+'</metadata>',
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{left}' y='25' font-size='16' fill='#203040'>{method.title()} correlation</text>",
        f"<text x='{left}' y='45' font-size='11' fill='#65758c'>{size//1000:,} kb · {bins_used:,} bins · {'nonzero union' if skip_zeros else 'whole genome'}</text>"]
    xpos={index:left+(j+.5)*cell for j,index in enumerate(order)}
    # Dendrogram uses a fixed 0–2 distance range, consistent across scales/methods.
    roots={(i,):(xpos[i],top-26) for i in order}
    for merge in data['merges']:
        a,b=tuple(merge['left']),tuple(merge['right']);xa,ya=roots[a];xb,yb=roots[b]
        y=top-26-merge['distance']*23
        parts.append(f"<path d='M{xa} {ya} V{y} H{xb} V{yb}' fill='none' stroke='#8c9bad' stroke-width='1.1'/>")
        roots[a+b]=((xa+xb)/2,y)
    for row,i in enumerate(order):
        sid=result['scopes'][i]['id'];label=esc(result['scopes'][i]['label']);y=top+row*cell
        parts.append(f"<text x='{left-12}' y='{y+cell/2+4}' text-anchor='end' font-size='13' fill='#34495e'><title>{label}</title>{esc(sid)}</text>")
        parts.append(f"<text x='{xpos[i]}' y='{top-8}' text-anchor='middle' font-size='12' fill='#34495e'>{esc(sid)}</text>")
        for col,j in enumerate(order):
            value=data['matrix'][i][j];text='—' if value is None else f'{value:.2f}';x=left+col*cell
            tip=esc(f"{result['scopes'][i]['label']} / {result['scopes'][j]['label']}: {text}")
            parts.append(f"<g><title>{tip}</title><rect x='{x}' y='{y}' width='{cell}' height='{cell}' fill='{color(value)}' stroke='white' stroke-width='2'/><text x='{x+cell/2}' y='{y+cell/2+4}' text-anchor='middle' font-size='13' fill='{('#fff' if value is not None and abs(value)>.65 else '#34495e')}'>{text}</text></g>")
    y=top+n*cell+22
    for i in range(100):parts.append(f"<rect x='{left+i*2}' y='{y}' width='2.1' height='8' fill='{color(-1+2*i/99)}'/>")
    for x,t in [(left,'−1'),(left+100,'0'),(left+200,'1')]:parts.append(f"<text x='{x}' y='{y+25}' text-anchor='middle' font-size='11' fill='#65758c'>{t}</text>")
    parts.append('</svg>');return ''.join(parts)


def write_assets(sample, scopes, records, sizes, primary, mapq, asset_dir):
    """Persist numeric data so reports remain re-renderable after BAM cleanup."""
    if len(scopes)<2:return None
    result,tables=calculate(scopes,records,sizes,primary,mapq)
    root=Path(asset_dir);root.mkdir(parents=True,exist_ok=True);prefix=f'{sample}.sublibraries.similarity'
    (root/f'{prefix}.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for size,table in tables.items():
        (root/f'{prefix}.{size}.counts.tsv.gz').write_bytes(gzip.compress(table.encode(),mtime=0))
        for method in ('spearman','pearson'):
            for skip_zeros in (False,True):
                suffix=method+('.nonzero' if skip_zeros else '')
                (root/f'{prefix}.{size}.{suffix}.svg').write_text(heatmap(result,size,method,skip_zeros))
                stream=io.StringIO();writer=csv.writer(stream,delimiter='\t',lineterminator='\n')
                ids=[s['id'] for s in result['scopes']];writer.writerow(['sublibrary']+ids)
                key='nonzero_correlations' if skip_zeros else 'correlations'
                for sid,row in zip(ids,result['scales'][str(size)][key][method]['matrix']):writer.writerow([sid]+['NA' if v is None else f'{v:.8g}' for v in row])
                (root/f'{prefix}.{size}.{suffix}.tsv').write_text(stream.getvalue())
    return result


CSS='''
.sublibrary-similarity{margin:22px 0 14px;padding-top:18px;border-top:1px solid #dce5f0}
.similarity-heading{display:flex;align-items:center;justify-content:space-between;gap:12px;flex-wrap:wrap;margin-bottom:12px}
.similarity-heading h3{margin:0!important;font-size:19px}.similarity-controls{display:flex;gap:12px;flex-wrap:wrap}
.similarity-controls label{font-size:12px;color:#65758c}.similarity-controls select{display:inline-block;width:auto;padding:5px 8px;margin-left:5px;font-size:13px;background:white;color:#34495e;border:1px solid #dce5f0;border-radius:5px}
.similarity-body{display:grid;grid-template-columns:minmax(280px,1fr) minmax(220px,.7fr);gap:24px;align-items:center;max-width:980px;margin:auto}
.similarity-panel[hidden]{display:none!important}.similarity-panel figure{margin:0}.similarity-panel img{display:block;width:100%;max-width:620px;margin:auto}
.similarity-info table{width:100%;font-size:13px}.similarity-info th,.similarity-info td{padding:7px 9px!important}.similarity-info p,.similarity-downloads{font-size:13px;line-height:1.6}
.similarity-warning{background:#fff2f1;border:1px solid #efcfcb;border-radius:6px;padding:9px 12px;font-size:13px;color:#93473f}
.similarity-method{font-size:13px;margin-top:10px;color:#65758c}.similarity-method summary{cursor:pointer}.similarity-method p{line-height:1.6}
@media(max-width:700px){.similarity-body{grid-template-columns:minmax(0,1fr);gap:8px}.similarity-controls{gap:8px}.similarity-heading h3{font-size:17px}.similarity-panel img{max-width:520px}}
@media print{.similarity-controls{display:none}.similarity-body{grid-template-columns:minmax(0,1fr) minmax(0,.65fr);gap:12px}.sublibrary-similarity{break-inside:avoid}.similarity-downloads{display:none}}
'''


def report_html(sample, asset_dir, asset_prefix='assets/'):
    prefix=f'{sample}.sublibraries.similarity';path=Path(asset_dir)/f'{prefix}.json'
    if not path.exists():return ''
    result=json.loads(path.read_text())
    if len(result['scopes'])<2:return ''
    esc=lambda v:html.escape(str(v),quote=True)
    parts=["<section class='sublibrary-similarity' aria-label='Sublibrary similarity'><style>"+CSS+'</style>',
        "<div class='similarity-heading'><h3>Sublibrary similarity</h3><div class='similarity-controls'><label>Correlation<select data-similarity-method aria-label='Similarity correlation method'><option value='spearman'>Spearman</option><option value='pearson'>Pearson</option></select></label><label>Window<select data-similarity-bin aria-label='Similarity genomic window'>"]
    for size in result['scales']:parts.append(f"<option value='{size}'{' selected' if int(size)==100000 else ''}>{'1 Mb' if int(size)==1000000 else str(int(size)//1000)+' kb'}</option>")
    parts.append("</select></label><label>Bins<select data-similarity-zeros aria-label='Similarity zero-bin policy'><option value='all'>All genomic bins</option><option value='nonzero'>Exclude shared zeros</option></select></label></div></div>")
    for size,scale in result['scales'].items():
        for method in ('spearman','pearson'):
            for skip_zeros in (False,True):
                stem=f'{asset_prefix}{prefix}.{size}';suffix=method+('.nonzero' if skip_zeros else '');url=stem+'.'+suffix+'.svg';hidden='' if size=='100000' and method=='spearman' and not skip_zeros else ' hidden'
                key='nonzero_correlations' if skip_zeros else 'correlations'
                bins_used=scale['signal_bins'] if skip_zeros else scale['bins_used']
                parts.append(f"<div class='similarity-panel' data-bin='{size}' data-method='{method}' data-zeros='{'nonzero' if skip_zeros else 'all'}'{hidden}><div class='similarity-body'><figure><a href='{esc(url)}' target='_blank' rel='noopener'><img loading='lazy' src='{esc(url)}' alt='{method.title()} sublibrary correlation, {int(size)//1000} kb bins'></a></figure><div class='similarity-info'><table><tr><th>Sublibrary</th><th>Alignments</th><th>Covered bins</th></tr>")
                for s,count,covered in zip(result['scopes'],result['eligible_alignments'],scale['covered_bins']):parts.append(f"<tr><td title='{esc(s['label'])}'>{esc(s['id'])}</td><td>{count:,}</td><td>{covered:,}</td></tr>")
                parts.append(f"</table><p>{bins_used:,} bins compared · MAPQ ≥{result['mapq']}<br>Primary nuclear chromosomes</p>")
                if skip_zeros:
                    parts.append("<p class='similarity-warning'>Excluding shared-zero bins can produce negative correlations at low depth.</p>")
                if any(n<200 for n in result['eligible_alignments']) or scale['signal_bins']<20:
                    parts.append("<p class='similarity-warning'>Low depth: correlations may be unstable.</p>")
                if any(scale[key][method]['matrix'][i][i] is None for i in range(len(result['scopes']))):parts.append("<p class='similarity-warning'>— indicates too few bins or no variation.</p>")
                parts.append(f"<p class='similarity-downloads'><a href='{esc(url)}' target='_blank' rel='noopener'>Open figure ↗</a> · <a href='{esc(stem+'.'+suffix+'.tsv')}'>Correlation matrix</a><br><a href='{esc(stem+'.counts.tsv.gz')}'>Bin counts</a></p></div></div></div>")
    parts.append("<details class='similarity-method'><summary>Method</summary><p>Each eligible alignment counts once in every window it overlaps. All genomic windows are retained by default. The optional zero filter removes only windows with zero counts in every sublibrary. Downloaded bin counts are sparse: omitted windows have zero counts in all sublibraries. Spearman uses average ranks for ties; Pearson uses raw counts (equivalent to per-library CPM scaling). Order uses average-linkage clustering on 1 − correlation. No downsampling or outlier removal.</p><p>Sparse sampling can lower correlations. This compares genomic distributions and does not establish replicate identity or experimental success.</p></details>")
    parts.append("""<script>(function(){var section=document.currentScript.closest('.sublibrary-similarity');function update(){var method=section.querySelector('[data-similarity-method]').value,bin=section.querySelector('[data-similarity-bin]').value,zeros=section.querySelector('[data-similarity-zeros]').value;section.querySelectorAll('.similarity-panel').forEach(function(panel){panel.hidden=panel.dataset.method!==method||panel.dataset.bin!==bin||panel.dataset.zeros!==zeros;});}section.querySelectorAll('select').forEach(function(select){select.addEventListener('change',update);});update();})();</script></section>""")
    return ''.join(parts)
