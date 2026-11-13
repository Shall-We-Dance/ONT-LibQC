#!/usr/bin/env python3
"""Prepare assembly-checked, TE-only BED from an operator-downloaded RepeatMasker file."""
import argparse,collections,gzip,hashlib,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'workflow/scripts'))
from report_regions import merge

TE_CLASSES={'DNA','LINE','SINE','LTR','RC','Retroposon','PLE','Helitron'}

def read_sizes(path):
    return {v[0]:int(v[1]) for line in Path(path).read_text().splitlines() if (v:=line.split())}

def prepare(source,fmt,sizes_path,source_sizes_path,output,assembly,source_url,aliases_path=None,primary=None):
    sizes=read_sizes(sizes_path);source_sizes=read_sizes(source_sizes_path);aliases={c:c for c in sizes}
    for c in sizes:aliases.setdefault(c[3:] if c.startswith('chr') else 'chr'+c,c)
    if aliases_path:
        for line in Path(aliases_path).read_text().splitlines():
            if line.startswith('#'):continue
            fields=line.split();matches={aliases[x] for x in fields if x in aliases}
            if len(matches)==1:
                for name in fields:aliases.setdefault(name,next(iter(matches)))
    matching={}
    for name,size in source_sizes.items():
        c=aliases.get(name)
        if c in sizes:
            if sizes[c]!=size:raise ValueError(f'Assembly chromosome length mismatch: {c}')
            matching[name]=c
    missing=set(primary or sizes)-set(matching.values())
    if missing:raise ValueError(f'Repeat source lacks required chromosomes: {sorted(missing)}')
    opener=gzip.open if str(source).endswith('.gz') else open
    intervals=collections.defaultdict(list);classes=collections.Counter();excluded=collections.Counter();rows=0
    with opener(source,'rt') as handle:
        for line in handle:
            fields=line.split()
            if not fields or not fields[0].isdigit():continue
            rows+=1
            if fmt=='rmsk':chrom,start,end,kind=fields[5],int(fields[6]),int(fields[7]),fields[11]
            else:chrom,start,end,kind=fields[4],int(fields[5])-1,int(fields[6]),fields[10]
            category=kind.split('/')[0].rstrip('?')
            if category not in TE_CLASSES:excluded[category]+=1;continue
            if chrom not in matching:excluded['nonshared_contig']+=1;continue
            c=matching[chrom]
            if 0<=start==end<=sizes[c]:excluded['zero_length_fragment']+=1;continue
            if not 0<=start<end<=sizes[c]:raise ValueError('Repeat interval exceeds validated reference bounds')
            intervals[c].append((start,end));classes[category]+=1
    if not rows:raise ValueError('No RepeatMasker records parsed')
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True);count=bases=0
    with gzip.open(output,'wt') as out:
        for c in sizes:
            for s,e in merge(intervals.get(c,[])):
                out.write(f'{c}\t{s}\t{e}\n');count+=1;bases+=e-s
    metadata=dict(assembly=assembly,source=source_url,format=fmt,source_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest(),
        bed_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),merged_intervals=count,covered_bases=bases,
        included_classes=dict(classes),excluded_classes=dict(excluded),matching_chromosomes=len(matching),
        definition='Union of classified transposable elements; excludes simple repeats, low complexity, satellites, RNA and unclassified repeats.')
    output.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return metadata

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ['source','format','sizes','source-sizes','output','assembly','source-url']:p.add_argument('--'+flag,required=True)
    p.add_argument('--aliases');p.add_argument('--primary',help='Comma-separated primary chromosome names')
    a=p.parse_args()
    result=prepare(a.source,a.format,a.sizes,a.source_sizes,a.output,a.assembly,a.source_url,a.aliases,a.primary.split(',') if a.primary else None)
    print(json.dumps(result,indent=2))
