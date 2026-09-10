#!/usr/bin/env python3
"""Host-competitive screening for a curated panel of candidate spike-in species."""
import argparse,collections,csv,gzip,json,math,re,subprocess
from pathlib import Path

def wilson(k,n):
 if not n:return [None,None]
 z=1.959963984540054;p=k/n;d=1+z*z/n;c=(p+z*z/(2*n))/d;w=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d;return [max(0,c-w),min(1,c+w)]
def signal_status(n,loci,denominator):
 if n>=3 and loci>=2:return 'CANDIDATE_SIGNAL'
 if n:return 'TRACE'
 return 'NOT_DETECTED' if denominator else 'NO_SCREENABLE_INSERTS'
def classify_alignment(f,contigs,minimum_bp=80,minimum_coverage=.8,minimum_identity=.85,mapq=20,score_margin=10):
 flag=int(f[1]);mq=int(f[4]);species=contigs.get(f[2],'unknown')
 if flag&(4|256|2048):return None
 if species=='unknown':raise ValueError('Competitive index and contig manifest disagree: '+f[2])
 ops=[(int(n),o) for n,o in re.findall(r'(\d+)([MIDNSHP=X])',f[5])];aligned=sum(n for n,o in ops if o in 'M=X');tags={x[:2]:x[5:] for x in f[11:] if len(x)>5}
 cov=aligned/max(1,len(f[9]));identity=1-int(tags.get('NM',999999))/max(1,aligned);margin=int(tags.get('AS',-999999))-int(tags.get('XS',-999999))
 if mq<mapq or mq==255 or aligned<minimum_bp or cov<minimum_coverage or identity<minimum_identity or margin<score_margin:return None
 return dict(species=species,contig=f[2],position=int(f[3])-1,locus=f[2]+':'+str((int(f[3])-1)//100),mapq=mq,aligned_bp=aligned,coverage=cov,identity=identity,score_margin=margin)
def scoped_summary(result,allowed=None):
 assignments=result.get('candidate_read_species',{});evidence=result.get('qualifying_read_loci',{});screened=set(result.get('screened_read_ids',[]))
 if allowed is not None:screened &= set(allowed)
 n=len(screened);selected={r:s for r,s in assignments.items() if r in screened};counts=collections.Counter(selected.values());loci={evidence[r] for r in selected};total=len(selected);status=signal_status(total,len(loci),n)
 return dict(spike_in_status=status,spike_in_label=', '.join(result.get('candidates',{}).get(k,{}).get('label',k) for k in sorted(counts)) or 'No candidate signal',spike_in_screened_reads=n,spike_in_reads=total,spike_in_fraction=total/n if n else None,spike_in_species_counts=dict(counts),spike_in_fraction_ci95=wilson(total,n),spike_in_declared=result.get('declared_spike_in',False))
def main():
 p=argparse.ArgumentParser();p.add_argument('--fastq',required=True);p.add_argument('--output',required=True);p.add_argument('--database',required=True);p.add_argument('--host',required=True);p.add_argument('--bowtie2',required=True);p.add_argument('--threads',type=int,default=4);p.add_argument('--mapq',type=int,default=20);p.add_argument('--minimum-aligned-bp',type=int,default=80);p.add_argument('--minimum-query-coverage',type=float,default=.8);p.add_argument('--minimum-identity',type=float,default=.85);p.add_argument('--minimum-score-margin',type=int,default=10);p.add_argument('--declared-spike-in',action='store_true');a=p.parse_args()
 db=json.loads(Path(a.database).read_text());host=db['host_indexes'][a.host];contigs=json.loads(Path(host['contigs']).read_text());entries={x['id']:x for x in db['entries']};output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
 manifest_dir=output.parent/'database_manifests';manifest_dir.mkdir(exist_ok=True);snapshot=manifest_dir/(output.stem+'.json')
 snapshot.write_text(json.dumps(db,indent=2)+'\n')
 command=[a.bowtie2,'--very-sensitive-local','-p',str(a.threads),'-x',host['index'],'-U',a.fastq]
 screened=set();hits={};loci={};rows=[];all_counts=collections.Counter()
 with output.with_suffix('.alignment.log').open('w') as log:
  process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=log,text=True)
  for line in process.stdout:
   if line.startswith('@'):continue
   f=line.rstrip().split('\t')
   if int(f[1])&(256|2048):continue
   if f[0] in screened:raise ValueError('Duplicate primary query in competitive screen')
   screened.add(f[0]);hit=classify_alignment(f,contigs,a.minimum_aligned_bp,a.minimum_query_coverage,a.minimum_identity,a.mapq,a.minimum_score_margin)
   if not hit:all_counts['ambiguous_or_unassigned']+=1;continue
   group=hit['species'];all_counts[group]+=1
   rows.append(dict(unit_id=f[0],**hit))
   if group in entries and entries[group]['role']=='spikein_candidate':hits[f[0]]=group;loci[f[0]]=hit['locus']
  if process.wait():raise RuntimeError('Competitive spike-in alignment failed; see '+str(output.with_suffix('.alignment.log')))
 n=len(screened);candidate_stats={}
 for key,entry in entries.items():
  evidence=[r for r in rows if r['species']==key];count=len(evidence);nl=len({r['locus'] for r in evidence});excluded=key in host['excluded_candidates']
  candidate_stats[key]=dict(label=entry['label'],role=entry['role'],reads=None if excluded else count,fraction=None if excluded or not n else count/n,loci=None if excluded else nl,status='HOST_NOT_DISTINGUISHABLE' if excluded else signal_status(count,nl,n),fraction_ci95=[None,None] if excluded else wilson(count,n))
 result=dict(mode='auto',label='Common spike-in candidates',host_reference=a.host,declared_spike_in=a.declared_spike_in,screened_reads=n,screened_read_ids=sorted(screened),qualifying_reads=len(hits),qualifying_fraction=len(hits)/n if n else None,qualifying_read_ids=sorted(hits),candidate_read_species=hits,qualifying_read_loci=loci,status=signal_status(len(hits),len(set(loci.values())),n),candidates=candidate_stats,assignment_counts=dict(all_counts),fraction_ci95=wilson(len(hits),n),zero_hit_upper95=1-math.pow(.05,1/n) if n else None,command=command,database_source=a.database,database_manifest=str(snapshot.resolve()),input_path=str(Path(a.fastq).resolve()),parameters={k:getattr(a,k) for k in ['mapq','minimum_aligned_bp','minimum_query_coverage','minimum_identity','minimum_score_margin']},interpretation='Candidate foreign-sequence signal, not proof that spike-in was intentionally added. Experimental spike-in intent must be established from metadata; positive signals require investigation. Zero hits do not establish absence, and the binomial bound assumes independent inserts and perfect detectability. Host and candidate genomes plus Mollicutes/bacterial decoys compete in one alignment.')
 with gzip.open(output.with_suffix('.assignments.tsv.gz'),'wt') as f:
  writer=csv.DictWriter(f,fieldnames=['unit_id','species','contig','position','locus','mapq','aligned_bp','coverage','identity','score_margin'],delimiter='\t');writer.writeheader();writer.writerows(rows)
 output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['host_reference','screened_reads','qualifying_reads','status','assignment_counts']},indent=2))
if __name__=='__main__':main()
