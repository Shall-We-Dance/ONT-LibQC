#!/usr/bin/env python3
"""Build reproducible host-competitive references for common spike-in candidates."""
import argparse,concurrent.futures,datetime,gzip,json,shutil,subprocess,urllib.request
from pathlib import Path
import yaml

def fasta(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        name=None;seq=[]
        for line in f:
            if line.startswith('>'):
                if name:yield name,''.join(seq)
                name=line[1:].split()[0];seq=[]
            else:seq.append(line.strip())
        if name:yield name,''.join(seq)
def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    cfg=yaml.safe_load(Path(a.config).read_text());root=Path(a.out).resolve();root.mkdir(parents=True,exist_ok=True);sources=root/'sources';sources.mkdir(exist_ok=True)
    downloads={
      's_cerevisiae':('https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz','sacCer3.fa.gz','Saccharomyces cerevisiae','sacCer3 / S288C','spikein_candidate','https://www.cellsignal.com/products/cut-run-kits-reagents/cut-run-pag-mnase-and-spike-in-dna/40366'),
      's_pombe':('https://www.pombase.org/data/genome_sequence_and_features/genome_sequence/Schizosaccharomyces_pombe_all_chromosomes.fa.gz','s_pombe.fa.gz','Schizosaccharomyces pombe','PomBase retrieved 2026-09-08','spikein_candidate','https://www.pombase.org/'),
      'lambda':('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001416.1&rettype=fasta&retmode=text','lambda.fa','Lambda phage','NC_001416.1','sequencing_control','https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1'),
      'phix':('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text','phix.fa','PhiX174','NC_001422.1','sequencing_control','https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1')}
    entries=[]
    for key,(url,name,species,assembly,role,evidence) in downloads.items():
        path=sources/name
        if not path.exists():
            print('Downloading',key,flush=True)
            with urllib.request.urlopen(url,timeout=90) as src,open(str(path)+'.tmp','wb') as dest:shutil.copyfileobj(src,dest)
            Path(str(path)+'.tmp').rename(path)
        records=list(fasta(path));assert records and sum(len(s) for _,s in records)>1000,(key,'invalid FASTA')
        entries.append(dict(id=key,label=species,assembly=assembly,role=role,source_url=url,evidence_url=evidence,input=str(path),size_bytes=path.stat().st_size))
    source_config=cfg.get('spike_in',{}).get('database_sources',{})
    bacterial=Path(source_config['bacterial_decoys'])
    ecoli=sources/'e_coli.fa';decoys=sources/'bacterial_decoys.fa'
    with ecoli.open('w') as ef,decoys.open('w') as df:
        for name,seq in fasta(bacterial):
            dest=ef if 'escherichia_coli' in name.lower() else df;dest.write('>'+name+'\n'+seq+'\n')
    assert sum(len(s) for _,s in fasta(ecoli))==4641652
    entries.append(dict(id='e_coli',label='Escherichia coli',assembly='GCF_000005845.2 / K-12 MG1655 representative',role='spikein_candidate',source_url='https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/',evidence_url='https://www.epicypher.com/content/documents/tds/18-1401.pdf',input=str(ecoli),size_bytes=ecoli.stat().st_size))
    fly=Path(source_config.get('fly_fasta',cfg.get('references',{}).get('dm6',{}).get('fasta','')))
    if not fly.is_file():raise ValueError('Set spike_in.database_sources.fly_fasta or references.dm6.fasta')
    entries.append(dict(id='d_melanogaster',label='Drosophila melanogaster',assembly='dm6',role='spikein_candidate',source_url='https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/',evidence_url='https://www.cellsignal.com/applications/cut-and-run',input=str(fly),size_bytes=fly.stat().st_size,host_reference='dm6'))
    decoy_paths=[decoys,Path(source_config['mollicutes'])]
    manifest=dict(version='spikein_candidates_v1',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),entries=entries,host_indexes={},decoys=[dict(path=str(x),size_bytes=x.stat().st_size) for x in decoy_paths],limitations=['Sequence matches identify candidate foreign material, not the intention to add spike-in.','Same-species spike-in cannot be separated from host without specific barcodes.','E. coli uses a K-12 representative; strain-specific regions and unlisted species can be missed.','Synthetic standards such as ERCC and barcoded nucleosomes require exact product sequences; absent from this default DNA panel.'])
    bowtie=cfg.get('tools',{}).get('bowtie2','bowtie2')
    build=cfg.get('tools',{}).get('bowtie2_build') or str(Path(shutil.which(bowtie) or bowtie).with_name('bowtie2-build'))
    def host_build(host):
        folder=root/host;folder.mkdir(exist_ok=True);fa=folder/'combined.fa';prefix=folder/'index';contigs={}
        inputs=[('host',Path(cfg['references'][host]['fasta']))]+[(x['id'],Path(x['input'])) for x in entries if x.get('host_reference')!=host]+[('decoy',x) for x in decoy_paths]
        source_files=[dict(role=group,path=str(path),size_bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns) for group,path in inputs]
        complete=folder/'complete.json';cached=json.loads(complete.read_text()) if complete.exists() else {}
        index_ready=any(all(Path(str(prefix)+'.'+part+'.'+suffix).exists() for part in ['1','2','3','4','rev.1','rev.2']) for suffix in ['bt2','bt2l'])
        reuse=cached.get('source_files')==source_files and index_ready and fa.exists() and (folder/'contigs.json').exists()
        print(('Reusing reference with matching file metadata' if reuse else 'Building reference'),host,flush=True)
        if not reuse:
            with fa.open('w') as out:
                for group,path in inputs:
                    for i,(name,seq) in enumerate(fasta(path)):
                        cname=f'{group}__{path.stem}__{i}__{name}';contigs[cname]=group
                        out.write('>'+cname+'\n')
                        for start in range(0,len(seq),80):out.write(seq[start:start+80]+'\n')
            (folder/'contigs.json').write_text(json.dumps(contigs)+'\n')
            command=[str(build),'--threads',str(a.threads),str(fa),str(prefix)]
            with (folder/'build.log').open('w') as log:subprocess.run(command,stdout=log,stderr=log,check=True)
            info=dict(source_files=source_files,index=str(prefix),contigs=str(folder/'contigs.json'),combined_fasta=str(fa),fasta_size_bytes=fa.stat().st_size,command=command,host=host,excluded_candidates=[x['id'] for x in entries if x.get('host_reference')==host])
            (folder/'complete.json').write_text(json.dumps(info,indent=2)+'\n')
        return host,json.loads((folder/'complete.json').read_text())
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for host,info in pool.map(host_build,sorted(cfg['references'])):manifest['host_indexes'][host]=info
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print('Database complete',root,flush=True)
if __name__=='__main__':main()
