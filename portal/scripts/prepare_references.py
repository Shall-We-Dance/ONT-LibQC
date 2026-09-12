#!/usr/bin/env python3
"""Prepare administrator-approved references; this is never called by an upload request."""
import argparse, fcntl, gzip, json, os, re, shutil, subprocess, sys, tempfile, time, urllib.request
from pathlib import Path
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ont_portal.catalog import public_catalog, bt2_ready, bwa_ready
from ont_portal.minimap_prep import prepare_minimap_indices
from ont_portal.reference_prep import GIB, memory_bytes, tool_memory, check_headroom, normalize_names

DEFAULT_SOURCES=Path(__file__).resolve().parents[1]/'deploy/references.common.yaml'

SAFE=re.compile(r'^[A-Za-z][A-Za-z0-9_-]{0,31}$')

def save(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)

def absolute(value,base):
    p=Path(value).expanduser();return str(p.absolute() if p.is_absolute() else (base/p).absolute())

def adopt(config,out):
    origin=Path(config).resolve();cfg=yaml.safe_load(origin.read_text());catalog={'references':{},'mycoplasma_regex':cfg.get('mycoplasma_screen',{}).get('target_contig_regex','mycoplasma_target')}
    for key,ref in cfg['references'].items():
        if not SAFE.fullmatch(key):raise ValueError('Unsafe reference ID')
        row={k:absolute(ref[k],origin.parent) for k in ('fasta','chrom_sizes','annotation','bowtie2_index','bwa_mem2_index')}
        if ref.get('minimap2_index'):row['minimap2_index']=absolute(ref['minimap2_index'],origin.parent)
        if ref.get('mycoplasma_index'):row['mycoplasma_index']=absolute(ref['mycoplasma_index'],origin.parent)
        row.update(assembly=ref.get('assembly',key),primary_chromosomes=ref['primary_chromosomes'],label=ref.get('label',{'hg38':'Human · hg38 / GRCh38','dm6':'Fruit fly · dm6'}.get(key,key)))
        spike=ref.get('spike_in',cfg.get('spike_in',{}))
        if spike.get('enabled',True) and (spike.get('database_manifest') or spike.get('combined_bowtie2_index')):
            row['spike_in']={k:spike[k] for k in ('mode','label','target_contig_regex') if k in spike}
            for k in ('database_manifest','combined_bowtie2_index'):
                if spike.get(k):row['spike_in'][k]=absolute(spike[k],origin.parent)
        catalog['references'][key]=row
    ready=public_catalog(catalog)
    if not all(x['ready'] for x in ready['references']):raise ValueError('One or more references lack FASTA, GTF, chromosome sizes, Bowtie2 or BWA-MEM2 indexes')
    save(out,catalog);save(Path(out).with_name('catalog.public.json'),ready)
    print(json.dumps({'catalog':str(Path(out).resolve()),'references':[{k:r[k] for k in ('id','ready','mappers','mycoplasma','spike_in')} for r in ready['references']]},indent=2))

def materialize(value,destination,base):
    """Download only explicit administrator HTTPS sources, or copy local files."""
    dest=Path(destination);partial=dest.with_name(dest.name+'.download')
    if str(value).startswith('https://'):
        for attempt in range(3):
            try:
                with urllib.request.urlopen(str(value),timeout=120) as src,partial.open('wb') as out:
                    if not src.url.startswith('https://'):raise ValueError('Reference download redirected away from HTTPS')
                    shutil.copyfileobj(src,out,length=1024**2)
                break
            except (OSError, TimeoutError):
                partial.unlink(missing_ok=True)
                if attempt==2:raise
                time.sleep(2*(attempt+1))
    elif '://' in str(value):raise ValueError('Only HTTPS reference downloads are supported')
    else:shutil.copyfile(absolute(value,base),partial)
    with partial.open('rb') as f:compressed=f.read(2)==b'\x1f\x8b'
    if compressed:
        with gzip.open(partial,'rb') as src,dest.open('wb') as out:shutil.copyfileobj(src,out,length=1024**2)
        partial.unlink()
    else:partial.replace(dest)
    if not dest.stat().st_size:raise ValueError('Empty reference source')

def fasta_records(path):
    name=None;chunks=[]
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if name:yield name,''.join(chunks)
                name=line[1:].split()[0];chunks=[]
            else:chunks.append(line.strip())
    if name:yield name,''.join(chunks)

def prefix_fasta(path,out,group,contigs):
    # Stream sequence lines; a human chromosome need not fit in RAM.
    count=0;seen=False
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                count+=1;seen=True;original=line[1:].split()[0]
                name=group+'__'+str(count)+'__'+re.sub(r'[^A-Za-z0-9_.-]','_',original)
                contigs[name]=group;out.write('>'+name+'\n')
            else:
                if not seen and line.strip():raise ValueError('Malformed reference FASTA')
                out.write(line)
    if not count:raise ValueError('Reference contains no FASTA records')

def build(args):
    source=Path(args.sources).resolve();spec=yaml.safe_load(source.read_text());root=Path(args.output).resolve();root.mkdir(parents=True,exist_ok=True)
    # Keep independent invocations from mutating the same prepared-reference tree.
    lock=(root/'.prepare.lock').open('a')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise ValueError('Another reference build is already using this output directory')
    selected=set(args.only or spec['references'])
    if selected-set(spec['references']):raise ValueError('Unknown reference IDs: '+', '.join(sorted(selected-set(spec['references']))))
    catalog={'references':{},'mycoplasma_regex':'^mycoplasma_target_'}
    if (root/'catalog.admin.json').is_file():catalog=json.loads((root/'catalog.admin.json').read_text())
    reused=json.loads(Path(args.reuse_catalog).read_text()) if args.reuse_catalog else {'references':{}}
    if reused.get('mycoplasma_regex'):catalog['mycoplasma_regex']=reused['mycoplasma_regex']
    def publish():
        catalog['references']={key:catalog['references'][key] for key in spec['references'] if key in catalog['references']}
        public=public_catalog(catalog)
        if not all(row['ready'] for row in public['references']):raise ValueError('Reference readiness check failed')
        save(root/'catalog.admin.json',catalog);save(root/'catalog.public.json',public)
    def tool(stage,argv,bases=0):
        budget=tool_memory(argv[0],bases,args.memory)
        check_headroom(stage,budget,memory_bytes(args.reserve_memory))
        name='ont-qc-reference-'+stage.name.removeprefix('.building-')+'-'+str(os.getpid())
        print('  '+argv[0]+': '+str(args.threads)+' CPU, '+str(budget//GIB)+' GiB cap',flush=True)
        cmd=['docker','run','--rm','--name',name,'--label','ont-qc.role=reference-build','--user',f'{os.getuid()}:{os.getgid()}','--network','none','--read-only','--cap-drop','ALL','--security-opt','no-new-privileges',
            '--cpus',str(args.threads),'--cpu-shares','128','--memory',str(budget),'--memory-swap',str(budget),'--pids-limit','256','--tmpfs','/tmp:rw,nosuid,nodev,size=512m',
            '-e','OMP_NUM_THREADS='+str(args.threads),'--mount',f'type=bind,src={stage},dst=/build','--entrypoint',argv[0],args.image,*argv[1:]]
        try:
            with (stage/'build.log').open('ab') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True,timeout=86400)
        except BaseException:
            subprocess.run(['docker','rm','-f',name],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=30)
            raise
    for key,r in sorted(spec['references'].items(),key=lambda item:int(item[1].get('sequence_bases',0))):
        if key not in selected:continue
        if not SAFE.fullmatch(key):raise ValueError('Reference IDs must be simple letters/numbers')
        folder=root/key;done=folder/'ready.json'
        if done.is_file():
            saved=json.loads(done.read_text())
            if saved['source_settings']!=r:raise ValueError('Reference sources changed. Build into a new output directory.')
            row=saved['reference']
            if not public_catalog({'references':{key:row}})['references'][0]['ready']:raise ValueError('Incomplete existing reference; use a new output directory')
            catalog['references'][key]=row;publish();print('Reusing',key,flush=True);continue
        if key in reused['references']:
            row=reused['references'][key]
            if row.get('assembly',key)!=r.get('assembly',key):raise ValueError('Reused reference has a different assembly: '+key)
            if not public_catalog({'references':{key:row}})['references'][0]['ready']:raise ValueError('Reused reference is incomplete: '+key)
            catalog['references'][key]=row;publish();print('Reusing prepared operator reference',key,flush=True);continue
        if folder.exists():raise ValueError('Incomplete reference directory exists; inspect it and choose a fresh output directory')
        bases=int(r.get('sequence_bases',0))
        if bases:tool_memory('bwa-mem2',bases,args.memory)
        if shutil.disk_usage(root).free < 10*GIB+16*bases:raise ValueError('Insufficient disk space for reference build plus reserved space')
        stage=root/('.building-'+key)
        if stage.exists():
            if not (stage/'source_settings.json').is_file() or json.loads((stage/'source_settings.json').read_text())!=r:
                raise ValueError('Interrupted build has different or unknown sources: '+str(stage))
            print('Resuming',key,flush=True)
        else:
            stage.mkdir();save(stage/'source_settings.json',r);print('Preparing',key,flush=True)
        progress_file=stage/'progress.json'
        progress=json.loads(progress_file.read_text()) if progress_file.is_file() else {}
        def step(name,action,exists):
            if progress.get(name) and exists():return
            action()
            if not exists():raise ValueError('Incomplete preparation step: '+name)
            progress[name]=True;save(progress_file,progress)
        for kind,filename in [('fasta','genome.fa'),('annotation','genes.gtf')]+([('assembly_report','assembly_report.txt')] if r.get('rename_from_assembly_report') else []):
            step('download_'+kind,lambda kind=kind,filename=filename:materialize(r[kind],stage/filename,source.parent),lambda filename=filename:(stage/filename).is_file() and (stage/filename).stat().st_size>0)
        if r.get('rename_from_assembly_report'):
            step('normalize',lambda:normalize_names(stage),lambda:(stage/'genome.fa').is_file() and (stage/'genes.gtf').is_file())
        step('faidx',lambda:tool(stage,['samtools','faidx','/build/genome.fa']),lambda:(stage/'genome.fa.fai').is_file())
        sizes=''.join('\t'.join(line.split('\t')[:2])+'\n' for line in (stage/'genome.fa.fai').read_text().splitlines())
        (stage/'chrom.sizes').write_text(sizes)
        chroms={x.split()[0] for x in sizes.splitlines()}
        if not set(r['primary_chromosomes'])<=chroms:raise ValueError('Primary chromosome names do not match FASTA')
        bases=sum(int(line.split()[1]) for line in sizes.splitlines())
        (stage/'bowtie2').mkdir(exist_ok=True);(stage/'bwa').mkdir(exist_ok=True)
        # Larger blocks reduce repeated genome scans within the existing per-tool memory cap.
        step('bowtie2',lambda:tool(stage,['bowtie2-build','--threads',str(args.threads),'--bmaxdivn','2','/build/genome.fa','/build/bowtie2/index'],bases),lambda:bt2_ready(stage/'bowtie2/index'))
        step('bwa',lambda:tool(stage,['bwa-mem2','index','-p','/build/bwa/genome','/build/genome.fa'],bases),lambda:bwa_ready(stage/'bwa/genome'))
        row={'assembly':r.get('assembly',key),'label':r.get('label',key),'primary_chromosomes':r['primary_chromosomes'],
             'fasta':str(folder/'genome.fa'),'annotation':str(folder/'genes.gtf'),'chrom_sizes':str(folder/'chrom.sizes'),
             'bowtie2_index':str(folder/'bowtie2/index'),'bwa_mem2_index':str(folder/'bwa/genome')}
        row.update({name:r[name] for name in ('assembly_accession','source_assembly','annotation_release','source_page') if name in r})
        screens=r.get('screens',{})
        for screen in ('mycoplasma','spike_in'):
            panel=screens.get(screen)
            if not panel:continue
            targetlist=panel.get('targets',[])
            if not targetlist:raise ValueError(screen+' needs target FASTA sources')
            screen_dir=stage/screen;screen_dir.mkdir(exist_ok=True);contigs={};entries=[]
            with (screen_dir/'combined.fa').open('w') as out:
                prefix_fasta(stage/'genome.fa',out,'host',contigs)
                for i,entry in enumerate(targetlist):
                    identifier=entry.get('id','target'+str(i))
                    if not SAFE.fullmatch(identifier):raise ValueError('Unsafe panel target ID')
                    if entry.get('host_reference')==key:continue
                    raw=screen_dir/(identifier+'.fa');materialize(entry['fasta'],raw,source.parent)
                    group='mycoplasma_target' if screen=='mycoplasma' else identifier
                    # Add the file index to ensure different targets cannot duplicate contig IDs.
                    temp={}
                    prefix_fasta(raw,out,group+'_'+str(i) if screen=='mycoplasma' else group,temp)
                    contigs.update({k:('mycoplasma_target' if screen=='mycoplasma' else group) for k in temp})
                    entries.append({'id':identifier,'label':entry.get('label',identifier),'role':entry.get('role','spikein_candidate'),'assembly':entry.get('assembly','operator-supplied')})
                    raw.unlink()
                for i,entry in enumerate(panel.get('decoys',[])):
                    raw=screen_dir/('decoy'+str(i)+'.fa');materialize(entry,raw,source.parent)
                    temp={};prefix_fasta(raw,out,'decoy'+str(i),temp);contigs.update({k:'decoy' for k in temp});raw.unlink()
            tool(stage,['bowtie2-build','--threads',str(args.threads),'/build/'+screen+'/combined.fa','/build/'+screen+'/index'],(screen_dir/'combined.fa').stat().st_size)
            (screen_dir/'combined.fa').unlink()
            if screen=='mycoplasma':row['mycoplasma_index']=str(folder/screen/'index');catalog['mycoplasma_regex']='^mycoplasma_target_'
            else:
                save(screen_dir/'contigs.json',contigs)
                save(screen_dir/'manifest.json',{'version':'operator-panel-v1','entries':entries,
                    'host_indexes':{key:{'index':str(folder/screen/'index'),'contigs':str(folder/screen/'contigs.json'),'excluded_candidates':[]}}})
                row['spike_in']={'mode':'auto','database_manifest':str(folder/screen/'manifest.json')}
        stage.replace(folder);save(done,{'source_settings':r,'reference':row});catalog['references'][key]=row
        publish();print('Ready',key,flush=True)
    catalog=prepare_minimap_indices(catalog,root/'minimap2',threads=args.threads,memory=args.memory,reserve_memory=args.reserve_memory,image=args.image,catalog_output=root)
    ready=public_catalog(catalog)
    if not all(r['ready'] for r in ready['references']):raise ValueError('Final reference readiness check failed')
    save(root/'catalog.admin.json',catalog);save(root/'catalog.public.json',ready)
    print('Catalog:',root/'catalog.admin.json')

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('adopt');a.add_argument('--config',required=True);a.add_argument('--output',required=True)
    b=sub.add_parser('build');b.add_argument('--sources',default=str(DEFAULT_SOURCES));b.add_argument('--output',required=True);b.add_argument('--threads',type=int,default=2);b.add_argument('--memory',default='96g',help='Maximum build memory; each tool gets a smaller estimated cap when possible');b.add_argument('--reserve-memory',default='16g');b.add_argument('--only',nargs='+',help='Reference IDs to prepare; default is the ten common species');b.add_argument('--reuse-catalog',help='Explicitly reuse prepared references of the same assembly');b.add_argument('--image',default='ont-qc-runner:1')
    m=sub.add_parser('minimap2',help='Add prepared Minimap2 indexes to an existing administrator catalog')
    m.add_argument('--catalog',required=True);m.add_argument('--output',required=True);m.add_argument('--threads',type=int,default=2)
    m.add_argument('--memory',default='96g');m.add_argument('--reserve-memory',default='16g');m.add_argument('--image',default='ont-qc-runner:1')
    args=p.parse_args()
    if args.command=='adopt':adopt(args.config,args.output)
    elif args.command=='minimap2':
        prepare_minimap_indices(json.loads(Path(args.catalog).read_text()),args.output,threads=args.threads,memory=args.memory,reserve_memory=args.reserve_memory,image=args.image)
        print('Catalog:',Path(args.output).resolve()/'catalog.admin.json')
    else:
        if not 1<=args.threads<=64:raise ValueError('Use between 1 and 64 threads')
        build(args)
if __name__=='__main__':main()
