"""Prepare reusable Minimap2 indexes with the workflow's short-unit settings."""
import copy
import fcntl
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from .catalog import ReferenceMounts, minimap_ready, public_catalog
from .reference_prep import GIB, check_headroom, memory_bytes, tool_memory


def save(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def prepare_minimap_indices(catalog, output, *, threads=2, memory='96g', reserve_memory='16g', image='ont-qc-runner:1', catalog_output=None):
    if not 1<=threads<=64:raise ValueError('Use between 1 and 64 threads')
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=True)
    if ',' in str(root):raise ValueError('Commas are not supported in reference directory paths')
    destination=Path(catalog_output).resolve() if catalog_output else root
    lock=(root/'.prepare.lock').open('a')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise ValueError('Another Minimap2 preparation is using this output directory')
    result=copy.deepcopy(catalog)
    if not all(row['ready'] for row in public_catalog(result)['references']):raise ValueError('Base reference preparation is incomplete')
    def publish():
        public=public_catalog(result)
        if not all(row['ready'] for row in public['references']):raise ValueError('Base reference preparation is incomplete')
        save(destination/'catalog.admin.json',result);save(destination/'catalog.public.json',public)
    for key,row in sorted(result['references'].items(),key=lambda item:Path(item[1]['fasta']).stat().st_size):
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]{0,31}',key):raise ValueError('Unsafe reference ID')
        if minimap_ready(row.get('minimap2_index')):
            print('Reusing prepared Minimap2 reference',key,flush=True);continue
        fasta=Path(row['fasta']).resolve();stat=fasta.stat()
        with fasta.open('rb') as source:
            if not source.read(1024).lstrip().startswith(b'>'):raise ValueError('Minimap2 preparation requires uncompressed reference FASTA')
        settings={'fasta':str(fasta),'size':stat.st_size,'mtime_ns':stat.st_mtime_ns,'k':9,'w':5,'flags':0}
        folder=root/key;folder.mkdir(exist_ok=True)
        metadata=folder/'source_settings.json'
        if not metadata.exists() and (folder/'index.mmi').exists():raise ValueError('Existing Minimap2 index has no source record; use a new output directory')
        if metadata.exists() and json.loads(metadata.read_text())!=settings:
            raise ValueError('Minimap2 source changed; use a new output directory')
        save(metadata,settings)
        index=folder/'index.mmi'
        if not minimap_ready(index):
            temporary=folder/'index.mmi.tmp';temporary.unlink(missing_ok=True)
            if shutil.disk_usage(root).free<12*stat.st_size+10*GIB:raise ValueError('Insufficient disk space for Minimap2 preparation')
            budget=tool_memory('minimap2',stat.st_size,memory)
            check_headroom(root,budget,memory_bytes(reserve_memory))
            mounts=ReferenceMounts();target=mounts.path(str(fasta))
            name='ont-qc-minimap-index-'+key+'-'+str(os.getpid())
            print('Preparing Minimap2',key,':',threads,'CPU,',budget//GIB,'GiB cap',flush=True)
            command=['docker','run','--rm','--name',name,'--label','ont-qc.role=minimap-build','--user',f'{os.getuid()}:{os.getgid()}',
                '--network','none','--read-only','--cap-drop','ALL','--security-opt','no-new-privileges','--cpus',str(threads),'--cpu-shares','128',
                '--memory',str(budget),'--memory-swap',str(budget),'--pids-limit','128','--tmpfs','/tmp:rw,nosuid,nodev,size=128m',
                *mounts.arguments(),'--mount',f'type=bind,src={folder},dst=/build','--entrypoint','minimap2',image,
                '-x','sr','-k','9','-w','5','-I',str(max(8_000_000_000,stat.st_size+1)),'-t',str(threads),'-d','/build/index.mmi.tmp',target]
            try:
                with (folder/'build.log').open('ab') as log:subprocess.run(command,stdout=log,stderr=log,check=True,timeout=3600)
            except BaseException:
                subprocess.run(['docker','rm','-f',name],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=30)
                raise
            after=fasta.stat()
            if (after.st_size,after.st_mtime_ns)!=(stat.st_size,stat.st_mtime_ns):raise ValueError('Reference FASTA changed during indexing')
            if not minimap_ready(temporary):raise ValueError('Incomplete Minimap2 index or mismatched indexing settings')
            temporary.replace(index)
        row['minimap2_index']=str(index)
        publish();print('Ready Minimap2',key,flush=True)
    publish()
    return result
