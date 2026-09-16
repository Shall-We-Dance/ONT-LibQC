"""Trusted host dispatcher. Launches a fresh restricted container for every job.

The public web container has no Docker socket. Job containers have no network,
credentials, database, other jobs, or writable references.
"""
import argparse, multiprocessing, fcntl, json, math, resource, logging, os, re, shutil, signal, stat, subprocess, time, zipfile
from pathlib import Path
import yaml
from .limits import DEFAULTS, load_limits
from .schema import Options
from .report_files import REPORT_FILE_SUFFIXES
from .store import Store, ACTIVE, TERMINAL, token_for
from .catalog import job_config, public_catalog

LOG=logging.getLogger('ont-worker')
STOP=False

def command(args,**kwargs):return subprocess.run(args,check=True,**kwargs)

def remove_tree(path):
    path=Path(path)
    if path.is_symlink():path.unlink()
    elif path.exists():shutil.rmtree(path)

def valid_id(value):
    return isinstance(value,str) and re.fullmatch(r'[0-9a-f]{32}',value) is not None

def copy_uploads(source,dest,expected,bed):
    dest.mkdir(mode=0o700)
    directory=os.open(source,os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW)
    try:
        for name,limit in [('input.upload',expected)]+([('regions.upload',4*1024**2)] if bed else []):
            fd=os.open(name,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK,dir_fd=directory)
            with os.fdopen(fd,'rb') as src:
                info=os.fstat(src.fileno())
                if not stat.S_ISREG(info.st_mode) or not 0<info.st_size<=limit:raise ValueError('Invalid uploaded file')
                if name=='input.upload' and info.st_size!=expected:raise ValueError('Incomplete upload')
                with (dest/name).open('xb') as out:
                    copied=0
                    while data:=src.read(1024**2):
                        copied+=len(data)
                        if copied>limit:raise ValueError('Upload changed during processing')
                        out.write(data)
    finally:os.close(directory)

def tree_size(path):
    total=0
    for directory,dirs,files in os.walk(path,followlinks=False):
        for name in files:
            p=Path(directory)/name
            try:total+=p.lstat().st_size
            except FileNotFoundError:pass
    return total

def publish_report(source,dest,max_bytes=DEFAULTS['max_report_bytes'],max_html_bytes=DEFAULTS['max_html_bytes']):
    """Treat the sandbox output as untrusted: no symlinks, devices or arbitrary files."""
    if source.is_symlink() or not source.is_dir(): raise ValueError('Report must be a real directory')
    total=0;files=[]
    for path in source.rglob('*'):
        info=path.lstat()
        if stat.S_ISDIR(info.st_mode):continue
        if path.relative_to(source).parts==('.snakemake_timestamp',) and stat.S_ISREG(info.st_mode):continue
        if not stat.S_ISREG(info.st_mode) or path.is_symlink() or path.suffix not in REPORT_FILE_SUFFIXES or any(p.startswith('.') for p in path.relative_to(source).parts):
            raise ValueError('Unsupported file in report: '+str(path.relative_to(source)))
        if re.search(r'\.(?:bam|bai|sam|cram|fastq|fq|paf)(?:\.|$)',path.name,re.I):raise ValueError('Intermediate in report')
        if path.suffix=='.html' and info.st_size>max_html_bytes:raise ValueError('Report page exceeds the display limit')
        total+=info.st_size;files.append(path)
        if total>max_bytes or len(files)>2000:raise ValueError('Report exceeds publication limits')
    if not (source/'index.html').is_file() or not files:raise ValueError('Report has no index')
    if dest.exists():raise ValueError('Report destination already exists')
    staging=dest.with_name(dest.name+'.publishing');remove_tree(staging);staging.mkdir(mode=0o700)
    try:
        for path in files:
            target=staging/path.relative_to(source);target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(path,target,follow_symlinks=False)
        with zipfile.ZipFile(staging/'report.zip','w',zipfile.ZIP_DEFLATED,compresslevel=3) as archive:
            for path in files:
                rel=path.relative_to(source)
                archive.write(staging/rel,'report/'+str(rel))
        staging.replace(dest)
    except BaseException:remove_tree(staging);raise

def limit_log_output(max_bytes):
    # This applies to the host Docker client only, not the daemon-created job.
    resource.setrlimit(resource.RLIMIT_FSIZE,(max_bytes,max_bytes))


def memory_bytes(value):
    match=re.fullmatch(r'([1-9][0-9]*)([kmg]?)',str(value).lower())
    if not match:raise ValueError('Use a positive byte count or k/m/g memory size')
    return int(match[1])*1024**({'':0,'k':1,'m':2,'g':3}[match[2]])


def available_memory():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):return int(line.split()[1])*1024
    return 0


class Worker:
    def __init__(self,config):
        # Older deployments may have a flat memory limit above the new default cap.
        if not config.get('limits_file') and 'max_memory' not in config and 'memory' in config:
            config=dict(config,max_memory=config['memory'])
        config=load_limits(config)
        self.cfg=config; self.state=Path(config['state']);self.store=Store(self.state)
        for name in ('work','logs'):(self.state/name).mkdir(exist_ok=True,mode=0o700)
        self.catalog=json.loads(Path(config['catalog_file']).read_text());self.public=public_catalog(self.catalog)
        self.pipeline=Path(config['pipeline_dir']).resolve()
        if self.state.resolve()==self.pipeline or self.pipeline in self.state.resolve().parents:raise ValueError('Private runtime must be outside the workflow tree')
        self.image=config.get('runner_image','ont-qc-runner:1')
        self.cores=int(config['cores']);self.memory=config['memory']
        self.max_jobs=int(config['max_concurrent_jobs'])
        if not 1<=self.max_jobs<=8:raise ValueError('Invalid analysis concurrency limit')
        if not 1<=self.cores<=64:raise ValueError('Invalid worker CPU limit')
        if memory_bytes(config['max_memory'])<memory_bytes(self.memory):raise ValueError('max_memory is smaller than memory')
        self.timeout=int(config['timeout_seconds']);self.uid=os.getuid();self.gid=os.getgid()
        self.secret=json.loads(Path(config['secrets_file']).read_text())['signing_key']

    def docker_stop(self,jid):
        subprocess.run(['docker','rm','-f','ont-qc-job-'+jid],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=30)

    def memory_for(self, options):
        required=memory_bytes(self.memory)
        mapper=options.mapper if options.mapper!='auto' else ('bwa-mem2' if options.assay=='microc' else 'bowtie2')
        if mapper=='bwa-mem2' or options.assay=='microc':
            prefix=self.catalog['references'][options.reference]['bwa_mem2_index']
            index_size=sum(Path(prefix+suffix).stat().st_size for suffix in ('.bwt.2bit.64','.0123','.pac'))
            required=max(required,math.ceil((index_size+3*1024**3)/1024**3)*1024**3)
        if mapper=='minimap2':
            prefix=self.catalog['references'][options.reference].get('minimap2_index')
            if not prefix:raise ValueError('Minimap2 needs an administrator-prepared index')
            required=max(required,math.ceil((Path(prefix).stat().st_size+3*1024**3)/1024**3)*1024**3)
        if required>memory_bytes(self.cfg['max_memory']):
            raise ValueError('This reference index exceeds the online memory budget; run the pipeline locally')
        return required

    def claim(self):
        with self.store.connect(immediate=True) as c:
            active=c.execute("SELECT options FROM jobs WHERE status IN ('running','cancelling')").fetchall()
            if len(active)>=self.max_jobs:return None
            reserved=0
            for running in active:
                try:reserved+=self.memory_for(Options.model_validate_json(running['options']).check_catalog(self.public))
                except (ValueError,KeyError,OSError):reserved+=memory_bytes(self.cfg['max_memory'])
            row=c.execute("SELECT * FROM jobs WHERE status='queued' ORDER BY created,rowid LIMIT 1").fetchone()
            if row and not valid_id(row['id']):
                c.execute('DELETE FROM jobs WHERE id=?',(row['id'],));row=None
            if row:
                try: budget=self.memory_for(Options.model_validate_json(row['options']).check_catalog(self.public))
                except (ValueError,KeyError,OSError): budget=memory_bytes(self.memory)
                if available_memory()<reserved+budget+self.cfg['min_free_memory_bytes']:return None
                if shutil.disk_usage(self.state).free<self.cfg['min_free_bytes']+(len(active)+1)*self.cfg['max_work_bytes']:return None
                if tree_size(self.state/'reports')+2*(len(active)+1)*self.cfg['max_report_bytes']>self.cfg['max_report_total_bytes']:return None
                c.execute("UPDATE jobs SET status='running',updated=?,progress=5 WHERE id=?",(time.time(),row['id']))
        return dict(row) if row else None

    def cleanup(self,jid):
        if not valid_id(jid):raise ValueError('Invalid task ID')
        for kind in ('uploads','work'):remove_tree(self.state/kind/jid)

    def maintain(self,skip=None):
        now=time.time();skip=set() if skip is None else {skip} if isinstance(skip,str) else set(skip)
        with self.store.connect() as c:
            c.execute('DELETE FROM admission_events WHERE created<=?',(now-86400,))
            rows=[dict(r) for r in c.execute('SELECT * FROM jobs WHERE status!=?',('expired',))]
            remote={r[0] for r in c.execute('SELECT id FROM remote_jobs')} if c.execute("SELECT 1 FROM sqlite_master WHERE name='remote_jobs'").fetchone() else set()
        for row in rows:
            jid=row['id'];status=row['status']
            if jid in skip or not valid_id(jid):continue
            if jid in remote and status in {'created','succeeded'}:continue  # Retain interrupted transfers and unacknowledged reports.
            if status=='cancelling':
                self.docker_stop(jid)
                # Wait until any upload handler has closed its file descriptor.
                if now-row['updated']<5:continue
                self.cleanup(jid);self.store.update(jid,status='cancelled',error=None)
            elif status in {'created','uploading','uploaded'} and now-row['created']>self.cfg.get('upload_hours',DEFAULTS['upload_hours'])*3600:
                self.store.update(jid,status='cancelling')
            elif status in TERMINAL and row['expires']<=now:
                self.cleanup(jid);remove_tree(self.state/'reports'/jid)
                (self.state/'logs'/(jid+'.log')).unlink(missing_ok=True)
                self.store.update(jid,status='expired',email='',options='{}',client='',filename='',expected_bytes=0,received_bytes=0,error=None,mail_status='expired')
            elif status in {'succeeded','failed','cancelled'}:self.cleanup(jid)

    def recover(self):
        with self.store.connect() as c:
            rows=[dict(r) for r in c.execute("SELECT * FROM jobs WHERE status='running'")]
        for row in rows:
            if not valid_id(row['id']):continue
            self.docker_stop(row['id']);self.cleanup(row['id'])
            self.store.update(row['id'],status='failed',error='The service restarted while processing this task. Please submit the file again.')

    def process(self,row):
        jid=row['id']; work=self.state/'work'/jid; start=time.monotonic();proc=None
        try:
            options=Options.model_validate_json(row['options']).check_catalog(self.public)
            if not valid_id(jid) or not 0<row['expected_bytes']<=self.cfg['max_upload_bytes']:raise ValueError('Invalid task envelope')
            budget=self.memory_for(options)
            work.mkdir(mode=0o700,exist_ok=False)
            copy_uploads(self.state/'uploads'/jid,work/'incoming',row['expected_bytes'],options.custom_regions)
            cfg,mounts=job_config(options,self.catalog,work,self.cores)
            for directory in mounts.parents:
                for private in [self.state.resolve(),Path(self.cfg['secrets_file']).resolve().parent]:
                    if directory==private or directory in private.parents:raise ValueError('Reference directory overlaps private service data')
            (work/'config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
            (work/'request.json').write_text(json.dumps({'compressed':row['filename'].lower().endswith('.gz'),
                'max_expanded_bytes':self.cfg['max_expanded_bytes'],'max_reads':self.cfg['max_reads'],
                'custom_regions':options.custom_regions,'cores':self.cores,'options':options.model_dump()}))
            args=['docker','run','--rm','--name','ont-qc-job-'+jid,'--network','none','--read-only',
                '--user',f'{self.uid}:{self.gid}','--cap-drop','ALL','--security-opt','no-new-privileges',
                '--cpus',str(self.cores),'--memory',str(budget),'--memory-swap',str(budget),'--pids-limit',str(self.cfg['pids_limit']),
                '--ulimit','nofile=1024:1024','--ulimit',f"fsize={self.cfg['max_file_bytes']}:{self.cfg['max_file_bytes']}",
                '--log-driver','none','--tmpfs',f"/tmp:rw,nosuid,nodev,size={self.cfg['tmp_bytes']}",
                '--env','USER=ont-qc','--env','LOGNAME=ont-qc','--env','OMP_NUM_THREADS=1','--env','OPENBLAS_NUM_THREADS=1','--env','MKL_NUM_THREADS=1',
                '--mount',f'type=bind,src={self.pipeline},dst=/app,readonly',
                '--mount',f'type=bind,src={work / "incoming"},dst=/input,readonly',
                '--mount',f'type=bind,src={work},dst=/job',*mounts.arguments(),self.image]
            logfile=self.state/'logs'/(jid+'.log')
            with logfile.open('wb') as log:
                proc=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT,preexec_fn=lambda:limit_log_output(self.cfg['max_log_bytes']))
                while proc.poll() is None:
                    self.store.heartbeat();latest=self.store.get(jid)
                    if STOP or latest['status']=='cancelling':
                        self.docker_stop(jid);proc.wait(timeout=30)
                        self.store.update(jid,status='cancelled' if not STOP else 'failed',error=None if not STOP else 'The service stopped during analysis. Please submit again.')
                        return
                    if time.monotonic()-start>self.timeout:
                        raise TimeoutError('The analysis exceeded the runtime limit. Use a smaller ONT file.')
                    if tree_size(work)>self.cfg['max_work_bytes'] or logfile.stat().st_size>self.cfg['max_log_bytes'] or shutil.disk_usage(self.state).free<self.cfg['min_free_bytes']:
                        raise TimeoutError('The analysis exceeded the storage limit. Use a smaller ONT file.')
                    # Progress represents completed workflow steps when Snakemake reports them.
                    with logfile.open('rb') as progress_log:
                        progress_log.seek(max(0,logfile.stat().st_size-32768));tail=progress_log.read().decode('utf8','replace')
                    matches=re.findall(r'(\d+) of (\d+) steps \((\d+)%\)',tail)
                    if matches:self.store.update(jid,progress=min(95,max(5,int(matches[-1][2]))))
                    time.sleep(3)
            result_file=work/'result.json'
            if result_file.is_symlink() or not result_file.is_file() or result_file.stat().st_size>16384:raise ValueError('Missing bounded job result')
            result=json.loads(result_file.read_text())
            if proc.returncode or not result.get('ok'):
                error=result.get('error','The analysis failed. Contact the service operator.')
                self.store.update(jid,status='failed',error=str(error)[:400]);return
            latest=self.store.get(jid)
            if latest['status']=='cancelling':self.store.update(jid,status='cancelled');return
            if tree_size(self.state/'reports')+2*tree_size(work/'report')>self.cfg['max_report_total_bytes']:
                raise TimeoutError('Report storage is full. Please retry later.')
            publish_report(work/'report',self.state/'reports'/jid,self.cfg['max_report_bytes'],self.cfg['max_html_bytes'])
            self.store.update(jid,status='succeeded',progress=100,error=None,expires=time.time()+self.cfg.get('retention_days',7)*86400)
        except (TimeoutError,ValueError) as exc:
            self.docker_stop(jid);self.store.update(jid,status='failed',error=str(exc))
        except Exception:
            LOG.exception('Job %s failed',jid)
            self.docker_stop(jid);self.store.update(jid,status='failed',error='The service could not complete this analysis. Please contact the operator with your task ID.')
        finally:
            if proc and proc.poll() is None:
                self.docker_stop(jid)
                try:proc.wait(timeout=30)
                except subprocess.TimeoutExpired:proc.kill();proc.wait()
            self.cleanup(jid)

    def notify(self):
        from .notifications import notify
        notify(self.cfg,limit=1)

    def run(self):
        self.recover();children={};context=multiprocessing.get_context('spawn');maintained=0
        try:
            while not STOP:
                self.store.heartbeat()
                for jid,proc in list(children.items()):
                    if proc.is_alive():continue
                    proc.join();del children[jid]
                    row=self.store.get(jid)
                    if row and row['status'] in {'running','cancelling'}:
                        self.docker_stop(jid);self.cleanup(jid)
                        self.store.update(jid,status='cancelled' if row['status']=='cancelling' else 'failed',error='The analysis process stopped unexpectedly.')
                if time.monotonic()-maintained>10:
                    self.maintain(skip=children);self.notify();maintained=time.monotonic()
                while len(children)<self.max_jobs and not STOP:
                    row=self.claim()
                    if not row:break
                    proc=context.Process(target=process_job,args=(self.cfg,row),daemon=True)
                    proc.start();children[row['id']]=proc
                time.sleep(2)
        finally:
            for proc in children.values():
                if proc.is_alive():proc.terminate()
            deadline=time.monotonic()+35
            for jid,proc in children.items():
                proc.join(timeout=max(0,deadline-time.monotonic()))
                if proc.is_alive():self.docker_stop(jid);proc.kill();proc.join();self.cleanup(jid)


def stop_worker(signum,frame):
    global STOP
    STOP=True


def process_job(config,row):
    signal.signal(signal.SIGTERM,stop_worker);signal.signal(signal.SIGINT,stop_worker)
    Worker(config).process(row)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True);args=parser.parse_args()
    config=yaml.safe_load(Path(args.config).read_text())
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    os.umask(0o077)
    worker=Worker(config)
    with (worker.state/'worker.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        signal.signal(signal.SIGTERM,stop_worker);signal.signal(signal.SIGINT,stop_worker)
        worker.run()

if __name__=='__main__':main()
