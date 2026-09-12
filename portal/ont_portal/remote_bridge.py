"""Pull frontend uploads into the shared compute queue; return reports over SSH."""
import argparse, fcntl, json, logging, os, shutil, signal, subprocess, threading, time
from pathlib import Path
import yaml
from .limits import load_limits
from .schema import Options
from .store import Store, ACTIVE, TERMINAL
from .catalog import public_catalog
from .worker import valid_id, remove_tree

LOG=logging.getLogger('ont-bridge');STOP=False


class Bridge:
    def __init__(self,config):
        self.cfg=load_limits(config);self.state=Path(config['state']);self.store=Store(self.state)
        self.public=public_catalog(json.loads(Path(config['catalog_file']).read_text()))
        self.ssh=config['ssh_command']
        if not isinstance(self.ssh,list) or not self.ssh or self.ssh[0]!='ssh' or not all(isinstance(x,str) for x in self.ssh):raise ValueError('Configure an SSH argument list')
        with self.store.connect() as c:
            c.execute('CREATE TABLE IF NOT EXISTS remote_jobs (id TEXT PRIMARY KEY)')

    def start(self,request):
        proc=subprocess.Popen(self.ssh,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
        timer=threading.Timer(120,proc.kill);timer.daemon=True;timer.start()
        proc.stdin.write(json.dumps(request,separators=(',',':')).encode()+b'\n');proc.stdin.flush()
        return proc,timer

    @staticmethod
    def response(proc):
        line=proc.stdout.readline(4*1024**2+1)
        if len(line)>4*1024**2 or not line.endswith(b'\n'):raise ValueError('Invalid remote response')
        data=json.loads(line)
        if data.get('ok') is not True:raise ValueError('Remote operation was declined')
        return data

    @staticmethod
    def finish(proc,timer):
        timer.cancel()
        if proc.poll() is None:proc.kill()
        proc.wait()
        for f in (proc.stdin,proc.stdout):
            if f:
                try:f.close()
                except BrokenPipeError:pass

    def rpc(self,request):
        proc,timer=self.start(request)
        try:
            proc.stdin.close();data=self.response(proc)
            if proc.wait(timeout=10):raise ValueError('Remote control operation failed')
            return data
        finally:self.finish(proc,timer)

    def discard_fetch(self,jid):
        if not valid_id(jid):raise ValueError('Invalid task ID')
        with self.store.connect(immediate=True) as c:
            if not c.execute('SELECT 1 FROM remote_jobs WHERE id=?',(jid,)).fetchone():return
            row=c.execute('SELECT status FROM jobs WHERE id=?',(jid,)).fetchone()
            if row and row['status']!='created':return
            c.execute('DELETE FROM remote_jobs WHERE id=?',(jid,));c.execute('DELETE FROM jobs WHERE id=?',(jid,))
        remove_tree(self.state/'uploads'/jid)

    def reserve(self,row):
        jid=row['id']
        if not valid_id(jid) or type(row['expected_bytes']) is not int or not 0<row['expected_bytes']<=self.cfg['max_upload_bytes']:raise ValueError('Invalid remote task')
        Options.model_validate_json(row['options']).check_catalog(self.public)
        # Reuse the public schema to exclude filename/email and configuration injection.
        from .schema import Submission
        Submission.model_validate({'email':row['email'],'filename':row['filename'],'size':row['expected_bytes'],'options':json.loads(row['options']),'consent':True})
        with self.store.connect(immediate=True) as c:
            existing=c.execute('SELECT status FROM jobs WHERE id=?',(jid,)).fetchone()
            if existing:return existing['status']=='created' and bool(c.execute('SELECT 1 FROM remote_jobs WHERE id=?',(jid,)).fetchone())
            marks=','.join('?' for _ in ACTIVE)
            count,reserved=c.execute(f'SELECT COUNT(*),COALESCE(SUM(expected_bytes+4194304),0) FROM jobs WHERE status IN ({marks})',ACTIVE).fetchone()
            required=row['expected_bytes']+4*1024**2
            if count>=self.cfg['max_pending'] or reserved+required>self.cfg['max_reserved_bytes']:return False
            if shutil.disk_usage(self.state).free<self.cfg['min_free_bytes']+reserved+required:return False
            columns=['id','created','updated','status','email','filename','expected_bytes','received_bytes','options','client','expires','bed_ready','mail_status']
            data={**row,'created':time.time(),'updated':time.time(),'status':'created','received_bytes':0,'client':'remote:'+row['client'],'mail_status':'disabled'}
            c.execute('INSERT INTO jobs ('+','.join(columns)+') VALUES ('+','.join('?' for _ in columns)+')',tuple(data[k] for k in columns))
            c.execute('INSERT INTO remote_jobs VALUES (?)',(jid,))
        return True

    def fetch(self,row):
        jid=row['id']
        if not self.reserve(row):return
        proc=timer=None
        try:
            folder=self.state/'uploads'/jid;folder.mkdir(mode=0o700,exist_ok=True)
            offsets={}
            for name in ['input.upload']+(['regions.upload'] if json.loads(row['options'])['custom_regions'] else []):
                part=folder/(name+'.part');final=folder/name
                if part.is_symlink() or final.is_symlink():raise ValueError('Invalid partial upload')
                offsets[name]=part.stat().st_size if part.exists() else final.stat().st_size if final.exists() else 0
            proc,timer=self.start({'action':'fetch','id':jid,'offsets':offsets});proc.stdin.close();response=self.response(proc)
            remote=response['job']
            if remote['id']!=jid or remote['options']!=row['options'] or remote['expected_bytes']!=row['expected_bytes']:raise ValueError('Task changed during transfer')
            expected={'input.upload':row['expected_bytes']}
            if json.loads(row['options'])['custom_regions']:expected['regions.upload']=4*1024**2
            files=response['files']
            if len(files)!=len(expected) or {f['name'] for f in files}!=set(expected):raise ValueError('Unexpected upload members')
            for item in files:
                name,size=item['name'],item['size']
                if type(size) is not int or not 0<size<=expected[name] or (name=='input.upload' and size!=expected[name]):raise ValueError('Invalid transfer size')
                offset=item['offset']
                if offset!=offsets[name]:raise ValueError('Unexpected resume offset')
                part=folder/(name+'.part');final=folder/name
                if final.exists() and not part.exists():final.rename(part)
                with part.open('ab') as dest:
                    left=size-offset
                    while left:
                        block=proc.stdout.read(min(left,1024**2))
                        if not block:raise ValueError('Upload transfer interrupted')
                        dest.write(block);dest.flush();left-=len(block)
                    os.fsync(dest.fileno())
            if proc.stdout.read(1) or proc.wait(timeout=10):raise ValueError('Invalid upload stream')
            for item in files:(folder/(item['name']+'.part')).replace(folder/item['name'])
            with self.store.connect(immediate=True) as c:
                c.execute("UPDATE jobs SET status='queued',received_bytes=expected_bytes,updated=? WHERE id=? AND status='created'",(time.time(),jid))
        finally:
            if proc:self.finish(proc,timer)

    def publish(self,row):
        jid=row['id'];path=self.state/'reports'/jid/'report.zip'
        if path.is_symlink() or not path.is_file():raise ValueError('Missing report ZIP')
        size=path.stat().st_size;resume=self.rpc({'action':'publish_offset','id':jid,'size':size})
        if resume.get('complete'):self.forget(jid);return
        offset=resume['offset']
        if type(offset) is not int or not 0<=offset<=size:raise ValueError('Invalid report resume offset')
        proc,timer=self.start({'action':'publish','id':jid,'size':size,'offset':offset})
        try:
            with path.open('rb') as source:
                source.seek(offset);shutil.copyfileobj(source,proc.stdin,length=1024**2)
            proc.stdin.close();self.response(proc)
            if proc.wait(timeout=10):raise ValueError('Report transfer failed')
            self.forget(jid)
        finally:self.finish(proc,timer)

    def forget(self,jid):
        # The frontend acknowledged its own durable copy; remove the compute copy.
        with self.store.connect(immediate=True) as c:
            c.execute('UPDATE jobs SET expires=0 WHERE id=?',(jid,));c.execute('DELETE FROM remote_jobs WHERE id=?',(jid,))

    def cycle(self):
        with self.store.connect() as c:rows=[dict(r) for r in c.execute('SELECT jobs.* FROM jobs JOIN remote_jobs USING(id)')]
        statuses=[{'id':r['id'],'status':'queued' if r['status']=='created' else 'failed' if r['status']=='expired' else r['status'],
                   'progress':r['progress'],'error':r['error'] or ('The report could not be delivered before expiry.' if r['status']=='expired' else None)} for r in rows]
        online=self.store.worker_alive()
        response=self.rpc({'action':'sync','online':online,'jobs':statuses})
        for row in rows:
            state=response['states'].get(row['id'],'expired')
            if state in {'cancelling','cancelled','expired'} and row['status'] in ACTIVE:
                if row['status']!='cancelling':self.store.transition(row['id'],row['status'],'cancelling')
            elif state in TERMINAL and row['status'] in TERMINAL:self.forget(row['id'])
            elif row['status']=='succeeded':self.publish(row)
        known={r['id'] for r in rows if r['status']!='created'}
        for row in response['jobs']:
            if row['id'] in known:continue
            if row['status']=='cancelling':
                self.rpc({'action':'sync','online':online,'jobs':[{'id':row['id'],'status':'cancelled'}]})
            elif online:
                self.fetch(row)
                break  # Keep heartbeat/status latency bounded during large transfers.

    def run(self):
        while not STOP:
            try:self.cycle()
            except Exception:LOG.warning('Frontend synchronization interrupted; pending transfers will retry.')
            for _ in range(5):
                if STOP:break
                time.sleep(1)


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);args=p.parse_args()
    os.umask(0o077);logging.basicConfig(level=logging.INFO)
    bridge=Bridge(yaml.safe_load(Path(args.config).read_text()))
    with (bridge.state/'remote-bridge.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        def stop(signum,frame):
            global STOP
            STOP=True
        signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
        bridge.run()

if __name__=='__main__':main()
