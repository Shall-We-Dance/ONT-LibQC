"""Restricted SSH RPC for a frontend host. Never executes analysis or shell commands."""
import argparse, fcntl, json, os, re, shutil, stat, sys, tempfile, time, zipfile
from pathlib import Path, PurePosixPath
import yaml
from .limits import load_limits
from .report_files import REPORT_FILE_SUFFIXES
from .store import Store, TERMINAL
from .worker import valid_id, remove_tree


def reply(data):
    sys.stdout.buffer.write(json.dumps(data,separators=(',',':')).encode()+b'\n');sys.stdout.buffer.flush()


class Frontend:
    def __init__(self,config):
        self.cfg=load_limits(config);self.state=Path(config['state']);self.store=Store(self.state)
        (self.state/'transfers').mkdir(exist_ok=True,mode=0o700)
        with self.store.connect() as c:
            c.execute('CREATE TABLE IF NOT EXISTS dispatch (job_id TEXT PRIMARY KEY)')

    def clean(self,jid,reports=False):
        remove_tree(self.state/'uploads'/jid)
        remove_tree(self.state/'transfers'/jid)
        if reports:remove_tree(self.state/'reports'/jid)

    def maintain(self):
        now=time.time()
        with self.store.connect() as c:
            c.execute('DELETE FROM admission_events WHERE created<=?',(now-86400,))
            rows=[dict(r) for r in c.execute("SELECT * FROM jobs WHERE status!='expired'")]
            dispatched={r[0] for r in c.execute('SELECT job_id FROM dispatch')}
        for row in rows:
            jid=row['id']
            if not valid_id(jid):continue
            if row['status'] in {'created','uploading','uploaded'} and now-row['created']>self.cfg['upload_hours']*3600:
                self.store.update(jid,status='cancelling')
            elif row['status']=='cancelling' and jid not in dispatched and now-row['updated']>self.cfg['upload_timeout_seconds']+5:
                self.clean(jid);self.store.update(jid,status='cancelled')
            elif row['status'] in TERMINAL and row['expires']<=now:
                self.clean(jid,reports=True)
                self.store.update(jid,status='expired',email='',options='{}',client='',filename='',expected_bytes=0,received_bytes=0,error=None,mail_status='expired')
                with self.store.connect() as c:c.execute('DELETE FROM dispatch WHERE job_id=?',(jid,))
            elif row['status'] in {'succeeded','failed','cancelled'}:self.clean(jid)

    def sync(self,request):
        updates=request.get('jobs',[])
        if not isinstance(updates,list) or len(updates)>64:raise ValueError('Too many status updates')
        with self.store.connect(immediate=True) as c:
            for update in updates:
                jid=update.get('id');status=update.get('status')
                if not valid_id(jid) or status not in {'queued','running','succeeded','failed','cancelled','expired','cancelling'}:raise ValueError('Invalid task update')
                row=c.execute('SELECT status FROM jobs WHERE id=?',(jid,)).fetchone()
                if not row:continue
                current=row['status']
                if current not in {'queued','running','cancelling'}:continue
                if status in {'failed','cancelled'}:
                    if current=='cancelling':status='cancelled'
                    error=str(update.get('error') or '')[:400] or None
                    c.execute('UPDATE jobs SET status=?,error=?,updated=?,expires=? WHERE id=?',(status,error,time.time(),time.time()+self.cfg['retention_days']*86400,jid))
                elif status=='running' and current!='cancelling':
                    progress=max(0,min(99,int(update.get('progress',0))))
                    c.execute("UPDATE jobs SET status='running',progress=?,updated=? WHERE id=?",(progress,time.time(),jid))
        if request.get('online') is True:self.store.heartbeat()
        else:
            with self.store.connect() as c:c.execute("INSERT OR REPLACE INTO service VALUES ('worker_heartbeat','0')")
        self.maintain()
        with self.store.connect() as c:
            rows=[dict(r) for r in c.execute("SELECT * FROM jobs WHERE status IN ('queued','running','cancelling') ORDER BY created LIMIT 32")]
            states={u['id']:r['status'] if (r:=c.execute('SELECT status FROM jobs WHERE id=?',(u['id'],)).fetchone()) else 'expired' for u in updates}
        return {'ok':True,'jobs':rows,'states':states}

    def fetch(self,jid,offsets=None):
        offsets=offsets or {}
        if not valid_id(jid):raise ValueError('Invalid task ID')
        with self.store.connect(immediate=True) as c:
            row=c.execute('SELECT * FROM jobs WHERE id=?',(jid,)).fetchone()
            if not row or row['status'] not in {'queued','running'}:raise ValueError('Task no longer accepts dispatch')
            row=dict(row);c.execute('INSERT OR IGNORE INTO dispatch VALUES (?)',(jid,))
        files=[];handles=[]
        try:
            for name,limit in [('input.upload',row['expected_bytes'])]+([('regions.upload',4*1024**2)] if json.loads(row['options'])['custom_regions'] else []):
                fd=os.open(self.state/'uploads'/jid/name,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK)
                f=os.fdopen(fd,'rb');handles.append(f);info=os.fstat(fd)
                if not stat.S_ISREG(info.st_mode) or not 0<info.st_size<=limit:raise ValueError('Invalid uploaded file')
                if name=='input.upload' and info.st_size!=limit:raise ValueError('Incomplete uploaded file')
                offset=offsets.get(name,0)
                if type(offset) is not int or not 0<=offset<=info.st_size:raise ValueError('Invalid resume offset')
                f.seek(offset);files.append({'name':name,'size':info.st_size,'offset':offset})
            reply({'ok':True,'job':row,'files':files})
            for f in handles:shutil.copyfileobj(f,sys.stdout.buffer,length=1024**2)
            sys.stdout.buffer.flush()
        finally:
            for f in handles:f.close()

    def transfer(self,jid,size):
        if not valid_id(jid) or type(size) is not int or not 0<size<=2*self.cfg['max_report_bytes']:raise ValueError('Invalid report transfer')
        row=self.store.get(jid)
        if not row or row['status'] not in {'queued','running','succeeded'}:raise ValueError('Task no longer accepts a report')
        if row['status']=='succeeded':return None,None
        folder=self.state/'transfers'/jid;folder.mkdir(mode=0o700,exist_ok=True)
        manifest=folder/'size.json';archive=folder/'report.zip.part'
        if manifest.exists() and json.loads(manifest.read_text())!=size:raise ValueError('Report changed during transfer')
        if not manifest.exists():manifest.write_text(json.dumps(size))
        if archive.is_symlink() or (archive.exists() and not archive.is_file()):raise ValueError('Invalid transfer file')
        return folder,archive

    def publish_offset(self,jid,size):
        folder,archive=self.transfer(jid,size)
        if folder is None:return {'ok':True,'complete':True,'offset':size}
        offset=archive.stat().st_size if archive.exists() else 0
        if offset>size:raise ValueError('Invalid saved transfer size')
        return {'ok':True,'complete':False,'offset':offset}

    def publish(self,jid,size,offset):
        folder,archive=self.transfer(jid,size)
        if folder is None:raise ValueError('Report already published; check transfer status')
        with (folder/'transfer.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            actual=archive.stat().st_size if archive.exists() else 0
            if type(offset) is not int or offset!=actual or not 0<=offset<=size:raise ValueError('Resume offset changed')
            if shutil.disk_usage(self.state).free<self.cfg['min_free_bytes']+size-offset+self.cfg['max_report_bytes']:raise ValueError('Insufficient frontend storage')
            with archive.open('ab') as f:
                left=size-offset
                while left:
                    data=sys.stdin.buffer.read(min(left,1024**2))
                    if not data:raise ValueError('Interrupted report transfer')
                    f.write(data);f.flush();left-=len(data)
                os.fsync(f.fileno())
            return self.install_report(jid,archive,size)

    def install_report(self,jid,archive,size):
        total=sum(p.stat().st_size for p in (self.state/'reports').rglob('*') if p.is_file())
        if total+size+self.cfg['max_report_bytes']>self.cfg['max_report_total_bytes']:raise ValueError('Frontend report storage is full')
        with tempfile.TemporaryDirectory(prefix='.report-transfer-',dir=self.state/'reports') as tmp:
            tmp=Path(tmp);stage=tmp/'report';stage.mkdir(mode=0o700)
            with zipfile.ZipFile(archive) as bundle:
                members=bundle.infolist();expanded=0;seen=set()
                if len(members)>2000:raise ValueError('Too many report files')
                for member in members:
                    parts=PurePosixPath(member.filename).parts
                    if not parts or parts[0]!='report' or len(parts)<2 or '\\' in member.filename or any(p.startswith('.') for p in parts):raise ValueError('Unsafe report path')
                    rel=PurePosixPath(*parts[1:]);mode=member.external_attr>>16
                    if member.is_dir() or stat.S_IFMT(mode) not in (0,stat.S_IFREG) or rel.suffix not in REPORT_FILE_SUFFIXES or rel in seen:raise ValueError('Unsupported report member')
                    if re.search(r'\.(?:bam|bai|sam|cram|fastq|fq|paf)(?:\.|$)',rel.name,re.I):raise ValueError('Intermediate in report')
                    if rel.suffix=='.html' and member.file_size>self.cfg['max_html_bytes']:raise ValueError('Report page is too large')
                    expanded+=member.file_size;seen.add(rel)
                    if expanded>self.cfg['max_report_bytes']:raise ValueError('Expanded report is too large')
                    dest=stage/rel;dest.parent.mkdir(parents=True,exist_ok=True)
                    with bundle.open(member) as src,dest.open('xb') as out:
                        copied=0
                        while data:=src.read(1024**2):
                            copied+=len(data)
                            if copied>member.file_size:raise ValueError('Report member size mismatch')
                            out.write(data)
            sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'workflow/scripts'))
            from report_bundle import validate
            validate(stage);shutil.copyfile(archive,stage/'report.zip')
            with self.store.connect(immediate=True) as c:
                current=c.execute('SELECT status FROM jobs WHERE id=?',(jid,)).fetchone()
                if not current or current['status'] not in {'queued','running'}:raise ValueError('Task was cancelled during transfer')
                dest=self.state/'reports'/jid
                # A process may have stopped after rename but before the database commit.
                # The task is still unpublished; replace that orphan with this validated copy.
                if dest.exists():remove_tree(dest)
                stage.rename(dest)
                c.execute("UPDATE jobs SET status='succeeded',progress=100,error=NULL,updated=?,expires=? WHERE id=?",(time.time(),time.time()+self.cfg['retention_days']*86400,jid))
        self.clean(jid)
        return {'ok':True}


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--maintain',action='store_true');args=p.parse_args()
    os.umask(0o077);front=Frontend(yaml.safe_load(Path(args.config).read_text()))
    if args.maintain:front.maintain();return
    # The SSH key has a forced command and no shell/forwarding/PTY privileges.
    if os.environ.get('SSH_ORIGINAL_COMMAND'):raise SystemExit('Interactive commands are not permitted')
    try:
        line=sys.stdin.buffer.readline(131073)
        if len(line)>131072 or not line.endswith(b'\n'):raise ValueError('Invalid request envelope')
        request=json.loads(line);action=request.get('action')
        if action=='sync':reply(front.sync(request))
        elif action=='fetch':front.fetch(request.get('id'),request.get('offsets'))
        elif action=='publish_offset':reply(front.publish_offset(request.get('id'),request.get('size')))
        elif action=='publish':reply(front.publish(request.get('id'),request.get('size'),request.get('offset')))
        else:raise ValueError('Unknown operation')
    except Exception:
        # Do not expose filesystem locations or submission contents to logs.
        reply({'ok':False,'error':'Frontend operation could not complete.'});raise SystemExit(1)

if __name__=='__main__':main()
