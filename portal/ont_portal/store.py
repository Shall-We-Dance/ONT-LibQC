import hashlib, hmac, json, sqlite3, time
from contextlib import contextmanager
from pathlib import Path

ACTIVE = ('created','uploading','uploaded','queued','running','cancelling')
TERMINAL = ('succeeded','failed','cancelled','expired')

def token_for(secret, job_id):
    return hmac.new(secret.encode(),('job:'+job_id).encode(),hashlib.sha256).hexdigest()

def report_session_token(secret, job_id):
    """Report-only browser credential; never accepted by task management APIs."""
    return hmac.new(secret.encode(),('report-session:'+job_id).encode(),hashlib.sha256).hexdigest()

def report_grant(secret, job_id, expires):
    expires=int(expires)
    signature=hmac.new(secret.encode(),f'report:{job_id}:{expires}'.encode(),hashlib.sha256).hexdigest()
    return str(expires)+'.'+signature

def valid_report_grant(secret, job_id, grant):
    try:
        expires=int(grant.split('.',1)[0])
        return time.time()<expires<=time.time()+3600 and hmac.compare_digest(grant,report_grant(secret,job_id,expires))
    except (ValueError,TypeError):return False

def client_key(secret, ip):
    return hmac.new(secret.encode(),('client:'+ip).encode(),hashlib.sha256).hexdigest()

class Store:
    def __init__(self,state):
        self.state=Path(state); self.state.mkdir(parents=True,exist_ok=True,mode=0o700)
        for child in ('uploads','reports','control'): (self.state/child).mkdir(exist_ok=True,mode=0o700)
        self.db=self.state/'control'/'jobs.sqlite3'
        with self.connect() as c:
            c.executescript('''
            PRAGMA journal_mode=WAL;
            CREATE TABLE IF NOT EXISTS jobs (
              id TEXT PRIMARY KEY, created REAL NOT NULL, updated REAL NOT NULL,
              status TEXT NOT NULL, email TEXT NOT NULL, filename TEXT NOT NULL,
              expected_bytes INTEGER NOT NULL, received_bytes INTEGER NOT NULL DEFAULT 0,
              options TEXT NOT NULL, client TEXT NOT NULL, error TEXT,
              progress INTEGER NOT NULL DEFAULT 0, mail_status TEXT NOT NULL DEFAULT 'pending',
              mail_attempts INTEGER NOT NULL DEFAULT 0, mail_updated REAL NOT NULL DEFAULT 0,
              expires REAL NOT NULL, bed_ready INTEGER NOT NULL DEFAULT 0
            );
            CREATE INDEX IF NOT EXISTS jobs_status_created ON jobs(status,created);
            CREATE INDEX IF NOT EXISTS jobs_client_created ON jobs(client,created);
            CREATE INDEX IF NOT EXISTS jobs_email_created ON jobs(email,created);
            CREATE TABLE IF NOT EXISTS admission_events (
              subject TEXT NOT NULL, created REAL NOT NULL, amount INTEGER NOT NULL
            );
            CREATE INDEX IF NOT EXISTS admission_subject_created ON admission_events(subject, created);
            CREATE TABLE IF NOT EXISTS service (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            ''')
        # Serialize the additive migration across web, worker and mail processes.
        with self.connect(immediate=True) as c:
            if 'batch_id' not in {r['name'] for r in c.execute('PRAGMA table_info(jobs)')}:
                c.execute('ALTER TABLE jobs ADD COLUMN batch_id TEXT')
            c.execute('CREATE INDEX IF NOT EXISTS jobs_batch ON jobs(batch_id)')
        self.db.chmod(0o600)

    @contextmanager
    def connect(self, immediate=False):
        c=sqlite3.connect(self.db, timeout=15);c.row_factory=sqlite3.Row
        try:
            c.execute('PRAGMA busy_timeout=15000')
            if immediate: c.execute('BEGIN IMMEDIATE')
            yield c
            c.commit()
        except BaseException: c.rollback();raise
        finally: c.close()

    def get(self,job_id):
        with self.connect() as c:
            row=c.execute('SELECT * FROM jobs WHERE id=?',(job_id,)).fetchone()
            return dict(row) if row else None

    def update(self,job_id,**fields):
        allowed={'status','updated','error','progress','received_bytes','mail_status','mail_attempts','mail_updated','bed_ready','expires','email','options','client','filename','expected_bytes'}
        if not fields.keys() <= allowed: raise ValueError('Invalid database column')
        fields['updated']=time.time()
        with self.connect() as c:
            c.execute('UPDATE jobs SET '+', '.join(k+'=?' for k in fields)+' WHERE id=?',(*fields.values(),job_id))

    def transition(self,job_id,old,new):
        with self.connect(immediate=True) as c:
            result=c.execute('UPDATE jobs SET status=?,updated=? WHERE id=? AND status=?',(new,time.time(),job_id,old))
            return result.rowcount==1

    def heartbeat(self):
        with self.connect() as c:
            c.execute("INSERT OR REPLACE INTO service VALUES ('worker_heartbeat',?)",(str(time.time()),))

    def worker_alive(self):
        with self.connect() as c:
            row=c.execute("SELECT value FROM service WHERE key='worker_heartbeat'").fetchone()
        return bool(row and time.time()-float(row[0])<90)
