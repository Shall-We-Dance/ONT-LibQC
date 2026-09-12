"""Private, bounded SMTP notifications after report publication.

The public container never receives SMTP credentials. A frontend-only deployment
runs this module from a small host timer, independently of the compute bridge.
"""
import argparse, fcntl, json, logging, os, smtplib, ssl, time
from pathlib import Path
import yaml
from .store import Store
from .mail_template import build_message

LOG=logging.getLogger('ont-libqc-mail')


def send(smtp,message):
    mode=smtp.get('tls_mode','starttls')
    if mode not in {'ssl','starttls'}:raise ValueError('Verified SMTP TLS is required')
    context=ssl.create_default_context()
    if mode=='ssl':client=smtplib.SMTP_SSL(smtp['host'],int(smtp.get('port',465)),timeout=15,context=context)
    else:client=smtplib.SMTP(smtp['host'],int(smtp.get('port',587)),timeout=15)
    try:
        if mode=='starttls':client.starttls(context=context)
        if smtp.get('username'):client.login(smtp['username'],smtp['password'])
        client.send_message(message)
        # A disconnect during QUIT cannot invalidate an accepted DATA response.
    finally:client.close()


def notify(config,limit=5):
    smtpfile=config.get('smtp_file')
    if not config.get('mail_enabled') or not smtpfile or not Path(smtpfile).is_file():return
    smtp=yaml.safe_load(Path(smtpfile).read_text()) or {}
    if not smtp.get('enabled'):return
    state=Path(config['state']);store=Store(state)
    with (state/'control/mail.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return
        secret=json.loads(Path(config['secrets_file']).read_text())['signing_key']
        with store.connect() as c:
            c.execute('CREATE TABLE IF NOT EXISTS notification_attempts (job_id TEXT PRIMARY KEY, first_attempt REAL NOT NULL)')
            c.execute('DELETE FROM notification_attempts WHERE job_id NOT IN (SELECT id FROM jobs WHERE status != ?)',('expired',))
            rows=[dict(r) for r in c.execute("SELECT * FROM jobs WHERE status IN ('succeeded','failed') AND expires>? AND mail_status IN ('pending','retry','sending') AND mail_attempts<5 AND mail_updated<? ORDER BY created LIMIT ?",(time.time(),time.time()-300,max(1,min(5,limit))))]
        for row in rows:
            jid=row['id']
            if row['status']=='succeeded' and not (state/'reports'/jid/'report.zip').is_file():continue
            with store.connect(immediate=True) as c:
                c.execute('INSERT OR IGNORE INTO notification_attempts VALUES (?,?)',(jid,time.time()))
                first=c.execute('SELECT first_attempt FROM notification_attempts WHERE job_id=?',(jid,)).fetchone()[0]
                # Resend keeps deduplication keys for 24 hours. Never replay an
                # ambiguous old delivery after that window (e.g. a long outage).
                if time.time()-first>23*3600:
                    c.execute("UPDATE jobs SET mail_status='failed' WHERE id=?",(jid,));continue
                changed=c.execute("UPDATE jobs SET mail_status='sending',mail_attempts=mail_attempts+1,mail_updated=? WHERE id=? AND status=? AND expires>? AND mail_status IN ('pending','retry','sending')",(time.time(),jid,row['status'],time.time())).rowcount
                if not changed:continue
            try:
                message=build_message(config,smtp,row,secret)
                send(smtp,message);status='sent'
            except Exception as exc:
                permanent=isinstance(exc,smtplib.SMTPRecipientsRefused) or (isinstance(exc,smtplib.SMTPResponseException) and 500<=exc.smtp_code<600)
                status='failed' if permanent or row['mail_attempts']+1>=5 else 'retry'
                # Never log provider responses, recipients, links or credentials.
                LOG.warning('Notification %s for task %s (%s)',status,jid,type(exc).__name__)
            with store.connect() as c:
                c.execute("UPDATE jobs SET mail_status=?,mail_updated=? WHERE id=? AND mail_status='sending'",(status,time.time(),jid))
        with store.connect() as c:
            c.execute("UPDATE jobs SET mail_status='failed' WHERE mail_status='sending' AND mail_attempts>=5 AND mail_updated<?",(time.time()-300,))


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);args=p.parse_args()
    os.umask(0o077);logging.basicConfig(level=logging.INFO)
    notify(yaml.safe_load(Path(args.config).read_text()))

if __name__=='__main__':main()
