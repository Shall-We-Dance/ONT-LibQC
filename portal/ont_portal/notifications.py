"""Private SMTP notifications, one durable delivery per completed submission."""
import argparse, fcntl, json, logging, os, smtplib, ssl, time
from email import policy
from email.parser import BytesParser
from pathlib import Path
import yaml
from .store import Store, TERMINAL
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
            c.execute('CREATE TABLE IF NOT EXISTS notification_messages (job_id TEXT PRIMARY KEY, payload BLOB NOT NULL)')
            # Also guard legacy/imported rows: no recipient means no send or retry.
            c.execute("UPDATE jobs SET mail_status='disabled' WHERE TRIM(email)='' AND mail_status IN ('pending','retry','sending')")
            # Frozen MIME makes retries identical even if a task expires meanwhile.
            c.execute("DELETE FROM notification_attempts WHERE job_id NOT IN (SELECT COALESCE(batch_id,id) FROM jobs WHERE status != 'expired')")
            c.execute('DELETE FROM notification_messages WHERE job_id NOT IN (SELECT job_id FROM notification_attempts)')
            candidates=[dict(r) for r in c.execute("SELECT * FROM jobs WHERE status IN ('succeeded','failed','cancelled') AND TRIM(email)!='' AND expires>? AND mail_status IN ('pending','retry','sending') AND mail_attempts<5 AND mail_updated<? ORDER BY created",(time.time(),time.time()-300))]
        seen=set();attempted=0
        for row in candidates:
            key=row['batch_id'] or row['id']
            if key in seen:continue
            seen.add(key)
            with store.connect(immediate=True) as c:
                rows=[dict(r) for r in c.execute('SELECT * FROM jobs WHERE batch_id=? ORDER BY created,id',(row['batch_id'],))] if row['batch_id'] else [dict(c.execute('SELECT * FROM jobs WHERE id=?',(row['id'],)).fetchone())]
                if any(r['status'] not in TERMINAL for r in rows):continue
                live=[r for r in rows if r['mail_status']!='expired']
                if not live or any(r['mail_status'] not in {'pending','retry','sending'} or r['mail_attempts']>=5 or r['mail_updated']>=time.time()-300 for r in live):continue
                if any(r['status']=='succeeded' and r['expires']>time.time() and not (state/'reports'/r['id']/'report.zip').is_file() for r in rows):continue
                if any(not r['email'].strip() for r in live) or len({r['email'] for r in live})!=1:continue
                # A batch keeps the same retry record even if its earliest task expires.
                row=live[0];jid=key;ids=[r['id'] for r in live]
                placeholders=','.join('?' for _ in ids)
                c.execute('INSERT OR IGNORE INTO notification_attempts VALUES (?,?)',(jid,time.time()))
                first=c.execute('SELECT first_attempt FROM notification_attempts WHERE job_id=?',(jid,)).fetchone()[0]
                if time.time()-first>23*3600:
                    c.execute(f"UPDATE jobs SET mail_status='failed' WHERE id IN ({placeholders})",ids);continue
                saved=c.execute('SELECT payload FROM notification_messages WHERE job_id=?',(jid,)).fetchone()
                if saved:message=BytesParser(policy=policy.default).parsebytes(saved[0])
                else:
                    message=build_message(config,smtp,row,secret,rows=rows if row['batch_id'] else None)
                    c.execute('INSERT INTO notification_messages VALUES (?,?)',(jid,message.as_bytes()))
                c.execute(f"UPDATE jobs SET mail_status='sending',mail_attempts=mail_attempts+1,mail_updated=? WHERE id IN ({placeholders})",(time.time(),*ids))
            try:
                send(smtp,message);status='sent'
            except Exception as exc:
                permanent=isinstance(exc,smtplib.SMTPRecipientsRefused) or (isinstance(exc,smtplib.SMTPResponseException) and 500<=exc.smtp_code<600)
                status='failed' if permanent or row['mail_attempts']+1>=5 else 'retry'
                LOG.warning('Notification %s for task %s (%s)',status,jid,type(exc).__name__)
            with store.connect() as c:
                c.execute(f"UPDATE jobs SET mail_status=?,mail_updated=? WHERE id IN ({placeholders}) AND mail_status='sending'",(status,time.time(),*ids))
                if status in {'sent','failed'}:c.execute('DELETE FROM notification_messages WHERE job_id=?',(jid,))
            attempted+=1
            if attempted>=max(1,min(5,limit)):break
        with store.connect() as c:
            c.execute("UPDATE jobs SET mail_status='failed' WHERE mail_status='sending' AND mail_attempts>=5 AND mail_updated<?",(time.time()-300,))


def main():
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);args=p.parse_args()
    os.umask(0o077);logging.basicConfig(level=logging.INFO)
    notify(yaml.safe_load(Path(args.config).read_text()))

if __name__=='__main__':main()
