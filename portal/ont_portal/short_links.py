"""Compact private capabilities; fragments never enter HTTP paths or referrers."""
import base64
import hashlib
import hmac
import json
import re

KEY = re.compile(r'^[A-Za-z0-9_-]{22}$')

def link_key(secret, kind, identity):
    digest = hmac.new(secret.encode(), ('short:' + kind + ':' + identity).encode(), hashlib.sha256).digest()
    return base64.urlsafe_b64encode(digest[:16]).decode().rstrip('=')

def key_hash(key):
    return hashlib.sha256(key.encode()).hexdigest()

def resolve_jobs(db, key):
    """Resolve a compact analysis link to exactly the jobs it covers."""
    digest = key_hash(key)
    rows = db.execute('SELECT id,short_link_hash FROM jobs WHERE short_link_hash=?', (digest,)).fetchall()
    if rows:
        return rows
    batch = db.execute('SELECT batch_id FROM batch_links WHERE key_hash=?', (digest,)).fetchone()
    if batch:
        return db.execute('SELECT id,short_link_hash FROM jobs WHERE batch_id=? ORDER BY created LIMIT 64',
                          (batch[0],)).fetchall()
    selection = db.execute('SELECT job_ids FROM saved_links WHERE key_hash=?', (digest,)).fetchone()
    if not selection:
        return []
    rows = [db.execute('SELECT id,short_link_hash FROM jobs WHERE id=?', (jid,)).fetchone()
            for jid in json.loads(selection[0])]
    return rows if all(rows) else []

def job_link(secret, row):
    if row.get('short_link_hash'):
        return '/j/#' + link_key(secret, 'job', row['id'])
    from .store import token_for
    return '/jobs/' + row['id'] + '#' + token_for(secret, row['id'])

def batch_link(secret, rows):
    batch = rows[0].get('batch_id') if rows else None
    if batch and all(row.get('short_link_hash') and row.get('batch_id') == batch for row in rows):
        return '/b/#' + link_key(secret, 'batch', batch)
    from .store import token_for
    return '/#jobs=' + ','.join(row['id'] + ':' + token_for(secret, row['id']) for row in rows)
