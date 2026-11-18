"""Links HTTP routes; registered by the application factory."""
import hmac, json, re
from fastapi import Request, HTTPException
from .store import token_for
from .short_links import KEY, link_key, key_hash, job_link, resolve_jobs
from .http_access import JOB_ID

def register(api, s, store, access, example_inputs):

    @api.post('/api/links/resolve')
    async def resolve_link(request: Request):
        try: data = await request.json()
        except ValueError: raise HTTPException(404, 'Private link not found')
        key = data.get('key') if isinstance(data, dict) else None
        if not isinstance(key, str) or not KEY.fullmatch(key):
            raise HTTPException(404, 'Private link not found')
        with store.connect() as c:
            rows = resolve_jobs(c, key)
        if not rows: raise HTTPException(404, 'Private link not found')
        return {'jobs': [{'id': row['id'], 'token': token_for(s.secret,row['id']),
                          'private_link': job_link(s.secret,dict(row))} for row in rows]}


    @api.post('/api/links/compact')
    async def compact_link(request: Request):
        try: data = await request.json()
        except ValueError: raise HTTPException(404, 'Private link not found')
        entries = data.get('jobs') if isinstance(data, dict) else None
        if not isinstance(entries, list) or not 1 <= len(entries) <= 64:
            raise HTTPException(404, 'Private link not found')
        ids = []
        for entry in entries:
            jid, token = (entry.get('id'), entry.get('token')) if isinstance(entry, dict) else (None, None)
            if (not isinstance(jid, str) or not JOB_ID.fullmatch(jid)
                    or not isinstance(token, str) or not re.fullmatch(r'[0-9a-f]{64}', token)
                    or not hmac.compare_digest(token, token_for(s.secret, jid))):
                raise HTTPException(404, 'Private link not found')
            if jid not in ids: ids.append(jid)
        with store.connect(immediate=True) as c:
            rows = [c.execute('SELECT batch_id FROM jobs WHERE id=?', (jid,)).fetchone() for jid in ids]
            if not all(rows):
                raise HTTPException(404, 'Private link not found')
            # Legacy batches have no registered short key. Register it without
            # disclosing it here; selection links still cover only supplied IDs.
            for batch in {row['batch_id'] for row in rows if row['batch_id']}:
                c.execute('INSERT OR IGNORE INTO batch_links VALUES (?,?)',
                          (batch, key_hash(link_key(s.secret, 'batch', batch))))
            jobs = []
            for jid in ids:
                key = link_key(s.secret, 'job', jid)
                c.execute('UPDATE jobs SET short_link_hash=? WHERE id=?', (key_hash(key), jid))
                jobs.append({'id': jid, 'token': token_for(s.secret, jid), 'private_link': '/j/#' + key})
            link = jobs[0]['private_link']
            if len(ids) > 1:
                # The capability covers exactly the authenticated selection, even
                # for subsets or links combining tasks from different batches.
                key = link_key(s.secret, 'selection', ','.join(ids))
                c.execute('INSERT OR IGNORE INTO saved_links VALUES (?,?)', (key_hash(key), json.dumps(ids)))
                link = '/b/#' + key
        return {'jobs': jobs, 'private_link': link}


