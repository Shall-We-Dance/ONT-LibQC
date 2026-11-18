"""The trusted editor is outside the sandbox used to display analysis reports."""
import hmac
import json
import time

from fastapi import HTTPException, Request
from fastapi.responses import FileResponse

from .short_links import KEY, name_edit_token, resolve_jobs
from .store import token_for
from .sublibrary_names import report_sample, scope_catalog, validate_names


def analysis_link_covers(db, secret, jid, link):
    if not isinstance(link, str) or len(link) > 8192 or not link.isascii():
        return False
    path, _, fragment = link.partition('#')
    if path in {'/j/', '/b/'} and KEY.fullmatch(fragment):
        return any(row['id'] == jid for row in resolve_jobs(db, fragment))
    # Older analysis links carry per-job credentials in the fragment.
    if path == '/jobs/' + jid:
        return hmac.compare_digest(fragment, token_for(secret, jid))
    if path == '/' and fragment.startswith('jobs='):
        entries = fragment[5:].split(',')
        expected = jid + ':' + token_for(secret, jid)
        return len(entries) <= 64 and any(hmac.compare_digest(entry, expected) for entry in entries)
    return False


def register(api, s, store, access):
    def report_record(key):
        record = access.short_report_record(key)
        if record is None:
            raise HTTPException(404, 'Report not found or expired.')
        return record

    def report_scopes(row, link):
        sample = report_sample(link['file'])
        scopes = scope_catalog(s.state/'reports'/row['id'], sample)
        if not scopes:
            raise HTTPException(404, 'This report has no sublibraries to name.')
        return row['id'], sample, scopes

    def editable(request, key):
        row, link = report_record(key)
        scheme, _, token = request.headers.get('authorization', '').partition(' ')
        if scheme != 'Bearer' or not KEY.fullmatch(token) or not hmac.compare_digest(token, name_edit_token(s.secret, row)):
            raise HTTPException(403, 'A private analysis link is required to edit names.')
        return report_scopes(row, link)

    @api.get('/sublibrary-names/{key}')
    async def editor(key: str):
        if not KEY.fullmatch(key):
            raise HTTPException(404)
        return FileResponse(s.static/'sublibrary-names.html')

    def current_names(jid, sample, scopes):
        with store.connect() as db:
            row = db.execute('SELECT revision,names FROM sublibrary_names WHERE job_id=? AND sample=?',
                             (jid, sample)).fetchone()
        return {'sample': sample, 'scopes': scopes, 'revision': row['revision'] if row else 0,
                'names': json.loads(row['names']) if row else {}}

    @api.get('/api/sublibrary-names/{key}')
    def get_names(request: Request, key: str):
        return current_names(*editable(request, key))

    @api.post('/api/sublibrary-names/{key}/unlock')
    async def unlock_names(request: Request, key: str):
        row, link = report_record(key)
        try:
            data = await request.json()
        except (ValueError, RecursionError):
            raise HTTPException(422, 'Enter a private analysis link.') from None
        private_link = data.get('link') if isinstance(data, dict) else None
        with store.connect() as db:
            if not analysis_link_covers(db, s.secret, row['id'], private_link):
                raise HTTPException(403, 'This analysis link does not give access to this report.')
        return current_names(*report_scopes(row, link)) | {'name_edit_token': name_edit_token(s.secret, row)}

    @api.put('/api/sublibrary-names/{key}')
    async def save_names(request: Request, key: str):
        jid, sample, scopes = editable(request, key)
        try:
            data = await request.json()
            if not isinstance(data, dict) or type(data.get('revision')) is not int:
                raise ValueError('Reload the editor and try again.')
            names = validate_names(data.get('names'), scopes)
        except (ValueError, RecursionError) as error:
            raise HTTPException(422, str(error) if type(error) is ValueError else 'Invalid names.') from None
        with store.connect(immediate=True) as db:
            if not db.execute("SELECT 1 FROM jobs WHERE id=? AND status='succeeded' AND expires>?",
                              (jid, time.time())).fetchone():
                raise HTTPException(410, 'This report has expired.')
            row = db.execute('SELECT revision FROM sublibrary_names WHERE job_id=? AND sample=?',
                             (jid, sample)).fetchone()
            revision = row[0] if row else 0
            if data['revision'] != revision:
                raise HTTPException(409, 'Names changed in another tab. Reload to see the saved names.')
            db.execute('INSERT OR REPLACE INTO sublibrary_names(job_id,sample,revision,names) VALUES(?,?,?,?)',
                       (jid, sample, revision + 1, json.dumps(names, ensure_ascii=False)))
        return {'revision': revision + 1}
