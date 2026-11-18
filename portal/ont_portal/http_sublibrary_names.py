"""The trusted editor is outside the sandbox used to display analysis reports."""
import hmac
import json
import re
import time

from fastapi import HTTPException, Request
from fastapi.responses import FileResponse

from .schema import SAFE_NAME
from .short_links import KEY, name_edit_token
from .sublibrary_names import scope_catalog, validate_names


def report_sample(file):
    match = re.fullmatch(r'samples/(' + SAFE_NAME + r')\.html', file)
    return match[1] if match else None


def register(api, s, store, access):
    def editable(request, key):
        record = access.short_report_record(key)
        if record is None:
            raise HTTPException(404, 'Report not found or expired.')
        row, link = record
        scheme, _, token = request.headers.get('authorization', '').partition(' ')
        if scheme != 'Bearer' or not KEY.fullmatch(token) or not hmac.compare_digest(token, name_edit_token(s.secret, row)):
            raise HTTPException(403, 'Open the report from your private analysis page to edit names.')
        sample = report_sample(link['file'])
        scopes = scope_catalog(s.state/'reports'/row['id'], sample)
        if not scopes:
            raise HTTPException(404, 'This report has no sublibraries to name.')
        return row['id'], sample, scopes

    @api.get('/sublibrary-names/{key}')
    async def editor(key: str):
        if not KEY.fullmatch(key):
            raise HTTPException(404)
        return FileResponse(s.static/'sublibrary-names.html')

    @api.get('/api/sublibrary-names/{key}')
    def get_names(request: Request, key: str):
        jid, sample, scopes = editable(request, key)
        with store.connect() as db:
            row = db.execute('SELECT revision,names FROM sublibrary_names WHERE job_id=? AND sample=?',
                             (jid, sample)).fetchone()
        return {'sample': sample, 'scopes': scopes, 'revision': row['revision'] if row else 0,
                'names': json.loads(row['names']) if row else {}}

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
