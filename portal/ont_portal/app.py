"""Public HTTP service. It cannot run tools or access the Docker socket/reference data."""
import asyncio, html, hmac, json, mimetypes, os, posixpath, re, shutil, time, uuid
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit, parse_qsl, urlencode
from fastapi import FastAPI, Request, HTTPException
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, JSONResponse, HTMLResponse, RedirectResponse
from starlette.middleware.trustedhost import TrustedHostMiddleware
from starlette.concurrency import run_in_threadpool
from .settings import Settings
from .report_files import REPORT_FILE_SUFFIXES
from .report_entry import online_report_links
from .schema import Submission, BatchSubmission, SAFE_NAME, schema_payload
from .security import RequestGate, RequestSecurity, client_address
from .store import Store, ACTIVE, token_for, client_key, report_grant, valid_report_grant
from .short_links import KEY, link_key, key_hash, job_link, batch_link

JOB_ID = re.compile(r'^[0-9a-f]{32}$')
APP_CSP = "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; font-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'"
# Native Print / Save PDF needs allow-modals; the report remains an opaque origin
# without cookies/storage, network requests, forms, or parent-page access.
REPORT_CSP = "sandbox allow-scripts allow-downloads allow-popups allow-modals; default-src 'none'; script-src 'unsafe-inline' 'self'; style-src 'unsafe-inline' 'self'; img-src 'self' data: blob:; object-src 'self'; font-src data:; connect-src 'none'; frame-ancestors 'none'; form-action 'none'; base-uri 'none'"
# Auto-injected Web Analytics posts to the existing same-origin /cdn-cgi/rum
# endpoint. Keep it out of private task pages and sandboxed report documents.
PUBLIC_CSP = APP_CSP.replace("script-src 'self'", "script-src 'self' https://static.cloudflareinsights.com")
EXAMPLE_INDEX_CSP = PUBLIC_CSP.replace("style-src 'self'", "style-src 'unsafe-inline' 'self'")


def create_app(settings=None):
    s=settings or Settings.load(); store=Store(s.state)
    api=FastAPI(docs_url=None,redoc_url=None,openapi_url=None)
    api.state.settings=s; api.state.store=store
    manifest=s.example_inputs/'manifest.json'
    example_inputs=json.loads(manifest.read_text())['examples'] if manifest.is_file() else []
    examples_by_file={entry['filename']: entry for entry in example_inputs}

    api.add_middleware(TrustedHostMiddleware,allowed_hosts=s.allowed_hosts,www_redirect=False)
    @api.exception_handler(RequestValidationError)
    async def bad_input(request,exc):
        return JSONResponse({'detail':'; '.join('.'.join(str(x) for x in e['loc'][1:])+': '+e['msg'] for e in exc.errors())},422)

    def task_authenticated(request,job_id):
        auth=request.headers.get('authorization','')
        token=auth[7:] if auth.startswith('Bearer ') else ''
        return bool(re.fullmatch(r'[0-9a-f]{64}',token) and hmac.compare_digest(token,token_for(s.secret,job_id)))

    def owned(request,job_id,report=False):
        if not JOB_ID.fullmatch(job_id): raise HTTPException(404,'Task not found')
        authenticated=task_authenticated(request,job_id)
        if report and not authenticated:
            authenticated=valid_report_grant(s.secret,job_id,request.query_params.get('view',''))
        if not authenticated:
            raise HTTPException(404,'Task not found or the private link is invalid')
        row=store.get(job_id)
        if not row: raise HTTPException(404,'Task not found')
        return row

    def public_job(row):
        opts=json.loads(row['options']) if row['options'] else {}
        # Only authenticated task responses contain these short-lived, read-only links.
        available=row['status']=='succeeded' and row['expires']>time.time()
        access_expires=int(row['expires']) if available else None
        with store.connect() as c:
            ahead=c.execute("SELECT COUNT(*) FROM jobs WHERE status IN ('queued','running') AND created<?",(row['created'],)).fetchone()[0]
        return {k:row[k] for k in ('id','status','created','updated','filename','expected_bytes','received_bytes','error','progress','mail_status','expires','bed_ready')} | {
            'sample':opts.get('sample','Expired task'),'assay':opts.get('assay'),'reference':opts.get('reference'),
            'queue_ahead':ahead if row['status']=='queued' else 0,'options':opts,
            # Manual deletion uses expires=0; cleanup preserves this marker and
            # natural expiry timestamps even after private task details are erased.
            'deletion_reason':('manual' if row['expires']==0 else 'retention') if row['status']=='expired' else None,
            'deletion_requested':row['status']!='expired' and row['expires']==0,
            'report_url':report_link(row['id'],f"samples/{opts['sample']}.html",access_expires) if available and re.fullmatch(SAFE_NAME,opts.get('sample','')) else None,
            'download_url':report_link(row['id'],'report.zip',access_expires) if available else None,
            'report_access_expires':access_expires,
            'private_link':job_link(s.secret,row),
            'stage':row.get('stage') if row['status']=='running' else None,
        }

    @api.get('/healthz')
    async def health(): return {'status':'ok','worker_online':store.worker_alive()}

    @api.get('/api/catalog')
    async def catalog():
        return s.public_catalog | schema_payload() | {'max_upload_bytes':s.max_upload_bytes,'max_expanded_bytes':s.max_expanded_bytes,'max_reads':s.max_reads,'max_files':s.max_files,'pipeline_url':s.pipeline_url,
            'retention_days':s.retention_days,'mail_enabled':s.mail_enabled,'submission_key_required':bool(s.submission_key),
            'worker_online':store.worker_alive(), 'examples':example_inputs}

    @api.get('/api/statistics')
    async def statistics():
        return await run_in_threadpool(store.statistics)

    @api.post('/api/links/resolve')
    async def resolve_link(request: Request):
        try: data = await request.json()
        except ValueError: raise HTTPException(404, 'Private link not found')
        key = data.get('key') if isinstance(data, dict) else None
        if not isinstance(key, str) or not KEY.fullmatch(key):
            raise HTTPException(404, 'Private link not found')
        digest = key_hash(key)
        with store.connect() as c:
            rows = c.execute('SELECT * FROM jobs WHERE short_link_hash=?', (digest,)).fetchall()
            if not rows:
                batch = c.execute('SELECT batch_id FROM batch_links WHERE key_hash=?', (digest,)).fetchone()
                if batch:
                    rows = c.execute('SELECT * FROM jobs WHERE batch_id=? ORDER BY created LIMIT 64', (batch[0],)).fetchall()
                else:
                    saved = c.execute('SELECT job_ids FROM saved_links WHERE key_hash=?', (digest,)).fetchone()
                    if saved:
                        rows = [c.execute('SELECT * FROM jobs WHERE id=?', (jid,)).fetchone() for jid in json.loads(saved[0])]
                        if not all(rows): rows = []
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

    def create_jobs(request, submissions):
        amount=len(submissions)
        if not 1<=amount<=s.max_files:
            raise HTTPException(422,f'Choose between 1 and {s.max_files} files per submission')
        if len({entry.email for entry in submissions})!=1:
            raise HTTPException(422,'Use the same email address for every file in a batch, or leave it blank for every file')
        for entry in submissions:
            if entry.website: raise HTTPException(400,'Submission rejected')
            if s.submission_key and not hmac.compare_digest(entry.submission_key,s.submission_key):
                raise HTTPException(403,'The submission code is incorrect')
            if entry.size>s.max_upload_bytes:
                raise HTTPException(413,f'Each FASTQ must be at most {s.max_upload_bytes/1_000_000:g} MB. Run larger files locally using the pipeline on GitHub.')
            try: entry.options.check_catalog(s.public_catalog,require_regions_confirmation=True)
            except ValueError as e: raise HTTPException(422,str(e))
        if not store.worker_alive(): raise HTTPException(503,'The analysis worker is offline. Please try again later.')
        client=client_key(s.secret,client_address(request.scope))
        email=submissions[0].email
        quotas=[('client',client,client)]
        # Missing addresses must not create a shared quota for unrelated users.
        if email:quotas.append(('email',email,client_key(s.secret,'email:'+email.casefold())))
        mail_status='pending' if s.mail_enabled and email else 'disabled'
        now=time.time()
        requested=sum(entry.size+(4*1024**2 if entry.options.custom_regions else 0) for entry in submissions)
        folders=[]; result=[]; batch_id=uuid.uuid4().hex
        try:
            # Reserve the whole batch in one transaction. Concurrent submissions cannot
            # pass the same capacity check or leave a partly accepted batch behind.
            with store.connect(immediate=True) as c:
                c.execute('DELETE FROM admission_events WHERE created<=?',(now-86400,))
                placeholders=','.join('?' for _ in ACTIVE)
                for column,value,subject in quotas:
                    match=f'{column}=?'+(' COLLATE NOCASE' if column=='email' else '')
                    used=c.execute('SELECT COALESCE(SUM(amount),0) FROM admission_events WHERE subject=?',(subject,)).fetchone()[0]
                    previous=c.execute(f'SELECT COUNT(*) FROM jobs WHERE {match} AND created>?',(value,now-86400)).fetchone()[0]
                    if max(used,previous)+amount>s.max_jobs_per_day:
                        raise HTTPException(429,'Daily submission limit reached. Please try tomorrow.')
                    pending=c.execute(f'SELECT COUNT(*) FROM jobs WHERE {match} AND status IN ({placeholders})',(value,*ACTIVE)).fetchone()[0]
                    if pending+amount>s.max_pending_per_client:
                        raise HTTPException(429,'Your active task limit is reached. Wait for earlier tasks to finish.')
                count,reserved=c.execute(f"SELECT COUNT(*), COALESCE(SUM(expected_bytes + CASE WHEN json_extract(options,'$.custom_regions') THEN 4194304 ELSE 0 END),0) FROM jobs WHERE status IN ({placeholders})",ACTIVE).fetchone()
                if count+amount>s.max_pending or reserved+requested>s.max_reserved_bytes:
                    raise HTTPException(503,'The queue is full. Please try again later.')
                if shutil.disk_usage(s.state).free-reserved-requested<s.min_free_bytes:
                    raise HTTPException(503,'Storage is temporarily full. Please try again later.')
                c.execute('INSERT INTO batch_links VALUES (?,?)',(batch_id,key_hash(link_key(s.secret,'batch',batch_id))))
                for entry in submissions:
                    jid=uuid.uuid4().hex
                    folder=s.state/'uploads'/jid; folder.mkdir(mode=0o700); folders.append(folder)
                    c.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires,mail_status,batch_id,short_link_hash) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',
                        (jid,now+len(result)/1_000_000,now,'created',entry.email,entry.filename,entry.size,entry.options.model_dump_json(),client,now+s.retention_days*86400,mail_status,batch_id,key_hash(link_key(s.secret,'job',jid))))
                    result.append({'id':jid,'token':token_for(s.secret,jid),'status':'created','mail_status':mail_status,
                                   'private_link':'/j/#'+link_key(s.secret,'job',jid),
                                   'batch_link':'/b/#'+link_key(s.secret,'batch',batch_id)})
                for _,_,subject in quotas:
                    c.execute('INSERT INTO admission_events VALUES (?,?,?)',(subject,now,amount))
        except BaseException:
            for folder in folders: shutil.rmtree(folder)
            raise
        return result

    @api.post('/api/jobs',status_code=201)
    async def create(request:Request, submission:Submission):
        return create_jobs(request,[submission])[0]

    @api.post('/api/batches',status_code=201)
    async def create_batch(request:Request, batch:BatchSubmission):
        return {'jobs':create_jobs(request,batch.files)}

    @api.put('/api/jobs/{job_id}/upload')
    async def upload(request:Request,job_id:str):
        row=owned(request,job_id)
        if request.headers.get('content-type','').split(';')[0]!='application/octet-stream': raise HTTPException(415,'Upload the FASTQ as binary data')
        if request.headers.get('content-length')!=str(row['expected_bytes']): raise HTTPException(400,'Upload size differs from the selected file')
        if not store.transition(job_id,'created','uploading'): raise HTTPException(409,'This file has already been uploaded or is being uploaded')
        folder=s.state/'uploads'/job_id; partial=folder/'input.part'; target=folder/'input.upload'
        count=0
        try:
            async with asyncio.timeout(s.upload_timeout_seconds):
                with partial.open('xb') as f:
                    async for chunk in request.stream():
                        count+=len(chunk)
                        if count>row['expected_bytes']: raise HTTPException(413,'Upload exceeds the expected size')
                        await run_in_threadpool(f.write,chunk)
                if count!=row['expected_bytes']: raise HTTPException(400,'The upload was interrupted; please retry')
                partial.replace(target)
            store.update(job_id,received_bytes=count)
            if not store.transition(job_id,'uploading','uploaded'):
                target.unlink(missing_ok=True)
                raise HTTPException(409,'The task was cancelled during upload')
        except BaseException:
            partial.unlink(missing_ok=True)
            store.transition(job_id,'uploading','created')
            raise
        return {'status':'uploaded','received_bytes':count}

    @api.put('/api/jobs/{job_id}/regions')
    async def upload_regions(request:Request,job_id:str):
        row=owned(request,job_id)
        if row['status'] not in {'created','uploaded'}: raise HTTPException(409,'This task no longer accepts BED uploads')
        if not json.loads(row['options'])['custom_regions']: raise HTTPException(400,'Custom regions were not selected')
        if request.headers.get('content-type','').split(';')[0]!='application/octet-stream': raise HTTPException(415,'Upload the BED as binary data')
        data=await request.body()
        if not data or len(data)>4*1024**2: raise HTTPException(413,'BED must be between 1 byte and 4 MiB')
        # Validation and coordinate normalization happen inside the job sandbox.
        dest=s.state/'uploads'/job_id/'regions.upload'
        with dest.open('wb') as f: f.write(data)
        store.update(job_id,bed_ready=1)
        return {'status':'uploaded'}

    @api.post('/api/jobs/{job_id}/submit')
    async def submit(request:Request,job_id:str):
        row=owned(request,job_id)
        if json.loads(row['options'])['custom_regions'] and not row['bed_ready']:
            raise HTTPException(400,'Upload the selected BED regions before submitting')
        if not store.transition(job_id,'uploaded','queued'): raise HTTPException(409,'The file is not ready for submission')
        return public_job(store.get(job_id))

    @api.get('/api/jobs/{job_id}')
    async def status(request:Request,job_id:str): return public_job(owned(request,job_id))

    @api.post('/api/jobs/{job_id}/cancel')
    async def cancel(request:Request,job_id:str):
        row=owned(request,job_id)
        if row['status'] not in ACTIVE: raise HTTPException(409,'This task has already finished')
        # The worker owns deletion; it first stops a running container or an interrupted upload.
        if not store.transition(job_id,row['status'],'cancelling'): raise HTTPException(409,'Task state changed; please retry')
        return {'status':'cancelling'}

    @api.delete('/api/jobs/{job_id}')
    async def delete(request:Request,job_id:str):
        row=owned(request,job_id)
        if row['status'] in ACTIVE: raise HTTPException(409,'Cancel the active task before deleting its files')
        # Repeated deletion must not relabel an already expired report as manual.
        if row['status']=='expired':return {'status':'expired'}
        if row['expires']>time.time():store.update(job_id,expires=0)
        return {'status':'deletion_requested'}

    def report_link(job_id,file,expires,connection=None):
        # A separate, read-only capability for one file; never resolves to a task
        # token or grants cancellation/deletion access. Only hashes are stored.
        expires=int(expires)
        key=link_key(s.secret,'report',f'{job_id}:{expires}:{file}')
        values=(key_hash(key),job_id,file,expires)
        if connection is None:
            with store.connect() as c:
                c.execute('DELETE FROM report_links WHERE expires<=?',(time.time(),))
                c.execute('INSERT OR IGNORE INTO report_links VALUES (?,?,?,?)',values)
        else:connection.execute('INSERT OR IGNORE INTO report_links VALUES (?,?,?,?)',values)
        return '/r/'+key

    def safe_report(root,file,links=None):
        if '\\' in file or any(part in {'.','..'} or part.startswith('.') for part in file.split('/')):
            raise HTTPException(404,'File not found')
        root=root.resolve(); path=root/file
        if not path.is_relative_to(root) or not path.is_file() or any(p.is_symlink() for p in [path,*path.parents] if p!=root.parent):
            raise HTTPException(404,'File not found')
        if not path.resolve().is_relative_to(root): raise HTTPException(404,'File not found')
        allowed=REPORT_FILE_SUFFIXES | {'.zip'}
        if path.suffix not in allowed: raise HTTPException(404,'File not found')
        media={'.svg':'image/svg+xml','.html':'text/html','.css':'text/css','.js':'text/javascript'}.get(path.suffix)
        if path.suffix in {'.zip','.gz','.tsv','.csv','.bed','.json','.txt','.md','.pdf'}:
            return FileResponse(path,media_type='application/pdf' if path.suffix=='.pdf' else 'application/octet-stream',filename=path.name)
        if path.suffix=='.html' and path.stat().st_size>s.max_html_bytes: raise HTTPException(413,'Report page exceeds the display limit; download the ZIP instead')
        if path.suffix=='.html' and links:
            # Give figures, tables and downloads short read-only links too. The
            # original relative links remain untouched in the downloaded bundle.
            def local_link(match):
                attr,quote,url=match.groups();parts=urlsplit(html.unescape(url))
                if parts.scheme or parts.netloc or not parts.path:return match.group(0)
                file_path=posixpath.normpath(posixpath.join(posixpath.dirname(file),parts.path))
                if file_path.startswith(('/', '../')) or '\\' in file_path:raise HTTPException(404,'File not found')
                query=dict(parse_qsl(parts.query));query.pop('view',None)
                target=urlunsplit(('', '', links(file_path), urlencode(query), parts.fragment))
                return attr+'='+quote+html.escape(target,quote=True)+quote
            body=re.sub(r"\b(href|src)\s*=\s*(['\"])(.*?)\2",local_link,online_report_links(path.read_text(),file),flags=re.I)
            return HTMLResponse(body)
        return FileResponse(path,media_type=media)

    @api.get('/reports/{job_id}/{file:path}')
    async def report(request:Request,job_id:str,file:str):
        try:
            row=owned(request,job_id,report=True)
        except HTTPException as error:
            document=request.headers.get('sec-fetch-dest')=='document' or 'text/html' in request.headers.get('accept','')
            if error.status_code==404 and JOB_ID.fullmatch(job_id) and document:
                return RedirectResponse('/r/',status_code=303)
            raise
        if row['expires']<time.time(): raise HTTPException(410,'This report has expired')
        if row['status']!='succeeded': raise HTTPException(404,'The report is not available')
        # Preserve the old grant's lifetime; conversion cannot turn a temporary
        # read-only URL into a permanent task-management credential.
        expires=int(row['expires'])
        if not task_authenticated(request,job_id):
            expires=min(expires,int(request.query_params['view'].split('.',1)[0]))
        if file=='index.html':
            sample=json.loads(row['options']).get('sample','')
            if re.fullmatch(SAFE_NAME,sample):file=f'samples/{sample}.html'
        safe_report(s.state/'reports'/job_id,file)  # Validate before registering an alias.
        return RedirectResponse(report_link(job_id,file,expires),status_code=303)

    @api.get('/r/')
    async def unavailable_report():
        return HTMLResponse("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Report link unavailable · ONT-LibQC</title><link rel='stylesheet' href='/static/style.css'></head><body><main><h1>Report link unavailable</h1><p>This report link has expired or is incomplete. Reopen the private analysis link from your email or saved links.</p><p><a href='/'>New analysis</a></p></main></body></html>")

    @api.get('/r/{key}')
    async def short_report(request:Request,key:str):
        if not KEY.fullmatch(key):raise HTTPException(404,'Report not found')
        with store.connect() as c:
            link=c.execute('SELECT * FROM report_links WHERE key_hash=?',(key_hash(key),)).fetchone()
        row=store.get(link['job_id']) if link else None
        if not link or link['expires']<=time.time() or not row or row['expires']<=time.time() or row['status']!='succeeded':
            document=request.headers.get('sec-fetch-dest')=='document' or 'text/html' in request.headers.get('accept','')
            if document:return RedirectResponse('/r/',status_code=303)
            raise HTTPException(404,'Report not found or expired')
        with store.connect() as c:
            return safe_report(s.state/'reports'/row['id'],link['file'],
                               lambda file:report_link(row['id'],file,min(link['expires'],int(row['expires'])),c))

    @api.get('/example/{file:path}')
    async def example(file:str): return safe_report(s.examples,file or 'index.html')

    @api.get('/example-inputs/{filename}')
    async def example_input(filename:str):
        if filename not in examples_by_file or not re.fullmatch(r'HL(?:120|240)\.fastq\.gz',filename):
            raise HTTPException(404,'Example input not found')
        path=s.example_inputs/filename
        if not path.is_file() or path.is_symlink(): raise HTTPException(404,'Example input not found')
        return FileResponse(path,media_type='application/gzip',filename=filename)

    @api.get('/help')
    async def help_page(): return FileResponse(s.static/'help.html',media_type='text/html')

    @api.get('/license')
    async def license_page(): return FileResponse(s.static/'license.html',media_type='text/html')

    @api.get('/LICENSE')
    async def license_text(): return FileResponse(s.static/'LICENSE.txt',media_type='text/plain')

    @api.get('/favicon.ico',include_in_schema=False)
    async def favicon():return FileResponse(s.static/'favicon.ico',media_type='image/vnd.microsoft.icon')

    @api.api_route('/static/{file}',methods=['GET','HEAD'])
    async def static(file:str):
        if file not in {'app.js','rna-consistency.js','library-structure.js','cite.js','ont-libqc.bib','style.css','mark.svg','favicon-32.png','apple-touch-icon.png','brand-logo.png','mail-logo.png','mail-logo.svg','social-card.png'}: raise HTTPException(404)
        if file=='ont-libqc.bib':
            return FileResponse(s.static/file,media_type='application/x-bibtex',filename='ont-libqc-biorxiv.bib')
        return FileResponse(s.static/file)

    @api.get('/')
    @api.get('/jobs/{job_id}')
    @api.get('/j/')
    @api.get('/b/')
    async def index(): return FileResponse(s.static/'index.html')
    api.add_middleware(RequestGate, settings=s)
    api.add_middleware(RequestSecurity, settings=s, app_csp=APP_CSP, report_csp=REPORT_CSP,
                       public_csp=PUBLIC_CSP, example_index_csp=EXAMPLE_INDEX_CSP)
    return api

app=create_app() if os.environ.get('PORTAL_CONFIG') else None
