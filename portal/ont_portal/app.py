"""Public HTTP service. It cannot run tools or access the Docker socket/reference data."""
import asyncio, html, hmac, json, mimetypes, os, re, shutil, time, uuid
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
from .store import Store, ACTIVE, token_for, client_key, report_grant, valid_report_grant, report_session_token

JOB_ID = re.compile(r'^[0-9a-f]{32}$')
APP_CSP = "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; font-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'"
REPORT_CSP = "sandbox allow-scripts allow-downloads allow-popups; default-src 'none'; script-src 'unsafe-inline' 'self'; style-src 'unsafe-inline' 'self'; img-src 'self' data: blob:; object-src 'self'; font-src data:; connect-src 'none'; frame-ancestors 'none'; form-action 'none'; base-uri 'none'"


def create_app(settings=None):
    s=settings or Settings.load(); store=Store(s.state)
    api=FastAPI(docs_url=None,redoc_url=None,openapi_url=None)
    api.state.settings=s; api.state.store=store
    api.add_middleware(TrustedHostMiddleware,allowed_hosts=s.allowed_hosts,www_redirect=False)
    @api.exception_handler(RequestValidationError)
    async def bad_input(request,exc):
        return JSONResponse({'detail':'; '.join('.'.join(str(x) for x in e['loc'][1:])+': '+e['msg'] for e in exc.errors())},422)

    def owned(request,job_id,report=False):
        if not JOB_ID.fullmatch(job_id): raise HTTPException(404,'Task not found')
        auth=request.headers.get('authorization','')
        token=auth[7:] if auth.startswith('Bearer ') else ''
        authenticated=bool(re.fullmatch(r'[0-9a-f]{64}',token) and hmac.compare_digest(token,token_for(s.secret,job_id)))
        if report and not authenticated:
            cookie=request.cookies.get('ont_report_'+job_id,'')
            authenticated=bool(re.fullmatch(r'[0-9a-f]{64}',cookie) and (
                hmac.compare_digest(cookie,report_session_token(s.secret,job_id)) or
                hmac.compare_digest(cookie,token_for(s.secret,job_id))))  # Existing browser sessions.
            if not authenticated:
                authenticated=valid_report_grant(s.secret,job_id,request.query_params.get('view',''))
        if not authenticated:
            raise HTTPException(404,'Task not found or the private link is invalid')
        row=store.get(job_id)
        if not row: raise HTTPException(404,'Task not found')
        return row

    def report_cookie(response, row):
        response.set_cookie('ont_report_'+row['id'],report_session_token(s.secret,row['id']),
                            httponly=True,secure=s.secure_cookies,samesite='lax',
                            path=f"/reports/{row['id']}/",max_age=max(0,int(row['expires']-time.time())))
        return response

    def public_job(row):
        opts=json.loads(row['options']) if row['options'] else {}
        with store.connect() as c:
            ahead=c.execute("SELECT COUNT(*) FROM jobs WHERE status IN ('queued','running') AND created<?",(row['created'],)).fetchone()[0]
        return {k:row[k] for k in ('id','status','created','updated','filename','expected_bytes','received_bytes','error','progress','mail_status','expires','bed_ready')} | {
            'sample':opts.get('sample','Expired task'),'assay':opts.get('assay'),'reference':opts.get('reference'),
            'queue_ahead':ahead if row['status']=='queued' else 0,'options':opts,
            'report_url':f"/reports/{row['id']}/samples/{opts['sample']}.html" if row['status']=='succeeded' and re.fullmatch(SAFE_NAME,opts.get('sample','')) else None,
            'download_url':f"/reports/{row['id']}/report.zip" if row['status']=='succeeded' else None,
        }

    @api.get('/healthz')
    async def health(): return {'status':'ok','worker_online':store.worker_alive()}

    @api.get('/api/catalog')
    async def catalog():
        return s.public_catalog | schema_payload() | {'max_upload_bytes':s.max_upload_bytes,'max_expanded_bytes':s.max_expanded_bytes,'max_reads':s.max_reads,'max_files':s.max_files,'pipeline_url':s.pipeline_url,
            'retention_days':s.retention_days,'mail_enabled':s.mail_enabled,'submission_key_required':bool(s.submission_key),
            'worker_online':store.worker_alive()}

    def create_jobs(request, submissions):
        amount=len(submissions)
        if not 1<=amount<=s.max_files:
            raise HTTPException(422,f'Choose between 1 and {s.max_files} files per submission')
        if len({entry.email for entry in submissions})!=1:
            raise HTTPException(422,'All files in a batch must use the same email address')
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
        subjects=[client,client_key(s.secret,'email:'+email.casefold())]
        now=time.time()
        requested=sum(entry.size+(4*1024**2 if entry.options.custom_regions else 0) for entry in submissions)
        folders=[]; result=[]; batch_id=uuid.uuid4().hex
        try:
            # Reserve the whole batch in one transaction. Concurrent submissions cannot
            # pass the same capacity check or leave a partly accepted batch behind.
            with store.connect(immediate=True) as c:
                c.execute('DELETE FROM admission_events WHERE created<=?',(now-86400,))
                placeholders=','.join('?' for _ in ACTIVE)
                for (column,value),subject in zip([('client',client),('email',email)],subjects):
                    used=c.execute('SELECT COALESCE(SUM(amount),0) FROM admission_events WHERE subject=?',(subject,)).fetchone()[0]
                    previous=c.execute(f'SELECT COUNT(*) FROM jobs WHERE {column}=? AND created>?',(value,now-86400)).fetchone()[0]
                    if max(used,previous)+amount>s.max_jobs_per_day:
                        raise HTTPException(429,'Daily submission limit reached. Please try tomorrow.')
                    pending=c.execute(f'SELECT COUNT(*) FROM jobs WHERE {column}=? AND status IN ({placeholders})',(value,*ACTIVE)).fetchone()[0]
                    if pending+amount>s.max_pending_per_client:
                        raise HTTPException(429,'Your active task limit is reached. Wait for earlier tasks to finish.')
                count,reserved=c.execute(f"SELECT COUNT(*), COALESCE(SUM(expected_bytes + CASE WHEN json_extract(options,'$.custom_regions') THEN 4194304 ELSE 0 END),0) FROM jobs WHERE status IN ({placeholders})",ACTIVE).fetchone()
                if count+amount>s.max_pending or reserved+requested>s.max_reserved_bytes:
                    raise HTTPException(503,'The queue is full. Please try again later.')
                if shutil.disk_usage(s.state).free-reserved-requested<s.min_free_bytes:
                    raise HTTPException(503,'Storage is temporarily full. Please try again later.')
                for entry in submissions:
                    jid=uuid.uuid4().hex
                    folder=s.state/'uploads'/jid; folder.mkdir(mode=0o700); folders.append(folder)
                    c.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires,mail_status,batch_id) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                        (jid,now+len(result)/1_000_000,now,'created',entry.email,entry.filename,entry.size,entry.options.model_dump_json(),client,now+s.retention_days*86400,'pending' if s.mail_enabled else 'disabled',batch_id))
                    result.append({'id':jid,'token':token_for(s.secret,jid),'status':'created'})
                for subject in subjects:
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

    @api.post('/api/jobs/{job_id}/session')
    async def session(request:Request,job_id:str):
        row=owned(request,job_id)
        return report_cookie(JSONResponse({'status':row['status']}),row)

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
        store.update(job_id,expires=0)
        return {'status':'deletion_requested'}

    def safe_report(root,file,grant=None):
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
        if path.suffix=='.html' and grant:
            # A sandboxed report has an opaque origin: SameSite cookies are not reliable
            # on its subresources. Use a short-lived, report-only grant on local links.
            def local_link(match):
                attr,quote,url=match.groups();parts=urlsplit(html.unescape(url))
                if parts.scheme or parts.netloc or not parts.path:return match.group(0)
                query=dict(parse_qsl(parts.query));query['view']=grant
                target=urlunsplit(('', '', parts.path, urlencode(query), parts.fragment))
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
                return RedirectResponse(f'/jobs/{job_id}?report=1',status_code=303)
            raise
        if row['expires']<time.time(): raise HTTPException(410,'This report has expired')
        if row['status']!='succeeded': raise HTTPException(404,'The report is not available')
        # Renew resource grants on every authorized document navigation, including
        # navigation from an older report page. The cookie carries read-only access.
        grant=report_grant(s.secret,job_id,min(int(time.time())+3590,int(row['expires'])))
        if file=='index.html':
            sample=json.loads(row['options']).get('sample','')
            if re.fullmatch(SAFE_NAME,sample):return report_cookie(RedirectResponse(f'/reports/{job_id}/samples/{sample}.html?view={grant}',status_code=307),row)
        response=safe_report(s.state/'reports'/job_id,file,grant)
        return report_cookie(response,row) if file.endswith('.html') else response

    @api.get('/example/{file:path}')
    async def example(file:str): return safe_report(s.examples,file or 'index.html')

    @api.get('/favicon.ico',include_in_schema=False)
    async def favicon():return FileResponse(s.static/'favicon.ico',media_type='image/vnd.microsoft.icon')

    @api.get('/static/{file}')
    async def static(file:str):
        if file not in {'app.js','style.css','mark.svg','favicon-32.png','apple-touch-icon.png','brand-logo.png','mail-logo.png','mail-logo.svg'}: raise HTTPException(404)
        return FileResponse(s.static/file)

    @api.get('/')
    @api.get('/jobs/{job_id}')
    async def index(): return FileResponse(s.static/'index.html')
    api.add_middleware(RequestGate, settings=s)
    api.add_middleware(RequestSecurity, settings=s, app_csp=APP_CSP, report_csp=REPORT_CSP)
    return api

app=create_app() if os.environ.get('PORTAL_CONFIG') else None
