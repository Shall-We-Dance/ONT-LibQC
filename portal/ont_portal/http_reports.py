"""Reports HTTP routes; registered by the application factory."""
import json, re, time
from pathlib import Path
from fastapi import Request, HTTPException
from fastapi.responses import HTMLResponse, RedirectResponse
from .report_download import archive_filename
from .schema import SAFE_NAME
from .short_links import KEY
from .http_access import JOB_ID

def register(api, s, store, access, example_inputs):
    owned = access.owned
    report_link = access.report_link
    safe_report = access.safe_report
    short_report_record = access.short_report_record
    task_authenticated = access.task_authenticated

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
        record=short_report_record(key)
        if record is None:
            document=request.headers.get('sec-fetch-dest')=='document' or 'text/html' in request.headers.get('accept','')
            if document:return RedirectResponse('/r/',status_code=303)
            raise HTTPException(404,'Report not found or expired')
        row,link=record
        with store.connect() as c:
            return safe_report(s.state/'reports'/row['id'],link['file'],
                               lambda file:report_link(row['id'],file,min(link['expires'],int(row['expires'])),c),
                               download_name=archive_filename(s.state/'reports'/row['id'],json.loads(row['options'] or '{}').get('sample'),row['created']) if link['file']=='report.zip' else None)


    def authorized_report_image(scope):
        """Exempt actual accessible images, never opaque links on trust alone."""
        path=scope['path']
        try:
            if path.startswith('/r/'):
                record=short_report_record(path[3:])
                if record is None:return False
                row,link=record;root=s.state/'reports'/row['id'];file=link['file']
            elif path.startswith('/reports/'):
                parts=path.split('/',3)
                if len(parts)!=4:return False
                _,_,job_id,file=parts
                if Path(file).suffix not in {'.svg','.png'}:return False
                row=owned(Request(scope),job_id,report=True)
                if row['status']!='succeeded' or row['expires']<=time.time():return False
                root=s.state/'reports'/job_id
            elif path.startswith('/example/'):
                root=s.examples;file=path[len('/example/'):]
            else:return False
            if Path(file).suffix not in {'.svg','.png'}:return False
            safe_report(root,file)
            return True
        except HTTPException:
            return False


    @api.get('/example/{file:path}')
    async def example(file:str): return safe_report(s.examples,file or 'index.html')


    return authorized_report_image
