"""Public HTTP routes; registered by the application factory."""
import re
from fastapi import HTTPException
from fastapi.responses import FileResponse
from starlette.concurrency import run_in_threadpool
from .schema import schema_payload

def register(api, s, store, access, example_inputs):
    examples_by_file = {entry["filename"]: entry for entry in example_inputs}

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
        if file not in {'app.js','sublibrary-names.js','task-views.js','server-status.js','library-structure.js','cite.js','ont-libqc.bib','style.css','mark.svg','favicon-32.png','apple-touch-icon.png','brand-logo.png','mail-logo.png','mail-logo.svg','social-card.png'}: raise HTTPException(404)
        if file=='ont-libqc.bib':
            return FileResponse(s.static/file,media_type='application/x-bibtex',filename='ont-libqc-biorxiv.bib')
        return FileResponse(s.static/file)


    @api.get('/')
    @api.get('/jobs/{job_id}')
    @api.get('/j/')
    @api.get('/b/')
    async def index(): return FileResponse(s.static/'index.html')


