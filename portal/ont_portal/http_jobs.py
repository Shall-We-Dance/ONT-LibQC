"""Jobs HTTP routes; registered by the application factory."""
import asyncio, hmac, json, shutil, time, uuid
from fastapi import Request, HTTPException
from starlette.concurrency import run_in_threadpool
from .schema import Submission, BatchSubmission
from .security import client_address
from .store import ACTIVE, token_for, client_key
from .short_links import link_key, key_hash

def register(api, s, store, access, example_inputs):
    owned = access.owned
    public_job = access.public_job

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


