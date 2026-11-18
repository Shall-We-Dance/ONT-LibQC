"""Private task authentication and report capability/file access."""
from types import SimpleNamespace
import html, hmac, json, posixpath, re, time
from urllib.parse import urlsplit, urlunsplit, parse_qsl, urlencode
from fastapi import HTTPException
from fastapi.responses import FileResponse, HTMLResponse
from .report_files import REPORT_FILE_SUFFIXES
from .report_entry import online_report_links
from .schema import SAFE_NAME
from .store import token_for, valid_report_grant
from .short_links import KEY, link_key, key_hash, job_link
JOB_ID = re.compile(r"^[0-9a-f]{32}$")

def access_services(s, store):
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


    def safe_report(root,file,links=None,download_name=None,page_transform=None):
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
            return FileResponse(path,media_type='application/pdf' if path.suffix=='.pdf' else 'application/octet-stream',filename=download_name if file == "report.zip" and download_name else path.name)
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
            return HTMLResponse(page_transform(body) if page_transform else body)
        return FileResponse(path,media_type=media)


    def short_report_record(key):
        if not KEY.fullmatch(key):return None
        with store.connect() as c:
            link=c.execute('SELECT * FROM report_links WHERE key_hash=?',(key_hash(key),)).fetchone()
        row=store.get(link['job_id']) if link else None
        if not link or link['expires']<=time.time() or not row or row['expires']<=time.time() or row['status']!='succeeded':
            return None
        return row,link


    return SimpleNamespace(task_authenticated=task_authenticated, owned=owned, public_job=public_job, report_link=report_link, safe_report=safe_report, short_report_record=short_report_record)
