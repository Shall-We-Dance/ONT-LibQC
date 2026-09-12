"""Bound HTTP work before body parsing, with bounded in-memory rate bookkeeping."""
import time
from collections import OrderedDict
from starlette.responses import JSONResponse

def client_address(scope):
    """Use the transport peer; proxy headers are accepted only on the isolated backend."""
    import ipaddress
    value = (scope.get('client') or ('unknown', 0))[0]
    try:
        address = ipaddress.ip_address(value)
        if address.version == 6:
            if address.ipv4_mapped:return str(address.ipv4_mapped)
            return str(ipaddress.ip_network(str(address)+'/64', strict=False))
        return str(address)
    except ValueError:return 'unknown'


class RequestGate:
    def __init__(self, app, settings):
        self.app, self.settings = app, settings
        self.active = self.uploads = 0
        self.clients = OrderedDict()
        self.inflight = {}
        self.global_bucket = (time.monotonic(), settings.global_request_burst)

    @staticmethod
    def consume(bucket, now, rate, burst):
        started, tokens = bucket
        tokens = min(burst, tokens + (now-started)*rate/60)
        allowed = tokens >= 1
        return (now, tokens-1 if allowed else tokens), allowed

    async def __call__(self, scope, receive, send):
        if scope['type'] != 'http':return await self.app(scope, receive, send)
        s = self.settings; now = time.monotonic(); client = client_address(scope)
        bucket = self.clients.pop(client, (now, s.request_burst))
        self.clients[client], client_ok = self.consume(bucket, now, s.requests_per_minute, s.request_burst)
        if len(self.clients)>4096:self.clients.popitem(last=False)
        self.global_bucket, global_ok = self.consume(self.global_bucket, now, s.global_requests_per_minute, s.global_request_burst)
        upload = scope['method']=='PUT' and scope['path'].endswith(('/upload','/regions'))
        active, uploads = self.inflight.get(client, (0,0))
        status, detail = 0, ''
        if not client_ok or active>=s.max_http_requests_per_client or (upload and uploads>=s.max_uploads_per_client):
            status, detail = 429, 'Too many requests from this network. Please retry shortly.'
        elif not global_ok or self.active>=s.max_http_requests or (upload and self.uploads>=s.max_concurrent_uploads):
            status, detail = 503, 'The service is busy. Please retry shortly.'
        if status:
            return await JSONResponse({'detail':detail}, status, headers={'Retry-After':'5','Cache-Control':'no-store','Connection':'close'})(scope, receive, send)
        # Reserve without awaiting; the deployment keeps one Uvicorn event loop.
        self.active += 1; self.uploads += int(upload)
        self.inflight[client] = (active+1, uploads+int(upload))
        try:await self.app(scope, receive, send)
        finally:
            self.active -= 1; self.uploads -= int(upload)
            active, uploads = self.inflight[client]
            if active==1:self.inflight.pop(client)
            else:self.inflight[client] = (active-1, uploads-int(upload))


class RequestSecurity:
    """ASGI body bounds preserve endpoint cleanup and avoid buffering uploads."""
    def __init__(self, app, settings, app_csp, report_csp):
        self.app, self.settings = app, settings
        self.app_csp, self.report_csp = app_csp, report_csp

    async def __call__(self, scope, receive, send):
        if scope['type'] != 'http':
            return await self.app(scope, receive, send)
        import asyncio, re
        from starlette.datastructures import Headers, MutableHeaders
        from starlette.exceptions import HTTPException
        s = self.settings
        path, method = scope['path'], scope['method']
        headers = Headers(scope=scope)
        started = False
        async def secured_send(message):
            nonlocal started
            if message['type'] == 'http.response.start':
                started = True
                response = MutableHeaders(scope=message)
                response['X-Content-Type-Options'] = 'nosniff'
                response['Referrer-Policy'] = 'no-referrer'
                response['X-Frame-Options'] = 'DENY'
                response['Permissions-Policy'] = 'camera=(), microphone=(), geolocation=()'
                response['Cache-Control'] = 'no-store'
                response['Content-Security-Policy'] = self.report_csp if path.startswith(('/reports/','/example/')) else self.app_csp
                if s.secure_cookies: response['Strict-Transport-Security'] = 'max-age=31536000'
            await send(message)
        async def deny(status, detail):
            await JSONResponse({'detail': detail}, status, headers={'Connection':'close'})(scope, receive, secured_send)
        if headers.get('content-encoding','identity').lower() != 'identity':
            return await deny(415, 'HTTP content encoding is not supported; send the file as binary data')
        if method in {'POST','PUT','PATCH','DELETE'}:
            if headers.get('origin') and headers['origin'] != s.public_url.rstrip('/'):
                return await deny(403, 'Origin is not allowed')
            if headers.get('x-portal-request') != '1':
                return await deny(403, 'Missing request protection header')
        limit = (s.max_upload_bytes if re.fullmatch(r'/api/jobs/[0-9a-f]{32}/upload',path)
                 else 4*1024**2 if path.endswith('/regions')
                 else s.max_files*64*1024 if path == '/api/batches' else 64*1024)
        length = headers.get('content-length')
        if length is None and method in {'POST','PUT','PATCH'}:
            return await deny(411, 'Content-Length is required')
        try:
            declared = int(length) if length is not None else None
            if declared is not None and declared < 0: raise ValueError()
        except ValueError:
            return await deny(400, 'Invalid content length')
        if declared is not None and declared > limit:
            return await deny(413, 'Request is too large')
        consumed, complete = 0, False
        deadline = time.monotonic() + (s.upload_timeout_seconds if method == 'PUT' else 30)
        async def bounded_receive():
            nonlocal consumed, complete
            if complete or method not in {'POST','PUT','PATCH'}: return await receive()
            remaining = min(s.body_idle_seconds, deadline - time.monotonic())
            if remaining <= 0: raise HTTPException(408, 'Request body timed out')
            try: message = await asyncio.wait_for(receive(), timeout=remaining)
            except TimeoutError: raise HTTPException(408, 'Request body timed out')
            if message['type'] == 'http.request':
                consumed += len(message.get('body', b''))
                if consumed > limit or (declared is not None and consumed > declared):
                    raise HTTPException(413, 'Request is too large')
                complete = not message.get('more_body', False)
                if complete and declared is not None and consumed != declared:
                    raise HTTPException(400, 'Request length does not match the body')
            return message
        try:
            await self.app(scope, bounded_receive, secured_send)
        except HTTPException as exc:
            if started: raise
            await deny(exc.status_code, exc.detail)
