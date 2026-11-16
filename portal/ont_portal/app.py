"""Public HTTP service. It cannot run tools or access the Docker socket/reference data."""

import json
import os
from fastapi import FastAPI
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from starlette.middleware.trustedhost import TrustedHostMiddleware
from .settings import Settings
from .security import RequestGate, RequestSecurity
from .store import Store
APP_CSP = "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; font-src 'self'; object-src 'none'; base-uri 'none'; frame-ancestors 'none'; form-action 'self'"


# Native Print / Save PDF requires allow-modals; private reports remain opaque
# origins without cookies, storage, network requests, forms or parent access.
REPORT_CSP = "sandbox allow-scripts allow-downloads allow-popups allow-modals; default-src 'none'; script-src 'unsafe-inline' 'self'; style-src 'unsafe-inline' 'self'; img-src 'self' data: blob:; object-src 'self'; font-src data:; connect-src 'none'; frame-ancestors 'none'; form-action 'none'; base-uri 'none'"

PUBLIC_CSP = APP_CSP.replace("script-src 'self'", "script-src 'self' https://static.cloudflareinsights.com")

EXAMPLE_INDEX_CSP = PUBLIC_CSP.replace("style-src 'self'", "style-src 'unsafe-inline' 'self'")

from .http_access import access_services
from . import http_jobs, http_links, http_reports, http_public
def create_app(settings=None):
    s = settings or Settings.load()
    store = Store(s.state)
    api = FastAPI(docs_url=None, redoc_url=None, openapi_url=None)
    api.state.settings = s
    api.state.store = store
    manifest = s.example_inputs / "manifest.json"
    example_inputs = json.loads(manifest.read_text())["examples"] if manifest.is_file() else []
    api.add_middleware(TrustedHostMiddleware, allowed_hosts=s.allowed_hosts, www_redirect=False)

    @api.exception_handler(RequestValidationError)
    async def bad_input(request,exc):
        return JSONResponse({'detail':'; '.join('.'.join(str(x) for x in e['loc'][1:])+': '+e['msg'] for e in exc.errors())},422)


    access = access_services(s, store)
    http_links.register(api, s, store, access, example_inputs)
    http_jobs.register(api, s, store, access, example_inputs)
    report_image = http_reports.register(api, s, store, access, example_inputs)
    http_public.register(api, s, store, access, example_inputs)
    api.add_middleware(RequestGate, settings=s, report_image=report_image)
    api.add_middleware(RequestSecurity, settings=s, app_csp=APP_CSP, report_csp=REPORT_CSP,
                       public_csp=PUBLIC_CSP, example_index_csp=EXAMPLE_INDEX_CSP)
    return api

app = create_app() if os.environ.get("PORTAL_CONFIG") else None
