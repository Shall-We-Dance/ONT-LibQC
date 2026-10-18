"""Whole report bursts must not consume the protected submission/API budget."""
import asyncio
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import httpx
from fastapi.testclient import TestClient
from ont_portal.app import create_app
from ont_portal.security import RequestGate
from ont_portal.settings import Settings


def settings(state, **overrides):
    return Settings(state=state, public_url='https://testserver', allowed_hosts=['testserver'],
                    secret='test-only-'*8, submission_key='', public_catalog={}, min_free_bytes=0, **overrides)


class ResourceBudgetTests(unittest.IsolatedAsyncioTestCase):
    async def asyncSetUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.release = asyncio.Event()
        self.started = []
        async def app(scope, receive, send):
            self.started.append(scope['path'])
            if scope['path'].endswith('/hold'):
                await self.release.wait()
            await send({'type':'http.response.start','status':200,'headers':[]})
            await send({'type':'http.response.body','body':b'ok'})
        self.app = app

    def client(self, gate, ip='192.0.2.1'):
        return httpx.AsyncClient(transport=httpx.ASGITransport(app=gate, client=(ip, 1234)), base_url='https://testserver')

    async def started_count(self, count):
        async with asyncio.timeout(3):
            while len(self.started) < count:
                await asyncio.sleep(.005)

    async def test_report_burst_and_api_pool_are_independent(self):
        gate = RequestGate(self.app, settings(Path(self.tmp.name)))
        async with self.client(gate) as client:
            pending = [asyncio.create_task(client.get('/api/hold')) for _ in range(8)]
            try:
                await self.started_count(8)
                pending += [asyncio.create_task(client.get('/example/assets/hold')) for _ in range(24)]
                await self.started_count(32)
                self.assertEqual((await client.get('/api/catalog')).status_code, 429)
                self.assertEqual((await client.get('/static/style.css')).status_code, 200)
            finally:
                self.release.set()
                results = await asyncio.gather(*pending)
            self.assertTrue(all(r.status_code == 200 for r in results))
            self.assertEqual(gate.api.active + gate.resources.active, 0)
            self.assertEqual(gate.api.inflight, {})
            self.assertEqual(gate.resources.inflight, {})

    async def test_resource_client_and_global_limits_still_apply(self):
        config = settings(Path(self.tmp.name), max_resource_requests=4, max_resource_requests_per_client=3)
        gate = RequestGate(self.app, config)
        async with self.client(gate) as first, self.client(gate, '192.0.2.2') as second:
            pending = [asyncio.create_task(first.get('/reports/task/assets/hold')) for _ in range(3)]
            try:
                await self.started_count(3)
                self.assertEqual((await first.get('/example/figure.svg')).status_code, 429)
                pending.append(asyncio.create_task(second.get('/static/hold')))
                await self.started_count(4)
                self.assertEqual((await second.get('/example/figure.svg')).status_code, 503)
                self.assertEqual((await first.get('/api/catalog')).status_code, 200)
            finally:
                self.release.set()
                await asyncio.gather(*pending)

    async def test_rate_limits_do_not_leak_between_pools_or_methods(self):
        gate = RequestGate(self.app, settings(Path(self.tmp.name), resource_request_burst=2, resource_requests_per_minute=1, request_burst=2, requests_per_minute=1))
        async with self.client(gate) as client:
            for path in ['/example/figure.svg', '/example-inputs/HL120.fastq.gz']:
                self.assertEqual((await client.get(path)).status_code, 200)
            self.assertEqual((await client.get('/static/style.css')).status_code, 429)
            self.assertEqual((await client.get('/api/catalog')).status_code, 200)
            self.assertEqual((await client.post('/example/fake')).status_code, 200)
            self.assertEqual((await client.post('/api/batches')).status_code, 429)

    async def test_upload_concurrency_remains_protected(self):
        async def upload_app(scope, receive, send):
            self.started.append(scope['path'])
            await self.release.wait()
            await send({'type':'http.response.start','status':200,'headers':[]})
            await send({'type':'http.response.body','body':b'ok'})
        gate = RequestGate(upload_app, settings(Path(self.tmp.name)))
        async with self.client(gate) as client:
            pending = asyncio.create_task(client.put('/api/jobs/a/upload', content=b'data'))
            try:
                await self.started_count(1)
                self.assertEqual((await client.put('/api/jobs/b/upload', content=b'data')).status_code, 429)
            finally:
                self.release.set()
                self.assertEqual((await pending).status_code, 200)


class PublicExamplesTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(); self.addCleanup(temp.cleanup)
        self.settings = settings(Path(temp.name))
        self.client = TestClient(create_app(self.settings), base_url='https://testserver')
        self.addCleanup(self.client.close)

    def test_examples_are_downloadable_valid_fastq_with_matching_metadata(self):
        entries = self.client.get('/api/catalog').json()['examples']
        self.assertEqual({e['id'] for e in entries}, {'HL120','HL240'})
        for entry in entries:
            response = self.client.get(entry['url'])
            self.assertEqual(response.status_code, 200)
            self.assertIn(entry['filename'], response.headers['content-disposition'])
            self.assertEqual(hashlib.sha256(response.content).hexdigest(), entry['sha256'])
            self.assertEqual(len(response.content), entry['bytes'])
            lines = gzip.decompress(response.content).splitlines()
            self.assertEqual(len(lines), entry['reads']*4)
            self.assertEqual(len({lines[i] for i in range(0,len(lines),4)}), entry['reads'])
            for i in range(0,len(lines),4):
                self.assertRegex(lines[i], rb'^@[A-Za-z0-9_.:-]+$')
                self.assertEqual(lines[i+2], b'+')
                self.assertEqual(len(lines[i+1]), len(lines[i+3]))
        for path in ['/example-inputs/manifest.json','/example-inputs/other.fastq.gz','/example-inputs/%2e%2e%2fLICENSE','/reports/'+'d'*32+'/assets/figure.svg']:
            self.assertEqual(self.client.get(path).status_code, 404, path)

    def test_public_help_and_license_do_not_require_credentials(self):
        for path in ['/help','/license','/LICENSE']:
            response = self.client.get(path)
            self.assertEqual(response.status_code, 200)
            self.assertNotIn('set-cookie', response.headers)
        self.assertEqual(self.client.get('/LICENSE').text, (self.settings.static.parents[1]/'LICENSE').read_text())
        help_text = self.client.get('/help').text
        for link in ['/?example=HL120','/?example=HL240','/example/samples/HL120.html','/example/samples/HL240.html']:
            self.assertIn(link, help_text)

    def test_cloudflare_analytics_is_limited_to_public_portal_pages(self):
        for path in ['/','/?example=HL120','/help','/license','/example/','/example/index.html']:
            with self.subTest(path=path):
                response=self.client.get(path)
                self.assertEqual(response.status_code,200)
                csp=response.headers['content-security-policy']
                directives={p.split()[0]:p.split()[1:] for p in csp.split(';') if p.strip()}
                self.assertEqual(directives['script-src'],["'self'",'https://static.cloudflareinsights.com'])
                self.assertEqual(directives['connect-src'],["'self'"])
                self.assertNotIn('sandbox',directives)
                self.assertNotIn('set-cookie',response.headers)
                self.assertEqual(response.headers['referrer-policy'],'no-referrer')

    def test_tasks_and_full_reports_exclude_proxy_injected_analytics(self):
        for path in ['/jobs/'+'e'*32,'/example/samples/HL120.html','/example/samples/HL240.html',
                     '/reports/'+'e'*32+'/samples/private.html']:
            with self.subTest(path=path):
                response=self.client.get(path)
                csp=response.headers['content-security-policy']
                self.assertNotIn('cloudflare',csp)
                self.assertEqual(response.headers['cache-control'],'no-store, no-transform')
                if not path.startswith('/jobs/'):
                    self.assertIn("connect-src 'none'",csp)
                    self.assertIn('sandbox ',csp)
                    self.assertNotIn('allow-same-origin',csp)

if __name__ == '__main__': unittest.main()
