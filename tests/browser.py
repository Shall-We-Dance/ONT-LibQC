"""Exercise the upload editor with synthetic files and a mocked catalog."""
import functools
import http.server
import json
from pathlib import Path
import sys
import threading
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
CATALOG = dict(worker_online=True, max_files=10, max_upload_bytes=10000000,
    max_expanded_bytes=20000000, max_reads=100000, retention_days=7,
    mail_enabled=True, submission_key_required=False, pipeline_url="https://example.org",
    params=[], ranges={}, references=[dict(id="dm6", label="Fruit fly · dm6",
    assembly="dm6", ready=True, primary_chromosomes=["chr2L"],
    mappers=["bowtie2", "bwa-mem2"], mycoplasma=True, spike_in=True)])

class Handler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass

server = http.server.ThreadingHTTPServer(("127.0.0.1", 0),
    functools.partial(Handler, directory=str(ROOT / "portal")))
threading.Thread(target=server.serve_forever, daemon=True).start()
origin = f"http://127.0.0.1:{server.server_port}"
scripts = {name: (ROOT / f"tests/portal/{name}.browser.js").read_text()
           for name in ("library_structure", "data_policy")}
try:
    with sync_playwright() as p:
        for engine in sys.argv[1:] or ["chromium", "webkit"]:
            browser = getattr(p, engine).launch()
            try:
                for width in (1100, 390):
                    for scenario, script in scripts.items():
                        page = browser.new_page(viewport={"width": width, "height": 900})
                        def route_request(route):
                            url = route.request.url
                            if url == origin + "/api/catalog":
                                route.fulfill(json=CATALOG)
                            elif url.startswith(origin + "/static/"):
                                route.continue_()
                            else:
                                route.abort()
                        page.route("**/*", route_request)
                        try:
                            page.goto(origin + "/static/index.html")
                            page.evaluate(script)
                            result = page.evaluate("window.uiVerification")
                            print(json.dumps({"browser": engine, "scenario": scenario, **result}), flush=True)
                            assert result["ok"], result
                            if scenario == "library_structure":
                                labels = page.locator(".structure-settings select").all_text_contents()
                                assert any("TruSeq" in label and "Nextera" in label for label in labels), labels
                        except Exception:
                            (ROOT / "test-results").mkdir(exist_ok=True)
                            page.screenshot(path=str(ROOT / f"test-results/{scenario}-{engine}-{width}.png"), full_page=True)
                            raise
                        finally:
                            page.close()
            finally:
                browser.close()
finally:
    server.shutdown()
