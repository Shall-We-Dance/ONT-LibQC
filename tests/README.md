# Regression checks

Run from the repository root:

```sh
python -m pip install -r tests/requirements.txt
python -m unittest discover -s tests -v
PYTHONPATH=portal python -m unittest discover -s tests/portal -v
python -m pip install -r tests/browser-requirements.txt
python -m playwright install --with-deps chromium webkit
python tests/browser.py
PYTHONPATH=portal python tests/portal/check_public_pages_browser.py
PYTHONPATH=portal python tests/portal/check_private_reports_browser.py
PYTHONPATH=portal python tests/portal/check_optional_email_browser.py
```

Tests use temporary synthetic reads, the two published research examples, mocked
services, and local test servers.
Browser tests cover desktop/mobile layouts, UMI editing and movement, custom
adapters, mapper defaults and conditional availability notes, required data
acknowledgments, and submission payloads. Route tests exercise example loading,
navigation, every report image, optional email, private links and ZIP downloads.
Uploads stay on temporary local servers; completion is simulated and mail is mocked.
Python regressions also check classic BWA index validation, memory budgeting,
and the aligner recorded in Micro-C reports, separate resource/API rate budgets,
report-grant isolation and the public-only analytics policy.
CI runs Python 3.11/3.13 and Chromium/WebKit; it does not run the full mapping
pipeline or require reference genomes.
