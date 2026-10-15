# Regression checks

Run from the repository root:

```sh
python -m pip install -r tests/requirements.txt
python -m unittest discover -s tests -v
PYTHONPATH=portal python -m unittest discover -s tests/portal -v
python -m pip install -r tests/browser-requirements.txt
python -m playwright install --with-deps chromium webkit
python tests/browser.py
```

Tests use temporary synthetic reads, mocked services, and a local static server.
Browser tests cover desktop/mobile layouts, UMI editing and movement, custom
adapters, mapper defaults and conditional availability notes, required data
acknowledgments, and submission payloads without uploading files or sending mail.
Python regressions also check classic BWA index validation, memory budgeting,
and the aligner recorded in Micro-C reports.
CI runs Python 3.11/3.13 and Chromium/WebKit; it does not run the full mapping
pipeline or require reference genomes.
