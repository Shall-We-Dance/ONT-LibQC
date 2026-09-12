#!/usr/bin/env python3
"""Build and validate a self-contained, report-only delivery directory."""
import argparse
import html
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import uuid
from urllib.parse import unquote, urlsplit

from report_brand import copy_mark, downloads_page

MARKER = '.ont-qc-work.json'
WORK_FOLDERS = {'bam', 'prepared', 'pairs', 'fragments', 'tracks', 'reference',
                'qc', 'plots', 'reports', 'summary', 'logs'}
FORBIDDEN = re.compile(r'\.(?:bam(?:\.bai)?|bai|sam|cram(?:\.crai)?|crai|fastq|fq|fa|fasta|bw|bigwig|bedgraph|bt2l?|mmi)(?:\.gz)?$', re.I)
ATTR = re.compile(r'\b(href|src)\s*=\s*([\'"])(.*?)\2', re.I)
RETIRED = ('collision-sensitivity', 'Collision hotspot', 'blacklist sensitivity')


def cleaned(value):
    """Omit retired sensitivity fields while retaining the ordinary O/E metric."""
    if isinstance(value, dict):
        return {k: cleaned(v) for k, v in value.items()
                if k not in {'collision_sensitivity', 'blacklist'}
                and not k.startswith(('collision_top', 'collision_blacklist'))}
    if isinstance(value, list):
        return [cleaned(v) for v in value]
    return value


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def check_locations(work, report, config=None):
    work, report = Path(work).resolve(), Path(report).resolve()
    if work in {Path('/'), Path.home(), Path.cwd().resolve()}:
        raise ValueError('output_dir must be a dedicated working directory')
    if work == report or work in report.parents or report in work.parents:
        raise ValueError('output_dir and report_dir must be separate, non-nested directories')
    if config:
        # Inputs are user-owned, even when a custom output_dir is configured.
        def check(value):
            if isinstance(value, dict):
                for item in value.values(): check(item)
            elif isinstance(value, list):
                for item in value: check(item)
            elif isinstance(value, str):
                path = Path(value).resolve()
                if path == work or work in path.parents or path == report or report in path.parents:
                    raise ValueError(f'Input/reference paths must be outside output directories: {value}')
        for key in ('samples', 'references', 'spike_in', 'mycoplasma_screen'):
            check(config.get(key, {}))
    return work, report


def initialize_work(work, report, config):
    work, report = check_locations(work, report, config)
    if work.exists() and any(work.iterdir()) and not (work / MARKER).is_file():
        raise ValueError(f'Refusing to adopt an existing unmarked directory: {work}. Choose a new output_dir.')
    work.mkdir(parents=True, exist_ok=True)
    write_json(work / MARKER, {'run_id': str(uuid.uuid4()), 'report': str(report)})


def cleanup_work(work, report):
    work, report = check_locations(work, report)
    marker, manifest = work / MARKER, report / 'manifest.json'
    if not marker.is_file() or not manifest.is_file():
        return False
    ownership = json.loads(marker.read_text())
    delivery = json.loads(manifest.read_text())
    if ownership.get('report') != str(report) or ownership['run_id'] != delivery.get('source_run_id'):
        return False  # A successful partial target does not authorize deleting other work.
    unknown = {p.name for p in work.iterdir()} - WORK_FOLDERS - {MARKER}
    if unknown:
        raise ValueError(f'Unexpected files in working directory; cleanup stopped: {sorted(unknown)}')
    validate(report)
    shutil.rmtree(work)
    return True


def destination(relative):
    parts = relative.parts
    if relative == Path('summary/index.html'):
        return Path('index.html')
    if parts[:2] == ('reports', 'assets'):
        return Path('assets', *parts[2:])
    if parts[0] == 'reports':
        return Path('samples', *parts[1:])
    if parts[:2] == ('plots', 'summary'):
        return Path('assets', *parts[2:])
    if parts[0] == 'qc' and len(parts) > 1 and parts[1] != 'collision':
        return Path('data', *parts[1:])
    if parts[0] == 'summary':
        return Path('data', 'summary', *parts[1:])
    raise ValueError(f'Not a report deliverable: {relative}')


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links, self.ids = [], set()
    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name in {'src', 'href'} and value:
                self.links.append(value)
            if name == 'id':
                self.ids.add(value)


def local_target(page, url, root):
    parsed = urlsplit(html.unescape(url))
    if parsed.scheme or parsed.netloc:
        if parsed.scheme in {'http', 'https', 'data', 'blob', 'mailto'}:
            return None, parsed
        raise ValueError(f'Non-portable URL in {page}: {url}')
    if parsed.path.startswith('/'):
        raise ValueError(f'Absolute link in {page}: {url}')
    target = (page.parent / unquote(parsed.path)).resolve() if parsed.path else page
    if target != root and root not in target.parents:
        raise ValueError(f'Link escapes report directory: {page}: {url}')
    return target, parsed


def validate(report):
    root = Path(report).resolve()
    if not (root / 'index.html').is_file():
        raise ValueError('Missing report/index.html')
    pages, count = {}, 0
    for path in root.rglob('*'):
        if path.is_symlink():
            raise ValueError(f'Symlink is not portable: {path}')
        if not path.is_file(): continue
        if FORBIDDEN.search(path.name) or '.matrix.' in path.name or path.name.endswith('.pyc'):
            raise ValueError(f'Intermediate file in report: {path}')
        if path.suffix in {'.html', '.svg'}:
            text = path.read_text()
            if path.suffix == '.html' and any(term.lower() in text.lower() for term in RETIRED):
                raise ValueError(f'Retired section in report: {path}')
            parser = Links(); parser.feed(text); pages[path] = parser
    for page, parser in pages.items():
        for url in parser.links:
            target, parsed = local_target(page, url, root)
            if target is None: continue
            if not target.is_file():
                raise ValueError(f'Broken link in {page.relative_to(root)}: {url}')
            if parsed.fragment and target in pages and unquote(parsed.fragment) not in pages[target].ids:
                # SVG fragment names are sometimes font glyph identifiers; check them too.
                raise ValueError(f'Missing anchor in {page.relative_to(root)}: {url}')
            count += 1
    return {'sample_reports': len(list((root / 'samples').glob('*.html'))),
            'html_pages': sum(p.suffix == '.html' for p in pages), 'local_links': count}


def build(work, report, config_path):
    import yaml
    work, report = check_locations(work, report)
    config = yaml.safe_load(Path(config_path).read_text())
    samples = list(config['samples'])
    report.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.report-build-', dir=report.parent) as staging:
        stage = Path(staging)
        queue = [work / 'summary/index.html'] + [work / 'reports' / f'{s}.html' for s in samples]
        for path in (work / 'summary').glob('*'):
            if path.suffix in {'.json', '.tsv'}: queue.append(path)
        for sample in samples:
            queue += [work / 'qc/metrics' / f'{sample}.json',
                      work / 'qc/adapters' / f'{sample}.summary.json',
                      work / 'qc/mycoplasma' / f'{sample}.json',
                      work / 'qc/spike_in' / f'{sample}.json']
            for category, pattern in [('gene_profile', f'{sample}*.profile.tsv'),
                                      ('custom_regions', f'{sample}.json'),
                                      ('ligation', f'{sample}.json')]:
                queue.extend((work / 'qc' / category).glob(pattern))
            snapshot = work / 'qc/spike_in/database_manifests' / f'{sample}.json'
            if snapshot.is_file(): queue.append(snapshot)
        done = set()
        while queue:
            source = queue.pop().resolve()
            if source in done: continue
            if work not in source.parents or not source.is_file():
                raise ValueError(f'Missing or external report dependency: {source}')
            relative = source.relative_to(work)
            if FORBIDDEN.search(source.name) or '.matrix.' in source.name:
                raise ValueError(f'Report links to an intermediate: {relative}')
            target = stage / destination(relative)
            target.parent.mkdir(parents=True, exist_ok=True)
            if source.suffix in {'.html', '.svg'}:
                def rewrite(match):
                    url = match[3]
                    linked, parsed = local_target(source, url, work)
                    if linked is None: return match[0]
                    queue.append(linked)
                    dest = stage / destination(linked.relative_to(work))
                    new = os.path.relpath(dest, target.parent) if parsed.path else ''
                    if parsed.query: new += '?' + parsed.query
                    if parsed.fragment: new += '#' + parsed.fragment
                    return match[1] + '=' + match[2] + html.escape(new, quote=True) + match[2]
                target.write_text(ATTR.sub(rewrite, source.read_text()))
            elif source.suffix == '.json':
                write_json(target, cleaned(json.loads(source.read_text())))
            else:
                shutil.copy2(source, target)
            done.add(source)
        write_json(stage / 'data/resolved_config.json', cleaned(config))
        shutil.copy2(Path(__file__).with_name('metric_definitions.json'), stage / 'data/metric_definitions.json')
        copy_mark(stage / "assets")
        # List every retained table/JSON/download, including evidence behind plots.
        downloads = sorted(p for p in stage.rglob('*') if p.is_file() and p.suffix not in {'.html', '.svg', '.png'})
        items = ''.join(f'<li><a href="{html.escape(str(p.relative_to(stage)), quote=True)}">{html.escape(str(p.relative_to(stage)))}</a></li>' for p in downloads)
        (stage / 'downloads.html').write_text(downloads_page(items))
        index = stage / 'index.html'
        index.write_text(index.read_text().replace('</body>', '<p><a href="downloads.html">All report data and provenance</a></p></body>'))
        (stage / 'README.md').write_text('<p><img src="assets/ont-libqc-mark.svg" width="40" height="40" alt="ONT-LibQC"></p>\n\n# ONT-LibQC report\n\nOpen `index.html` directly in a browser. Copy this entire folder to share it.\n\n- `samples/`: all sample reports, including sublibraries and local BED upload.\n- `assets/`: figures and figure source tables.\n- `data/`: report metrics, downloadable evidence and provenance.\n- `downloads.html`: complete download directory.\n\nNo BAM, SAM, FASTQ, coverage track or matrix is needed to view or interact with the report.\n')
        audit = validate(stage)
        if audit['sample_reports'] != len(samples): raise ValueError('Sample report count mismatch')
        marker = work / MARKER
        manifest = {'schema_version': 1, 'samples': samples, 'validation': audit,
                    'source_run_id': json.loads(marker.read_text())['run_id'] if marker.is_file() else None}
        write_json(stage / 'manifest.json', manifest)
        validate(stage)
        # Validate first; replacing a previously generated report is then safe.
        if report.exists():
            if not (report / 'manifest.json').is_file():
                raise ValueError(f'Refusing to replace an unrecognized report directory: {report}')
            shutil.rmtree(report)
        shutil.copytree(stage, report)
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work')
    parser.add_argument('--report', default='report')
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--validate', action='store_true')
    args = parser.parse_args()
    if args.validate:
        result = validate(args.report)
    else:
        if not args.work: parser.error('--work is required when building')
        result = build(args.work, args.report, args.config)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
