#!/usr/bin/env python3
"""Validate, preview and transactionally refresh retained report presentation.

Deployment inventory lives in an external manifest, never in the repository.
This command neither maps reads nor changes retained metric/analysis records.
"""
import argparse
import datetime
import fcntl
import hashlib
import html
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import uuid
import zipfile

from report_bundle import validate
from report_export import add_pdf_export
from report_metrics import COLLISION_ROWS, METRIC_META, REPORT_CSS, REPORT_TOOLTIP_SCRIPT, scoped_metric_rows
from report_profiles import render_profiles_card
from report_provenance import RELEASE_VERSION
from species_composition import render as render_species

VOID = {'area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link', 'meta', 'param', 'source', 'track', 'wbr'}


class Elements(HTMLParser):
    """Locate balanced source ranges without rewriting unrelated HTML or scripts."""
    def __init__(self, text):
        super().__init__(convert_charrefs=False)
        self.text, self.stack, self.elements = text, [], []
        self.lines = [0]
        for line in text.splitlines(True): self.lines.append(self.lines[-1]+len(line))
        self.feed(text)
        self.close()

    def position(self):
        line, col = self.getpos()
        return self.lines[line-1]+col

    def handle_starttag(self, tag, attrs):
        start = self.position()
        node = dict(tag=tag, attrs=dict(attrs), start=start, end=start+len(self.get_starttag_text()),
                    parent=self.stack[-1] if self.stack else None)
        self.elements.append(node)
        if tag not in VOID: self.stack.append(node)

    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs)
        if tag not in VOID: self.stack.pop()

    def handle_endtag(self, tag):
        # Only matched, explicit closing tags are eligible replacement boundaries.
        for i in range(len(self.stack)-1, -1, -1):
            if self.stack[i]['tag'] == tag:
                node = self.stack[i]
                node['end'] = self.text.index('>', self.position())+1
                node['closed'] = True
                del self.stack[i:]
                return

    def within(self, node, parent):
        return parent['start'] < node['start'] and node['end'] < parent['end']

    def fragment(self, node):
        if not node.get('closed'): raise ValueError('Unclosed presentation element')
        return self.text[node['start']:node['end']]


def replace_ranges(text, replacements):
    end = len(text)
    for node, replacement in sorted(replacements, key=lambda pair: pair[0]['start'], reverse=True):
        if not node.get('closed') or node['end'] > end: raise ValueError('Overlapping presentation replacements')
        text = text[:node['start']]+replacement+text[node['end']:]
        end = node['start']
    return text


def refresh_profile_layout(root, sample, text):
    parsed = Elements(text)
    sections = [n for n in parsed.elements if n['attrs'].get('id') == 'profiles']
    if len(sections) != 1: raise ValueError('Expected one profiles section: '+sample)
    section = sections[0]
    metric = json.loads((root/'data/metrics'/f'{sample}.json').read_text())
    scopes = [dict(metric, scope_id='total', scope_name='Total library')]+metric.get('sub_libraries', [])
    ids = [scope.get('scope_id') for scope in scopes]
    if any(not isinstance(sid, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sid) for sid in ids) or len(set(ids)) != len(ids):
        raise ValueError('Invalid retained scope identifiers')
    rna_path = root/'data/rna'/sample/'summary.json'
    rna = json.loads(rna_path.read_text()) if metric.get('assay') in {'rnaseq', 'clip'} and rna_path.is_file() else None
    scope_html = {}
    for node in parsed.elements:
        if 'profile-scope' not in node['attrs'].get('class', '').split() or not parsed.within(node, section): continue
        ancestor = node['parent']
        sid = 'total' if len(scopes) == 1 else None
        while ancestor and ancestor is not section:
            if ancestor['attrs'].get('data-scope'): sid = ancestor['attrs']['data-scope']; break
            if 'pooled-profile' in ancestor['attrs'].get('class', '').split(): sid = 'total'; break
            ancestor = ancestor['parent']
        if sid is None and node['parent'] is section:
            # Legacy reports put every scope directly under the profiles heading.
            images = [img['attrs'].get('src', '') for img in parsed.elements
                      if img['tag'] == 'img' and parsed.within(img, node)]
            found = [scope for scope in ids if any(Path(src).name.startswith(f'{sample}.{scope}.') for src in images)]
            if len(found) == 1:
                sid = found[0]
        if sid is None or sid in scope_html: raise ValueError('Cannot identify historical profile scope: '+sample)
        scope_html[sid] = parsed.fragment(node)
    if set(scope_html) != {s['scope_id'] for s in scopes}: raise ValueError('Retained profile scopes do not match metrics: '+sample)
    if len(scopes) == 1 and rna is None:
        # Single-library figures are already retained verbatim. Older reports
        # embed them inside Library quality; do not add a second card frame.
        new = parsed.fragment(section)
        ancestor = section['parent']
        while ancestor and ancestor['attrs'].get('id') != 'library-shape':
            ancestor = ancestor['parent']
        if ancestor and 'report-section' in section['attrs'].get('class', '').split():
            new = ("<div class='library-profiles' id='profiles'><h3>Size and genomic profiles</h3>"
                   + scope_html['total'] + '</div>')
    else:
        new = render_profiles_card(sample, scopes, [scope_html[s['scope_id']] for s in scopes],
                                   root/'assets', root/'data/gene_profile', '../assets/',
                                   rna_summary=rna, rna_prefix='../data/rna/'+sample)
    replacements = [(section, new)]
    if rna is not None:
        parent = section['parent']
        if parent and parent['attrs'].get('id') == 'library-shape':
            # Promote the old embedded profiles to their own report section.
            replacements = [(parent, parsed.fragment(parent).replace(parsed.fragment(section), '', 1) + new)]
        for node in parsed.elements:
            parent = node['parent']
            if parent and parent['attrs'].get('id') == 'rna-qc':
                if 'rna-comparison' in node['attrs'].get('class', '').split():
                    replacements.append((node, ''))
                elif node['tag'] == 'style' and '.rna-comparison{' in parsed.fragment(node):
                    replacements.append((node, ''))
    return replace_ranges(text, replacements)


def refresh_collision_rows(root, sample, text):
    """Update only the collision rows using the retained per-scope counts."""
    label = METRIC_META["collision_oe"][0]
    if label not in text:
        return text
    metric = json.loads((root/'data/metrics'/f'{sample}.json').read_text())
    scopes = [dict(metric, scope_name='Total library'), *metric.get('sub_libraries', [])]
    parsed = Elements(text)
    labels = [METRIC_META[k][0] for k in COLLISION_ROWS]
    rows = [n for n in parsed.elements if n['tag'] == 'tr' and
            any('>' + name + ' <b>' in parsed.fragment(n) for name in labels)]
    if not rows:
        return text
    parent = rows[0]['parent']
    headers = [html.unescape(parsed.fragment(n)[4:-5]) for n in parsed.elements
               if n['tag'] == 'th' and n['parent']['parent'] is parent]
    if (headers[1:-1] != [s.get('scope_name', s.get('scope_id')) for s in scopes]
            or any(n['parent'] is not parent for n in rows)
            or '>' + label + ' <b>' not in parsed.fragment(rows[0])):
        raise ValueError('Collision table columns do not match the retained metrics')
    rendered = scoped_metric_rows(scopes, COLLISION_ROWS, {})
    return replace_ranges(text, [(rows[0], rendered)] + [(n, '') for n in rows[1:]])


def refresh_species_layout(root, sample, text):
    data = root / 'data/species_composition' / f'{sample}.json'
    if not data.is_file():
        return text
    parsed = Elements(text)
    sections = {n['attrs'].get('id'): n for n in parsed.elements if n['tag'] == 'section'}
    links = {n['attrs'].get('href'): n for n in parsed.elements
             if n['tag'] == 'a' and n['parent'] and n['parent']['tag'] == 'nav'}
    section = render_species(json.loads(data.read_text()))
    mixture = sections.get('sublibrary-composition')
    anchor = mixture or sections['sequencing']
    original = parsed.fragment(anchor)
    nav_anchor = links['#sequencing']
    changes = [(anchor, original + section if mixture else section + original),
               (nav_anchor, "<a href='#species-composition'>Species</a>" + parsed.fragment(nav_anchor))]
    for previous in (sections.get('species-composition'), links.get('#species-composition')):
        if previous:
            changes.append((previous, ''))
    return replace_ranges(text, changes)


def refresh_tables(root, sample, text):
    text = refresh_collision_rows(root, sample, text)
    text = refresh_species_layout(root, sample, text)
    return refresh_presentation(text)


def refresh_page(root, sample):
    page = root/'samples'/f'{sample}.html'
    text = refresh_profile_layout(root, sample, page.read_text())
    page.write_text(refresh_tables(root, sample, text))


def refresh_presentation(text):
    """Refresh owned resources without replacing independent component styles.

    Legacy reports have an unnamed base style. Export/print styles also mention
    hero and tooltip selectors, so content alone cannot identify ownership.
    Reinstall the export component to repair reports affected by older refreshes.
    """
    parsed = Elements(text)
    base = [node for node in parsed.elements if node['tag'] == 'style' and
            (node['attrs'].get('id') == 'report-base-style' or
             (not node['attrs'].get('id') and '.report-hero' in parsed.fragment(node)
              and '.metric-tip' in parsed.fragment(node)))]
    if not base:
        raise ValueError('Expected one base report style')
    # Older refreshes replaced the named export style with a second, unnamed
    # copy of the base CSS. Only collapse byte-identical duplicates, never
    # discard a different unowned component based on selector overlap.
    contents = [parsed.fragment(node).split('>', 1)[1].rsplit('</style>', 1)[0] for node in base]
    if len(set(contents)) != 1:
        raise ValueError('Ambiguous base report styles')
    replacements = [(base[0], '<style id="report-base-style">'+REPORT_CSS+'</style>')]
    replacements.extend((node, '') for node in base[1:])
    for node in parsed.elements:
        parent = node['parent']
        if (node['tag'] == 'span' and parent and 'hero-meta' in parent['attrs'].get('class', '').split()):
            content = parsed.fragment(node).split('>', 1)[1].rsplit('</span>', 1)[0]
            match = re.fullmatch(r'v[0-9]+\.[0-9]+\.[0-9]+( · [^<>]+)', content)
            if match:
                replacements.append((node, '<span class="report-version" title="Report release and original analysis date">'
                                      + RELEASE_VERSION + match[1] + '</span>'))
        if node['tag'] == 'script' and "tip.id = 'metric-tooltip'" in parsed.fragment(node):
            replacements.append((node, REPORT_TOOLTIP_SCRIPT))
    return add_pdf_export(replace_ranges(text, replacements), refresh=True)



def remove_retired_cohort(root):
    page = root/'index.html'
    text = page.read_text(); parsed = Elements(text)
    nodes = [n for n in parsed.elements if n['attrs'].get('id') == 'rna-cohort']
    if nodes: page.write_text(replace_ranges(text, [(n, '') for n in nodes]))


def hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob('*')) if p.is_file()}


def build_archive(root):
    target = root/'report.zip'
    if target.exists():
        with zipfile.ZipFile(target) as old:
            if any(not name.startswith('report/') for name in old.namelist()):
                raise ValueError('Unrecognized report ZIP layout')
    with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED, compresslevel=3) as archive:
        for path in sorted(root.rglob('*')):
            if path.is_file() and path != target: archive.write(path, 'report/'+str(path.relative_to(root)))
    with zipfile.ZipFile(target) as archive:
        if archive.testzip(): raise ValueError('Report ZIP validation failed')


def refresh(root, apply=False, expected_manifest_sha256=None):
    original = Path(root)
    if original.is_symlink(): raise ValueError('Report root must not be a symlink')
    root = original.resolve()
    manifest = root/'manifest.json'
    if expected_manifest_sha256 and hashlib.sha256(manifest.read_bytes()).hexdigest() != expected_manifest_sha256:
        raise ValueError('Report manifest has changed')
    data = json.loads(manifest.read_text())
    samples = data.get('samples')
    if data.get('schema_version') != 1 or not isinstance(samples, list) or not samples or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', s) for s in samples):
        raise ValueError('Unrecognized report manifest')
    validate(root)
    with (root.parent/('.'+root.name+'.refresh.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        before = hashes(root)
        with tempfile.TemporaryDirectory(prefix='.'+root.name+'.refresh-', dir=root.parent) as tmp:
            stage = Path(tmp)/'report'
            shutil.copytree(root, stage)
            for sample in samples: refresh_page(stage, sample)
            remove_retired_cohort(stage)
            audit = validate(stage)
            if audit['sample_reports'] != len(samples): raise ValueError('Sample count changed')
            # Updating presentation cannot modify scientific inputs or provenance.
            for name, digest in before.items():
                protected = name.startswith('data/') or (name != 'manifest.json' and Path(name).suffix in {'.tsv', '.csv', '.json', '.gz'})
                if protected and hashlib.sha256((stage/name).read_bytes()).hexdigest() != digest:
                    raise ValueError('Scientific record changed: '+name)
            data['report_version'] = RELEASE_VERSION
            data['validation'] = audit
            data['presentation_refresh'] = {'schema_version': 1, 'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
            (stage/'manifest.json').write_text(json.dumps(data, indent=2)+'\n')
            if (root/'report.zip').is_file(): build_archive(stage)
            after = hashes(stage)
            changed = sorted(name for name in set(before)|set(after) if before.get(name) != after.get(name))
            # Preserve service-account ownership, including newly generated files.
            root_stat = root.stat()
            for target in [stage, *stage.rglob('*')]:
                source = root/target.relative_to(stage)
                stat = source.stat() if source.exists() else root_stat
                current = target.stat()
                if (current.st_uid, current.st_gid) != (stat.st_uid, stat.st_gid):
                    os.chown(target, stat.st_uid, stat.st_gid)
            if hashes(root) != before: raise ValueError('Report changed during refresh; retry when idle')
            result = dict(report=str(root), applied=apply, changed_files=changed, validation=audit)
            if apply:
                backup = root.parent/('.'+root.name+'.backup-'+uuid.uuid4().hex)
                root.rename(backup)
                try: stage.rename(root)
                except BaseException:
                    backup.rename(root)
                    raise
                result['backup'] = str(backup)
            return result


def rollback(root, backup):
    original, saved = Path(root), Path(backup)
    if original.is_symlink() or saved.is_symlink(): raise ValueError('Symlink report roots are not supported')
    root, backup = original.resolve(), saved.resolve()
    if backup.parent != root.parent or not re.fullmatch(re.escape('.'+root.name+'.backup-')+r'[0-9a-f]{32}', backup.name):
        raise ValueError('Backup does not belong to this report')
    validate(backup)
    with (root.parent/('.'+root.name+'.refresh.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        displaced = root.parent/('.'+root.name+'.backup-'+uuid.uuid4().hex)
        root.rename(displaced)
        try: backup.rename(root)
        except BaseException:
            displaced.rename(root)
            raise
    return dict(report=str(root), restored=True, backup=str(displaced))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, help='External JSON inventory: schema_version=1, reports=[{root: ...}]')
    parser.add_argument('--apply', action='store_true', help='Publish after validation; default is dry run')
    parser.add_argument('--rollback', nargs=2, metavar=('REPORT', 'BACKUP'))
    args = parser.parse_args()
    if args.rollback:
        if args.manifest or args.apply: parser.error('--rollback is a separate operation')
        print(json.dumps(rollback(*args.rollback), indent=2)); return
    if not args.manifest: parser.error('--manifest is required')
    config = json.loads(args.manifest.read_text())
    if config.get('schema_version') != 1 or not isinstance(config.get('reports'), list) or not config['reports']:
        raise ValueError('Invalid refresh inventory')
    for entry in config['reports']:
        print(json.dumps(refresh(entry['root'], args.apply, entry.get('manifest_sha256')), indent=2))


if __name__ == '__main__': main()
