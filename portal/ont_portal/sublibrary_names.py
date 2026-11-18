"""Display names are annotations on stable report scope IDs."""
import html
import json
import re
import unicodedata
from pathlib import Path

from .schema import SAFE_NAME

SCOPE = re.compile(r'lib[0-9]{2,6}')
MAX_NAME = 48
MAP_FILE = 'data/sublibrary_names.tsv'
TABLE = re.compile(r"(<h3>Physical (?:indexes|index pairs)</h3>)(<div class=['\"]pair-table-scroll['\"][^>]*><table class=['\"]pair-table['\"]>)(.*?)</table>", re.S)
LIBRARY_CELL = re.compile(r'<tr(?: data-sublibrary="(lib[0-9]{2,6})")?><td><strong>([^<]*)</strong></td>')


def scope_catalog(root, sample):
    if not isinstance(sample, str) or not re.fullmatch(SAFE_NAME, sample):
        return {}
    root = Path(root)
    path = root/'data'/'metrics'/(sample + '.json')
    if root.is_symlink() or not path.is_file() or path.resolve() != root.resolve()/'data'/'metrics'/path.name or path.stat().st_size > 8 * 1024**2:
        return {}
    try:
        scopes = json.loads(path.read_bytes()).get('sub_libraries') or []
        catalog = {}
        for row in scopes:
            sid = row['scope_id']
            if not SCOPE.fullmatch(sid) or sid in catalog:
                return {}
            fields = {key: row.get(key) for key in ('scope_name', 'i5', 'i7')}
            if any(value is not None and not isinstance(value, str) for value in fields.values()):
                return {}
            catalog[sid] = fields
        return catalog
    except (ValueError, KeyError, TypeError, AttributeError):
        return {}


def has_sublibraries(root, primary):
    # A retained bundle can contain more than the task's primary sample.
    if scope_catalog(root, primary):
        return True
    folder = root/'data'/'metrics'
    return (folder.is_dir() and folder.resolve() == root.resolve()/'data'/'metrics'
            and any(scope_catalog(root, p.stem) for p in folder.glob('*.json') if p.stem != primary))


def validate_names(values, catalog):
    if not isinstance(values, dict) or not values.keys() <= catalog.keys():
        raise ValueError('Choose a sublibrary listed in this report.')
    names = {}
    for sid, value in values.items():
        if not isinstance(value, str):
            raise ValueError('Enter a text name for ' + sid + '.')
        value = unicodedata.normalize('NFC', value).strip(' ')
        if not value:
            continue
        if (len(value) > MAX_NAME or unicodedata.category(value[0])[0] not in 'LN'
                or any(unicodedata.category(c)[0] not in 'LNM' and c not in ' ._()+-' for c in value)):
            raise ValueError(f'{sid}: use up to {MAX_NAME} characters, starting with a letter or number. Spaces and . _ ( ) + - are allowed.')
        names[sid] = value
    return names


def name_map_tsv(annotations):
    def index(value):
        return value if isinstance(value, str) and re.fullmatch('[ACGTNacgtn]+', value) else ''
    lines = ['sample\tscope_id\tname\ti5\ti7']
    for sample, (catalog, names) in annotations.items():
        for sid, row in catalog.items():
            lines.append('\t'.join([sample, sid, names.get(sid, ''), index(row['i5']), index(row['i7'])]))
    return '\n'.join(lines) + '\n'


def annotate_page(body, catalog, names, *, editor=None):
    """Render the index crosswalk and use it to label the report's tables and legends."""
    if not catalog:
        return body
    # Older published pages use the configured scope name in the Library cell.
    ids = {}
    for sid, row in catalog.items():
        label = row.get('scope_name') or sid
        ids[label] = sid if label not in ids else None
    def table(match):
        heading, opening, rows = match.groups()
        seen = set()
        def add_cell(cell):
            label = html.unescape(cell[2])
            sid = cell[1] or ids.get(label)
            if sid not in catalog or sid in seen:
                return cell[0]
            seen.add(sid)
            original = '' if label == sid else label
            value = html.escape(names.get(sid) or original or '—')
            return (f'<tr><td><strong>{sid}</strong></td><td class="scope-name" data-scope-id="{sid}" '
                    f'data-original="{html.escape(label, quote=True)}">{value}</td>')
        rows = LIBRARY_CELL.sub(add_cell, rows)
        if seen != set(catalog) or len(seen) != rows.count('<tr>') - 1 or '<th>Library</th>' not in rows:
            return match[0]
        rows = rows.replace('<th>Library</th>', '<th>Library</th><th>Name</th>', 1)
        if editor:
            heading = '<h3>Physical indexes<a class="edit-sublibraries" href="' + html.escape(editor, quote=True) + '">Edit names</a></h3>'
        return heading + opening + rows + '</table>'
    body, matched = TABLE.subn(table, body, count=1)
    if not matched or 'class="scope-name"' not in body:
        return body
    style = """<style>
@media screen and (min-width:901px){.identity-panel:has(.scope-name){grid-template-columns:minmax(235px,.61fr) minmax(0,1.4fr)}}
.scope-name,.scope-display-name{overflow-wrap:anywhere}.scope-name{min-width:6em}
.edit-sublibraries{float:right;font-size:inherit;font-weight:500}.scope-display-name{font-weight:normal}
@media print{.edit-sublibraries{display:none}.pair-table .scope-name{min-width:0}}
</style>"""
    body = body.replace('</head>', style + '</head>', 1)
    if names:
        body = body.replace('</body>', '<script>' + LABEL_SCRIPT + '</script></body>', 1)
    return body


LABEL_SCRIPT = Path(__file__).with_name('sublibrary_labels.js').read_text()
