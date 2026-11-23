"""Present physical sample-index calls alongside single-cell chemistry results."""
import hashlib
import html
import json
import re
from pathlib import Path

from report_render import identity_panel_html


def library_identity(summary):
    library = summary.get('physical_library') or {}
    pairing = library.get('index_pairing')
    if pairing is None:
        # Early reports retained the indexes but not the evidence for a mode or share.
        pairing = {'mode': 'not_recorded'}
        saved = summary.get('physical_indexes')
        pairs = None if saved is None else [dict(indexes, scope_id=scope)
                                             for scope, indexes in sorted(saved.items())]
    else:
        pairs = pairing.get('pairs', [])
    metric = dict(library, reference_assembly=summary.get('assembly', 'unspecified'), assay='scRNA-seq')
    if pairs and not metric.get('sub_libraries'):
        metric['sub_libraries'] = [{'scope_id': p['scope_id']} for p in pairs]
    return identity_panel_html(metric, pairing, pairs, show_support_counts=True)


def naming_assets(summary):
    """Stable IDs stay in the data files. Names annotate the matching report scopes."""
    pairing = (summary.get('physical_library') or {}).get('index_pairing')
    if pairing is None:
        catalog = summary.get('physical_indexes') or {}
    else:
        catalog = {p['scope_id']: p for p in pairing.get('pairs', [])}
    scopes = sorted(s for s in catalog if re.fullmatch(r'lib[0-9]{2,6}', s))
    if not scopes:
        return ''
    # A new assignment of the same library IDs to different indexes starts its own local names.
    identity = [summary.get('input_sha256'), summary.get('sample'), summary.get('assembly'),
                [(s, catalog[s].get('i5'), catalog[s].get('i7')) for s in scopes]]
    key = hashlib.sha256(json.dumps(identity, ensure_ascii=True).encode('ascii')).hexdigest()
    fields = ''.join(f'<label><strong>{s}</strong><input type="text" data-scope-id="{s}" '
                     f'aria-label="Name for {s}" autocomplete="off" spellcheck="false"></label>' for s in scopes)
    script = Path(__file__).with_name('scrna_names.js').read_text()
    style = Path(__file__).with_name('scrna_names.css').read_text()
    return (f'<style>{style}</style><dialog id="sc-name-editor" data-report-key="{key}" '
            'aria-labelledby="sc-name-title"><form><h2 id="sc-name-title">Sublibrary names</h2>'
            f'<p class="sc-name-context">{html.escape(summary.get("sample", ""))}</p>'
            f'<fieldset><legend class="sr-only">Names</legend>{fields}</fieldset>'
            '<p class="sc-name-storage"></p><p class="sc-name-status" role="status"></p>'
            '<div class="sc-name-actions"><button type="button">Cancel</button>'
            '<button type="submit">Save</button></div></form></dialog>'
            f'<script>{script}</script>')
