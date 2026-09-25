"""Sequence-aware library diagrams from saved reconstruction evidence."""
from collections import Counter
import html
from library_structure import normalize_structure, umi_lengths, handles, POSITIONS, adapter_families

FAMILY_LABELS = {'truseq': 'TruSeq', 'nextera': 'Nextera',
                 'custom': 'Custom adapters', 'custom_truseq': 'Custom + TruSeq',
                 'custom_nextera': 'Custom + Nextera'}


def family_label(value):
    return FAMILY_LABELS.get(value, str(value).replace('_', ' ').capitalize())


def _escape(value):
    return html.escape(str(value), quote=True)


def _positive(value):
    return value is True or str(value).lower() in {'1', 'true'}


def library_structure_html(metric, adapter_summary=None, unit_rows=None, *, adapter_templates):
    summary = adapter_summary or {}
    structure = normalize_structure(metric.get('library_structure') or summary.get('library_structure'))
    adapter_templates = adapter_families(structure, adapter_templates)
    lengths = umi_lengths(structure)
    p5, p7 = handles(structure)
    rows = list(unit_rows) if unit_rows is not None else None
    pairing = metric.get('index_pairing') or summary.get('index_pairing', {})
    pairs = pairing.get('pairs', [])
    if pairing.get('mode') == 'single_library' and pairs:
        pairs = [max(pairs, key=lambda p: p.get('informative_support', 0))]
    if not pairs and metric.get('sub_libraries'):
        pairs = [dict(scope_id=s.get('scope_id'), name=s.get('scope_name'),
                      i5=s.get('i5'), i7=s.get('i7'), support_units=s.get('library_units', 0))
                 for s in metric['sub_libraries']]
    has_subs = len(pairs) > 1 or bool(metric.get('sub_libraries'))

    def sequence(value):
        if value is None:
            return "<span class='library-unresolved'>Not detected</span>"
        if value == '':
            return "<span class='library-unresolved'>No index · 0 bp</span>"
        return f"<code class='library-sequence'>{_escape(value)}</code>"

    def index_content(name, pair=None):
        if pair is not None:
            values = [pair.get(name)]
        elif pairs:
            values = list(dict.fromkeys(p.get(name) for p in pairs))
        else:
            value = summary.get(name, metric.get(name))
            values = [value.get('sequence') if isinstance(value, dict) else value]
        if len(values) > 1:
            resolved = sum(v is not None for v in values)
            note = f'{resolved} resolved ' + ('sequence' if resolved == 1 else 'sequences')
            if None in values:
                note += ' · some unresolved'
            return f"<small>{note}</small><span class='library-unresolved'>See sublibraries below</span>"
        value = values[0] if values else None
        note = f"<small>{len(value)} bp</small>" if value else ''
        return note + sequence(value)

    sequence_groups = []

    def diagram(scope_rows, pair=None, scope_name='Total library'):
        records = []
        if scope_rows is not None:
            families = Counter(row.get('family') for row in scope_rows if row.get('family') in adapter_templates)
        elif pair is None:
            families = Counter(summary.get('adapter_families', {}))
        else:
            families = Counter()
        # Older metrics still identify a matched family even without its count.
        if not families and scope_rows is None and pair is None and metric.get('adapter_family') in adapter_templates:
            families[metric['adapter_family']] = metric.get('library_units', 0)
        families = {name: count for name, count in families.items() if count > 0}
        def compact(name, note=''):
            return f"<small class='library-type'>{_escape(name)}</small>" + (f"<span class='library-unresolved'>{_escape(note)}</span>" if note else '')
        def handle_content(name, motif):
            if scope_rows is not None:
                count = sum(_positive(row.get(name.lower() + '_detected')) for row in scope_rows)
            elif pair is None:
                observations = summary.get('i5' if name == 'P5' else 'i7', {})
                count = observations.get('all_flanked_observations') if isinstance(observations, dict) else None
            else:
                count = None
            configured = structure.get(name.lower() + '_sequence')
            label = 'Custom' if configured else 'Standard'
            if count:
                records.append((name, label, f'Matched · {count:,} units', motif))
                return compact(label)
            status = 'Configured · not detected' if configured else ('Not detected' if count == 0 else 'Detection unavailable')
            records.append((name, label, status, configured))
            return compact(label, 'Not detected' if configured or count == 0 else 'Unavailable')
        def adapter_content(side):
            configured = structure.get('adapter_i5_sequence' if side == 0 else 'adapter_i7_sequence')
            element = 'Adapter · ' + ('i5 side' if side == 0 else 'i7 side')
            if not families:
                label = 'Custom' if configured else 'Unresolved'
                records.append((element, label, 'Configured · not detected' if configured else 'Not detected', configured))
                return compact(label, 'Not detected')
            motifs = Counter()
            for name, count in families.items():
                motif = adapter_templates[name][side]
                if side == 1:
                    motif = motif.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                label = 'Custom' if configured else family_label(name.removeprefix('custom_'))
                motifs[(label, motif)] += count
            for (label, motif), count in sorted(motifs.items(), key=lambda item: -item[1]):
                records.append((element, label, f'Matched · {count:,} units', motif))
            return compact(' / '.join(dict.fromkeys(label for label, motif in motifs)))
        blocks = [('P5', 'handle', handle_content('P5', p5)), ('i5', 'index', index_content('i5', pair)),
                  ('Adapter', 'adapter', "<small>i5 side</small>" + adapter_content(0)),
                  ('Insert', 'insert', "<small>Sample DNA</small>"),
                  ('Adapter', 'adapter', "<small>i7 side</small>" + adapter_content(1)),
                  ('i7', 'index', index_content('i7', pair)), ('P7', 'handle', handle_content('P7', p7))]
        parts = ["<span class='library-end' style='grid-column:1'>5′</span>"]
        number = 0
        for i, (label, kind, content) in enumerate(blocks):
            parts.append(f"<span class='library-block {kind}' style='grid-column:{2 + 2 * i}'><strong>{label}</strong>{content}</span>")
            if i < len(POSITIONS):
                if POSITIONS[i] in lengths:
                    number += 1
                    parts.append(f"<span class='library-block umi' data-position='{POSITIONS[i]}' style='grid-column:{3 + 2 * i}'><strong>UMI {number}</strong><small>{lengths[POSITIONS[i]]} bp</small></span>")
                else:
                    parts.append(f"<span class='library-connector' aria-hidden='true' style='grid-column:{3 + 2 * i}'>—</span>")
        parts.append("<span class='library-end' style='grid-column:15'>3′</span>")
        start, end = (7 if 'before_insert' in lengths else 8), (10 if 'after_insert' in lengths else 9)
        first = 'the i5-side UMI' if 'before_insert' in lengths else 'the left end of the insert'
        second = 'the i7-side UMI' if 'after_insert' in lengths else 'the right end of the insert on the opposite strand'
        parts.append(f"<div class='library-read-directions' role='group' aria-label='NGS read starts and directions' style='grid-column:{start}/{end}'>"
                     f"<span class='library-read forward' title='Read 1 starts at {first}'>R1<i aria-hidden='true'></i></span>"
                     f"<span class='library-read reverse' title='Read 2 starts at {second}'>R2<i aria-hidden='true'></i></span></div>")
        sequence_groups.append(((pair or {}).get('scope_id') or 'total', scope_name, records))
        return "<div class='library-diagram-scroll' tabindex='0' aria-label='Library sequences; scroll horizontally to view the full diagram'><div class='library-diagram'>" + ''.join(parts) + "</div></div>"

    mode = 'Auto detection' if structure['mode'] == 'auto' else 'Customized'
    sections = ["<details class='metric-group library-structure-report' open><summary>Library structure "
                f"<span class='library-mode'>{mode}</span></summary><div class='library-structure-body'>",
                diagram(rows),
                "<p class='library-diagram-note'>i5/i7 show resolved index sequences. Short arrows mark read starts and directions, not read lengths. Full P5/P7 and adapter sequences are listed below.</p>"]
    if has_subs:
        sections.append("<div class='library-sublibraries'><h3>Resolved sublibraries</h3>")
        sub_metrics = {s.get('scope_id'): s for s in metric.get('sub_libraries', [])}
        for pair in pairs:
            scope = pair.get('scope_id')
            subset = [row for row in rows if row.get('scope_id') == scope and row.get('assignment_status') == 'assigned'] if rows is not None else None
            label = pair.get('name') or sub_metrics.get(scope, {}).get('scope_name') or scope or 'Sublibrary'
            count = pair.get('support_units', len(subset) if subset is not None else 0)
            share = count / metric['library_units'] if metric.get('library_units') else None
            support = f"{count:,} assigned units" + (f" · {share:.1%} of all units" if share is not None else '')
            sections.append(f"<section class='library-sublibrary'><header><h4>{_escape(label)}</h4><span>{support}</span></header>{diagram(subset, pair, label)}</section>")
        sections.append("</div>")
    if lengths:
        umi = metric.get('umi_summary') or summary.get('umi_summary', {})
        sections.append("<div class='library-umi-summary'><strong>UMI extraction</strong><div class='library-umi-counts'>" + ''.join(
            f"<span><b>{umi.get(field, 0):,}</b>{label}</span>" for field, label in
            [('complete_units', 'Complete units'), ('partial_units', 'Partial units'), ('rejected_short_units', 'Too short for layout')]) +
            "</div><p>UMIs are recorded in the per-unit table and excluded from mapped inserts. Full-library lengths include them. Counts and deduplication follow the assay settings; fixed-offset extraction can be affected by ONT indels. Unreconstructed-read dimer template screening is unavailable for random-UMI layouts.</p></div>")
    sections.append("<details class='library-sequence-details'><summary>P5/P7 &amp; adapter sequences</summary>"
                    "<p>Sequences follow the physical 5′ P5-to-P7 strand. Detected sequences are matched reference motifs, not de novo consensus calls. Custom sequences without supporting matches are marked as configured. Standard P5/P7 handles are shared across library types and do not identify a preparation kit.</p>")
    for scope_id, scope_name, records in sequence_groups:
        sections.append(f"<section class='library-sequence-group' data-scope='{_escape(scope_id)}'><h4>{_escape(scope_name)}</h4><div class='library-sequence-table-scroll'><table><thead><tr><th>Element</th><th>Type</th><th>Evidence</th><th>Sequence · 5′ → 3′</th></tr></thead><tbody>")
        for element, label, status, motif in records:
            bases = sequence(motif) if motif else "<span class='library-unresolved'>—</span>"
            sections.append(f"<tr><td>{_escape(element)}</td><td>{_escape(label)}</td><td>{_escape(status)}</td><td>{bases}</td></tr>")
        sections.append("</tbody></table></div></section>")
    sections.append("</details></div></details>")
    return ''.join(sections)
