"""Verify display geometry against the actual extraction rules."""
import sys
import unittest
from base64 import b64decode
from pathlib import Path
from xml.etree import ElementTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from scrna_chemistry import C4, C4_COMPACT, LABELS, PARSE, TENX, TENX_TSO
from scrna_chemistry_report import _segments, render_chemistry, structure
from scrna_report import table


def call(kit, status='identified'):
    return {'kit': kit, 'status': status, 'units': 60, 'evidence': {
        k: dict(supported_units=42 if k == kit else 0, parents=40 if k == kit else 0,
                exact_barcode_units=46 if k == kit else 0, opposite_end_exact=0,
                distinct_barcodes=21 if k == kit else 0) for k in LABELS}}


def block_labels(markup):
    return [node.text for node in ElementTree.fromstring(markup).iter('strong')]


class TestChemistryDisplay(unittest.TestCase):
    def test_protocol_position_lengths_and_orientation(self):
        blocks, _ = _segments(PARSE)
        self.assertEqual([k for k, _, _ in blocks],
                         ['UMI', 'BC3', 'Linker 1', 'BC2', 'Linker 2', 'BC1', 'cDNA'])
        self.assertEqual(sum(n for _, _, n in blocks[:-1]), 58)
        self.assertEqual(block_labels(structure(PARSE)),
                         ['cDNA', 'BC1', 'Linker 2', 'BC2', 'Linker 1', 'BC3', 'UMI'])
        for kit, (_, length, anchor) in TENX.items():
            blocks, _ = _segments(kit)
            self.assertEqual([v[2] for v in blocks[:-1]], [16, length, '≥6 bp' if anchor == 'polyT' else len(TENX_TSO)])
            self.assertEqual(block_labels(structure(kit)), [v[0] for v in blocks])
        self.assertEqual([n for _, _, n in _segments(C4)[0][:-2]], [6, 10, 6, 10, 5, 10])
        self.assertEqual([n for _, _, n in _segments(C4_COMPACT)[0][:-2]], [10, 10, 10])

    def test_one_diagram_per_supported_layout_and_escaped_evidence(self):
        confirmed = call(PARSE)
        confirmed['linkers'] = ['A'*12, 'C'*12]
        summary = dict(scopes={'lib1<': confirmed, 'Unassigned': call(PARSE), 'lib2': call('10x-3p-v2')},
                       physical_indexes={'lib1<': {'i5': None, 'i7': 'X<Y'}}, adapter_family_counts={'truseq': 80})
        fragment = render_chemistry(summary, table)
        self.assertEqual(fragment.count('class="sc-molecule"'), 2)
        self.assertEqual(fragment.count('class="library-diagram"'), 3)
        self.assertEqual(fragment.count('<td>42 / 60</td><td>70.0%</td>'), 3)
        self.assertNotIn('lib1<', fragment)
        for key in ('Matched / tested', 'all three cell barcodes and both linkers', 'X&lt;Y', 'AAAAAAAAAAAA',
                    'Read 2 linker 1', 'TruSeq (80 units)'):
            self.assertIn(key, fragment)

    def test_outer_indexes_are_scoped_to_the_physical_library(self):
        pair = {'i5': '', 'i7': None}
        outer = structure(PARSE, {'one': pair})
        self.assertEqual(block_labels(outer)[:7], ['P5', 'i5', 'Adapter', 'Insert', 'Adapter', 'i7', 'P7'])
        self.assertIn('No index', outer)
        self.assertIn('Unresolved', outer)
        mixed = structure(PARSE, {'one': pair, 'two': {'i5': 'ACGT', 'i7': 'TGCA'}})
        self.assertEqual(mixed.count('Varies by library'), 2)
        self.assertNotIn('Physical library', structure(C4))

    def test_no_diagram_for_unresolved_insufficient_or_mixed_call(self):
        for status in ('unresolved', 'insufficient_evidence', 'mixed_or_ambiguous'):
            fragment = render_chemistry({'scopes': {'lib1': call(None, status)}}, table)
            self.assertNotIn('class="sc-molecule"', fragment)
            self.assertIn('Detection result', fragment)
            self.assertIn('Distinct ONT reads', fragment)
        self.assertIn('No complete non-dimer units', render_chemistry({'scopes': {}}, table))

    def test_support_plot_counts_all_units_and_keeps_unassigned_details_closed(self):
        main = call(PARSE)
        main.update(units=2100, protocol_support={'parse': 1200, '10x': 1, 'mgi': 0, 'multiple': 3, 'none': 896})
        unassigned = call(PARSE)
        unassigned.update(units=100, protocol_support={'parse': 91, '10x': 0, 'mgi': 0, 'multiple': 0, 'none': 9})
        page = render_chemistry({'scopes': {'lib01': main, 'unassigned': unassigned}}, table)
        body = ElementTree.fromstring('<main>' + page + '</main>')
        figure = body.find("figure[@class='sc-protocol-figure']")
        self.assertIn('2,200 reconstructed non-dimer units', figure.findtext('figcaption/p'))
        svg = figure.find('{*}svg')
        self.assertEqual([x.text for x in svg.findall('{*}text')][1::2],
                         ['1,291 (58.7%)', '1 (<0.1%)', '0 (0.0%)', '3 (0.1%)', '905 (41.1%)'])
        download = figure.find('a')
        self.assertEqual(download.get('download'), 'chemistry_support.svg')
        exported = ElementTree.fromstring(b64decode(download.get('href').split(',', 1)[1]))
        self.assertEqual(svg.findtext('{*}desc'), exported.findtext('{*}desc'))
        hidden = body.find("details[@class='sc-unassigned']")
        self.assertNotIn('open', hidden.attrib)
        self.assertIn('Unassigned index', ''.join(hidden.itertext()))
        body.remove(hidden)
        self.assertNotIn('Unassigned index', ''.join(body.itertext()))
        self.assertEqual(page.count('class="sc-molecule"'), 1)

    def test_unassigned_only_plot_and_legacy_per_kit_evidence(self):
        legacy = call(PARSE)
        page = render_chemistry({'scopes': {'unassigned': legacy}}, table)
        self.assertNotIn('sc-protocol-figure', page)
        body = ElementTree.fromstring('<main>' + page + '</main>')
        self.assertEqual([child.tag for child in body], ['details'])
        self.assertIn('42 / 60', ''.join(body.itertext()))
        legacy['protocol_support'] = {'parse': 42, '10x': 0, 'mgi': 0, 'multiple': 0, 'none': 18}
        page = render_chemistry({'scopes': {'unassigned': legacy}}, table)
        body = ElementTree.fromstring('<main>' + page + '</main>')
        self.assertEqual([child.tag for child in body], ['figure', 'details'])
        self.assertNotIn('open', body.find('details').attrib)
        with self.assertRaises(ValueError):
            legacy['protocol_support']['multiple'] = 1
            render_chemistry({'scopes': {'unassigned': legacy}}, table)

    def test_unmatched_breakdown_pools_scopes_and_exposes_positional_evidence(self):
        def scope_with_diagnostics(units, assessed):
            scope = call(PARSE)
            scope.update(units=units + 2, protocol_support=dict(parse=2, multiple=0, mgi=0, none=units, **{'10x': 0}))
            scope['unmatched_parse'] = dict(units=units, assessed=assessed, partial=assessed,
                categories=dict(barcode=assessed, linker=0, shifted=0, short_end=0, mixed=0, unresolved=0),
                sites={key: dict(exact=0, near=assessed, mismatch=0, incomplete=0)
                       for key in ('bc1', 'bc2', 'bc3', 'linker1', 'linker2')},
                shorter_than_58=0, short_end_patterns=dict(short_insert=0, early_anchors=0), opposite_end_partial=0)
            return scope

        assigned = scope_with_diagnostics(24, 24)
        unassigned = scope_with_diagnostics(7, 7)
        unknown = scope_with_diagnostics(5, 0)
        assigned['unmatched_parse']['raw_r2_boundary'] = dict(measured=24, at_least_1=24, at_least_12=24, min_bp=32)
        unassigned['unmatched_parse']['raw_r2_boundary'] = dict(measured=7, at_least_1=7, at_least_12=7, min_bp=38)
        scopes = dict(lib01=assigned, lib02=unknown, unassigned=unassigned)
        body = ElementTree.fromstring('<main>' + render_chemistry(dict(scopes=scopes), table) + '</main>')
        self.assertEqual(len(body.findall("div[@id='chemistry-evidence']/figure")), 2)
        fig = body.find(".//figure[@class='sc-partial-figure']")
        self.assertIn('31 of 36 unmatched units (86.1%)', fig.findtext('figcaption/p'))
        svg = fig.find('{*}svg')
        labels = [n.text for n in svg.findall('{*}text')]
        self.assertEqual(labels[:2], ['Barcode mismatch', '31 (86.1%)'])
        self.assertEqual(labels[-2:], ['No linker consensus', '5 (13.9%)'])
        self.assertEqual(fig.find('a').get('download'), 'unmatched_parse.svg')
        detail = body.find("details[@class='sc-site-detail']")
        self.assertNotIn('open', detail.attrib)
        self.assertIn('31 units with a known linker pair', ''.join(detail.itertext()))
        self.assertIn('31 of 31 units checked', ''.join(detail.itertext()))
        self.assertIn('shortest extension was 32 bp', ''.join(detail.itertext()))
        self.assertEqual([node.text for node in detail.find('.//tbody/tr')], ['BC3', '0', '31', '0', '0'])
        unassigned.pop('unmatched_parse')
        self.assertNotIn('sc-partial-figure', render_chemistry(dict(scopes=scopes), table))
        self.assertNotIn('sc-partial-figure', render_chemistry(dict(scopes={'unknown': unknown}), table))
