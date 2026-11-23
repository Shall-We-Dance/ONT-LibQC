import copy
import unittest
from dataclasses import replace

from test_scrna_qc import C4, C4_COMPACT, PARSE, REF, detect, fixture


class ProtocolSupport(unittest.TestCase):
    def test_each_unit_votes_once_even_if_multiple_layouts_match(self):
        refs = copy.deepcopy(REF)
        tenx = fixture('10x-3p-v3', 0)
        sibling = replace(tenx, name='second-unit-from-same-ONT-read')
        hybrid = fixture(C4_COMPACT, 1)
        # Positions 30 and 28 of this unit begin overlapping poly(T) anchors
        # for MGI C4 and 10x. Its first 16 bases also appear on two 10x lists.
        hybrid.seq = hybrid.seq[:20] + 'ACGTACGT' + 'T' * 12 + hybrid.seq[39:]
        hybrid.qual = 'I' * len(hybrid.seq)
        refs['tenx']['3p-v3'].add(hybrid.seq[:16])
        refs['tenx']['3p-v4'].update((tenx.seq[:16], hybrid.seq[:16]))
        missing_anchor = replace(tenx, name='barcode-only', seq=tenx.seq[:16] + 'C' * 60,
                                 qual='I' * 76)
        group = [tenx, sibling, hybrid, fixture(C4, 2), missing_anchor]
        call, _ = detect(group, refs)

        self.assertEqual(call['protocol_support'], {'parse': 0, '10x': 2, 'mgi': 1,
                                                   'multiple': 1, 'none': 1})
        self.assertEqual(sum(call['protocol_support'].values()), call['units'])
        self.assertEqual(call['evidence']['10x-3p-v3']['supported_units'], 3)
        self.assertEqual(call['evidence']['10x-3p-v4']['supported_units'], 3)
        self.assertEqual(call['evidence']['10x-3p-v3']['exact_barcode_units'], 4)
        self.assertEqual(call['evidence'][C4_COMPACT]['supported_units'], 1)

    def test_parse_and_empty_scope_counts(self):
        group = [fixture(PARSE, i) for i in range(12)]
        call, _ = detect(group, REF)
        self.assertEqual(call['protocol_support']['parse'], 12)
        self.assertEqual(sum(call['protocol_support'].values()), len(group))
        empty, _ = detect([], REF)
        self.assertEqual(empty['protocol_support'], dict.fromkeys(call['protocol_support'], 0))
