import unittest
from dataclasses import replace

from test_scrna_qc import LINKS, PARSE, REF, detect, fixture, rc
from scrna_parse_diagnostics import CATEGORIES, Cassette, diagnose
from scrna_qc import raw_bases_beyond_r2


class ParseDiagnostics(unittest.TestCase):
    def setUp(self):
        self.unit = fixture(PARSE, 0)
        self.r2 = rc(self.unit.seq)
        self.model = Cassette(REF['parse'], LINKS)

    def classify(self, seq):
        return self.model.classify(seq, self.model.canonical(seq))

    def test_fixed_and_locally_realigned_cassettes(self):
        s = self.r2
        self.assertEqual(self.classify(s[:10]+'N'*8+s[18:])[0], 'barcode')
        self.assertEqual(self.classify(s[:22]+'NNN'+s[25:])[0], 'linker')
        self.assertEqual(self.classify(s[:24]+s[25:])[0], 'linker')  # one-base linker deletion
        self.assertEqual(self.classify(s[:32]+s[33:])[0], 'barcode')  # deletion between linkers
        self.assertEqual(self.classify('G'+s)[0], 'shifted')
        self.assertEqual(self.classify(s[5:]), ('short_end', 'early_anchors'))
        self.assertEqual(self.classify(s[22:]), ('short_end', 'early_anchors'))
        self.assertEqual(self.classify(s[:54]), ('short_end', 'short_insert'))
        self.assertEqual(self.classify('G'*54)[0], 'unresolved')
        self.assertEqual(self.classify('A'*18+LINKS[0]+'A'*36)[0], 'unresolved')

    def test_rounds_are_oriented_and_near_assignment_is_unique(self):
        onlist = ({'AAAAAAAA', 'AAAAAACC'}, {'CCCCCCCC'}, {'GGGGGGGG'})
        model = Cassette(onlist, LINKS)
        seq = 'C'*10+'AAAAAAAC'+LINKS[0]+'CCCCCCCT'+LINKS[1]+'GGGG'
        self.assertEqual(model.canonical(seq), ['mismatch', 'exact', 'near', 'exact', 'incomplete'])
        masked = seq[:18] + 'N' + seq[19:37] + 'N' + seq[38:]
        self.assertEqual(model.canonical(masked)[1:3], ['near', 'near'])
        self.assertEqual(self.classify(rc(self.r2))[0], 'unresolved')

    def test_raw_boundary_direction_and_missing_coordinates(self):
        row = {'insert_start': '13', 'insert_end': '118', 'raw_length': '123'}
        distances = [raw_bases_beyond_r2(row, False), raw_bases_beyond_r2(row, True),
                     raw_bases_beyond_r2({}, False), 0]
        self.assertEqual(distances, [5, 13, None, 0])
        with self.assertRaises(ValueError):
            raw_bases_beyond_r2(dict(row, insert_end='124'), False)
        units = [replace(self.unit, raw_bases_beyond_r2=n) for n in distances]
        d = diagnose(units, REF['parse'], None)['raw_r2_boundary']
        self.assertEqual(d, {'measured': 3, 'at_least_1': 2, 'at_least_12': 1, 'min_bp': 0})

    def test_all_unmatched_and_only_unmatched_are_counted(self):
        group = [fixture(PARSE, i) for i in range(15)]
        for i, seq in enumerate((self.r2[:54], self.r2[5:], 'A'*35)):
            group.append(replace(self.unit, name='q'+str(i), parent='p'+str(i), seq=rc(seq), qual='I'*len(seq)))
        call, _ = detect(group, REF)
        d = call['unmatched_parse']
        self.assertEqual(d['assessed'], call['protocol_support']['none'])
        self.assertEqual((d['units'], d['partial'], d['shorter_than_58']), (3, 2, 2))
        self.assertEqual(d['short_end_patterns'], {'short_insert': 1, 'early_anchors': 1})
        self.assertEqual(d['categories']['unresolved'], 1)
        for site in d['sites'].values():
            self.assertEqual(sum(site.values()), d['assessed'])
        missing = diagnose([self.unit], REF['parse'], None)
        self.assertEqual((missing['units'], missing['assessed']), (1, 0))
        self.assertEqual(missing['categories'], dict.fromkeys(CATEGORIES, 0))


if __name__ == '__main__':
    unittest.main()
