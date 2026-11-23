import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
import rna_qc
import scrna_regions
from scrna_report import region_section


def feature(gene, tx, typ, start, end, strand='+'):
    return f'chr1\ttest\t{typ}\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene}"; transcript_id "{tx}";\n'


def profile(features, reads):
    with tempfile.TemporaryDirectory() as t:
        path = Path(t) / 'genes.gtf'
        path.write_text(''.join(features))
        annotation = rna_qc.read_annotation(path, transcript_regions=True)
        rows, mappings = [], {}
        for i, (cigar, status, scope) in enumerate(reads):
            name = str(i)
            blocks, _, _, _ = rna_qc.parse_cigar(cigar, 0)
            rows.append({'unit_id': name, 'scope_id': scope, 'mapping_status': status})
            mappings[name] = [dict(chrom='chr1', blocks=blocks, secondary=False, supplementary=False)]
        return scrna_regions.summarize(rows, mappings, annotation)


class RegionComposition(unittest.TestCase):
    def test_reverse_strand_stop_codon_and_spliced_blocks(self):
        g = lambda kind, start, end: feature('g', 'tx', kind, start, end, '-')
        values = profile([g('exon', 1, 30), g('exon', 41, 80), g('CDS', 21, 30),
                          g('CDS', 41, 65), g('stop_codon', 18, 20)],
                         [('30M10N40M', 'unique', 'lib01'), ('80M', 'multimapped', 'lib01')])
        self.assertEqual(values['unique_units'], 1)
        self.assertEqual(values['aligned_bases'], 70)
        self.assertEqual([values['counts'][key] for key, _, _ in scrna_regions.REGIONS], [15, 38, 17])
        self.assertEqual(values['coding_region_bases'], 70)
        self.assertEqual(values['scopes']['lib01']['counts'], values['counts'])

    def test_aliases_and_isoform_boundaries_and_competing_utrs(self):
        g = lambda tx, kind, start, end: feature('g', tx, kind, start, end)
        features = [g('a', 'exon', 1, 50), g('a', 'CDS', 10, 20), g('a', '5UTR', 4, 9),
                    g('b', 'exon', 31, 60), g('b', 'CDS', 41, 55), g('b', '3UTR', 56, 59)]
        values = profile(features, [('60M', 'unique', 'lib01')])
        counts = values['counts']
        self.assertEqual(counts['five_prime_UTR'], 6)
        self.assertEqual(counts['CDS'], 26)
        self.assertEqual(counts['three_prime_UTR'], 14)
        self.assertEqual(counts['ambiguous_UTR'], 10)
        self.assertEqual(counts['other_exon'], 4)
        self.assertEqual(sum(counts.values()), values['aligned_bases'])

    def test_no_coding_or_unique_support_and_html(self):
        feats = [feature('g','a','exon',1,20)]
        counts = profile(feats, [('15M5D5M', 'unique', ''), ('25M', 'low_mapq', '')])
        self.assertEqual(counts['scopes']['unassigned']['aligned_bases'], 20)
        self.assertEqual(counts['coding_region_bases'], 0)
        self.assertIsNone(scrna_regions.figure(counts))
        self.assertIn('None of the 20', region_section(counts, None, 20))
        empty = profile(feats, [('20M', 'unmapped', 'lib01')])
        self.assertIn('No uniquely mapped', region_section(empty, None, 20))
        good = profile(feats + [feature('g','a','CDS',5,15)], [('20M','unique','lib01')])
        svg = scrna_regions.figure(good)
        self.assertIn(b'Coding transcript regions', svg)
        self.assertIn('55.0%', region_section(good, svg, 20))


if __name__ == '__main__':
    unittest.main()
