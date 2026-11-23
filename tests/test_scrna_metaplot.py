import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
import rna_qc
import scrna_metaplot


class TranscriptMetaplot(unittest.TestCase):
    def test_spliced_fragments_are_gene_weighted_stranded_and_scope_specific(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'genes.gtf'
            features = [('A', 'a', '+', 'exon', 1, 25), ('A', 'a', '+', 'exon', 76, 150),
                        ('A', 'a', '+', 'CDS', 11, 25), ('A', 'a', '+', 'CDS', 76, 130),
                        ('A', 'a', '+', 'stop_codon', 131, 133), ('A', 'a', '+', '5UTR', 3, 11),
                        ('B', 'b', '-', 'exon', 301, 400), ('B', 'alt', '-', 'exon', 401, 410)]
            path.write_text(''.join(f'chr1\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t'
                                    f'gene_id "{gene}"; transcript_id "{tx}";\n'
                                    for gene, tx, strand, kind, start, end in features))
            annotation = rna_qc.read_annotation(path, transcript_regions=True)
        units, alignments, families = {}, {}, []

        def unit(name, scope='lib01', kit='parse-wt-58', start=0, cigar='25M50N25M',
                 reverse=False, mapq=60, nh=1, sa=False, family='protein_coding'):
            units[name] = dict(scope_id=scope, kit=kit, orientation='p5_to_p7', assignment_status='assigned')
            blocks, junctions, end, _ = rna_qc.parse_cigar(cigar, start)
            alignments[name] = [dict(chrom='chr1', start=start, end=end, blocks=blocks, junctions=junctions,
                                     secondary=False, supplementary=False, nh=nh, sa=sa, mapq=mapq, reverse=reverse)]
            families.append({'unit_id': name, 'family': family})

        for i in range(4):
            unit(str(i))
        unit('unindexed', 'unassigned', '10x-3p-v3', start=300, cigar='50M')
        unit('wrong_strand', reverse=True)
        unit('multimapped', nh=2)
        unit('low_mapq', mapq=19)
        unit('chimera', sa=True)
        unit('panel_rrna', family='rRNA')
        unit('mito_rrna', family='mt-rRNA')
        unit('intronic', start=32, cigar='10M')
        unit('alternate_exon', 'lib02', '10x-5p-v3', start=400, cigar='10M', reverse=True)

        rows, _ = rna_qc.analyze_units(units, alignments, annotation, 20, collect_profiles=False)
        profile = scrna_metaplot.summarize(rows, alignments, annotation, units, families)
        self.assertEqual(profile['annotation_sha256'], annotation.metadata['sha256'])
        scopes = profile['scopes']
        whole = scopes['all']['gene_body']
        self.assertEqual((whole['supporting_units'], whole['contributing_genes']), (5, 2))
        self.assertEqual(whole['mean_relative_coverage'], [1.] * 100)
        for value in whole['sem_across_genes']:
            self.assertAlmostEqual(value, 1.)
        self.assertEqual(scopes['lib01']['gene_body']['mean_relative_coverage'], [2.] * 50 + [0.] * 50)
        self.assertEqual(scopes['unassigned']['gene_body']['mean_relative_coverage'], [0.] * 50 + [2.] * 50)
        self.assertEqual(scopes['lib02']['gene_body']['mean_relative_coverage'], [None] * 100)
        self.assertEqual(scopes['lib01']['gene_body']['sem_across_genes'], [None] * 100)
        self.assertEqual(scopes['all']['bias']['status'], 'insufficient_terminal_supported_genes')
        cds = scopes['all']['CDS']
        self.assertEqual((cds['contributing_genes'], cds['covered_genes'], cds['genes_missing_region']), (1, 1, 1))
        self.assertEqual(scopes['all']['three_prime_UTR']['mean_relative_coverage'], [0.] * 30)
        self.assertEqual(scopes['all']['three_prime_UTR']['supporting_units'], 0)
        representative = annotation.genes['A@chr1:+']
        self.assertEqual({k: rna_qc.interval_length(v) for k, v in representative.segments.items()},
                         {'CDS': 73, 'five_prime_UTR': 8, 'three_prime_UTR': 17})


if __name__ == '__main__':
    unittest.main()
