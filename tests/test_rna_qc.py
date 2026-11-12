"""Synthetic RNA regressions: biological coordinate semantics and safeguards."""
import collections
import csv
import gzip
import importlib.util
import itertools
import json
import random
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import rna_qc as rna
import rna_report


def gtf_line(chrom, feature, start, end, strand, gene, transcript=None, biotype='protein_coding', name=None):
    attrs = f'gene_id "{gene}"; gene_name "{name or gene}"; gene_type "{biotype}";'
    if transcript:
        attrs += f' transcript_id "{transcript}";'
    return f'{chrom}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n'


def fixture_gtf():
    # Gene G contains two 100 bp exons, separated by a 900 bp intron.
    return ''.join([
        gtf_line('chr1', 'gene', 101, 1200, '+', 'G'),
        gtf_line('chr1', 'exon', 101, 200, '+', 'G', 'G.1'),
        gtf_line('chr1', 'exon', 1101, 1200, '+', 'G', 'G.1'),
        gtf_line('chr1', 'CDS', 151, 200, '+', 'G', 'G.1'),
        gtf_line('chr1', 'CDS', 1101, 1150, '+', 'G', 'G.1'),
        gtf_line('chr1', 'gene', 2001, 2100, '-', 'R', biotype='rRNA'),
        gtf_line('chr1', 'exon', 2001, 2100, '-', 'R', 'R.1', biotype='rRNA'),
        gtf_line('chrM', 'gene', 1, 100, '+', 'MR', biotype='Mt_rRNA'),
        gtf_line('chrM', 'exon', 1, 100, '+', 'MR', 'MR.1', biotype='Mt_rRNA'),
        gtf_line('chrM', 'gene', 201, 300, '+', 'MC'),
        gtf_line('chrM', 'exon', 201, 300, '+', 'MC', 'MC.1'),
    ])


def sam_record(name, start=101, cigar='100M900N100M', flag=0, chrom='chr1', mapq=60, tags='NH:i:1'):
    _, _, _, length = rna.parse_cigar(cigar, start - 1)
    return f'{name}\t{flag}\t{chrom}\t{start}\t{mapq}\t{cigar}\t*\t0\t0\t' + 'A' * length + '\t' + 'I' * length + '\t' + tags + '\n'


class RNAUnitTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.gtf = self.root / 'test.gtf'
        self.gtf.write_text(fixture_gtf())
        self.annotation = rna.read_annotation(self.gtf)

    def tearDown(self):
        self.temp.cleanup()

    def record(self, cigar='100M900N100M', start=100, reverse=False, nh=1):
        blocks, junctions, end, length = rna.parse_cigar(cigar, start)
        return {'chrom': 'chr1', 'start': start, 'end': end, 'reverse': reverse,
            'secondary': False, 'supplementary': False, 'mapq': 60,
            'blocks': blocks, 'junctions': junctions, 'query_length': length, 'cigar': cigar, 'nh': nh, 'sa': False}

    def units(self, names=('u1',), orientation='p5_to_p7'):
        return {name: {'unit_id': name, 'read_id': name, 'scope_id': 'index01', 'assignment_status': 'assigned',
            'insert_length': '200', 'orientation': orientation, 'umi_status': 'not_configured'} for name in names}

    def test_spliced_span_does_not_cover_intron(self):
        record = self.record()
        self.assertEqual(record['end'] - record['start'], 1100)
        self.assertEqual(rna.interval_length(record['blocks']), 200)
        regions = rna.classify_alignment(record, self.annotation)['regions']
        self.assertEqual(regions['intron'], 0)
        self.assertEqual(regions['CDS'], 100)
        self.assertEqual(regions['five_prime_UTR'], 50)
        self.assertEqual(regions['three_prime_UTR'], 50)
        self.assertEqual(sum(regions.values()), 200)

    def test_deletions_and_introns_are_not_sequenced_bases(self):
        blocks, junctions, end, length = rna.parse_cigar('5S10M3D7M900N20M2I4M', 100)
        self.assertEqual(rna.interval_length(blocks), 41)
        self.assertEqual(length, 48)
        self.assertEqual(end, 1044)
        self.assertEqual(junctions, [(120, 1020)])

    def test_overlapping_annotation_counts_once(self):
        gene = self.annotation.genes['G@chr1:+']
        values = rna.region_bases([(100, 200)], [gene, gene])
        self.assertEqual(sum(values.values()), 100)

    def test_reversed_ont_read_recovers_same_physical_strand(self):
        units = self.units(('forward', 'reverse'))
        units['reverse']['orientation'] = 'p7_to_p5'
        rows, _ = rna.analyze_units(units, {'forward': [self.record()], 'reverse': [self.record(reverse=True)]}, self.annotation)
        self.assertEqual([row['physical_strand'] for row in rows], ['+', '+'])
        self.assertEqual([row['strand_evidence'] for row in rows], ['sense', 'sense'])
        self.assertEqual(rows[0]['coordinate_key'], rows[1]['coordinate_key'])

    def test_unknown_adapter_orientation_remains_unknown(self):
        rows, _ = rna.analyze_units(self.units(orientation=''), {'u1': [self.record()]}, self.annotation)
        self.assertEqual(rows[0]['strand_evidence'], 'unknown')

    def test_opposite_annotation_does_not_force_strand(self):
        self.gtf.write_text(fixture_gtf() + gtf_line('chr1', 'exon', 101, 1200, '-', 'opposite', 'opposite.1'))
        annotation = rna.read_annotation(self.gtf)
        rows, _ = rna.analyze_units(self.units(), {'u1': [self.record()]}, annotation)
        self.assertEqual(rows[0]['strand_evidence'], 'annotation_conflict')
        self.assertFalse(rows[0]['eligible_gene'])

    def test_multifamily_alternatives_are_retained(self):
        primary = self.record(nh=2)
        alternative = self.record(cigar='100M', start=2000)
        alternative['secondary'] = True
        rows, _ = rna.analyze_units(self.units(), {'u1': [primary, alternative]}, self.annotation)
        self.assertEqual(rows[0]['family'], 'ambiguous_family')
        self.assertEqual(rows[0]['families'], ['protein_coding', 'rRNA'])
        self.assertEqual(rows[0]['mapping_status'], 'multimapped')
        self.assertFalse(rows[0]['eligible_gene'])

    def test_known_junction_exact_match(self):
        rows, _ = rna.analyze_units(self.units(), {'u1': [self.record()]}, self.annotation)
        self.assertEqual(rows[0]['known_junctions'], 1)
        self.assertEqual(rows[0]['junctions'][0]['start'], 200)
        self.assertEqual(rows[0]['junctions'][0]['end'], 1100)

    def test_minus_gene_profile_reverses_exon_coordinates(self):
        forward = rna.profile_vector([(100, 120)], [(100, 200), (1100, 1200)], '+', 10)
        reverse = rna.profile_vector([(100, 120)], [(100, 200), (1100, 1200)], '-', 10)
        self.assertEqual(forward[0], 20)
        self.assertEqual(reverse[-1], 20)
        self.assertEqual(forward, list(reversed(reverse)))

    def test_shared_transcript_normalization_keeps_segment_enrichment(self):
        # Only the first 50 bp (5' UTR) are covered: UTR mean is 4, CDS/3' mean 0.
        rows, vectors = rna.analyze_units(self.units(), {'u1': [self.record(cigar='50M')]}, self.annotation)
        profiles = rna.summarize_profiles(rows, vectors, 1, 1)
        self.assertAlmostEqual(statistics_mean(profiles['gene_body']['mean_relative_coverage']), 1)
        self.assertAlmostEqual(statistics_mean(profiles['five_prime_UTR']['mean_relative_coverage']), 4)
        self.assertEqual(sum(profiles['CDS']['mean_relative_coverage']), 0)
        self.assertEqual(profiles['CDS']['contributing_genes'], 1)
        self.assertEqual(profiles['bias']['status'], 'available')
        self.assertEqual(profiles['bias']['five_to_three_ratio_status'], 'five_prime_only')
        self.assertEqual(profiles['bias']['median_end_balance'], 1)

    def test_metaplot_keeps_single_insert_with_default_thresholds(self):
        rows, vectors = rna.analyze_units(self.units(), {'u1': [self.record(cigar='50M')]}, self.annotation)
        profiles = rna.summarize_profiles(rows, vectors)
        self.assertEqual(profiles['five_prime_UTR']['contributing_genes'], 1)
        self.assertEqual(profiles['five_prime_UTR']['supporting_units'], 1)
        self.assertGreater(sum(profiles['five_prime_UTR']['mean_relative_coverage']), 0)
        self.assertEqual(profiles['CDS']['supporting_units'], 0)
        self.assertEqual(sum(profiles['CDS']['mean_relative_coverage']), 0)
        self.assertIsNone(profiles['bias']['median_end_balance'])
        plot = rna.profile_svg(profiles, 'Observed regions')
        ET.fromstring(plot)
        self.assertIn('5′ UTR', plot)
        self.assertIn('3′ UTR', plot)
        self.assertNotIn('Gene body', plot)
        self.assertNotIn('Insufficient depth', plot)

    def test_empty_metaplot_has_no_invented_zero_coverage(self):
        profiles = rna.summarize_profiles([], {})
        self.assertTrue(all(x is None for x in profiles['CDS']['mean_relative_coverage']))
        self.assertIn('No usable alignments', rna.profile_svg(profiles, 'Empty'))

    def test_bias_gate_requires_covered_genes(self):
        rows, vectors = rna.analyze_units(self.units(), {'u1': [self.record()]}, self.annotation)
        gated = rna.summarize_profiles(rows, vectors, 1, 2)
        self.assertIsNone(gated['bias']['five_to_three_ratio'])
        self.assertAlmostEqual(rna.summarize_profiles(rows, vectors, 1, 1)['bias']['five_to_three_ratio'], 1)

    def test_single_end_coverage_bias_is_retained_without_nonfinite_json(self):
        for start, balance, ratio in ((100, 1, None), (1160, -1, 0)):
            rows, vectors = rna.analyze_units(self.units(), {'u1': [self.record(cigar='40M', start=start)]}, self.annotation)
            result = rna.summarize_profiles(rows, vectors, 1, 1)['bias']
            self.assertEqual(result['status'], 'available')
            self.assertEqual(result['eligible_terminal_supported_genes'], 1)
            self.assertEqual(result['eligible_genes_with_both_ends'], 0)
            self.assertEqual(result['median_end_balance'], balance)
            self.assertEqual(result['five_to_three_ratio'], ratio)
            json.dumps(result, allow_nan=False)

    def test_middle_only_coverage_does_not_support_end_bias(self):
        rows, vectors = rna.analyze_units(self.units(), {'u1': [self.record(cigar='40M', start=150)]}, self.annotation)
        result = rna.summarize_profiles(rows, vectors, 1, 1)['bias']
        self.assertEqual(result['status'], 'insufficient_terminal_supported_genes')
        self.assertEqual(result['eligible_terminal_supported_genes'], 0)
        self.assertIsNone(result['median_end_balance'])

    def test_rna_families_and_missing_reference_evidence(self):
        self.assertEqual(self.annotation.genes['MR@chrM:+'].family, 'mt-rRNA')
        self.assertEqual(self.annotation.genes['MC@chrM:+'].family, 'mt-mRNA')
        self.assertIn('tRNA', self.annotation.metadata['missing_families'])

    def test_refflat_style_untyped_gtf_uses_cds_evidence_without_guessing_lncrna(self):
        text = fixture_gtf()
        text = re_sub_biotype(text)
        self.gtf.write_text(text + gtf_line('chrM', 'exon', 401, 450, '+', 'TRNP', 'TRNP.1', biotype=''))
        annotation = rna.read_annotation(self.gtf)
        self.assertEqual(annotation.genes['G@chr1:+'].family, 'protein_coding')
        self.assertEqual(annotation.genes['G@chr1:+'].family_evidence, 'CDS_annotation')
        self.assertEqual(annotation.genes['R@chr1:-'].family, 'unknown')
        self.assertEqual(annotation.genes['TRNP@chrM:+'].family, 'tRNA')
        self.assertEqual(annotation.genes['TRNP@chrM:+'].family_evidence, 'gene_symbol_heuristic')

    def test_pseudogenes_do_not_become_rrna_or_snrna_by_symbol_prefix(self):
        self.assertEqual(rna.classify_biotype('', 'RNA5SP123', 'chr1', {'chrM'}), 'unknown')
        self.assertEqual(rna.classify_biotype('rRNA_pseudogene', 'RNA5S1', 'chr1', {'chrM'}), 'other')
        self.assertEqual(rna.classify_biotype('', 'RNU6-1P', 'chr1', {'chrM'}), 'unknown')
        self.assertEqual(rna.classify_biotype('', 'RNR1', 'chrM', {'chrM'}), 'mt-rRNA')

    def test_rarefaction_matches_exact_enumeration(self):
        genes = ['a', 'a', 'b', 'c']
        for draw in range(5):
            samples = list(itertools.combinations(range(len(genes)), draw))
            expected = sum(len({genes[i] for i in sample}) for sample in samples) / len(samples)
            self.assertAlmostEqual(rna.expected_detected([2, 1, 1], draw), expected)

    def test_shuffle_preserves_length_composition_and_dinucleotides(self):
        sequence = 'AACGTTCGATACCGTTAAAAACCCGGGTTTNNAT'
        dinuc = lambda seq: collections.Counter(zip(seq, seq[1:]))
        variants = set()
        for seed in range(20):
            value = rna.dinucleotide_shuffle(sequence, random.Random(seed))
            self.assertEqual(len(value), len(sequence))
            self.assertEqual(collections.Counter(value), collections.Counter(sequence))
            self.assertEqual(dinuc(value), dinuc(sequence))
            variants.add(value)
        self.assertGreater(len(variants), 1)

    def test_motifs_depth_gate_counts_distinct_sequences(self):
        self.assertEqual(rna.discover_motifs([('a', 'ACGT' * 20)], 2)['status'], 'insufficient_units')
        result = rna.discover_motifs([(str(i), 'ACGT' * 20) for i in range(6)], 5, shuffles=2)
        self.assertEqual(result['status'], 'insufficient_units')
        self.assertEqual(result['distinct_sequences'], 1)
        self.assertEqual(result['engine'], 'HOMER')
        self.assertEqual(result['motifs'], [])

    def write_inputs(self, extra_empty=False):
        sam, table, fastq = self.root / 'reads.sam', self.root / 'adapters.tsv', self.root / 'reads.fq'
        sam.write_text('@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:3000\n@SQ\tSN:chrM\tLN:1000\n' + sam_record('u1') + sam_record('u2', flag=16))
        fields = ['unit_id', 'read_id', 'scope_id', 'assignment_status', 'orientation', 'insert_length', 'is_adapter_dimer', 'umi_status', 'umi_sequence']
        rows = [dict(zip(fields, ['u1', 'raw1', 'index01', 'assigned', 'p5_to_p7', 200, 0, 'complete', 'ACGT'])),
                dict(zip(fields, ['u2', 'raw2', 'index02', 'assigned', 'p7_to_p5', 200, 0, 'complete', 'ACGT']))]
        if extra_empty:
            rows.append(dict(zip(fields, ['empty', 'raw3', 'index01', 'assigned', 'p5_to_p7', 15, 0, 'not_configured', ''])))
        rna.write_tsv(table, rows, fields)
        fastq.write_text(''.join(f'@{name}\n' + 'ACGT' * 50 + '\n+\n' + 'I' * 200 + '\n' for name in ('u1', 'u2')))
        return sam, table, fastq

    def test_end_to_end_scopes_files_html_and_json_are_self_contained(self):
        sam, table, fastq = self.write_inputs()
        out = self.root / 'report'
        result = rna.run_analysis(sam, self.gtf, table, fastq, 'safe<script>', 'clip', out,
            min_gene_units=1, min_bias_genes=1, min_motif_units=2, assembly='test')
        self.assertEqual(result['scopes']['total']['units'], 2)
        self.assertEqual(result['scopes']['index01']['units'], 1)
        self.assertEqual(result['scopes']['index02']['units'], 1)
        self.assertEqual(result['scopes']['total']['umi']['observed_exact_coordinate_umi_groups'], 1)
        for name in result['files']:
            self.assertTrue((out / name).exists(), name)
            if name.endswith('.svg'):
                ET.parse(out / name)
        rendered = rna_report.render_rna_section(out / 'summary.json', '../data/rna/test')
        self.assertIn('id="rna-qc"', rendered)
        self.assertIn('scope_', rendered)
        self.assertNotIn('<script>', rendered)
        self.assertIn('summary.json', rendered)
        self.assertEqual(json.loads((out / 'summary.json').read_text())['assembly'], 'test')

    def test_missing_fastq_is_not_silently_accepted(self):
        sam, table, fastq = self.write_inputs(extra_empty=True)
        with self.assertRaisesRegex(ValueError, 'unexpected adapter units'):
            rna.run_analysis(sam, self.gtf, table, fastq, 'test', 'rnaseq', self.root / 'out')

    def test_explicit_trimmed_empty_remains_in_unit_denominator(self):
        sam, table, fastq = self.write_inputs(extra_empty=True)
        # Known P5/P7 technical removal, 20 bp total; 15 bp physical insert is empty.
        fastq.write_text(''.join(f'@{name}\n' + 'ACGT' * 45 + '\n+\n' + 'I' * 180 + '\n' for name in ('u1', 'u2')))
        prep = self.root / 'preparation.json'
        prep.write_text(json.dumps({'input_units': 3, 'retained_units': 2, 'empty_after_trim': 1, 'p5_trim_bp': 10, 'p7_trim_bp': 10}))
        result = rna.run_analysis(sam, self.gtf, table, fastq, 'test', 'rnaseq', self.root / 'out', preparation_json=prep)
        total = result['scopes']['total']
        self.assertEqual(total['units'], 3)
        self.assertEqual(total['mapping_counts']['empty_after_technical_trim'], 1)
        with gzip.open(self.root / 'out' / 'units.tsv.gz', 'rt') as handle:
            rows = {row['unit_id']: row for row in csv.DictReader(handle, delimiter='\t')}
        self.assertEqual(rows['u1']['insert_length'], '200')
        self.assertEqual(rows['u1']['analysis_length'], '180')
        self.assertEqual(rows['empty']['analysis_length'], '0')

    def test_preparation_sidecar_lengths_umi_and_empty_are_validated(self):
        import rna_prepare
        sam, table, fastq = self.write_inputs(extra_empty=True)
        with fastq.open('a') as handle:
            handle.write('@empty\n' + 'A' * 15 + '\n+\n' + 'I' * 15 + '\n')
        prep_json, prep_tsv = self.root / 'prep.json', self.root / 'prep.tsv.gz'
        prepared = self.root / 'prepared.fastq.gz'
        rna_prepare.prepare(fastq, table, prepared, prep_json, p5=10, p7=10, units_tsv=prep_tsv)
        result = rna.run_analysis(sam, self.gtf, table, prepared, 'test', 'clip', self.root / 'sidecar_out',
            preparation_json=prep_json, preparation_tsv=prep_tsv)
        self.assertEqual(result['scopes']['total']['units'], 3)
        self.assertEqual(result['scopes']['total']['mapping_counts']['empty_after_technical_trim'], 1)
        self.assertEqual(result['scopes']['total']['umi']['observed_exact_coordinate_umi_groups'], 1)
        self.assertIn('preparation_units.tsv.gz', result['files'])
        with self.assertRaisesRegex(ValueError, 'requires --preparation-tsv'):
            rna.run_analysis(sam, self.gtf, table, prepared, 'test', 'clip', self.root / 'missing_sidecar', preparation_json=prep_json)
        with gzip.open(prep_tsv, 'rt') as handle:
            contents = handle.read()
        with gzip.open(prep_tsv, 'wt') as handle:
            handle.write(contents.replace('manual_trim', 'modified_status'))
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            rna.run_analysis(sam, self.gtf, table, prepared, 'test', 'clip', self.root / 'tampered',
                preparation_json=prep_json, preparation_tsv=prep_tsv)

    def test_source_panel_is_independent_of_mapping_and_excludes_gene_motif_qc(self):
        sam, table, fastq = self.write_inputs()
        source = self.root / 'source.tsv'
        source.write_text('unit_id\tsource_class\treference_name\tmatched_fraction\tidentity\taligned_bases\ttaxon\ttaxon_ambiguous\tstatus\n'
            'u1\trRNA\tRRNA_REFERENCE\t1.0\t0.99\t198\tBacteria A,Bacteria B\t1\tassigned\n'
            'u2\t\t\t\t\t\t\t0\tunassigned\n')
        meta = self.root / 'sources.json'
        meta.write_text(json.dumps({'panel_sha256': 'example', 'limitations': ['Candidate references only']}))
        result = rna.run_analysis(sam, self.gtf, table, fastq, 'test', 'clip', self.root / 'source_out',
            min_gene_units=1, min_motif_units=1, source_assignments=source, source_summary=meta)
        total = result['scopes']['total']
        self.assertEqual(total['mapped_units'], 2)
        self.assertEqual(total['source_evidence']['assigned_units'], 1)
        self.assertEqual(total['source_evidence']['ambiguous_taxon_units'], 1)
        self.assertEqual(total['family_counts']['rRNA'], 1)
        self.assertEqual(total['gene_assigned_units'], 1)
        self.assertEqual(total['motifs']['eligible_units'], 1)
        self.assertEqual(result['source_screen']['panel_sha256'], 'example')
        self.assertIn('source_assignments.tsv.gz', result['files'])
        with gzip.open(self.root / 'source_out' / 'units.tsv.gz', 'rt') as handle:
            rows = {row['unit_id']: row for row in csv.DictReader(handle, delimiter='\t')}
        self.assertEqual(rows['u1']['genome_family'], 'protein_coding')
        self.assertEqual(rows['u1']['excluded_from_gene_qc_reason'], 'rRNA_source_panel')
        source.write_text(source.read_text().split('u2')[0])
        with self.assertRaisesRegex(ValueError, 'exactly the analysis FASTQ'):
            rna.run_analysis(sam, self.gtf, table, fastq, 'test', 'clip', self.root / 'missing_source', source_assignments=source)

    def test_umi_coordinate_groups_include_intergenic_without_gene_assignment(self):
        units = self.units(('one', 'two'))
        for unit in units.values():
            unit.update(umi_status='complete', umi_sequence='ACGTACGTAC')
        alignments = {name: [self.record(cigar='50M', start=1500)] for name in units}
        rows, vectors = rna.analyze_units(units, alignments, self.annotation)
        self.assertTrue(all(not row['gene_id'] and row['coordinate_key'] for row in rows))
        summary = rna.summarize_scope('total', rows, vectors, self.annotation, {}, 'rnaseq')
        self.assertEqual(summary['gene_assigned_units'], 0)
        self.assertEqual(summary['umi']['eligible_coordinate_umi_units'], 2)
        self.assertEqual(summary['umi']['observed_exact_coordinate_umi_groups'], 1)
        self.assertEqual(summary['umi']['exact_repeat_units'], 1)
        units['one']['orientation'] = 'unknown'
        rows, _ = rna.analyze_units(units, alignments, self.annotation)
        self.assertEqual(rows[0]['coordinate_key'], '')

    def test_candidate_reference_comparison_never_calls_host_mixture_foreign(self):
        rows = [{'source_status': 'assigned', 'source_class': 'rRNA', 'source_taxon': taxon} for taxon in
                ('Homo sapiens', 'Homo sapiens,Mus musculus', 'Escherichia coli', 'Unknown candidate')]
        result = rna.reference_source_comparison(rows, 'hg38')
        self.assertEqual(result['host_only_candidate_units'], 1)
        self.assertEqual(result['host_and_other_candidate_units'], 1)
        self.assertEqual(result['non_host_only_candidate_units'], 1)
        self.assertEqual(result['unresolved_candidate_units'], 1)
        self.assertEqual(result['non_host_only_fraction_of_all_units'], .25)
        unknown = rna.reference_source_comparison(rows, 'custom_reference')
        self.assertEqual(unknown['status'], 'assembly_taxon_unknown')
        self.assertIsNone(unknown['non_host_only_candidate_units'])
        unavailable = rna.reference_source_comparison([{'source_status': 'not_available'}], 'hg38')
        self.assertEqual(unavailable['status'], 'source_screen_not_available')
        self.assertIsNone(unavailable['non_host_only_candidate_units'])

    def test_cohort_compatibility_gates_and_constant_profile(self):
        paths = []
        for index, (assembly, assay, offset) in enumerate([('a', 'rnaseq', 0), ('a', 'rnaseq', 1), ('b', 'rnaseq', 1), ('a', 'clip', 1)]):
            counts = {f'g{i}': 10 + i + offset for i in range(25)}
            value = {'sample': f's{index}', 'assembly': assembly, 'assay': assay, 'annotation': {'sha256': 'same'},
                     'scopes': {'total': {'gene_counts': counts, 'gene_assigned_units': sum(counts.values())}}}
            path = self.root / f's{index}.json'
            path.write_text(json.dumps(value))
            paths.append(path)
        result = rna.cohort_analysis(paths, self.root / 'cohort')
        self.assertEqual(result['pairs'][0]['status'], 'available')
        self.assertEqual(result['pairs'][1]['status'], 'different_or_unknown_assembly')
        self.assertEqual(result['pairs'][2]['status'], 'different_assay')
        for path in paths[:2]:
            value = json.loads(path.read_text())
            value['scopes']['total']['gene_counts'] = {f'g{i}': 10 for i in range(25)}
            value['scopes']['total']['gene_assigned_units'] = 250
            path.write_text(json.dumps(value))
        constant = rna.cohort_analysis(paths[:2], self.root / 'constant')['pairs'][0]
        self.assertIsNone(constant['log_cpm_pearson'])


def re_sub_biotype(text):
    import re
    return re.sub(r'gene_type "[^"]*";', '', text)


def statistics_mean(values):
    return sum(values) / len(values)


if __name__ == '__main__':
    unittest.main()
