"""Equally weighted gene-body profiles for strand-compatible, non-ribosomal cDNA."""
import rna_qc
from scrna_chemistry import expected_cdna_strand


def summarize(rows, alignments, annotation, tagged, families):
    pooled = rna_qc.GeneProfileCounts()
    scopes = {scope: rna_qc.GeneProfileCounts() for scope in sorted({unit['scope_id'] for unit in tagged.values()})}
    for row, family in zip(rows, families, strict=True):
        name = row['unit_id']
        if family['unit_id'] != name:
            raise ValueError('RNA family and transcript profile input order differs')
        if row['mapping_status'] != 'unique' or not row['eligible_gene'] or family['family'] in ('rRNA', 'mt-rRNA'):
            continue
        unit = tagged[name]
        if row['strand_evidence'] != expected_cdna_strand(unit['kit']):
            continue
        record, = (r for r in alignments[name] if not r['secondary'] and not r['supplementary'])
        gene = annotation.genes[row['gene_id']]
        vectors = rna_qc.transcript_profile_vectors(record['blocks'], gene)
        pooled.add(gene.gene_id, vectors)
        scopes[unit['scope_id']].add(gene.gene_id, vectors)
    return {'annotation_sha256': annotation.metadata['sha256'],
            'scopes': {'all': rna_qc.summarize_gene_profiles(pooled),
                       **{scope: rna_qc.summarize_gene_profiles(counts) for scope, counts in scopes.items()}}}
