"""Shared mapper selection for admission, resource budgeting and job configuration."""
RNA_ASSAYS = frozenset({'rnaseq', 'clip'})
BWA_INDEX_KEYS = {'bwa-mem': 'bwa_index', 'bwa-mem2': 'bwa_mem2_index'}


def microc_default(reference):
    # Index availability must not silently change the default aligner.
    return reference.get('microc_mapper') or 'bwa-mem2'


def resolve_mapper(requested, assay, reference):
    if requested != 'auto':
        return requested
    if assay in RNA_ASSAYS:
        return 'starlong'
    return microc_default(reference) if assay == 'microc' else 'bowtie2'


def microc_bwa_mapper(mapper, reference):
    return mapper if mapper in BWA_INDEX_KEYS else microc_default(reference)
