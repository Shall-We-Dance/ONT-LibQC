"""Describe trusted workflow milestones without returning log contents."""
import json
from pathlib import Path
import re

RULE_LABELS = {
    'prepare_library_units': 'Recovering inserts and resolving sublibraries',
    'rna_prepare': 'Preparing RNA inserts and removing technical sequences',
    'align': 'Mapping inserts to the reference genome',
    'mycoplasma_screen': 'Checking for Mycoplasma sequences',
    'spikein_screen': 'Checking spike-in and other reference sequences',
    'rna_sources': 'Identifying RNA families and source genomes',
    'rna_qc': 'Calculating RNA regions, strand bias and splice evidence',
    'microc_pairtools': 'Classifying Micro-C contact pairs',
    'microc_full_insert_ligation': 'Analyzing Micro-C ligation fragments',
    'sample_metrics': 'Calculating library QC metrics',
    'reference_gene_bed': 'Preparing gene annotations',
    'chromatin_bigwig': 'Calculating genome coverage',
    'chromatin_gene_profile': 'Plotting total-library gene coverage',
    'resolved_report_config': 'Preparing the report',
    'reports': 'Preparing plots and the report',
    'report_bundle': 'Packaging report files',
    'all': 'Packaging report files',
}
RULE = re.compile(r'^(?:localrule|rule) ([a-z_]+):\s*$', re.MULTILINE)
SAFE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]{0,99}$')

def report_stage(work):
    root = Path(work) / '.work'
    expected_gene, expected_chromosome = [], []
    try:
        for path in (root / 'qc/metrics').glob('*.json'):
            if path.is_symlink() or path.stat().st_size > 8 * 1024**2:
                continue
            metric = json.loads(path.read_text())
            sample = metric.get('sample', '')
            if not isinstance(sample, str) or not SAFE.fullmatch(sample):
                continue
            scopes = [row.get('scope_id') for row in metric.get('sub_libraries', [])]
            scopes = [scope for scope in scopes if isinstance(scope, str) and SAFE.fullmatch(scope)]
            if metric.get('assay') == 'chromatin':
                expected_gene.extend(root / 'reports/assets' / f'{sample}.{scope}.TSS_TES.svg' for scope in scopes)
            expected_chromosome.extend(root / 'reports/assets' / f'{sample}.{scope}.karyotype.svg' for scope in ['total', *scopes])
        for paths, label in ((expected_gene, 'Plotting sublibrary gene coverage'),
                             (expected_chromosome, 'Plotting chromosome coverage')):
            done = sum(path.is_file() and path.stat().st_size > 0 for path in paths)
            if done < len(paths):
                return f'{label} ({done} of {len(paths)} complete)'
    except (OSError, ValueError, TypeError, AttributeError):
        # A metric/asset may still be in the process of being written.
        pass
    return RULE_LABELS['reports']

def analysis_stage(work, log_tail):
    rules = [name for name in RULE.findall(log_tail) if name in RULE_LABELS]
    if not rules:
        return 'Validating input reads'
    rule = rules[-1]
    return report_stage(work) if rule == 'reports' else RULE_LABELS[rule]
