"""Run the offline single-cell CLI and package its public deliverables for the portal."""
import json
from pathlib import Path
import shutil
import sys


def scrna_command(root, request, resources, pipeline):
    options = request['options']
    command = [sys.executable, str(pipeline / 'workflow/scripts/ont_qc.py'), 'scrna',
               '--fastq', str(root / 'input.fastq.gz'), '--out', str(root / 'scrna'),
               '--sample', options['sample'], '--threads', str(request['cores']),
               '--mapq', str(options.get('params', {}).get('mapq', 20))]
    for key in ('assembly', 'star_index', 'gtf', 'barcodes'):
        command.extend(['--' + key.replace('_', '-'), resources[key]])
    for key in ('rna_family_gtf', 'rna_source_panel', 'repeatmasker', 'species_panel'):
        if resources.get(key):
            command.extend(['--' + key.replace('_', '-'), resources[key]])
    return command


def package_scrna_report(source, report, sample, raw_reads):
    from report_brand import copy_mark, downloads_page
    from report_bundle import write_json
    from report_export import add_pdf_export
    from report_provenance import RELEASE_VERSION

    if source.is_symlink() or not source.is_dir() or Path(sample).name != sample:
        raise ValueError('Invalid single-cell report source or sample')
    # Preserve the publication boundary even when copying from the preceding phase.
    for path in source.iterdir():
        if path.is_symlink() or not (path.is_file() or path.is_dir()):
            raise ValueError('Unexpected file type in single-cell report output')
    summary = json.loads((source / 'summary.json').read_text())
    if summary['sample'] != sample or summary['raw_reads'] != raw_reads:
        raise ValueError('Single-cell report does not match the submitted sample and read count')
    samples = report / 'samples'
    samples.mkdir(parents=True)
    data = report / 'data/scrna'
    data.mkdir(parents=True)
    assets = report / 'assets'
    copy_mark(assets)

    # The CLI HTML embeds its figures, including per-library plots. Keep its local
    # editor usable offline; the server annotates physical index rows when viewed online.
    footer = ('<footer class="report-footer"><p class="report-download-links">'
              '<a href="../downloads.html">Browse report files</a></p></footer>')
    page = (source / 'report.html').read_text().replace('</body>', footer + '</body>', 1)
    (samples / (sample + '.html')).write_text(add_pdf_export(page))

    # Only these CLI outputs are report artifacts. Never collect reads or work files.
    items = []
    for filename, label, folder in (
        ('summary.json', 'Single-cell QC summary (JSON)', data),
        ('barcode_gene_counts.tsv.gz', 'Barcode and gene UMI counts (TSV.gz)', data),
        ('rna_families.svg', 'RNA family distribution (SVG)', assets),
        ('transcript_regions.svg', 'Coding transcript regions (SVG)', assets),
        ('transcript_metaplot.svg', 'Transcript metaplot (SVG)', assets),
    ):
        path = source / filename
        if folder == assets and not path.is_file():
            continue
        shutil.copyfile(path, folder / filename)
        link = (folder / filename).relative_to(report).as_posix()
        items.append(f'<li><a href="{link}">{label}</a></li>')
    (report / 'downloads.html').write_text(downloads_page(''.join(items)))

    pairs = summary['physical_library']['index_pairing']['pairs']
    catalog = [dict(scope_id=p['scope_id'], scope_name=p['name'], i5=p['i5'], i7=p['i7']) for p in pairs]
    write_json(report / 'data/metrics' / (sample + '.json'), {'sample': sample, 'assay': 'scrna', 'sub_libraries': catalog})
    write_json(report / 'data/summary/provenance.json', {'analyses': {sample: summary['analysis_provenance']}})
    write_json(report / 'manifest.json', {'schema_version': 1, 'report_version': RELEASE_VERSION, 'samples': [sample]})
