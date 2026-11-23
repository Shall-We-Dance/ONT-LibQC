"""Operator resources for offline shallow single-cell analyses."""
import hashlib
import json
from functools import lru_cache
from pathlib import Path
from .star_resources import STAR_FILES


@lru_cache(maxsize=256)
def _digest(path, identity):
    # The stat identity invalidates the cache on replacement or an in-place edit.
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def file_digest(path):
    path = Path(path).resolve(strict=True)
    stat = path.stat()
    return _digest(path, (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns))


def barcode_hashes():
    resource = Path(__file__).resolve().parents[2] / 'workflow/resources/scrna_barcodes.json'
    return json.loads(resource.read_text())


def barcodes_ready(directory):
    if not directory:
        return False
    try:
        return all(file_digest(Path(directory) / name) == digest for name, digest in barcode_hashes().items())
    except (OSError, ValueError, TypeError):
        return False


def reference_ready(key, reference):
    """The CLI compares the assembly and the exact source GTF with the index manifest."""
    try:
        metadata = json.loads((Path(reference['star_index']) / 'ontlibqc-star-index.json').read_text())
        if metadata['reference_id'] != key or metadata['source_sha256']['annotation'] != file_digest(reference['annotation']):
            return False
        for field in ('rna_family_gtf', 'repeatmasker'):
            if reference.get(field) and not Path(reference[field]).is_file():
                return False
        if reference.get('rna_source_panel'):
            panel = Path(reference['rna_source_panel'])
            record = json.loads((panel / 'panel.json').read_text())
            if record['fasta_sha256'] != file_digest(panel / 'panel.fa'):
                return False
            entries = record['references'].values()
            if not entries or any(entry.get('source_class') != 'rRNA' or not entry.get('taxon') for entry in entries):
                return False
        return True
    except (AttributeError, KeyError, OSError, TypeError, ValueError):
        return False


def mount_species_panel(key, catalog, jobdir, mounts, index_files):
    panel = catalog.get('species_composition') or {}
    if not panel.get('indexes'):
        return None
    indexes = dict(panel['indexes'])
    aliases = panel.get('assembly_aliases', {})
    for ref_id, row in catalog['references'].items():
        species = aliases.get(ref_id, ref_id)
        if species not in indexes and index_files(row.get('bowtie2_index')):
            indexes[species] = row['bowtie2_index']
    files = {species: paths for species, prefix in indexes.items() if (paths := index_files(prefix))}
    if len(files) < 2 or not (key in files or aliases.get(key) in files):
        return None
    mapped = {}
    for i, (species, parts) in enumerate(files.items()):
        destination = f'/refs/scrna/species/{i}/index'
        basename = Path(indexes[species]).name
        for part in parts:
            mounts.bind(part, destination + part.name[len(basename):])
        mapped[species] = destination
    (jobdir / 'scrna-species.json').write_text(json.dumps({'indexes': mapped, 'assembly_aliases': aliases}) + '\n')
    return '/job/scrna-species.json'


def job_resources(key, reference, catalog, jobdir, mounts, index_files):
    if not barcodes_ready(catalog.get('scrna_barcodes')) or not reference_ready(key, reference):
        raise ValueError('Single-cell RNA reference or barcode resources are unavailable')
    base = '/refs/scrna'
    barcodes = Path(catalog['scrna_barcodes'])
    for name in barcode_hashes():
        mounts.bind(barcodes / name, f'{base}/barcodes/{name}')
    index = Path(reference['star_index'])
    star = base + '/star'
    optional_index_files = ('sjdbInfo.txt', 'sjdbList.out.tab', 'sjdbList.fromGTF.out.tab',
                            'exonGeTrInfo.tab', 'exonInfo.tab', 'geneInfo.tab', 'transcriptInfo.tab')
    for name in (*STAR_FILES, 'ontlibqc-star-index.json', *optional_index_files):
        if (index / name).is_file():
            mounts.bind(index / name, f'{star}/{name}')
    gtf = reference['annotation']
    gtf_suffix = '.gtf.gz' if str(gtf).endswith('.gz') else '.gtf'
    cfg = {
        'assembly': key,
        'star_index': star,
        'gtf': mounts.bind(gtf, base + '/annotation' + gtf_suffix),
        'barcodes': base + '/barcodes',
    }
    for field in ('rna_family_gtf', 'repeatmasker'):
        if reference.get(field):
            suffix = '.gz' if str(reference[field]).endswith('.gz') else ''
            cfg[field] = mounts.bind(reference[field], f'{base}/{field}{suffix}')
    if reference.get('rna_source_panel'):
        cfg['rna_source_panel'] = base + '/rna-source-panel'
        for name in ('panel.fa', 'panel.json'):
            mounts.bind(Path(reference['rna_source_panel']) / name, cfg['rna_source_panel'] + '/' + name)

    species_panel = mount_species_panel(key, catalog, jobdir, mounts, index_files)
    if species_panel:
        cfg['species_panel'] = species_panel
    return cfg
