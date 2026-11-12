"""Validate operator-prepared STAR indexes and reserve independent job memory."""
import json
import math
from pathlib import Path

GIB = 1024 ** 3
STAR_VERSION = '2.7.11b'
# STAR records the index format version separately from its executable version.
STAR_GENOME_VERSION = '2.7.4a'
STAR_FILES = ('Genome', 'SA', 'SAindex', 'genomeParameters.txt', 'chrNameLength.txt',
              'chrName.txt', 'chrLength.txt', 'chrStart.txt')
STAR_ARRAYS = ('Genome', 'SA', 'SAindex')


def star_ready(directory):
    if not directory:
        return False
    try:
        root = Path(directory)
        metadata = json.loads((root / 'ontlibqc-star-index.json').read_text())
        if metadata.get('star_version') != STAR_VERSION:
            return False
        parameters = metadata.get('parameters', {})
        if parameters.get('genomeSAsparseD') != 4 or not 1 <= parameters.get('genomeSAindexNbases', 0) <= 12:
            return False
        sizes = metadata.get('index_files', {})
        for name in STAR_FILES:
            path = root / name
            if not path.is_file() or path.stat().st_size <= 0 or sizes.get(name) != path.stat().st_size:
                return False
        parameters_on_disk = {}
        for line in (root / 'genomeParameters.txt').read_text().splitlines():
            fields = line.split()
            if len(fields) > 1:
                parameters_on_disk[fields[0]] = fields[1]
        return (parameters_on_disk.get('versionGenome') == STAR_GENOME_VERSION
                and parameters_on_disk.get('genomeSAsparseD') == '4'
                and int(parameters_on_disk.get('genomeSAindexNbases', 0)) == parameters['genomeSAindexNbases'])
    except (OSError, ValueError, TypeError, AttributeError):
        return False


def star_memory_bytes(reference):
    directory = reference.get('star_index')
    if not star_ready(directory):
        raise ValueError('RNA/CLIP requires a validated administrator-prepared STAR index')
    root = Path(directory)
    metadata = json.loads((root / 'ontlibqc-star-index.json').read_text())
    required = sum((root / name).stat().st_size for name in STAR_ARRAYS) + 3 * GIB
    for declared in (reference.get('star_alignment_memory_bytes'), metadata.get('star_alignment_memory_bytes')):
        if declared is not None:
            if type(declared) is not int or declared <= 0:
                raise ValueError('Invalid STAR alignment memory reservation')
            required = max(required, declared)
    return math.ceil(required / GIB) * GIB
