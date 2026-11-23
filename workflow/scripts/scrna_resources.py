"""Versioned single-cell barcode resources. A run reads only verified local files."""
import gzip
import hashlib
import json
import re
import tempfile
from pathlib import Path
from urllib.request import urlopen

ONT_DATA = 'https://raw.githubusercontent.com/epi2me-labs/wf-single-cell/cb39932e3dfca08ad6dae1e6d08579c957874fcf/data/'
TENX = {
    'v2': '737K-august-2016.txt.gz',
    '3p-v3': '3M-february-2018.txt.gz',
    '3p-v4': '3M-3pgex-may-2023.txt.gz',
    '5p-v3': '3M-5pgex-jan-2023.txt.gz',
}
HASHES = json.loads((Path(__file__).resolve().parents[1] / 'resources' / 'scrna_barcodes.json').read_text())
URLS = {name: ONT_DATA + name for name in TENX.values()}
URLS.update({
    'C4.json': 'https://raw.githubusercontent.com/MGI-tech-bioinformatics/DNBelab_C_Series_scRNA-analysis-software/2b5c226c189ae24c1dc3ebab242e8cdd65a112f9/config/DNBelabC4_scRNA_readStructure.json',
    'v3_onlist.txt': 'https://caltech-bioinformatics-resource-center.github.io/Tutorials/files/v3_onlist.txt',
})


def digest(path):
    sha = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()


def prepare(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name, expected in HASHES.items():
        target = directory / name
        if target.is_file() and digest(target) == expected:
            continue
        with tempfile.NamedTemporaryFile(dir=directory, delete=False) as tmp:
            partial = Path(tmp.name)
        try:
            with urlopen(URLS[name], timeout=120) as response, partial.open('wb') as out:
                while block := response.read(1024 * 1024):
                    out.write(block)
            if digest(partial) != expected:
                raise ValueError('Barcode resource checksum mismatch: ' + name)
            partial.replace(target)
        finally:
            partial.unlink(missing_ok=True)


def load(directory, observed_10x):
    directory = Path(directory)
    for name, expected in HASHES.items():
        if not (directory / name).is_file() or digest(directory / name) != expected:
            raise ValueError(f'Missing or mismatched {name}. Run prepare-barcodes first.')
    rows = [line.split() for line in (directory / 'v3_onlist.txt').read_text().splitlines() if line.strip()]
    parse = tuple({row[i] for row in rows if row[i] != '-'} for i in range(3))
    config = json.loads(re.sub(r',\s*([}\]])', r'\1', (directory / 'C4.json').read_text()))
    mgi = tuple(set(r['white list']) for r in config['cell barcode'])
    # Only exact sequences encountered at the R1 boundary are used for 10x calls.
    # Intersecting while streaming also bounds memory for the largest on-lists.
    needles = {seq.encode('ascii') for seq in observed_10x}
    tenx = {}
    for label, name in TENX.items():
        with gzip.open(directory / name, 'rb') as f:
            tenx[label] = {seq.decode('ascii') for line in f if (seq := line.strip()) in needles}
    return {'parse': parse, 'mgi': mgi, 'tenx': tenx}
