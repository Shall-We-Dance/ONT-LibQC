"""Versioned numerical plot data. SVG parsing is a historical-report fallback only."""
import csv
import gzip
import json
import math
from pathlib import Path


def read_gene_profile(path):
    if not Path(path).is_file(): return None
    with Path(path).open() as handle: rows = list(csv.reader(handle, delimiter='\t'))
    if len(rows) == 1 and rows[0][0].startswith('No '): return None
    if len(rows) != 3 or rows[0][:1] != ['bin labels'] or rows[1][:1] != ['bins']:
        raise ValueError('Unsupported gene profile table')
    labels, xs, ys = rows[0][2:], [float(v) for v in rows[1][2:]], [float(v) for v in rows[2][2:]]
    if len(xs) < 2 or len(xs) != len(ys) or len(xs) != len(labels) or any(not math.isfinite(x) for x in xs) or any(b <= a for a, b in zip(xs, xs[1:])):
        raise ValueError('Invalid gene profile axis')
    if any(math.isinf(y) or y < 0 for y in ys): raise ValueError('Invalid gene coverage')
    return labels, xs, ys


def write_data(path, kind, data):
    Path(path).write_text(json.dumps({'schema_version': 1, 'kind': kind, 'data': data}, allow_nan=False, indent=2) + '\n')


def read_data(path, kind):
    value = json.loads(Path(path).read_text())
    if value.get('schema_version') != 1 or value.get('kind') != kind or not isinstance(value.get('data'), dict):
        raise ValueError('Unsupported numerical plot data: ' + str(path))
    return value['data']


def save_profile(matrix, table, output):
    """Stream the same masked mean and population SD / sqrt(rows) as deepTools plot_single.

    Each pipeline matrix contains one sample and one gene group. Reject other
    layouts instead of silently pooling them. Missing bins remain null in JSON.
    """
    gene = read_gene_profile(table)
    if gene is None:
        Path(output).unlink(missing_ok=True)
        return
    labels, xs, expected = gene
    counts, means, sums = [0]*len(xs), [0.0]*len(xs), [0.0]*len(xs)
    rows = 0
    with gzip.open(matrix, 'rt') as handle:
        header = handle.readline()
        if not header.startswith('@'): raise ValueError('Missing computeMatrix header')
        meta = json.loads(header[1:])
        if meta.get('sample_boundaries') != [0, len(xs)] or len(meta.get('group_boundaries', [])) != 2:
            raise ValueError('Expected one profile sample and one gene group')
        for line in handle:
            values = [float(v) for v in line.rstrip().split('\t')[6:]]
            if len(values) != len(xs): raise ValueError('Matrix profile width mismatch')
            rows += 1
            for i, value in enumerate(values):
                if math.isnan(value): continue
                if not math.isfinite(value) or value < 0: raise ValueError('Invalid matrix coverage')
                counts[i] += 1
                delta = value-means[i]
                means[i] += delta/counts[i]
                sums[i] += delta*(value-means[i])
    if not rows or meta['group_boundaries'] != [0, rows]: raise ValueError('Matrix row count mismatch')
    lower, upper, values = [], [], []
    for i, count in enumerate(counts):
        mean = means[i] if count else math.nan
        if not ((math.isnan(mean) and math.isnan(expected[i])) or math.isclose(mean, expected[i], rel_tol=1e-5, abs_tol=1e-6)):
            raise ValueError('Matrix and retained profile means disagree')
        sem = math.sqrt(max(0, sums[i])/count)/math.sqrt(rows) if count else math.nan
        values.append(mean if count else None)
        lower.append(mean-sem if count else None)
        upper.append(mean+sem if count else None)
    write_data(output, 'gene_profile', dict(labels=labels, x=xs, mean=values, lower=lower, upper=upper,
               rows=rows, valid_counts=counts, band='population_sd_over_sqrt_total_rows'))


def read_numeric_profile(path, gene):
    data = read_data(path, 'gene_profile')
    if gene is None: raise ValueError('Numeric profile has no matching TSV')
    labels, xs, means = gene
    if data.get('labels') != labels or data.get('x') != xs: raise ValueError('Numeric profile axis mismatch')
    arrays = [data.get(k) for k in ('mean', 'lower', 'upper')]
    if any(not isinstance(a, list) or len(a) != len(xs) for a in arrays): raise ValueError('Numeric profile width mismatch')
    for expected, mean, low, high in zip(means, *arrays):
        if mean is None:
            if not math.isnan(expected) or low is not None or high is not None: raise ValueError('Invalid missing profile bin')
        elif (not all(isinstance(v, (float, int)) and math.isfinite(v) for v in (mean, low, high))
              or not low <= mean <= high or not math.isclose(mean, expected, rel_tol=1e-5, abs_tol=1e-6)):
            raise ValueError('Invalid profile band or mean')
    convert = lambda a: [math.nan if v is None else v for v in a]
    return (labels, xs, convert(arrays[0])), (convert(arrays[1]), convert(arrays[2]))
