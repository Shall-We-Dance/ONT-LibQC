"""Stable, path-free pipeline identity and analysis timestamps for reports."""
import datetime
import hashlib
from pathlib import Path


def capture_analysis(root=None, now=None):
    root = Path(root) if root is not None else Path(__file__).resolve().parents[2]
    paths = [root / 'Snakefile']
    paths += sorted(p for p in (root / 'workflow').rglob('*')
                    if p.is_file() and p.suffix in {'.py', '.R', '.smk', '.json', '.yaml', '.yml', '.css', '.js', '.html'})
    digest = hashlib.sha256()
    for path in paths:
        if not path.is_file():
            continue
        name = path.relative_to(root).as_posix().encode()
        blob = path.read_bytes()
        digest.update(len(name).to_bytes(8, 'big') + name + len(blob).to_bytes(8, 'big') + blob)
    code_hash = digest.hexdigest()
    timestamp = now or datetime.datetime.now(datetime.timezone.utc)
    return {'pipeline_version': 'build.' + code_hash[:12], 'pipeline_sha256': code_hash,
            'analysis_completed_utc': timestamp.astimezone(datetime.timezone.utc).isoformat()}


def analysis_record(metric, provenance=None):
    provenance = provenance or {}
    recorded = metric.get('analysis_provenance') or provenance.get('analyses', {}).get(metric.get('sample'), {})
    if recorded:
        return dict(recorded)
    # Historical data does not identify the source revision. Never label it as
    # the current build merely because the report presentation was refreshed.
    return {'pipeline_version': None, 'pipeline_sha256': None,
            'analysis_completed_utc': provenance.get('created_utc')}


def analysis_label(record):
    version = record.get('pipeline_version') or 'Version unrecorded'
    timestamp = record.get('analysis_completed_utc')
    try:
        parsed = datetime.datetime.fromisoformat(str(timestamp).replace('Z', '+00:00'))
        if parsed.tzinfo is None:
            parsed = parsed.replace(tzinfo=datetime.timezone.utc)
        date = parsed.astimezone(datetime.timezone.utc).strftime('%Y-%m-%d UTC')
    except (ValueError, TypeError):
        date = 'Date unrecorded'
    return f'{version} · {date}'
