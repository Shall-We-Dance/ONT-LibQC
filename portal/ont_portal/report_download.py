"""Human-readable archive names without changing private report URLs."""
import datetime
import json
import re
from pathlib import Path
from .schema import SAFE_NAME


def archive_filename(root, sample, created):
    sample = sample if isinstance(sample, str) and re.fullmatch(SAFE_NAME, sample) else 'sample'
    date = datetime.datetime.fromtimestamp(created, datetime.timezone.utc).date()
    root = Path(root).resolve()
    path = root/'data/summary/provenance.json'
    if path.is_file() and path.resolve().is_relative_to(root) and path.stat().st_size <= 8 * 1024 * 1024:
        try:
            provenance = json.loads(path.read_text())
            stamp = (provenance.get('analyses', {}).get(sample, {}).get('analysis_completed_utc')
                     or provenance.get('created_utc'))
            parsed = datetime.datetime.fromisoformat(str(stamp).replace('Z', '+00:00'))
            if parsed.tzinfo is None:
                parsed = parsed.replace(tzinfo=datetime.timezone.utc)
            date = parsed.astimezone(datetime.timezone.utc).date()
        except (ValueError, TypeError, AttributeError, OSError):
            pass
    return f'report_{sample}_{date.isoformat()}.zip'
