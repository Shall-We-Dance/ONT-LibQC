from dataclasses import dataclass
from pathlib import Path
import json, os
import yaml
from urllib.parse import urlsplit
from .limits import DEFAULTS, load_limits

@dataclass
class Settings:
    state: Path
    public_url: str
    allowed_hosts: list[str]
    secret: str
    submission_key: str
    public_catalog: dict
    max_upload_bytes: int = DEFAULTS['max_upload_bytes']
    max_expanded_bytes: int = DEFAULTS['max_expanded_bytes']
    max_reads: int = DEFAULTS['max_reads']
    max_files: int = DEFAULTS['max_files']
    max_pending: int = DEFAULTS['max_pending']
    max_pending_per_client: int = DEFAULTS['max_pending_per_client']
    max_jobs_per_day: int = DEFAULTS['max_jobs_per_day']
    max_reserved_bytes: int = DEFAULTS['max_reserved_bytes']
    min_free_bytes: int = DEFAULTS['min_free_bytes']
    retention_days: int = DEFAULTS['retention_days']
    upload_hours: int = DEFAULTS['upload_hours']
    max_concurrent_uploads: int = DEFAULTS['max_concurrent_uploads']
    max_http_requests: int = DEFAULTS['max_http_requests']
    max_http_requests_per_client: int = DEFAULTS['max_http_requests_per_client']
    max_resource_requests: int = DEFAULTS['max_resource_requests']
    max_resource_requests_per_client: int = DEFAULTS['max_resource_requests_per_client']
    resource_request_burst: int = DEFAULTS['resource_request_burst']
    resource_requests_per_minute: int = DEFAULTS['resource_requests_per_minute']
    global_resource_request_burst: int = DEFAULTS['global_resource_request_burst']
    global_resource_requests_per_minute: int = DEFAULTS['global_resource_requests_per_minute']
    max_uploads_per_client: int = DEFAULTS['max_uploads_per_client']
    request_burst: int = DEFAULTS['request_burst']
    global_requests_per_minute: int = DEFAULTS['global_requests_per_minute']
    global_request_burst: int = DEFAULTS['global_request_burst']
    requests_per_minute: int = DEFAULTS['requests_per_minute']
    upload_timeout_seconds: int = DEFAULTS['upload_timeout_seconds']
    body_idle_seconds: int = DEFAULTS['body_idle_seconds']
    max_html_bytes: int = DEFAULTS['max_html_bytes']
    pipeline_url: str = DEFAULTS['pipeline_url']
    require_https: bool = True
    mail_enabled: bool = False
    static: Path = Path(__file__).resolve().parents[1] / 'static'
    examples: Path = Path(__file__).resolve().parents[2] / 'example_report'

    example_inputs: Path = Path(__file__).resolve().parents[1] / 'example_inputs'

    def __post_init__(self):
        for name in ['max_upload_bytes', 'max_expanded_bytes', 'max_reads', 'max_files', 'max_pending', 'max_pending_per_client', 'max_jobs_per_day', 'max_reserved_bytes', 'min_free_bytes', 'retention_days', 'upload_hours', 'max_concurrent_uploads', 'max_http_requests', 'max_http_requests_per_client', 'max_uploads_per_client', 'request_burst', 'global_requests_per_minute', 'global_request_burst', 'requests_per_minute', 'upload_timeout_seconds', 'body_idle_seconds', 'max_html_bytes', 'max_resource_requests', 'max_resource_requests_per_client', 'resource_request_burst', 'resource_requests_per_minute', 'global_resource_request_burst', 'global_resource_requests_per_minute']:
            value = getattr(self, name)
            if type(value) is not int or value < (0 if name == 'min_free_bytes' else 1):
                raise ValueError('Invalid positive limit: ' + name)
        if self.max_resource_requests > 256 or self.max_resource_requests_per_client > self.max_resource_requests:
            raise ValueError('Invalid resource request concurrency limit')
        if self.max_files > 64 or self.max_http_requests > 128 or self.max_concurrent_uploads > self.max_http_requests:
            raise ValueError('Invalid batch or request concurrency limit')
        if self.max_http_requests_per_client > self.max_http_requests or self.max_uploads_per_client > min(self.max_concurrent_uploads,self.max_http_requests_per_client):
            raise ValueError('Per-client concurrency must fit global limits')
        if self.max_files > min(self.max_pending, self.max_pending_per_client, self.max_jobs_per_day):
            raise ValueError('Batch size must fit the queue, per-client queue and daily allowance')
        url = urlsplit(self.pipeline_url)
        if url.scheme != 'https' or url.hostname != 'github.com' or url.username or url.password or url.query or url.fragment:
            raise ValueError('pipeline_url must be an HTTPS GitHub repository URL')

    @classmethod
    def load(cls):
        cfg = load_limits(yaml.safe_load(Path(os.environ['PORTAL_CONFIG']).read_text()))
        keys = json.loads(Path(cfg['secrets_file']).read_text())
        catalog = json.loads(Path(cfg['public_catalog_file']).read_text())
        values = {k: v for k, v in cfg.items() if k in cls.__dataclass_fields__}
        # Preserve HTTPS enforcement for existing deployments using the former name.
        if 'require_https' not in cfg and 'secure_cookies' in cfg:
            values['require_https'] = cfg['secure_cookies']
        values.update(state=Path(cfg['state']), secret=keys['signing_key'],
                      submission_key=keys.get('submission_key', ''), public_catalog=catalog)
        for k in ('static', 'examples', 'example_inputs'):
            if k in values: values[k] = Path(values[k])
        result = cls(**values)
        if not result.public_url.startswith('https://') and result.require_https:
            raise ValueError('HTTPS is required for secure public deployments')
        if len(result.secret) < 40 or '*' in result.allowed_hosts:
            raise ValueError('Configure a strong signing key and explicit allowed hosts')
        return result
