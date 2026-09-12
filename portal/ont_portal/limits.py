"""Shared effective defaults and the runtime configuration template. Sizes are bytes."""
DEFAULTS = {
    'max_upload_bytes': 50_000_000,
    'max_expanded_bytes': 250_000_000,
    'max_reads': 250_000,
    'max_files': 16,
    'max_pending': 32,
    'max_pending_per_client': 16,
    'max_jobs_per_day': 128,
    'max_reserved_bytes': 2_000_000_000,
    'min_free_bytes': 10 * 1024**3,
    'retention_days': 7,
    'upload_hours': 3,
    'max_concurrent_uploads': 10,
    'max_http_requests': 16,
    'max_http_requests_per_client': 8,
    'max_uploads_per_client': 1,
    'request_burst': 120,
    'global_requests_per_minute': 1200,
    'global_request_burst': 160,
    'requests_per_minute': 360,
    'upload_timeout_seconds': 600,
    'body_idle_seconds': 20,
    'pipeline_url': 'https://github.com/Shall-We-Dance/ONT-LibQC',
    'max_html_bytes': 4 * 1024**2,
    'max_report_bytes': 100 * 1024**2,
    'max_report_total_bytes': 2 * 1024**3,
    'max_concurrent_jobs': 2,
    'cores': 2,
    'memory': '8g',
    'max_memory': '20g',
    'min_free_memory_bytes': 2 * 1024**3,
    'timeout_seconds': 3600,
    'max_work_bytes': 4 * 1024**3,
    'max_file_bytes': 1024**3,
    'tmp_bytes': 512 * 1024**2,
    'pids_limit': 128,
    'max_log_bytes': 10 * 1024**2,
}


# Keep routine operator settings visible in new deployments. Advanced HTTP,
# temporary-filesystem and logging limits retain the defaults above; operators
# can still add any DEFAULTS key to limits.yaml when an override is needed.
RUNTIME_DEFAULTS = {name: DEFAULTS[name] for name in (
    'max_upload_bytes', 'max_expanded_bytes', 'max_reads', 'max_files',
    'max_pending', 'max_pending_per_client', 'max_jobs_per_day',
    'max_reserved_bytes', 'min_free_bytes', 'retention_days', 'upload_hours',
    'max_concurrent_uploads', 'max_uploads_per_client', 'pipeline_url',
    'max_report_bytes', 'max_report_total_bytes', 'max_concurrent_jobs',
    'cores', 'memory', 'max_memory', 'min_free_memory_bytes', 'timeout_seconds',
    'max_work_bytes', 'max_file_bytes',
)}


def load_limits(config):
    """One optional operator file controls limits for both web and worker."""
    from pathlib import Path
    import yaml
    overrides = {}
    if config.get('limits_file'):
        overrides = yaml.safe_load(Path(config['limits_file']).read_text())
        if not isinstance(overrides, dict) or overrides.keys() - DEFAULTS.keys():
            raise ValueError('limits_file contains unknown resource or quota settings')
    return DEFAULTS | config | overrides
