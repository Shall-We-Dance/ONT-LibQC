"""Bound public HTTP connections before parsing requests in the application."""
from pathlib import Path
from urllib.parse import urlsplit
import ipaddress, re, ssl, subprocess


def configure(compose, runtime, public_url, limits):
    runtime = Path(runtime)
    hostname = urlsplit(public_url).hostname or ''
    try:ipaddress.ip_address(hostname)
    except ValueError:
        if not re.fullmatch(r'[A-Za-z0-9.-]+',hostname):raise ValueError('Invalid HTTPS hostname')
    tls=runtime/'tls'
    if not (tls/'backend.crt').exists() and not (tls/'backend.key').exists():
        subprocess.run(['openssl','req','-x509','-newkey','rsa:2048','-nodes','-days','3650','-subj','/CN=ont-qc-backend',
            '-addext','subjectAltName=DNS:ont-qc-backend,IP:127.0.0.1','-keyout',str(tls/'backend.key'),'-out',str(tls/'backend.crt')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        (tls/'backend.key').chmod(0o600)
    if not all((tls/name).is_file() for name in ('backend.crt','backend.key')):raise ValueError('Incomplete backend TLS keypair')
    web = compose['services']['web']
    web['volumes']=[v for v in web['volumes'] if not any(v.endswith(':'+mount+':ro') for mount in ('/tls','/tls/backend.crt','/tls/backend.key'))]
    web['volumes'] += [f'{tls / name}:/tls/{name}:ro' for name in ('backend.crt','backend.key')]
    ports = web.pop('ports',None) or compose['services']['edge']['ports']
    web['networks'] = ['backend']
    # No host port or external network reaches this listener. The sole backend
    # peer is our edge, which overwrites forwarded headers using its socket peer.
    web['command'] = ['python','-m','uvicorn','ont_portal.app:app','--host','0.0.0.0','--port','8443',
        '--proxy-headers','--forwarded-allow-ips','*','--no-access-log','--workers','1',
        '--limit-concurrency','32','--backlog','32','--h11-max-incomplete-event-size','16384',
        '--timeout-keep-alive','3','--timeout-graceful-shutdown','15',
        '--ssl-keyfile','/tls/backend.key','--ssl-certfile','/tls/backend.crt']
    compose['networks']['backend'] = {'internal':True}
    compose['services']['edge'] = {
        'image':'nginx:stable-alpine','entrypoint':['nginx'],'command':['-g','daemon off;'],
        'init':True,'restart':'unless-stopped','user':web['user'],'read_only':True,
        'cap_drop':['ALL'],'security_opt':['no-new-privileges:true'],
        'cpus':0.25,'mem_limit':'64m','memswap_limit':'64m','pids_limit':32,
        'ulimits':{'nofile':{'soft':256,'hard':256}},'tmpfs':['/tmp:rw,nosuid,nodev,size=32m'],
        'ports':ports,'networks':['frontend','backend'],
        'depends_on':{'web':{'condition':'service_healthy'}},
        'volumes':[f'{runtime / "edge.conf"}:/etc/nginx/nginx.conf:ro',f'{runtime / "tls"}:/tls:ro',
                   f'{runtime / "edge-ca.crt"}:/etc/nginx/edge-ca.crt:ro'],
        'logging':{'driver':'json-file','options':{'max-size':'2m','max-file':'2'}},
    }
    certificate=(runtime/'tls/backend.crt').read_text()
    trusted=''.join(ssl.DER_cert_to_PEM_cert(cert) for cert in ssl.create_default_context().get_ca_certs(binary_form=True))
    (runtime/'edge-ca.crt').write_text(certificate+'\n'+trusted)
    (runtime/'edge-ca.crt').chmod(0o600)
    config = r'''worker_processes 1;
pid /tmp/nginx.pid;
# Do not log private report grants, credentials, query strings or uploads.
error_log /dev/stderr crit;
events { worker_connections 128; }
http {
    access_log off;
    server_tokens off;
    client_body_temp_path /tmp/client;
    proxy_temp_path /tmp/proxy;
    fastcgi_temp_path /tmp/fastcgi;
    uwsgi_temp_path /tmp/uwsgi;
    scgi_temp_path /tmp/scgi;
    client_header_timeout 10s;
    client_body_timeout @IDLE@s;
    send_timeout 20s;
    keepalive_timeout 3s;
    keepalive_requests 100;
    reset_timedout_connection on;
    client_header_buffer_size 1k;
    large_client_header_buffers 2 8k;
    client_body_buffer_size 64k;
    client_max_body_size @UPLOAD@;
    limit_conn_zone $binary_remote_addr zone=client_connections:1m;
    limit_conn_zone $server_name zone=all_connections:64k;
    limit_req_zone $binary_remote_addr zone=client_requests:1m rate=@CLIENT_RATE@r/m;
    limit_req_zone $server_name zone=all_requests:64k rate=@GLOBAL_RATE@r/m;
    map "$request_method:$uri" $submission_client {
        default "";
        "POST:/api/jobs" $binary_remote_addr;
        "POST:/api/batches" $binary_remote_addr;
    }
    limit_req_zone $submission_client zone=submissions:1m rate=6r/m;
    server {
        listen 8443 ssl;
        server_name _;
        ssl_certificate /tls/server.crt;
        ssl_certificate_key /tls/server.key;
        ssl_protocols TLSv1.2 TLSv1.3;
        ssl_session_cache shared:TLS:1m;
        ssl_session_timeout 10m;
        ssl_session_tickets off;
        limit_conn client_connections @CLIENT_CONNECTIONS@;
        limit_conn all_connections 32;
        limit_conn_status 429;
        limit_req zone=client_requests burst=@CLIENT_BURST@ nodelay;
        limit_req zone=all_requests burst=@GLOBAL_BURST@ nodelay;
        limit_req zone=submissions burst=4 nodelay;
        limit_req_status 429;
        # Error bodies do not reflect request input. Application errors pass through.
        error_page 429 = @busy;
        error_page 502 503 504 = @unavailable;
        location @busy {
            default_type application/json;
            add_header Retry-After 10 always;
            add_header Cache-Control no-store always;
            add_header X-Content-Type-Options nosniff always;
            add_header Referrer-Policy no-referrer always;
            return 429 '{"detail":"Too many requests. Please retry shortly."}';
        }
        location @unavailable {
            default_type application/json;
            add_header Retry-After 10 always;
            add_header Cache-Control no-store always;
            add_header X-Content-Type-Options nosniff always;
            add_header Referrer-Policy no-referrer always;
            return 503 '{"detail":"The service is busy. Please retry shortly."}';
        }
        location / {
            proxy_pass https://web:8443;
            proxy_http_version 1.1;
            proxy_set_header Host $http_host;
            proxy_set_header Connection "";
            proxy_set_header X-Forwarded-For $remote_addr;
            proxy_set_header X-Forwarded-Proto https;
            proxy_set_header Forwarded "";
            proxy_set_header X-Real-IP "";
            proxy_ssl_server_name on;
            proxy_ssl_name @HOST@;
            proxy_ssl_verify on;
            proxy_ssl_verify_depth 5;
            proxy_ssl_trusted_certificate /etc/nginx/edge-ca.crt;
            proxy_connect_timeout 5s;
            proxy_send_timeout @IDLE@s;
            proxy_read_timeout 30s;
            proxy_request_buffering off;
            proxy_buffering off;
            proxy_max_temp_file_size 0;
            proxy_next_upstream off;
            # Bound parallel report downloads without writing bodies to proxy disk.
            limit_rate_after 1m;
            limit_rate 2m;
        }
    }
}
'''
    values={'HOST':'ont-qc-backend','UPLOAD':max(limits['max_upload_bytes'],4*1024**2,limits['max_files']*65536),
        'IDLE':limits['body_idle_seconds'],'CLIENT_RATE':limits['requests_per_minute'],
        'GLOBAL_RATE':limits['global_requests_per_minute'],'CLIENT_BURST':limits['request_burst'],
        'GLOBAL_BURST':limits['global_request_burst'],'CLIENT_CONNECTIONS':limits['max_http_requests_per_client']}
    for key,value in values.items():config=config.replace('@'+key+'@',str(value))
    from .domain_edge import configure as configure_domain
    config=configure_domain(config,compose,runtime)
    (runtime/'edge.conf').write_text(config)
    (runtime/'edge.conf').chmod(0o600)
    return compose
