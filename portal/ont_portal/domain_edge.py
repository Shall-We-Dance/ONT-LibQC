"""Optional domain HTTPS and Cloudflare boundary for the existing edge."""
import ipaddress, json, re
from pathlib import Path
import yaml


def configure(config, compose, runtime):
    runtime=Path(runtime);path=runtime/'https.yaml'
    if not path.is_file():return config
    settings=yaml.safe_load(path.read_text())
    domain=settings['domain'];aliases=settings.get('aliases',[])
    for name in [domain,*aliases]:
        if not re.fullmatch(r'[a-z0-9](?:[a-z0-9.-]*[a-z0-9])?',name):raise ValueError('Invalid domain')
    real_ip=restriction=''
    if settings.get('cloudflare'):
        networks=json.loads((runtime/'cloudflare-ips.json').read_text())
        if not networks or len(networks)>100:raise ValueError('Invalid Cloudflare networks')
        networks=[str(ipaddress.ip_network(n)) for n in networks]
        # Never trust this header from an arbitrary Internet client.
        real_ip='real_ip_header CF-Connecting-IP;\n'+''.join('set_real_ip_from '+n+';\n' for n in networks)
        real_ip+='geo $realip_remote_addr $trusted_edge { default 0; 127.0.0.1/32 1; ::1/128 1; '+''.join(n+' 1; ' for n in networks)+'}\n'
        restriction='if ($trusted_edge = 0) { return 403; }'
    cert='domain' if all((runtime/'tls'/('domain.'+ext)).is_file() for ext in ('crt','key')) else 'server'
    config=config.replace('    access_log off;','    '+real_ip+'\n    access_log off;',1)
    config=config.replace('server_name _;',f'server_name {domain};\n        {restriction}',1)
    config=config.replace('/tls/server.crt',f'/tls/{cert}.crt').replace('/tls/server.key',f'/tls/{cert}.key')
    extra=f'''
    server {{
        listen 8080 default_server;
        server_name {domain} {' '.join(aliases)};
        limit_conn client_connections 4;
        limit_req zone=client_requests burst=10 nodelay;
        location ~ ^/\\.well-known/acme-challenge/[A-Za-z0-9_-]+$ {{
            root /acme; default_type text/plain; try_files $uri =404;
        }}
        location / {{ return 308 https://{domain}$request_uri; }}
    }}
    server {{
        listen 8443 ssl default_server;
        server_name _;
        ssl_certificate /tls/{cert}.crt;
        ssl_certificate_key /tls/{cert}.key;
        ssl_protocols TLSv1.2 TLSv1.3;
        return 444;
    }}
'''
    if aliases:
        extra+=f'''
    server {{
        listen 8443 ssl;
        server_name {' '.join(aliases)};
        ssl_certificate /tls/{cert}.crt;
        ssl_certificate_key /tls/{cert}.key;
        ssl_protocols TLSv1.2 TLSv1.3;
        {restriction}
        return 308 https://{domain}$request_uri;
    }}
'''
    ports=['0.0.0.0:443:8443','0.0.0.0:80:8080']
    legacy=settings.get('legacy_port')
    if legacy:
        if type(legacy) is not int or not 1024<=legacy<=65535:raise ValueError('Invalid legacy port')
        ports.append(f'0.0.0.0:{legacy}:9443')
        extra+=f'''
    server {{
        listen 9443 ssl;
        server_name _;
        ssl_certificate /tls/legacy.crt;
        ssl_certificate_key /tls/legacy.key;
        ssl_protocols TLSv1.2 TLSv1.3;
        return 308 https://{domain}$request_uri;
    }}
'''
    config=config.rsplit('}',1)[0]+extra+'}\n'
    compose['services']['edge']['ports']=ports
    compose['services']['edge']['volumes'].append(f'{runtime / "acme-webroot"}:/acme:ro')
    return config
