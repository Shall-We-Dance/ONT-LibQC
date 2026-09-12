#!/usr/bin/env python3
"""Prepare domain HTTPS, install an ACME certificate, and automate renewal.

Run as the deployment owner, after installing Certbot on the host. Runtime
certificates and ACME account credentials always remain outside this source tree.
"""
import argparse, ipaddress, json, os, re, shlex, shutil, subprocess, sys, urllib.request
from pathlib import Path
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ont_portal.edge import configure
from ont_portal.limits import load_limits


def write(path,data):
    path.write_text(yaml.safe_dump(data,sort_keys=False));path.chmod(0o600)


def regenerate(root):
    cfg=yaml.safe_load((root/('frontend.yaml' if (root/'frontend.yaml').is_file() else 'worker.yaml')).read_text())
    compose=yaml.safe_load((root/'compose.yaml').read_text())
    configure(compose,root,cfg['public_url'],load_limits(cfg));write(root/'compose.yaml',compose)


def certbot_args(root):
    return ['--config-dir',str(root/'acme'),'--work-dir',str(root/'acme-work'),'--logs-dir',str(root/'acme-logs')]


def install(root):
    cfg=yaml.safe_load((root/'https.yaml').read_text());source=root/'acme/live'/cfg['domain'];tls=root/'tls'
    # Validate the key pair and all requested hostnames before replacing anything.
    cert=source/'fullchain.pem';key=source/'privkey.pem'
    import ssl
    context=ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER);context.load_cert_chain(cert,key)
    for name in [cfg['domain'],*cfg.get('aliases',[])]:
        subprocess.run(['openssl','x509','-in',str(cert),'-noout','-checkhost',name],check=True,stdout=subprocess.DEVNULL)
    subprocess.run(['openssl','x509','-in',str(cert),'-noout','-checkend','86400'],check=True,stdout=subprocess.DEVNULL)
    for src,dest in [(key,tls/'domain.key'),(cert,tls/'domain.crt')]:
        tmp=dest.with_suffix(dest.suffix+'.new');shutil.copyfile(src,tmp);tmp.chmod(0o600);tmp.replace(dest)
    regenerate(root)
    docker=cfg.get('docker_command',['docker'])
    for command in (['nginx','-t'],['nginx','-s','reload']):
        subprocess.run(docker+['compose','-f',str(root/'compose.yaml'),'exec','-T','edge',*command],check=True)
    print('Certificate installed and HTTPS edge reloaded.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','issue','install','renew-timer']);p.add_argument('--runtime',required=True)
    p.add_argument('--domain');p.add_argument('--alias',action='append',default=[]);p.add_argument('--cloudflare',action='store_true');p.add_argument('--legacy-port',type=int);p.add_argument('--docker-sudo',action='store_true');p.add_argument('--email')
    args=p.parse_args();root=Path(args.runtime).resolve();os.umask(0o077)
    if args.action=='prepare':
        if not args.domain:raise SystemExit('--domain is required')
        for name in [args.domain,*args.alias]:
            if not re.fullmatch(r'[a-z0-9](?:[a-z0-9.-]*[a-z0-9])?',name):raise SystemExit('Invalid domain')
        if args.legacy_port and not 1024<=args.legacy_port<=65535:raise SystemExit('Invalid legacy port')
        cfg={'domain':args.domain,'aliases':args.alias,'cloudflare':args.cloudflare,'legacy_port':args.legacy_port,'docker_command':['sudo','-n','docker'] if args.docker_sudo else ['docker']}
        if args.cloudflare:
            request=urllib.request.Request('https://api.cloudflare.com/client/v4/ips',headers={'User-Agent':'ONT-LibQC HTTPS setup'})
            with urllib.request.urlopen(request,timeout=20) as r:result=json.loads(r.read(32768))
            if result.get('success') is not True:raise SystemExit('Cannot retrieve official Cloudflare networks')
            networks=[str(ipaddress.ip_network(n)) for n in result['result']['ipv4_cidrs']+result['result']['ipv6_cidrs']]
            if not 5<=len(networks)<=100:raise SystemExit('Unexpected Cloudflare network list')
            (root/'cloudflare-ips.json').write_text(json.dumps(networks,indent=2)+'\n')
        if args.legacy_port:
            for ext in ('crt','key'):
                dest=root/'tls'/('legacy.'+ext)
                if not dest.exists():shutil.copyfile(root/'tls'/('server.'+ext),dest);dest.chmod(0o600)
        (root/'acme-webroot/.well-known/acme-challenge').mkdir(parents=True,exist_ok=True)
        write(root/'https.yaml',cfg);regenerate(root)
        print('Prepared ports 80/443. Recreate the edge, then run issue. Keep port 80 reachable for renewal.')
    elif args.action=='install':install(root)
    elif args.action=='issue':
        cfg=yaml.safe_load((root/'https.yaml').read_text())
        command=['certbot','certonly','--non-interactive','--agree-tos','--webroot','-w',str(root/'acme-webroot'),'--cert-name',cfg['domain'],*certbot_args(root)]
        command+=['--email',args.email] if args.email else ['--register-unsafely-without-email']
        for name in [cfg['domain'],*cfg.get('aliases',[])]:command+=['-d',name]
        subprocess.run(command,check=True);install(root)
    else:
        unit=Path.home()/'.config/systemd/user';unit.mkdir(parents=True,exist_ok=True)
        # systemd command quoting differs from shell quoting; no shell runs here.
        def q(s):return '"'+str(s).replace('\\','\\\\').replace('"','\\"').replace('%','%%')+'"'
        hook=shlex.join([sys.executable,str(Path(__file__).resolve()),'install','--runtime',str(root)])
        command=[shutil.which('certbot') or 'certbot','renew',*certbot_args(root),'--quiet','--deploy-hook',hook]
        (unit/'ont-libqc-certificate.service').write_text('[Unit]\nDescription=ONT-LibQC automatic certificate renewal\n[Service]\nType=oneshot\nExecStart='+' '.join(q(x) for x in command)+'\nUMask=0077\nMemoryMax=192M\nCPUQuota=25%\nTasksMax=48\nTimeoutStartSec=15min\n')
        (unit/'ont-libqc-certificate.timer').write_text('[Unit]\nDescription=ONT-LibQC twice-daily certificate check\n[Timer]\nOnCalendar=*-*-* 00,12:00:00\nRandomizedDelaySec=1h\nPersistent=true\n[Install]\nWantedBy=timers.target\n')
        subprocess.run(['systemctl','--user','daemon-reload'],check=True)
        subprocess.run(['systemctl','--user','enable','--now','ont-libqc-certificate.timer'],check=True)
        print('Automatic certificate renewal enabled; keep user lingering enabled.')

if __name__=='__main__':main()
