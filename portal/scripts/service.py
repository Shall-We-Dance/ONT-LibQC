#!/usr/bin/env python3
"""Initialize and operate a portable server deployment without changing other services."""
import argparse, ipaddress, json, os, secrets, shutil, subprocess, sys
from pathlib import Path
from urllib.parse import urlsplit
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ont_portal.catalog import public_catalog
from ont_portal.limits import DEFAULTS, RUNTIME_DEFAULTS, load_limits
from ont_portal.edge import configure as configure_edge

ROOT=Path(__file__).resolve().parents[2]

def write(path,data):
    path.write_text(yaml.safe_dump(data,sort_keys=False));path.chmod(0o600)

def init(args):
    root=Path(args.runtime).resolve()
    if root==ROOT or ROOT in root.parents:raise SystemExit('Runtime data must be outside the workflow/source tree')
    if (root/'worker.yaml').exists() or (root/'web.yaml').exists():raise SystemExit('Deployment already initialized. Edit the existing configuration instead.')
    os.umask(0o077);root.mkdir(parents=True,exist_ok=True);root.chmod(0o700)
    url=urlsplit(args.public_url)
    if url.scheme!='https' or not url.hostname or url.path not in ('','/') or url.username or url.query or url.fragment:
        raise SystemExit('Use a root HTTPS origin such as https://server.example:8443')
    if (url.port or 443)!=args.port:raise SystemExit('The public URL port must match --port')
    if args.port!=443 and not 1024<=args.port<=65535:raise SystemExit('Choose port 443 or a port between 1024 and 65535')
    state=root/'state'
    for name in ('control','uploads','reports','work','logs'): (state/name).mkdir(parents=True,exist_ok=True,mode=0o700)
    catalog=json.loads(Path(args.catalog).resolve().read_text());public=catalog if args.frontend_only else public_catalog(catalog)
    if not all(r['ready'] for r in public['references']):raise SystemExit('Prepare references before starting the service')
    if not args.frontend_only:(root/'catalog.admin.json').write_text(json.dumps(catalog,indent=2)+'\n')
    (root/'catalog.public.json').write_text(json.dumps(public,indent=2)+'\n')
    (root/'secrets.json').write_text(json.dumps({'signing_key':secrets.token_hex(32),'submission_key':'' if args.open_submissions else secrets.token_urlsafe(18)},indent=2)+'\n')
    write(root/'limits.yaml',RUNTIME_DEFAULTS)
    common={}
    common.update(public_url=args.public_url.rstrip('/'),allowed_hosts=list(dict.fromkeys([url.hostname,'localhost','127.0.0.1'])),mail_enabled=False,secure_cookies=True)
    write(root/'web.yaml',common|{'limits_file':'/config/limits.yaml','state':'/state','secrets_file':'/config/secrets.json','public_catalog_file':'/config/catalog.public.json','examples':'/examples','static':'/app/static'})
    write(root/('frontend.yaml' if args.frontend_only else 'worker.yaml'),common|{'limits_file':str(root/'limits.yaml'),'state':str(state),'secrets_file':str(root/'secrets.json'),'public_catalog_file':str(root/'catalog.public.json'),
        'catalog_file':str(root/'catalog.admin.json'),'pipeline_dir':str(ROOT),'runner_image':'ont-qc-runner:1',
        'smtp_file':str(root/'smtp.yaml')})
    write(root/'smtp.yaml',{'enabled':False,'host':'smtp.example.org','port':587,'tls_mode':'starttls','from_address':'ONT-LibQC <ont-qc@example.org>','username':'','password':''})
    tls=root/'tls';tls.mkdir(mode=0o700)
    if args.cert and args.key:
        shutil.copyfile(args.cert,tls/'server.crt');shutil.copyfile(args.key,tls/'server.key')
    elif args.cert or args.key:raise SystemExit('Provide both --cert and --key')
    else:
        try:ipaddress.ip_address(url.hostname);san='IP:'+url.hostname
        except ValueError:san='DNS:'+url.hostname
        subprocess.run(['openssl','req','-x509','-newkey','rsa:3072','-nodes','-days','30','-subj','/CN='+url.hostname,
            '-addext','subjectAltName='+san+',DNS:localhost,IP:127.0.0.1','-keyout',str(tls/'server.key'),'-out',str(tls/'server.crt')],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    (tls/'server.key').chmod(0o600)
    volumes=[f'{root / "web.yaml"}:/config/web.yaml:ro',f'{root / "limits.yaml"}:/config/limits.yaml:ro',f'{root / "secrets.json"}:/config/secrets.json:ro',
        f'{root / "catalog.public.json"}:/config/catalog.public.json:ro',f'{tls}:/tls:ro',f'{ROOT / "example_report"}:/examples:ro',
        f'{state / "control"}:/state/control',f'{state / "uploads"}:/state/uploads',f'{state / "reports"}:/state/reports:ro']
    compose={'name':'ont-qc-portal','services':{'web':{'image':'ont-qc-web:1','init':True,'restart':'unless-stopped',
        'user':f'{os.getuid()}:{os.getgid()}','read_only':True,'cap_drop':['ALL'],'security_opt':['no-new-privileges:true'],
        'cpus':0.5,'mem_limit':'256m','memswap_limit':'256m','pids_limit':80,'ulimits':{'nofile':{'soft':512,'hard':512}},'tmpfs':['/tmp:rw,nosuid,nodev,size=64m'],
        'ports':[f'{args.bind}:{args.port}:8443'],'environment':{'PORTAL_CONFIG':'/config/web.yaml'},'volumes':volumes,
        'logging':{'driver':'json-file','options':{'max-size':'5m','max-file':'2'}},
        'healthcheck':{'test':['CMD','python','-c',"import urllib.request,ssl; urllib.request.urlopen('https://127.0.0.1:8443/healthz',context=ssl._create_unverified_context(),timeout=4)"],'interval':'30s','timeout':'5s','retries':3},
        'networks':['frontend']}},'networks':{'frontend':{}}}
    configure_edge(compose,root,args.public_url,DEFAULTS)
    write(root/'compose.yaml',compose)
    if args.frontend_only:
        print('Frontend initialized:',root)
        print('Public URL:',args.public_url)
        print('Connect the remote compute bridge before accepting submissions.')
        return
    python=Path(args.python).absolute()
    if not python.is_file():raise SystemExit('Worker Python interpreter not found')
    def quoted(value):return '"'+str(value).replace('\\','\\\\').replace('"','\\"').replace('%','%%')+'"'
    unit='\n'.join(['[Unit]','Description=ONT-LibQC isolated task dispatcher','After=network.target','',
        '[Service]','Type=simple','WorkingDirectory='+str(ROOT/'portal').replace('%','%%'),
        'Environment='+quoted('PYTHONPATH='+str(ROOT/'portal')),'Environment=PYTHONDONTWRITEBYTECODE=1',
        'ExecStart='+quoted(python)+' -m ont_portal.worker --config '+quoted(root/'worker.yaml'),
        'CPUQuota=50%','MemoryMax=512M','TasksMax=128','Restart=always','RestartSec=5','TimeoutStopSec=45','UMask=0077','Nice=10','NoNewPrivileges=true','',
        '[Install]','WantedBy=default.target',''])
    (root/'ont-qc-portal-worker.service').write_text(unit)
    print('Deployment initialized:',root)
    print('Public URL:',args.public_url)
    print('SMTP is disabled. Reports are collected through private task links.')
    print('Submission code: '+('not required (bounded public submissions)' if args.open_submissions else 'stored in '+str(root/'secrets.json')))

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('init');a.add_argument('--runtime',required=True);a.add_argument('--catalog',required=True);a.add_argument('--public-url',required=True)
    a.add_argument('--port',type=int,default=8443);a.add_argument('--bind',default='0.0.0.0');a.add_argument('--python',default=sys.executable)
    a.add_argument('--frontend-only',action='store_true',help='Host uploads/reports; analysis is dispatched by an SSH bridge');a.add_argument('--open-submissions',action='store_true');a.add_argument('--cert');a.add_argument('--key')
    for name in ('up','down','status','harden'):
        b=sub.add_parser(name);b.add_argument('--runtime',required=True)
    args=p.parse_args();root=Path(args.runtime).resolve()
    if args.action=='init':init(args);return
    if args.action=='harden':
        config=yaml.safe_load((root/('worker.yaml' if (root/'worker.yaml').exists() else 'frontend.yaml')).read_text())
        compose=yaml.safe_load((root/'compose.yaml').read_text())
        configure_edge(compose,root,config['public_url'],load_limits(config))
        write(root/'compose.yaml',compose)
        print('Edge configuration prepared. Recreate both web and edge services to apply it.')
        return
    if args.action=='up' and (root/'frontend.yaml').exists():
        subprocess.run(['docker','compose','-f',str(root/'compose.yaml'),'up','-d'],check=True)
        return
    if args.action=='up':
        subprocess.run(['docker','image','inspect','ont-qc-runner:1','ont-qc-web:1'],check=True,stdout=subprocess.DEVNULL)
        subprocess.run(['systemctl','--user','link',str(root/'ont-qc-portal-worker.service')],check=True)
        subprocess.run(['systemctl','--user','daemon-reload'],check=True)
        subprocess.run(['systemctl','--user','enable','--now','ont-qc-portal-worker.service'],check=True)
        subprocess.run(['systemctl','--user','is-active','--quiet','ont-qc-portal-worker.service'],check=True)
        subprocess.run(['docker','compose','-f',str(root/'compose.yaml'),'up','-d'],check=True)
    elif args.action=='down':
        subprocess.run(['docker','compose','-f',str(root/'compose.yaml'),'down'],check=True)
        if (root/'worker.yaml').exists():subprocess.run(['systemctl','--user','stop','ont-qc-portal-worker.service'],check=True)
    else:
        subprocess.run(['docker','compose','-f',str(root/'compose.yaml'),'ps'],check=True)
        if (root/'worker.yaml').exists():subprocess.run(['systemctl','--user','status','ont-qc-portal-worker.service','--no-pager'],check=False)
if __name__=='__main__':main()
