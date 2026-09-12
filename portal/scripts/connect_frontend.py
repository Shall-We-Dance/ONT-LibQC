#!/usr/bin/env python3
"""Connect a frontend-only host to this host's existing compute dispatcher."""
import argparse,json,os,shlex,subprocess,sys
from pathlib import Path
import yaml


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--frontend',required=True,help='Existing trusted SSH host alias');p.add_argument('--runtime',required=True);p.add_argument('--remote-root',default='/srv/ont-libqc');args=p.parse_args()
    os.umask(0o077);root=Path(args.runtime).resolve();remote=Path(args.remote_root)
    if not remote.is_absolute() or any(c in str(remote) for c in '\n\r'):raise SystemExit('Use an absolute remote deployment directory')
    cfg=yaml.safe_load((root/'worker.yaml').read_text());directory=root/'remote-bridge';directory.mkdir(mode=0o700,exist_ok=True)
    key=directory/'id_ed25519'
    if not key.exists():subprocess.run(['ssh-keygen','-q','-t','ed25519','-N','','-C','ont-qc-report-bridge','-f',str(key)],check=True)
    effective=subprocess.check_output(['ssh','-G',args.frontend],text=True,stderr=subprocess.DEVNULL)
    connection=dict(line.split(' ',1) for line in effective.splitlines() if ' ' in line)
    host,user,port=connection['hostname'],connection['user'],connection['port']
    hostkey=subprocess.check_output(['ssh','-T','-o','BatchMode=yes',args.frontend,'cat /etc/ssh/ssh_host_ed25519_key.pub'],text=True).split()
    known=directory/'known_hosts';known.write_text((host if port=='22' else '['+host+']:'+port)+' '+' '.join(hostkey[:2])+'\n')
    wrapper=remote/'run-bridge'
    script='#!/bin/sh\nexec env PYTHONDONTWRITEBYTECODE=1 '+shlex.quote('PYTHONPATH='+str(remote/'source/portal'))+' '+shlex.quote(str(remote/'venv/bin/python'))+' -m ont_portal.frontend_agent --config '+shlex.quote(str(remote/'runtime/frontend.yaml'))+'\n'
    public=key.with_suffix('.pub').read_text().strip()
    authorized='restrict,command="'+str(wrapper).replace('\\','\\\\').replace('"','\\"')+'" '+public
    timer='[Unit]\nDescription=ONT-LibQC frontend retention timer\n[Timer]\nOnBootSec=2min\nOnUnitActiveSec=5min\nPersistent=true\n[Install]\nWantedBy=timers.target\n'
    service='[Unit]\nDescription=ONT-LibQC frontend retention cleanup\n[Service]\nType=oneshot\nEnvironment=PYTHONDONTWRITEBYTECODE=1\nEnvironment=PYTHONPATH='+str(remote/'source/portal')+'\nExecStart='+str(remote/'venv/bin/python')+' -m ont_portal.frontend_agent --config '+str(remote/'runtime/frontend.yaml')+' --maintain\nMemoryMax=128M\nCPUQuota=25%\nNoNewPrivileges=true\nUMask=0077\n'
    install="""import json,os,sys,subprocess
from pathlib import Path
os.umask(0o077);d=json.load(sys.stdin)
p=Path(d['wrapper']);p.write_text(d['script']);p.chmod(0o700)
a=Path.home()/'.ssh/authorized_keys';a.parent.mkdir(mode=0o700,exist_ok=True)
old=a.read_text() if a.exists() else ''
if d['authorized'] not in old:a.write_text(old.rstrip()+'\\n'+d['authorized']+'\\n')
a.chmod(0o600)
u=Path.home()/'.config/systemd/user';u.mkdir(parents=True,exist_ok=True)
(u/'ont-qc-frontend-cleanup.service').write_text(d['service']);(u/'ont-qc-frontend-cleanup.timer').write_text(d['timer'])
subprocess.run(['systemctl','--user','daemon-reload'],check=True)
subprocess.run(['systemctl','--user','enable','--now','ont-qc-frontend-cleanup.timer'],check=True)
"""
    subprocess.run(['ssh','-T','-o','BatchMode=yes',args.frontend,'python3 -c '+shlex.quote(install)],input=json.dumps({'wrapper':str(wrapper),'script':script,'authorized':authorized,'timer':timer,'service':service}),text=True,check=True)
    cfg['ssh_command']=['ssh','-F','/dev/null','-T','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes','-o','IdentitiesOnly=yes','-o','ConnectTimeout=10',
        '-o','ServerAliveInterval=10','-o','ServerAliveCountMax=3','-o','UserKnownHostsFile='+str(known),'-i',str(key),'-p',port,user+'@'+host]
    config=directory/'config.yaml';config.write_text(yaml.safe_dump(cfg,sort_keys=False))
    unit=root/'ont-qc-remote-bridge.service'
    def q(x):return '"'+str(x).replace('\\','\\\\').replace('"','\\"').replace('%','%%')+'"'
    unit.write_text('[Unit]\nDescription=ONT-LibQC frontend transfer bridge\nAfter=network.target ont-qc-portal-worker.service\n[Service]\nType=simple\nEnvironment=PYTHONDONTWRITEBYTECODE=1\nEnvironment='+q('PYTHONPATH='+str(Path(__file__).resolve().parents[1]))+'\nExecStart='+q(sys.executable)+' -m ont_portal.remote_bridge --config '+q(config)+'\nRestart=always\nRestartSec=5\nTimeoutStopSec=130\nMemoryMax=256M\nCPUQuota=50%\nTasksMax=64\nNoNewPrivileges=true\nUMask=0077\n[Install]\nWantedBy=default.target\n')
    print('Restricted SSH connection and resumable bridge configured:',config)
    print('Start after validating both hosts: systemctl --user link '+str(unit))
    print('systemctl --user daemon-reload && systemctl --user enable --now ont-qc-remote-bridge.service')

if __name__=='__main__':main()
