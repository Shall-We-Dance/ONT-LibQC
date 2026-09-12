#!/usr/bin/env python3
"""Install the host-side notification timer; SMTP credentials stay in the runtime."""
import argparse, os, subprocess, sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',required=True);args=p.parse_args()
    os.umask(0o077);root=Path(args.runtime).resolve();portal=Path(__file__).resolve().parents[1]
    config=root/('frontend.yaml' if (root/'frontend.yaml').exists() else 'worker.yaml')
    if not config.is_file():raise SystemExit('Initialize the deployment first')
    def q(s):return '"'+str(s).replace('\\','\\\\').replace('"','\\"').replace('%','%%')+'"'
    unit=Path.home()/'.config/systemd/user';unit.mkdir(parents=True,exist_ok=True)
    (unit/'ont-libqc-mail.service').write_text('[Unit]\nDescription=ONT-LibQC report notifications\n[Service]\nType=oneshot\nEnvironment=PYTHONDONTWRITEBYTECODE=1\nEnvironment='+q('PYTHONPATH='+str(portal))+'\nExecStart='+q(sys.executable)+' -m ont_portal.notifications --config '+q(config)+'\nUMask=0077\nMemoryMax=96M\nCPUQuota=20%\nTasksMax=16\nTimeoutStartSec=180\nNoNewPrivileges=true\n')
    (unit/'ont-libqc-mail.timer').write_text('[Unit]\nDescription=ONT-LibQC report notification timer\n[Timer]\nOnBootSec=30s\nOnUnitInactiveSec=30s\nAccuracySec=5s\n[Install]\nWantedBy=timers.target\n')
    subprocess.run(['systemctl','--user','daemon-reload'],check=True)
    subprocess.run(['systemctl','--user','enable','--now','ont-libqc-mail.timer'],check=True)
    print('Private host notification timer enabled.')

if __name__=='__main__':main()
