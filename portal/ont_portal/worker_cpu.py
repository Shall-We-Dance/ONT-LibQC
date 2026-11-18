"""Keep one CPU per running job and lend idle CPUs to every active container."""
import logging
import os
import subprocess


class CpuPool:
    def __init__(self, cpu_ids):
        if (not isinstance(cpu_ids, list) or not 1 <= len(cpu_ids) <= 64
                or any(type(cpu) is not int for cpu in cpu_ids)
                or len(set(cpu_ids)) != len(cpu_ids) or not set(cpu_ids) <= os.sched_getaffinity(0)):
            raise ValueError('shared_cpus must contain unique CPUs available to the worker')
        self.cpus = set(cpu_ids)
        self.owners = {}
        self.applied = {}

    def reserve(self, jid):
        cpu = min(self.cpus.difference(self.owners.values()))
        self.owners[jid] = cpu
        try:
            # Reclaim the CPU from its borrowers before starting the next job.
            for other in self.owners:
                if other != jid:
                    self.update(other)
        except (OSError, subprocess.SubprocessError):
            self.release(jid)
            raise
        # Containers always launch on their reserved CPU, even when setup is slow.
        self.applied[jid] = {cpu}
        return cpu

    def release(self, jid):
        del self.owners[jid]
        self.applied.pop(jid, None)

    def update(self, jid):
        desired = self.cpus.difference(self.owners.values()) | {self.owners[jid]}
        if desired == self.applied.get(jid):
            return
        # An interrupted Docker client may have already applied the change.
        self.applied.pop(jid, None)
        mask = ','.join(map(str, sorted(desired)))
        result = subprocess.run(['docker', 'update', '--cpuset-cpus', mask, 'ont-qc-job-' + jid],
                                capture_output=True, text=True, timeout=5)
        if result.returncode:
            if 'no such container' in result.stderr.lower():
                # The child may still be staging its inputs or may already have finished.
                return
            result.check_returncode()
        self.applied[jid] = desired

    def rebalance(self):
        for jid in self.owners:
            try:
                self.update(jid)
            except (OSError, subprocess.SubprocessError):
                logging.getLogger('ont-worker').warning('CPU allocation will be retried', exc_info=True)


def cpu_limits(config):
    max_jobs = int(config['max_concurrent_jobs'])
    if not 1 <= max_jobs <= 8:
        raise ValueError('Invalid analysis concurrency limit')
    if config['shared_cpus'] != []:
        pool = CpuPool(config['shared_cpus'])
        return len(pool.cpus), min(max_jobs, len(pool.cpus)), pool
    cores = int(config['cores'])
    if not 1 <= cores <= 64:
        raise ValueError('Invalid worker CPU limit')
    return cores, max_jobs, None
