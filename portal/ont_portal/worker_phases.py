"""Trusted dispatcher accounting between serial workflow containers."""
import logging
import subprocess

LOG = logging.getLogger('ont-worker')


class PhasePool:
    def __init__(self, cpus):
        self.cpus = cpus
        self.reserved = {}
        self.channels = {}
        self.waiting = {}

    def start(self, jid, budget, channel):
        cpu = self.cpus.reserve(jid) if self.cpus else None
        self.reserved[jid] = budget
        self.channels[jid] = channel
        return cpu

    def release(self, jid):
        self.reserved.pop(jid, None)
        self.waiting.pop(jid, None)
        channel = self.channels.pop(jid, None)
        if channel is not None:
            channel.close()
        if self.cpus and jid in self.cpus.owners:
            self.cpus.release(jid)

    def advance(self, fits):
        # The host child sends only after docker run --rm returned. Collect
        # every release first, including when all three preparations finish.
        for jid, channel in self.channels.items():
            if channel is None or jid in self.waiting or not channel.poll():
                continue
            try:
                budget = channel.recv()
            except EOFError:
                continue
            self.reserved[jid] = 0
            self.waiting[jid] = budget
            if self.cpus:
                self.cpus.release(jid)
        for jid, budget in list(self.waiting.items()):
            if not fits(self.reserved, budget):
                continue
            try:
                cpu = self.cpus.reserve(jid) if self.cpus else None
            except (OSError, subprocess.SubprocessError):
                LOG.warning('Will retry CPU reservation for the next phase of %s', jid)
                continue
            self.reserved[jid] = budget
            try:
                self.channels[jid].send(cpu)
            except (BrokenPipeError, EOFError, OSError):
                # Keep the reservation until the supervisor reaps this child.
                pass
            del self.waiting[jid]
