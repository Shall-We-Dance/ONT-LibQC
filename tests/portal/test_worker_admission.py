"""Keep concurrent admission safe across startup, host pressure, and job exits."""
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from ont_portal import worker

GIB=1024**3


class AdmissionTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.worker=worker.Worker.__new__(worker.Worker)
        self.worker.state=Path(self.directory.name)
        self.worker.store=worker.Store(self.worker.state)
        self.worker.max_jobs=3
        self.worker.memory='8g'
        self.worker.public={}
        self.worker.cfg={'max_memory':'24g','min_free_memory_bytes':4*GIB,
            'min_free_bytes':10*GIB,'max_work_bytes':4*GIB,
            'max_report_bytes':100*1024**2,'max_report_total_bytes':2*GIB}
        self.worker.memory_for=lambda option:option.budget*GIB
        options=SimpleNamespace(model_validate_json=lambda raw:SimpleNamespace(
            check_catalog=lambda public:SimpleNamespace(**json.loads(raw))))
        self.options_patch=patch.object(worker,'Options',options)
        self.options_patch.start();self.addCleanup(self.options_patch.stop)
        self.disk_patch=patch.object(worker.shutil,'disk_usage',return_value=SimpleNamespace(free=300*GIB))
        self.disk_patch.start();self.addCleanup(self.disk_patch.stop)
        self.total_patch=patch.object(worker,'total_memory',return_value=31*GIB)
        self.total_patch.start();self.addCleanup(self.total_patch.stop)

    def job(self,number,status='queued',budget=8):
        jid=f'{number:032x}'
        with self.worker.store.connect() as conn:
            conn.execute('INSERT INTO jobs (id,created,updated,status,email,filename,expected_bytes,options,client,expires) VALUES (?,?,?,?,?,?,?,?,?,?)',
                (jid,number,number,status,'','fixture.fastq',100,json.dumps({'budget':budget}),'fixture',9999999999))
        return jid

    def claim(self,available,allocated=0):
        with patch.object(worker,'available_memory',return_value=available*GIB),patch.object(worker,'allocated_job_memory',return_value=allocated*GIB):
            return self.worker.claim()

    def test_three_light_jobs_can_start_before_containers_allocate(self):
        ids=[self.job(i) for i in range(1,5)]
        self.assertEqual([self.claim(30)['id'] for _ in range(3)],ids[:3])
        self.assertIsNone(self.claim(30))
        self.assertEqual(self.worker.store.get(ids[3])['status'],'queued')

    def test_third_job_can_start_after_first_two_charge_anonymous_ram(self):
        self.job(1,'running');self.job(2,'running');third=self.job(3)
        # 16 GiB is already included in MemAvailable's reduction.
        self.assertEqual(self.claim(14,allocated=8)['id'],third)

    def test_unallocated_running_caps_remain_reserved(self):
        self.job(1,'running');self.job(2,'running');self.job(3)
        self.assertIsNone(self.claim(25,allocated=0))

    def test_file_cache_is_not_treated_as_anonymous_credit(self):
        self.job(1,'running');self.job(2,'running');self.job(3)
        # A RAM/cache report of 8 GiB per job cannot be used when only 2 GiB is anon.
        self.assertIsNone(self.claim(20,allocated=2))

    def test_host_pressure_blocks_even_when_three_budgets_fit_total(self):
        self.job(1,'running');self.job(2,'running');self.job(3)
        self.assertIsNone(self.claim(11,allocated=8))

    def test_large_human_and_another_job_do_not_fit_total_headroom(self):
        self.job(1,'running',20);self.job(2)
        self.assertIsNone(self.claim(30,allocated=20))

    def test_cattle_and_one_light_job_fit_but_third_waits(self):
        self.job(1,'running',18);second=self.job(2);self.job(3)
        self.assertEqual(self.claim(30)['id'],second)
        self.assertIsNone(self.claim(30))

    def test_cancelling_job_reservation_is_kept(self):
        self.job(1,'cancelling',20);self.job(2)
        self.assertIsNone(self.claim(30))

    def test_budget_is_released_when_job_finishes(self):
        previous=self.job(1,'running',20);next_id=self.job(2)
        self.assertIsNone(self.claim(30))
        self.worker.store.update(previous,status='succeeded')
        self.assertEqual(self.claim(30)['id'],next_id)

    def test_unavailable_memtotal_fails_closed(self):
        self.job(1)
        with patch.object(worker,'total_memory',return_value=0):
            self.assertIsNone(self.claim(30))


class CgroupCreditTests(unittest.TestCase):
    def read(self,cgroup,stat_text):
        container='a'*64
        result=SimpleNamespace(stdout=container+' 123\n')
        texts={'/proc/123/cgroup':cgroup,
            '/sys/fs/cgroup/system.slice/docker-'+container+'.scope/memory.stat':stat_text,
            '/sys/fs/cgroup/docker/'+container+'/memory.stat':stat_text}
        class FakePath:
            def __init__(self,path):self.path=Path(path)
            @property
            def parts(self):return self.path.parts
            def __truediv__(self,other):return FakePath(self.path/(other.path if isinstance(other,FakePath) else other))
            def read_text(self):return texts[str(self.path)]
        with patch.object(worker.subprocess,'run',return_value=result),patch.object(worker,'Path',FakePath):
            return worker.allocated_job_memory('b'*32)

    def test_only_anonymous_pages_credited_from_systemd_cgroup(self):
        self.assertEqual(self.read('0::/system.slice/docker-'+'a'*64+'.scope\n','anon 1048576\nfile 9000000000\nkernel 100000\n'),1048576)

    def test_cgroupfs_driver_supported(self):
        self.assertEqual(self.read('0::/docker/'+'a'*64+'\n','anon 123\nfile 999\n'),123)

    def test_other_process_cgroup_cannot_receive_credit(self):
        self.assertEqual(self.read('0::/user.slice/unrelated\n','anon 9000000000\n'),0)

    def test_not_started_container_has_no_credit(self):
        with patch.object(worker.subprocess,'run',return_value=SimpleNamespace(stdout='a'*64+' 0\n')):
            self.assertEqual(worker.allocated_job_memory('b'*32),0)

    def test_failed_inspect_and_timeouts_have_no_credit(self):
        for error in (subprocess.CalledProcessError(1,'docker'),subprocess.TimeoutExpired('docker',2),OSError('unavailable')):
            with self.subTest(error=type(error).__name__),patch.object(worker.subprocess,'run',side_effect=error):
                self.assertEqual(worker.allocated_job_memory('b'*32),0)

    def test_invalid_job_id_never_calls_docker(self):
        with patch.object(worker.subprocess,'run') as run:
            self.assertEqual(worker.allocated_job_memory('../other'),0)
            run.assert_not_called()


if __name__=='__main__':unittest.main(verbosity=2)
