"""Edge configuration must honor operator budgets without dropping isolation."""
import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from ont_portal.edge import configure
from ont_portal.limits import DEFAULTS

class EdgeCapacityTests(unittest.TestCase):
    def test_configured_api_resource_and_submission_limits_are_used(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'tls').mkdir()
            (root/'tls/backend.key').write_text('test-key');(root/'tls/backend.crt').write_text('test-cert')
            compose={'services':{'web':{'user':'1000:1000','volumes':[], 'ports':['127.0.0.1:8443:8443']}},'networks':{}}
            limits=DEFAULTS | {'max_http_requests':32,'max_http_requests_per_client':16,
                'max_resource_requests':192,'max_resource_requests_per_client':96,
                'requests_per_minute':600,'global_requests_per_minute':2400,
                'resource_requests_per_minute':1800,'global_resource_requests_per_minute':9000,
                'submissions_per_minute':12,'submission_burst':6}
            with patch('ont_portal.domain_edge.configure',side_effect=lambda config,*args:config):
                configure(compose,root,'https://test.example.org',limits)
            config=(root/'edge.conf').read_text()
            for expected in ['limit_conn client_connections 16;', 'limit_conn all_connections 32;',
                'limit_conn resource_client_connections 96;', 'limit_conn resource_all_connections 192;',
                'zone=client_requests:1m rate=600r/m;', 'zone=all_requests:64k rate=2400r/m;',
                'zone=resource_client_requests:1m rate=1800r/m;', 'zone=resource_all_requests:64k rate=9000r/m;',
                'zone=submissions:1m rate=12r/m;', 'zone=submissions burst=6 nodelay;',
                'proxy_ssl_verify on;', 'proxy_request_buffering off;', 'limit_rate 2m;']:
                self.assertIn(expected,config)
            self.assertNotIn('@SUBMISSION_',config)
            command=compose['services']['web']['command']
            self.assertGreater(int(command[command.index('--limit-concurrency')+1]),32+192)
            self.assertEqual(compose['networks']['backend'],{'internal':True})
            self.assertNotIn('ports',compose['services']['web'])
            self.assertTrue(compose['services']['edge']['read_only'])

if __name__=='__main__':unittest.main()
