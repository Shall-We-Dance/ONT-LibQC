"""Edge configuration must honor operator budgets without dropping isolation."""
import copy
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch
from ont_portal.edge import configure, exempt_report_image_limits
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
            body=re.search(r'map "\$request_method:\$uri:\$http_sec_fetch_dest" \$resource_request \{\n(.*?)\n    \}',config,re.S).group(1)
            def pool(request):
                for line in body.splitlines():
                    expression,value=line.strip().rsplit(' ',1)
                    expression=expression.strip('"')
                    if expression.startswith('~') and re.search(expression[1:],request):
                        return int(value.rstrip(';'))
                return 0
            for request,expected in [
                ('GET:/r/'+'a'*22+':image',2),
                ('HEAD:/r/'+'a'*22+':image',2),
                ('GET:/reports/job/assets/plot.svg:image',2),
                ('GET:/example/assets/plot.png:document',2),
                ('GET:/r/'+'a'*22+':document',1),
                ('GET:/r/'+'a'*22+':',1),
                ('GET:/reports/job/report.zip:document',1),
                ('GET:/example/samples/report.html:document',1),
                ('GET:/static/style.css:style',1),
                ('GET:/example-inputs/HL120.fastq.gz:',1),
                ('POST:/r/'+'a'*22+':image',0),
                ('POST:/reports/job/assets/plot.svg:image',0),
                ('GET:/api/jobs/anything:image',0),
                ('POST:/api/batches:image',0),
            ]:
                with self.subTest(request=request):self.assertEqual(pool(request),expected)
            self.assertIn('map $resource_request $api_client { default $binary_remote_addr; 1 ""; 2 ""; }',config)
            self.assertIn('map $resource_request $api_server { default $server_name; 1 ""; 2 ""; }',config)
            self.assertIn('map $resource_request $resource_client { default ""; 1 $binary_remote_addr; }',config)
            self.assertEqual(exempt_report_image_limits(config),config)

    def test_image_map_migration_preserves_operator_configuration(self):
        prefix='events { worker_connections 1024; }\nhttp {\n'
        suffix='    server { server_name private.example; limit_conn all_connections 32; }\n}\n'
        old=prefix+'''    map "$request_method:$uri" $resource_request {
        default 0;
        ~^(GET|HEAD):/(reports|example|static|example-inputs)/ 1;
    }
    map $resource_request $api_client { default $binary_remote_addr; 1 ""; }
    map $resource_request $api_server { default $server_name; 1 ""; }
    map $resource_request $resource_client { default ""; 1 $binary_remote_addr; }
    map $resource_request $resource_server { default ""; 1 $server_name; }
'''+suffix
        migrated=exempt_report_image_limits(old)
        self.assertTrue(migrated.startswith(prefix));self.assertTrue(migrated.endswith(suffix))
        self.assertIn('~^(GET|HEAD):/(reports|r|example|static|example-inputs)/ 1;',migrated)
        self.assertEqual(exempt_report_image_limits(migrated),migrated)
        with self.assertRaises(ValueError):exempt_report_image_limits(prefix+suffix)

if __name__=='__main__':unittest.main()
