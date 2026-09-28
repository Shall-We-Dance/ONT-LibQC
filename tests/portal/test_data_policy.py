"""No task or upload reservation without an explicit data acknowledgment."""
from pathlib import Path
import tempfile
import unittest
from fastapi.testclient import TestClient
from ont_portal.app import create_app
from ont_portal.settings import Settings

PUBLIC = {'references': [{'id':'dm6','assembly':'dm6','ready':True,'mappers':['bowtie2'],'mycoplasma':False,'spike_in':False,'chromosomes':['chr2L']}]}
BASE = {'email':'test@example.org','filename':'example.fastq','size':20,'consent':True,'options':{'sample':'example','assay':'chromatin','reference':'dm6'}}

class DataPolicyTests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup);self.root=Path(tmp.name)
        self.app=create_app(Settings(state=self.root,public_url='https://testserver',allowed_hosts=['testserver'],secret='test-only-'*8,submission_key='',public_catalog=PUBLIC,min_free_bytes=0))
        self.app.state.store.heartbeat()
        self.client=TestClient(self.app,base_url='https://testserver',headers={'X-Portal-Request':'1'})
        self.addCleanup(self.client.close)

    def test_single_and_mixed_batch_reject_missing_false_or_coerced_acknowledgments(self):
        for value in ('missing',False,None,0,1,'true','false',[],{}):
            bad=dict(BASE)
            if value!='missing':bad['data_policy_acknowledged']=value
            for endpoint in ('/api/jobs','/api/batches'):
                with self.subTest(value=value,endpoint=endpoint):
                    body={'files':[BASE|{'data_policy_acknowledged':True},bad]} if endpoint.endswith('batches') else bad
                    response=self.client.post(endpoint,json=body)
                    self.assertEqual(response.status_code,422,response.text)
                    self.assertIn('data_policy_acknowledged',response.json()['detail'])
                    with self.app.state.store.connect() as db:
                        self.assertEqual(db.execute('SELECT COUNT(*) FROM jobs').fetchone()[0],0)
                        self.assertEqual(db.execute('SELECT COUNT(*) FROM admission_events').fetchone()[0],0)
                    self.assertFalse(list((self.root/'uploads').iterdir()))

    def test_acknowledged_single_and_batch_are_accepted_but_permission_is_still_required(self):
        good=BASE|{'data_policy_acknowledged':True}
        for endpoint in ('/api/jobs','/api/batches'):
            body={'files':[good,good]} if endpoint.endswith('batches') else good
            self.assertEqual(self.client.post(endpoint,json=body).status_code,201)
        response=self.client.post('/api/jobs',json=good|{'consent':False})
        self.assertEqual(response.status_code,422)
        with self.app.state.store.connect() as db:
            self.assertEqual(db.execute('SELECT COUNT(*) FROM jobs').fetchone()[0],3)

if __name__=='__main__':unittest.main()
