"""Small synthetic checkpoint/worker integration; no project responses read."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import run_formal_evaluation as runner

class FormalWorkerIntegration(unittest.TestCase):
    def setUp(self):
        self.config=json.loads((ROOT/'config/formal_model_v1.json').read_text(encoding='utf-8'))
        self.config['feature'].update(expression_genes=6,matched_expression_genes=4,reo_genes=4)
        self.config['elastic_net'].update(alpha_ratios=[1.0],l1_ratios=[1.0])
        self.config['random_forest'].update(n_estimators=3,max_features=['sqrt'],min_samples_leaf=[1],reo_max_pairs=4)
        self.config['ridge_alphas']=[0.0]
        self.ids=pd.DataFrame({'canonical_model_id':[f'm{i:02d}' for i in range(30)],
             'biological_group_id':[f'g{i:02d}' for i in range(30)],'tissue':['A']*15+['B']*15})
        self.folds,_=runner.make_global_split_records(self.ids,self.config)
        self.temp=tempfile.TemporaryDirectory()
        self.out=Path(self.temp.name);(self.out/'checkpoints').mkdir()
        rng=np.random.default_rng(7)
        self.x=np.column_stack([rng.uniform(.5,20,(30,8)),np.repeat([0,1],15),np.arange(30)])
        self.data=self.ids.copy();self.data['expression_row_index']=np.arange(30)
        self.data['drug_name']='synthetic';self.data['response_auc']=rng.uniform(.1,.9,30)
        runner.WORKER.update(run=self.out,config=self.config,freeze_sha='frozen',drug=1931,
              data=self.data,x=self.x,folds=self.folds,stable=tuple(f'gene{i}' for i in range(8)),
              blocked=[],access=set())

    def tearDown(self):
        runner.WORKER.clear();self.temp.cleanup()

    def test_every_model_writes_valid_outer_checkpoint(self):
        for model in self.config['models']:
            task={'drug_id':1931,'repeat':0,'outer_fold':0,'model':model}
            result=runner.execute_task(task)
            self.assertEqual(result['status'],'completed',result.get('error'))
            path=self.out/'checkpoints'/runner.checkpoint_name(task)
            self.assertTrue(runner.validate_resume_checkpoint(path,task,'frozen'))
            saved=json.loads(path.read_text(encoding='utf-8'))
            train=set(saved['expected_training_expression_ids']);test=set(saved['outer_test_expression_ids'])
            self.assertFalse(train & test)
            for fit in saved['fit_events']:
                actual=set(fit['fit_sample_ids'])
                expected=train if fit['stage']=='outer_training_refit' else set(saved['inner_training_expression_ids'][fit['inner_fold']])
                self.assertEqual(actual,expected)
                self.assertFalse(actual & test)
            if model in ['expression','expression_matched','reo']:
                self.assertEqual(saved['result']['outer_fit_audit']['n_nonzero_molecular_coeffs'],0)
            self.assertEqual(runner.find_pending_tasks([task],self.out/'checkpoints','frozen'),[])

    def test_same_repeat_outer_and_inner_full_coverage(self):
        for repeat in range(3):
            block=self.folds[self.folds.repeat==repeat]
            tests=block[(block.inner_fold==-1)&(block.role=='test')]
            self.assertTrue((tests.groupby('canonical_model_id').size()==1).all())
            for fold in range(5):
                task={'repeat':repeat,'outer_fold':fold}
                train,test,inner=runner.local_split_indices(self.data,self.folds,task,3)
                self.assertEqual(sorted(np.r_[train,test].tolist()),list(range(30)))
                covered=np.zeros(len(train),int)
                for it,iv in inner:
                    self.assertFalse(set(it)&set(iv));covered[iv]+=1
                np.testing.assert_array_equal(covered,np.ones(len(train),int))

    def test_worker_error_is_preserved_and_never_successful(self):
        task={'drug_id':1931,'repeat':0,'outer_fold':0,'model':'unsupported'}
        result=runner.execute_task(task)
        self.assertEqual(result['status'],'failed')
        path=self.out/'checkpoints'/runner.checkpoint_name(task)
        self.assertFalse(runner.validate_resume_checkpoint(path,task,'frozen'))
        self.assertEqual(runner.find_pending_tasks([task],self.out/'checkpoints','frozen'),[task])
        with self.assertRaises(ValueError):runner.validate_resume_checkpoint(path,task,'other_freeze')

if __name__=='__main__':unittest.main()
