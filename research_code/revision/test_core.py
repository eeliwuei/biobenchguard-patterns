import numpy as np
from core import seed,generate,corr_rank,scope_stat,split_stat,threshold,wilson,provenance_status

def test_finite_rank_boundary_and_ties():
 scores=np.arange(20.)
 rejects=sum(scores[i]>threshold(np.delete(scores,i),.05) for i in range(20))
 assert rejects==1
 assert not 0>threshold(np.zeros(20),.05)
 assert np.isinf(threshold(np.ones(2),.05))

def test_generation_axes_are_separate():
 X=np.random.default_rng(11).normal(size=(100,30));g=np.repeat(np.arange(5),20);c={'task':'regression','signal':.5}
 a,y,m=generate(X,g,c,1,0,0,False);b,y2,m2=generate(X,g,c,1,400,0,True)
 assert np.array_equal(y,y2) and np.array_equal(a,b[:,:30])
 _,y3,m3=generate(X,g,c,1,0,.8,False)
 assert np.isclose(m3['group_signal_sd'],.8) and not np.array_equal(y,y3)
 assert seed('calibration',1)!=seed('evaluation',1)

def test_fold_internal_selector_does_not_read_test_labels():
 rng=np.random.default_rng(2);X=rng.normal(size=(100,30));y=rng.normal(size=100);tr=np.arange(70);yy=y.copy();yy[70:]+=1000
 assert np.array_equal(corr_rank(X[tr],y[tr],10),corr_rank(X[tr],yy[tr],10))

def test_observed_scope_different_from_xy_contrast():
 rng=np.random.default_rng(3);X=rng.normal(size=(60,20));y=rng.normal(size=60)
 a=scope_stat(X,y,return_predictions=True);b=scope_stat(X,y)
 assert a['T_scope']==b['T_scope']
 good=[{'fit_ids':r['B_selector_fit_ids'],'test_ids':r['test_ids']} for r in a['fitlog']]
 bad=[{'fit_ids':r['C_selector_fit_ids'],'test_ids':r['test_ids']} for r in a['fitlog']]
 assert provenance_status(good)=='NO_PROHIBITED_OVERLAP_IN_RECORDED_SCOPE'
 assert provenance_status(bad)=='PROHIBITED_FIT_SCOPE_OBSERVED'
 assert provenance_status(None)=='INSUFFICIENT_PIPELINE_PROVENANCE'
 for r in a['fitlog']:assert not set(r['transform_fit_ids'])&set(r['test_ids'])

def test_group_folds_equal_train_sizes_and_valid_predictions():
 rng=np.random.default_rng(4);X=rng.normal(size=(100,20));y=rng.normal(size=100);g=np.repeat(np.arange(5),20)
 r=split_stat(X,y,g,'regression','linear',10,32)
 assert len(r['fold_sizes'])==5 and np.isfinite(r['T_split'])

def test_binomial_intervals_contain_observed_rate():
 for k in [0,1,5,50,100]:
  lo,hi=wilson(k,100);assert lo-1e-15<=k/100<=hi+1e-15


def test_candidate_nesting_and_model_pairing():
 X=np.random.default_rng(12).normal(size=(50,20));g=np.repeat(np.arange(5),10);c={'task':'regression','signal':.5}
 a,y,_=generate(X,g,c,1,400,0,False);b,yy,_=generate(X,g,c,1,2500,0,False)
 assert np.array_equal(a,b[:,:420]) and np.array_equal(y,yy)
 from experiments import jobs
 config={'conditions':[],'mixed':{'replicates':1,'L':[0],'eta':[0]}}
 j=jobs(config);aa=[x for x in j if x['condition']['base']=='wheat' and not x['encoding']]
 assert aa[0]['instance_seed']==aa[1]['instance_seed'] and aa[0]['split_seed']==aa[1]['split_seed'] and aa[0]['instance_id']!=aa[1]['instance_id']
