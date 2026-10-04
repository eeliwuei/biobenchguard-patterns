from pathlib import Path
import importlib.util,json,numpy as np
p=Path(__file__).with_name('run_e1.py');sp=importlib.util.spec_from_file_location('e1',p);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
d,X,y,g,features=m.getdata()
assert len(set(g))==32 and 'name' not in features and 'status' not in features
for seed in [73101,110003,111012]:
 for scheme in ['record_wise','subject_wise']:
  sf=m.splits(y,g,seed,scheme);assert all(np.array_equal(a,c) and np.array_equal(b,e) for (a,b),(c,e) in zip(sf,m.splits(y,g,seed,scheme)))
  seen=np.concatenate([b for a,b in sf]);assert sorted(seen.tolist())==list(range(len(y)))
  for tr,te in sf:
   result=m.audit_manifest({'split_id':'test','target_population':'previously_unseen_subjects','train_subject_ids':g[tr].tolist(),'test_subject_ids':g[te].tolist()})
   if scheme=='subject_wise':assert not set(g[tr])&set(g[te]) and result['overlap_subjects']==0
   else:assert result['overlap_subjects']>0
s=[a for scheme in ['record_wise','subject_wise'] for a in m.splits(y,g,110003,scheme)]
for k in ['subject_label','record_label']:
 yp,_=m.permutation_labels(y,g,510001,k,s)
 if k=='subject_label':assert all(len(set(yp[g==z]))==1 for z in set(g))
 else:assert sorted(yp.tolist())==sorted(y.tolist())
# Mechanism negative and violated positive are independently constructed expectations; no labels needed.
assert m.audit_manifest({'split_id':'independent_unit_fixture','target_population':'previously_unseen_subjects','train_subject_ids':['a','b'],'test_subject_ids':['c','d']})['overlap_subjects']==0
assert m.audit_manifest({'split_id':'independent_unit_fixture','target_population':'previously_unseen_subjects','train_subject_ids':['a','b'],'test_subject_ids':['a','c']})['decision']=='repair_subject_disjoint_split'
print(json.dumps({'status':'PASS','tests':['authentic_data_schema','no_label_or_ID_features','deterministic_splits','complete_disjoint_record_coverage','subject_repair_no_overlap','label_control_invariants','target_contract_unit_fixtures'],'note':'Unit fixtures are engineering tests, not independently scrutinized real negative controls.'}))
