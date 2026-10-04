#!/usr/bin/env python3
"""Controlled replay of a published identity-confounding workflow, not exact R replication."""
from pathlib import Path
import argparse, json, hashlib, datetime, time, platform, sys
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

BASE=Path(__file__).resolve().parents[3]
CODE=Path(__file__).resolve()

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,obj): Path(p).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def getdata():
    p=BASE/'data_manifest/raw/e1/parkinsons.data';d=pd.read_csv(p)
    g=d.name.str.split('_').str[2].to_numpy();y=d.status.to_numpy(dtype=int)
    features=[c for c in d.columns if c not in ['name','status']]
    X=d[features].to_numpy(float)
    assert X.shape==(195,22) and len(set(g))==32 and not np.isnan(X).any()
    assert d.name.is_unique and all(len(set(y[g==z]))==1 for z in set(g))
    assert sum(y)==147 and sum([y[g==z][0] for z in set(g)])==24
    return d,X,y,g,features

def models(name,seed):
    if name=='random_forest':return RandomForestClassifier(n_estimators=500,max_features='sqrt',min_samples_leaf=1,bootstrap=True,n_jobs=1,random_state=seed)
    return make_pipeline(StandardScaler(),LogisticRegression(C=1.0,max_iter=5000,solver='lbfgs',random_state=seed))

def splits(y,g,seed,scheme):
    rng=np.random.default_rng(seed)
    if scheme=='record_wise':tr=np.sort(rng.choice(len(y),round(len(y)/2),replace=False))
    else:
        # Original author logic: sample half of cases and half of controls at SUBJECT level.
        ids=np.unique(g);trainids=[]
        for c in [0,1]:
            ids_c=np.array([z for z in ids if y[g==z][0]==c]);trainids.extend(rng.choice(ids_c,round(len(ids_c)/2),replace=False))
        tr=np.flatnonzero(np.isin(g,trainids))
    te=np.setdiff1d(np.arange(len(y)),tr)
    assert len(np.intersect1d(tr,te))==0
    return [(tr,te),(te,tr)]

def audit_manifest(case):
    """Only declared pipeline records and subject IDs enter. No expected status or outcome labels."""
    assert 'expected_status' not in case and 'reference_truth' not in case
    tr=set(case['train_subject_ids']);te=case['test_subject_ids'];overlap=tr.intersection(te)
    frac=float(np.mean([z in tr for z in te]))
    breach=case['target_population']=='previously_unseen_subjects' and len(overlap)>0
    return {'split_id':case['split_id'],'overlap_subjects':len(overlap),'test_record_identity_overlap_fraction':frac,
            'decision':'repair_subject_disjoint_split' if breach else 'no_cross_subject_overlap_for_this_mechanism',
            'calibration':'deterministic zero-overlap target contract; no statistical calibration',
            'evidence_kind':'observed_pipeline_information_path'}

def permutation_labels(y,g,seed,kind,all_splits):
    rng=np.random.default_rng(seed);ids=np.unique(g);gy=np.array([y[g==z][0] for z in ids]);attempts=0
    while True:
        attempts+=1
        if kind=='subject_label':
            yp=gy[rng.permutation(len(ids))][np.searchsorted(ids,g)]
        else:yp=rng.permutation(y)
        if all(len(set(yp[a]))==2 and len(set(yp[b]))==2 for a,b in all_splits):return yp,attempts
        assert attempts<10000

def fit_pair(X,y,g,scheme,seed,name,rep):
    pred=np.empty(len(y));folds=[]
    for k,(tr,te) in enumerate(splits(y,g,seed,scheme)):
        m=models(name,seed+17*k);m.fit(X[tr],y[tr]);p=m.predict_proba(X[te])[:,1];pred[te]=p
        case={'split_id':f'{scheme}_r{rep:02d}_f{k}','target_population':'previously_unseen_subjects',
              'train_subject_ids':g[tr].tolist(),'test_subject_ids':g[te].tolist()}
        a=audit_manifest(case)
        folds.append({'model':name,'repeat':rep,'seed':seed,'scheme':scheme,'fold':k,'split_id':case['split_id'],
          'n_train':len(tr),'n_test':len(te),'n_train_subjects':len(set(g[tr])),'n_test_subjects':len(set(g[te])),
          'auc':float(roc_auc_score(y[te],p)),'balanced_accuracy':float(balanced_accuracy_score(y[te],p>=.5)),
          **a,'train_rows':tr.tolist(),'test_rows':te.tolist(),'case':case})
    return pred,folds

def bootstrap_summary(predictions,y,g,B,seed):
    # Cluster bootstrap is conditional on fitted OOF predictions; it does not refit models.
    ids=np.unique(g);gy=np.array([y[g==z][0] for z in ids]);pos=ids[gy==1];neg=ids[gy==0]
    rng=np.random.default_rng(seed)
    wp=np.array([np.bincount(rng.integers(len(pos),size=len(pos)),minlength=len(pos)) for _ in range(B)])
    wn=np.array([np.bincount(rng.integers(len(neg),size=len(neg)),minlength=len(neg)) for _ in range(B)])
    prow=np.searchsorted(pos,g[y==1]);nrow=np.searchsorted(neg,g[y==0]);Wpos=wp[:,prow];Wneg=wn[:,nrow]
    result={};draws={}
    for (name,scheme),P in predictions.items():
        wins=np.mean((P[:,y==1,None]>P[:,None,y==0])+.5*(P[:,y==1,None]==P[:,None,y==0]),axis=0)
        aucdraw=np.sum((Wpos@wins)*Wneg,axis=1)/(Wpos.sum(axis=1)*Wneg.sum(axis=1))
        aucs=np.array([roc_auc_score(y,p) for p in P]);ba=np.array([balanced_accuracy_score(y,p>=.5) for p in P])
        subjectP=np.stack([P[:,g==z].mean(axis=1) for z in ids],axis=1)
        subjectaucs=[roc_auc_score(gy,p) for p in subjectP]
        key=name+'__'+scheme;draws[key]=aucdraw
        result[key]={'model':name,'scheme':scheme,'record_weighted_oof_auc_mean':float(aucs.mean()),
          'cluster_ci_lower':float(np.quantile(aucdraw,.025)),'cluster_ci_upper':float(np.quantile(aucdraw,.975)),
          'split_auc_sd':float(aucs.std(ddof=1)) if len(aucs)>1 else 0.,
          'split_auc_min':float(aucs.min()),'split_auc_max':float(aucs.max()),
          'balanced_accuracy_mean':float(ba.mean()),'subject_aggregated_auc_mean':float(np.mean(subjectaucs)),
          'n_rows':len(y),'n_subjects':len(ids),'n_repeats':P.shape[0],'n_bootstrap':B,
          'ci_scope':'stratified subject bootstrap of fixed OOF predictions; conditional, no model refitting',
          'within_group_permutation':'not identifiable: status constant within every subject'}
    for name in ['random_forest','logistic_regression']:
        a=name+'__record_wise';b=name+'__subject_wise';d=draws[a]-draws[b]
        result[name+'__gap']={'model':name,'statistic':'record_minus_subject_oof_auc',
            'value':result[a]['record_weighted_oof_auc_mean']-result[b]['record_weighted_oof_auc_mean'],
            'ci_lower':float(np.quantile(d,.025)),'ci_upper':float(np.quantile(d,.975)),
            'interpretation':'predictive split sensitivity, not a verified error rate or sole proof of misuse'}
    return result,draws

def run(stage):
    start=time.time();d,X,y,g,features=getdata();protocol=json.loads((BASE/'protocol/e1/protocol.json').read_text())
    if stage=='main':
        frozen=json.loads((BASE/'protocol/e1/freeze.json').read_text())
        assert frozen['script_sha256']==sha(CODE) and frozen['protocol_sha256']==sha(BASE/'protocol/e1/protocol.json')
    R=1 if stage=='pilot' else protocol['repeats'];B=100 if stage=='pilot' else protocol['cluster_bootstrap']
    nperm=2 if stage=='pilot' else protocol['permutation_replicates_per_scheme']
    raw=BASE/'results_raw/e1'/stage;raw.mkdir(parents=True,exist_ok=True)
    summary=BASE/'results_summary/e1'/stage;summary.mkdir(parents=True,exist_ok=True)
    rowlist=[];foldlist=[];case_inputs=[];pipeline_provenance=[];detector=[];truth=[];allpred={}
    seeds=[73101] if stage=='pilot' else protocol['main_split_seeds']
    for name in ['random_forest','logistic_regression']:
        for scheme in ['record_wise','subject_wise']:
            P=[]
            for rep,seed in enumerate(seeds):
                pred,folds=fit_pair(X,y,g,scheme,seed,name,rep);P.append(pred)
                for f in folds:
                    tr=f.pop('train_rows');te=f.pop('test_rows');case=f.pop('case');foldlist.append(f)
                    if name=='random_forest':
                        case_inputs.append(case);detector.append(audit_manifest(case))
                        pipeline_provenance.append({'split_id':case['split_id'],'train_row_ids':d.name.iloc[tr].tolist(),'test_row_ids':d.name.iloc[te].tolist(),
                          'scaling_scope':'training rows only for logistic; no scaling for forest','feature_selection':'none; 22 original voice features',
                          'hyperparameter_scope':'fixed before main, no tuning','group_ID_feature':False})
                        truth.append({'split_id':case['split_id'],'reference_expected_path':'overlap' if scheme=='record_wise' else 'subject_disjoint',
                          'basis':'published record-wise mechanism and independently inspectable saved train/test assignments; not inferred from AUC'})
                    for i in te:rowlist.append({'model':name,'scheme':scheme,'repeat':rep,'fold':f['fold'],'seed':seed,'row_index':int(i),'record_id':d.name.iloc[i],
                        'subject_id':g[i],'y_true':int(y[i]),'pred_probability':float(pred[i])})
                print(json.dumps({'stage':stage,'part':'real_labels','model':name,'scheme':scheme,'repeat':rep,'elapsed_s':round(time.time()-start,1)}),flush=True)
            allpred[(name,scheme)]=np.array(P)
    pd.DataFrame(rowlist).to_csv(raw/'oof_predictions.csv',index=False)
    pd.DataFrame(foldlist).to_csv(raw/'fold_metrics.csv',index=False)
    dump(raw/'case_inputs.json',case_inputs);dump(raw/'pipeline_provenance.json',pipeline_provenance)
    dump(raw/'reference_truth_private.json',truth);dump(raw/'detector_predictions.json',detector)
    stats,draws=bootstrap_summary(allpred,y,g,B,protocol['bootstrap_seed'])
    dump(summary/'performance_summary.json',stats)
    pd.DataFrame(draws).to_csv(raw/'subject_bootstrap_auc_draws.csv',index=False)
    pd.DataFrame([v for k,v in stats.items() if not k.endswith('__gap')]).to_csv(summary/'real_failure_recovery.csv',index=False)
    dump(summary/'evaluation_comparison.json',{'mechanism_audit':{'flagged_record_splits':sum(x['decision']=='repair_subject_disjoint_split' for x in detector),'total_splits':len(detector)},
      'predictive_gap':{k:v for k,v in stats.items() if k.endswith('__gap')},
      'legacy_data_only_statistic':'X,y,g are identical for record and repaired cases. Recomputing a random-minus-group diagnostic from those inputs yields the same value for both case descriptors and cannot establish which actual pipeline ran.',
      'independent_real_negative':'NOT RECOVERED; repaired subject-disjoint control is local and mechanism-specific only'})
    # Fixed first split, fresh held-out seeds. Outcomes below are diagnostic controls, never substituted for real-label results.
    nullrows=[];nullpred=[];seed=seeds[0]
    fixed_splits={s:splits(y,g,seed,s) for s in ['record_wise','subject_wise']};both=[v for ss in fixed_splits.values() for v in ss]
    for kind in ['subject_label','record_label']:
        for b in range(nperm):
            pseed=protocol['permutation_seed']+10000*(kind=='record_label')+b
            yp,attempts=permutation_labels(y,g,pseed,kind,both)
            for name in ['random_forest','logistic_regression']:
                for scheme in ['record_wise','subject_wise']:
                    pred=np.empty(len(y));fas=[]
                    for fold,(tr,te) in enumerate(fixed_splits[scheme]):
                        m=models(name,seed+17*fold);m.fit(X[tr],yp[tr]);pred[te]=m.predict_proba(X[te])[:,1];fas.append(float(roc_auc_score(yp[te],pred[te])))
                    nullrows.append({'null_scheme':kind,'permutation':b,'permutation_seed':pseed,'label_generation_attempts':attempts,'model':name,'scheme':scheme,
                      'oof_auc':float(roc_auc_score(yp,pred)),'mean_fold_auc':float(np.mean(fas)),
                      'preserves_subject_label_consistency':kind=='subject_label','original_split_seed':seed,'label_vector_sha256':hashlib.sha256(yp.tobytes()).hexdigest()})
                    for i in range(len(y)):nullpred.append({'null_scheme':kind,'permutation':b,'model':name,'scheme':scheme,'row_index':i,'y_permuted':int(yp[i]),'pred_probability':float(pred[i])})
            if b%5==0: print(json.dumps({'stage':stage,'part':kind,'permutation':b,'elapsed_s':round(time.time()-start,1)}),flush=True)
    N=pd.DataFrame(nullrows);N.to_csv(raw/'permutation_controls.csv',index=False)
    pd.DataFrame(nullpred).to_csv(raw/'permutation_predictions.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    s=[]
    for keys,ns in N.groupby(['null_scheme','model','scheme']):
        v=ns.oof_auc.to_numpy();s.append(dict(zip(['null_scheme','model','scheme'],keys))|{'n_permutations':len(v),'mean_auc':float(v.mean()),'median_auc':float(np.median(v)),
         'mc_q025':float(np.quantile(v,.025)),'mc_q975':float(np.quantile(v,.975)),
         'interpretation':'conditional randomized-label distribution, not real-case FPR and not a confidence interval'})
    pd.DataFrame(s).to_csv(summary/'permutation_control_summary.csv',index=False)
    dump(raw/'run_manifest.json',{'stage':stage,'started_utc':datetime.datetime.fromtimestamp(start,datetime.timezone.utc).isoformat(),
       'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.time()-start,
       'script_sha256':sha(CODE),'protocol_sha256':sha(BASE/'protocol/e1/protocol.json'),'data_sha256':sha(BASE/'data_manifest/raw/e1/parkinsons.data'),
       'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__,
       'feature_columns':features,'n_rows':len(y),'n_subjects':len(set(g)),'real_label_models_fit':R*2*2*2,'control_models_fit':nperm*2*2*2*2,
       'status':'COMPLETED','no_hyperparameter_tuning':True,'independent_negative_available':False})
    print('COMPLETED',stage,'elapsed',time.time()-start,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--stage',choices=['pilot','main'],required=True);run(a.parse_args().stage)
