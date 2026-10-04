"""Independent checks from saved predictions, fold memberships and manifests."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()


ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);args=ap.parse_args()
root=Path(args.root).resolve();raw=root/'results_raw/e2';summ=root/'results_summary/e2';checks={}
X=np.load(raw/'crc_strict_inputs.npz');p=pd.read_csv(raw/'loso_predictions.csv');q=pd.read_csv(raw/'random_matched_size_predictions.csv')
folds=pd.read_csv(raw/'fold_membership.csv');metrics=pd.read_csv(summ/'loso_per_study.csv');fits=pd.read_csv(raw/'fit_execution.csv')
meta=json.loads((raw/'input_manifest.json').read_text());freeze=json.loads((root/'protocol/e2/execution_freeze.json').read_text())
checks['input_hash']=sha(raw/'crc_strict_inputs.npz')==freeze['input_npz_sha256']
checks['script_hash']=sha(root/'research_code/revision/e2/run_e2.py')==freeze['script_sha256']
checks['protocol_hash']=sha(root/'protocol/e2/analysis_plan.md')==freeze['plan_sha256']
checks['raw_hashes']=all(sha(root/'data_manifest/raw'/f)==h for f,h in meta['source_files_sha256'].items())
checks['fit_budget_matched']=fits.groupby('model').size().to_dict()=={'logistic_C1':72,'random_forest_120':72}
checks['no_fit_warnings']=bool((fits.warnings=='[]').all())
checks['finite_predictions']=bool(np.isfinite(p.prediction).all() and np.isfinite(q.prediction).all())
checks['primary_once_per_subject_model_arm']=bool((p.groupby(['sample_id','model','arm']).size()==1).all())
checks['primary_full_coverage']=bool((p.groupby(['model','arm']).size()==len(X['y'])).all())
labelmap=dict(zip(X['ids'],X['y']))
checks['labels_match_original_input']=all(int(r.y)==int(labelmap[r.sample_id]) for r in pd.concat([p,q]).itertuples())
checks['all_folds_disjoint']=all(not(set(z[z.role=='train'].sample_id)&set(z[z.role=='test'].sample_id)) for _,z in folds.groupby(['protocol','split_id']))
checks['all_folds_full_population']=all(len(z)==len(X['y']) and z.sample_id.nunique()==len(X['y']) for _,z in folds.groupby(['protocol','split_id']))
checks['LOSO_study_separation']=all((z[z.role=='test'].study==sid).all() and (z[z.role=='train'].study!=sid).all() for (prot,sid),z in folds.groupby(['protocol','split_id']) if prot=='LOSO')
errors=[]
for (study,model,arm),z in p.groupby(['split_id','model','arm']):
    expected=metrics[(metrics.study==study)&(metrics.model==model)&(metrics.arm==arm)].auroc.iloc[0]
    errors.append(abs(roc_auc_score(z.y,z.prediction)-expected))
checks['independent_sklearn_auc_recompute']=max(errors)<1e-12
random_scores=pd.read_csv(summ/'random_matched_size_scores.csv')
errors=[]
for (sid,model,arm),z in q.groupby(['split_id','model','arm']):
    st=sid.split('_repeat')[0].replace('matched_','');rep=int(sid.split('_repeat')[1]);row=random_scores[(random_scores.matched_loso_study==st)&(random_scores.repeat==rep)&(random_scores.model==model)&(random_scores.arm==arm)].iloc[0]
    errors.append(abs(roc_auc_score(z.y,z.prediction)-row.auroc))
    expected_test=int(metrics[metrics.study==st].n_test.iloc[0])
    assert len(z)==expected_test and row.n_train==len(X['y'])-expected_test
checks['random_auc_recompute_and_exact_sizes']=max(errors)<1e-12
index={x:i for i,x in enumerate(X['ids'])};feature_index={x:i for i,x in enumerate(X['features'])}
selections=pd.read_csv(raw/'selected_features.csv')
checks['feature_selection_uses_only_train_variable_universe']=all(
    np.var(X['X'][[index[x] for x in folds[(folds.protocol=='LOSO')&(folds.split_id==r.split_id)&(folds.role=='train')].sample_id],feature_index[r.feature]])>1e-12
    for r in selections.itertuples())
checks={k:bool(v) for k,v in checks.items()}
out={'status':'PASS'if all(checks.values())else'FAIL','checks':checks,'main_fits':len(fits),'primary_prediction_rows':len(p),'random_control_prediction_rows':len(q),
     'analysis_scope':'Saved outputs independently validated; no assertion that cross-study population calibration or clinical validity has been established'}
(raw/'independent_output_QA.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
assert all(checks.values())
