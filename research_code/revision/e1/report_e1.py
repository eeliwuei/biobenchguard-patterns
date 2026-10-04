from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path(__file__).resolve().parents[3];r=B/'results_raw/e1/main';s=B/'results_summary/e1/main'
assert json.loads((r/'run_manifest.json').read_text())['status']=='COMPLETED'
stats=json.loads((s/'performance_summary.json').read_text());f=pd.read_csv(r/'fold_metrics.csv');N=pd.read_csv(r/'permutation_controls.csv');P=pd.read_csv(r/'oof_predictions.csv')
# Every row has exactly one held-out prediction in each repeat/model/scheme.
counts=P.groupby(['model','scheme','repeat','row_index']).size();assert counts.eq(1).all() and len(P)==195*30*4
assert f.query("scheme == 'subject_wise'").overlap_subjects.eq(0).all()
assert f.query("scheme == 'record_wise'").overlap_subjects.gt(0).all()
assert len(N)==99*2*2*2
labels={'random_forest':'Random forest','logistic_regression':'Logistic regression'}
colors={'record_wise':'#A04B2D','subject_wise':'#236B70'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
fig,axs=plt.subplots(1,2,figsize=(8.0,3.8),gridspec_kw={'width_ratios':[1,1.15]})
ax=axs[0]
for j,name in enumerate(labels):
 for k,scheme in enumerate(colors):
  v=stats[name+'__'+scheme];x=j+(-.13 if k==0 else .13);a=v['record_weighted_oof_auc_mean']
  ax.errorbar(x,a,yerr=[[a-v['cluster_ci_lower']],[v['cluster_ci_upper']-a]],fmt='o',markersize=5,capsize=3,color=colors[scheme],label=scheme.replace('_',' ') if j==0 else None)
  ax.annotate(f'{a:.3f}',(x,a),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
ax.axhline(.5,color='#7C8588',lw=.8,ls='--');ax.set(xticks=[0,1],xticklabels=['Random\nforest','Logistic\nregression'],ylim=(.45,1.025),xlim=(-.45,1.45),ylabel='ROC AUC',title='A  Original disease labels')
ax.legend(frameon=False,loc='lower right',fontsize=8);ax.text(0,-.24,'30 repeated two-fold splits; conditional\n95% subject-bootstrap intervals',transform=ax.transAxes,fontsize=8)
ax=axs[1]
positions=[];data=[];c=[]
for j,name in enumerate(labels):
 for k,scheme in enumerate(colors):
  q=N.query("null_scheme == 'subject_label' and model == @name and scheme == @scheme").oof_auc.to_numpy();data.append(q);positions.append(j*3+k);c.append(colors[scheme])
bp=ax.boxplot(data,positions=positions,widths=.6,patch_artist=True,showfliers=True,medianprops={'color':'white','lw':1.2},flierprops={'markersize':2,'marker':'.'})
for patch,col in zip(bp['boxes'],c):patch.set(facecolor=col,alpha=.85)
ax.axhline(.5,color='#7C8588',lw=.8,ls='--');ax.set(xticks=[.5,3.5],xticklabels=['Random forest','Logistic regression'],ylim=(.15,1.02),ylabel='ROC AUC',title='B  Subject-label permutation control')
ax.text(0,-.24,'99 conditional randomizations on a fixed split;\nlabels shuffled across subjects, not within them',transform=ax.transAxes,fontsize=8)
fig.subplots_adjust(left=.08,right=.98,bottom=.26,top=.9,wspace=.37)
fig.savefig(s/'e1_real_case.pdf');fig.savefig(s/'e1_real_case.png',dpi=220);plt.close(fig)
rows=[]
for name in labels:
 for scheme in colors:
  v=stats[name+'__'+scheme];ff=f.query('model == @name and scheme == @scheme')
  rows.append({'model':name,'workflow':scheme,'auc_mean':v['record_weighted_oof_auc_mean'],'auc_ci_low':v['cluster_ci_lower'],'auc_ci_high':v['cluster_ci_upper'],
  'balanced_accuracy':v['balanced_accuracy_mean'],'overlap_test_record_fraction_min':float(ff.test_record_identity_overlap_fraction.min()),'overlap_test_record_fraction_max':float(ff.test_record_identity_overlap_fraction.max()),
  'n_folds':len(ff),'feature_scope':'22 native acoustic features; no identity/label feature','tuning':'none',
  'action':'replace record-wise evaluation for unseen-subject claim' if scheme=='record_wise' else 'retain subject-disjoint evaluation; external validation still needed'})
pd.DataFrame(rows).to_csv(s/'action_decision_table.csv',index=False)
lines=['# E1: real-label identity-confounding replay','',
'This is a controlled replay of the mechanism documented by [Chaibub Neto et al. (2019)](https://doi.org/10.1038/s41746-019-0178-x), using authentic status labels from the [UCI Parkinsons dataset](https://archive.ics.uci.edu/dataset/174/parkinsons). The exact R implementation was not rerun; all departures are recorded in protocol/e1/original_to_replay_mapping.csv.','',
'The unchanged file contains 195 voice records and 32 parsed subjects (24 Parkinson cases and 8 controls). The current UCI prose reports 31 subjects while the author paper and actual file agree on 32; this discrepancy is explicitly retained in the source audit.','',
'| Model | Record-wise AUC (95% conditional CI) | Subject-wise AUC (95% conditional CI) | Paired difference (95% conditional CI) |','|---|---|---|---|']
for name,label in labels.items():
 a=stats[name+'__record_wise'];b=stats[name+'__subject_wise'];d=stats[name+'__gap']
 lines.append(f"| {label} | {a['record_weighted_oof_auc_mean']:.3f} ({a['cluster_ci_lower']:.3f}, {a['cluster_ci_upper']:.3f}) | {b['record_weighted_oof_auc_mean']:.3f} ({b['cluster_ci_lower']:.3f}, {b['cluster_ci_upper']:.3f}) | {d['value']:.3f} ({d['ci_lower']:.3f}, {d['ci_upper']:.3f}) |")
lines += ['', 'Intervals resample subjects within status classes and keep all repeated predictions for each resampled subject together. They are conditional on the fitted models and these 32 subjects; they do not capture full training or population uncertainty. Split-to-split variability is reported separately. Repeated folds share data and are not independent real cases; their counts must not be used as binomial false-positive-rate denominators.', '',
'The saved manifests expose identity overlap in all 60 record-wise folds and none of the 60 repaired subject-wise folds. This is an observed information-path check for the stated unseen-subject target. It is also what ordinary group-aware validation can already reveal. Performance reduction by itself is not used as reference truth or evidence of methodological novelty.', '',
'The same X, y and subject vector underlie both workflows. A data-only random-minus-group performance diagnostic therefore gives the same dataset-level sensitivity result regardless of which actual workflow is described; the observed pipeline manifest is necessary to distinguish the two. The new report adds effect magnitude, conditional uncertainty, model sensitivity and explicit limits to a conventional split audit. It does not establish superior detection to existing identity-confounding methods.', '',
'## Randomized-label controls','',
'Label permutation across subjects preserves repeated-record identity structure while removing the original subject-to-status allocation; row permutation destroys that allocation and repeated-record consistency. Neither is a collection of independent real clean cases, and neither estimates real-case FPR. Within-subject status permutation is a no-op because every subject has a constant status label, so that proposed null is not identifiable here. These controls use a single frozen split and 99 draws each; their quantiles are conditional randomized-label distributions, not confidence intervals.', '',
'## Remaining evidence boundary','',
'An independently sourced, scrutinized real negative workflow has not been recovered. The repaired subject-disjoint workflow is a local mechanism-specific control, not globally clean and not independently certified. R1-M3’s independent-real-negative branch remains open. No clinical validity, causal identification, broad biomarker transportability, or universal detection calibration is claimed.', '',
'## Reproduction','',
'Run `python research_code/revision/e1/test_e1.py`, then `python research_code/revision/e1/run_e1.py --stage main`, then `python research_code/revision/e1/report_e1.py` with the supplied scientific environment. The main run checks the frozen script and protocol hashes before executing. Raw predictions, fold assignments, randomized labels, bootstrap draws, source receipts and runtime versions are saved under their E1 directories.']
(s/'E1_results.md').write_text('\n'.join(lines)+'\n')
(s/'positive_negative_control_audit.md').write_text('''# Control audit

| Case | Reference evidence | Status | What is not established |
|---|---|---|---|
| Record-wise workflow on real Parkinson status | Published identity-confounding mechanism; saved row assignment and explicit subject overlap | Documented mechanism replay on genuine labels | Exact numerical R-paper replication; dataset intrinsically erroneous |
| Subject-wise repair on same real labels | Subject-disjoint sampling; train/test intersections independently enumerable from provenance file | Local negative for the specified cross-subject overlap pathway | Independence from case-construction team; absence of any other leakage/bias |
| Independent real negative | No external, independently scrutinized negative workflow recovered in this bounded case analysis | NOT AVAILABLE | Cannot close R1-M3’s real dual-control branch |
| Subject-label randomizations | Preserve labels within a subject while destroying authentic subject-status assignment | Mechanism control | Genuine-label clinical performance; real clean-case FPR |
| Row-label randomizations | Destroy repeated-label structure and outcome allocation | Randomization baseline | A real negative workflow |

The manifest auditor receives only subject IDs, train/test assignments and the declared target population. The private reference file is generated separately and never read by `audit_manifest`. Expected statuses are forbidden auditor inputs. The classifiers use only 22 native acoustic features; subject ID is not a feature. Engineering unit fixtures are not counted as real negative controls.
''')
# Invariance is an explicit input-interface check, separate from performance claims.
datahash=hashlib.sha256((B/'data_manifest/raw/e1/parkinsons.data').read_bytes()).hexdigest()
invariance=[]
for case in ['record_wise','subject_wise']:
 for name in labels:invariance.append({'case_pipeline':case,'model':name,'data_sha256':datahash,'diagnostic_input_interface':['X','y','group'],'data_only_gap':stats[name+'__gap']['value'],'actual_pipeline_manifest_read_by_data_only_gap':False,'interpretation':'same input statistics reused transparently; a dataset-level diagnostic cannot identify the actual workflow'})
(s/'same_data_different_pipeline.json').write_text(json.dumps(invariance,indent=2)+'\n')
manifest={str(p.relative_to(B)):hashlib.sha256(p.read_bytes()).hexdigest() for p in list(r.glob('*'))+list(s.glob('*')) if p.is_file() and p.name!='output_checksums.json'}
(s/'output_checksums.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'status':'PASS','real_oof_predictions':len(P),'randomization_rows':len(N),'figure':str(s/'e1_real_case.pdf'),'outputs':len(manifest)}))
