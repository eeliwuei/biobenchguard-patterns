"""Summaries are calculated solely from actual reference-level outputs."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd
from core import threshold,wilson,seed

def summarize(run,out):
 run=Path(run);out=Path(out);out.mkdir(parents=True,exist_ok=True)
 m=json.loads((run/'run_manifest.json').read_text())
 if m['status']!='COMPLETE' or m['pilot']:raise ValueError('Only complete confirmatory runs may be summarized for manuscript')
 assert hashlib.sha256((run/'raw_results.jsonl').read_bytes()).hexdigest()==m['raw_sha256'],'raw hash mismatch'
 rows=[json.loads(x) for x in (run/'raw_results.jsonl').read_text().splitlines()];assert len(rows)==m['n_attempted'],'attempt count mismatch'
 from experiments import jobs
 assert {j['instance_id'] for j in jobs(m['config'])}=={r['instance_id'] for r in rows},'planned instance set mismatch'
 df=pd.DataFrame(rows).sort_values('instance_id').reset_index(drop=True);e=df[df.experiment=='E3'];rates=[];trans=[];cals=[];lods=[]
 assert df.instance_id.is_unique and (df.status=='COMPLETE').all()
 for cid,d in e.groupby('condition_id'):
  cal=d[d.role=='calibration'];ev=d[d.role=='evaluation'];assert not set(cal.reference_id)&set(ev.reference_id)
  tau=threshold(cal.T_scope);calid='cal_'+cid
  cals.append({'calibration_id':calid,'condition_id':cid,'n_calibration':len(cal),'cutoff':tau,'alpha':.05,'quantile_rule':'ceil((B+1)*.95),strict greater-than','null':'no additional candidate features at this outcome law and protocol; intentionally full-label C arm retained'})
  for L,dd in d[d.role!='calibration'].groupby('L'):
   n=len(dd);k=int((dd.T_scope>tau).sum());lo,hi=wilson(k,n);kb=int((dd.T_scope>.05).sum());bl,bh=wilson(kb,n)
   rates.append({'condition_id':cid,'base':dd.base.iloc[0],'model':dd.model.iloc[0],'L':int(L),'calibration_id':calid,'cutoff':tau,'n_calibration':len(cal),'n_attempted':n,'n_valid':n,'n_failed':0,'n_positive':k,'rate':k/n,'ci_lower':lo,'ci_upper':hi,'fixed_cutoff':.05,'fixed_n_positive':kb,'fixed_rate':kb/n,'fixed_ci_lower':bl,'fixed_ci_upper':bh,'ci_scope':'Wilson95 conditional on fixed calibration threshold; conditional-on-X references'})
  # No monotonic smoothing: discrete first crossing and censoring, bootstrap BOTH threshold & paired alternative instances.
  al=d[d.role=='alternative'];levels=sorted(int(x) for x in al.L.unique())
  if len(levels)>1:
   a=al.pivot(index='reference_id',columns='L',values='T_scope').sort_index();v=a[levels].to_numpy();rng=np.random.default_rng(seed('lod_bootstrap',cid));boot=[]
   for i in range(1000):
    tb=threshold(rng.choice(cal.T_scope.to_numpy(),len(cal),replace=True));vb=v[rng.integers(0,len(v),len(v))];ps=(vb>tb).mean(0)
    boot.append([next((int(l) for l,p in zip(levels,ps) if p>=target),float('inf')) for target in [.8,.95]])
   ps=(v>tau).mean(0)
   for j,target in enumerate([.8,.95]):
    vals=np.array(boot)[:,j];cross=next((int(l) for l,p in zip(levels,ps) if p>=target),None)
    con=next((int(l) for l,p in zip(levels,ps) if wilson(int(round(p*len(v))),len(v))[0]>=target),None)
    def q(p):
     vv=np.sort(vals);v=vv[min(len(vv)-1,int(np.ceil(p*len(vv)))-1)];return '>'+str(max(levels)) if np.isinf(v) else str(int(v))
    lods.append({'condition_id':cid,'target_exceedance':target,'grid_crossing':cross if cross is not None else '>'+str(max(levels)),'pointwise_lower_bound_grid_crossing':con if con is not None else '>'+str(max(levels)),'n_bootstrap':1000,'n_finite':int(np.isfinite(vals).sum()),'n_right_censored':int(np.isinf(vals).sum()),'n_evaluation_per_severity':len(v),'calibration_id':calid,'tested_grid':str(levels),'fraction_right_censored':float(np.isinf(vals).mean()),'unconditional_median':q(.5),'ci_lower':q(.025),'ci_upper':q(.975),'grid_max':max(levels),'lower_bound_scope':'pointwise Wilson conditional on fixed cutoff; not simultaneous after grid search','smoothing':'none','interpretation':'grid crossing of extra-candidate reference exceedance; not a universal leakage detection limit'})
 for target,d in e[e.role=='evaluation'].groupby('condition_id'):
  task=d.task.iloc[0];sources=[x for x in cals if e[e.condition_id==x['condition_id']].task.iloc[0]==task]
  for src in sources:
   source_ev=e[(e.role=='evaluation')&(e.condition_id==src['condition_id'])];source_k=int((source_ev.T_scope>src['cutoff']).sum());
   n=len(d);k=int((d.T_scope>src['cutoff']).sum());lo,hi=wilson(k,n);one_lo,_=wilson(k,n,z=1.6448536269514722)
   source_lo,source_hi=wilson(source_k,len(source_ev));p_target=k/n;p_source=source_k/len(source_ev);diff=p_target-p_source
   diff_lo=diff-np.sqrt((p_target-lo)**2+(source_hi-p_source)**2);diff_hi=diff+np.sqrt((hi-p_target)**2+(p_source-source_lo)**2)
   if src['condition_id']==target:diff_lo=diff_hi=0.
   trans.append({'source_theta':src['condition_id'],'target_theta':target,'calibration_id':src['calibration_id'],'n_calibration':src['n_calibration'],'cutoff':src['cutoff'],'nominal_alpha':.05,'observed_k':k,'observed_n':n,'observed_rate':k/n,'source_eval_k':source_k,'source_eval_n':len(source_ev),'source_eval_rate':source_k/len(source_ev),'target_minus_source_rate':diff,'difference_ci_lower':diff_lo,'difference_ci_upper':diff_hi,'matrix_inference':'cellwise exploratory; Newcombe95 difference conditional on cutoff; no multiplicity-controlled global discovery claim','ci_lower':lo,'ci_upper':hi,'one_sided_lower95':one_lo,'transport_status':'CUTOFF_EXCEEDS_NOMINAL_RECALIBRATE' if one_lo>.05 else 'INFLATION_NOT_DEMONSTRATED_NOT_VALIDATION','scheme':'task-matched' if src['condition_id']==target else ('geometry-only' if src['condition_id']==d.base.iloc[0]+'_baseline' else 'transported')})
 for name,data in [('operating_rates',rates),('calibration_registry',cals),('calibration_transfer_matrix',trans),('grid_crossing',lods)]:pd.DataFrame(data).to_csv(out/f'{name}.csv',index=False)
 mixed=df[df.experiment=='E4'];summ=[]
 for keys,d in mixed.groupby(['base','model','L','eta','encoding']):
  for stat in ['T_scope','T_split']:
   v=d[stat].to_numpy();rng=np.random.default_rng(seed('mixedCI',*keys,stat));b=np.mean(rng.choice(v,(2000,len(v)),replace=True),axis=1);lo,hi=np.quantile(b,[.025,.975]);summ.append(dict(zip(['base','model','L','eta','encoding'],keys),statistic=stat,n=len(v),mean=v.mean(),ci_lower=lo,ci_upper=hi))
 pd.DataFrame(summ).to_csv(out/'mixed_response.csv',index=False)
 # Paired difference-in-differences as descriptive interaction at highest L and eta.
 inter=[]
 for keys,d in mixed.groupby(['base','model','encoding']):
  for stat in ['T_scope','T_split']:
   p=d.pivot(index='reference_id',columns=['L','eta'],values=stat);v=(p[(2500,.8)]-p[(2500,0.)]-p[(0,.8)]+p[(0,0.)]).to_numpy();rng=np.random.default_rng(seed('interaction',*keys,stat));b=rng.choice(v,(2000,len(v)),replace=True).mean(1);lo,hi=np.quantile(b,[.025,.975]);inter.append(dict(zip(['base','model','encoding'],keys),statistic=stat,n=len(v),interaction_mean=v.mean(),ci_lower=lo,ci_upper=hi))
 pd.DataFrame(inter).to_csv(out/'mixed_interactions.csv',index=False)
 # Reference-level values exported without truth annotations entering computations.
 columns=[c for c in df if c not in ['condition','fold_sizes','pred_B','pred_C','pred_group','observed_y']];df[columns].sort_values('instance_id').to_csv(out/'reference_results.csv',index=False)
 (out/'summary_manifest.json').write_text(json.dumps({'source_raw_sha256':hashlib.sha256((run/'raw_results.jsonl').read_bytes()).hexdigest(),'config_sha256':m['config_sha256'],'summarizer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'core_sha256':hashlib.sha256((Path(__file__).parent/'core.py').read_bytes()).hexdigest(),'n_rows':len(df),'n_E3':len(e),'n_E4':len(mixed),'n_failures':int((df.status!='COMPLETE').sum()),'inference':'No universal FPR, absence, joint alarm, or real-task threshold validity claim. Fixed cutoff .05 is analyst-defined.'},indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--run',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();summarize(a.run,a.out)
