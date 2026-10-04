"""Publication vectors generated from observed source tables; no image synthesis."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from matplotlib import colors
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':600})
C=['#176A75','#C57A22','#525F8A','#797979']
def save(fig,out,name):
 for ext in ['pdf','svg','png']:fig.savefig(out/f'{name}.{ext}',bbox_inches='tight',pad_inches=.09)
 plt.close(fig)
def box(ax,x,y,w,h,title,body,color='#EAF2F3'):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.015',ec='#31434B',fc=color,lw=1));ax.text(x+w/2,y+h*.72,title,ha='center',va='center',fontweight='bold',fontsize=9);ax.text(x+w/2,y+h*.35,body,ha='center',va='center',fontsize=8)
def arrow(ax,a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=11,lw=1,color='#526168'))
def build(root):
 root=Path(root);out=root/'figures';out.mkdir(exist_ok=True);src=root/'results_summary';manifest=[]
 fig,ax=plt.subplots(figsize=(7.1,4.0));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 for x,t,b in [(.025,'1  Specify target','New subject, study or group\nRecord executed fit scopes'),(.355,'2  Define contrast','Selection-scope response\nGroup-split response'),(.685,'3  Define reference','Geometry and outcome law\nModel, selector and splits')]:box(ax,x,.65,.29,.27,t,b)
 arrow(ax,(.325,.785),(.34,.785));arrow(ax,(.655,.785),(.67,.785))
 for x,t,b in [(.685,'4  Calibrate / evaluate','Separate random instances\nCutoffs, counts and intervals'),(.355,'5  Check transfer','Fixed source cutoff\nSource and target rates'),(.025,'6  Interpret real case','Observed workflow evidence\nNext validation action')]:box(ax,x,.23,.29,.27,t,b,'#F4F0E7')
 arrow(ax,(.83,.63),(.83,.52));arrow(ax,(.67,.365),(.655,.365));arrow(ax,(.34,.365),(.325,.365))
 ax.text(.5,.075,'Missing provenance or incompatible calibration: withhold a cleanliness conclusion',ha='center',color='#8C3B2C',fontsize=8.6)
 save(fig,out,'Fig1_workflow');manifest.append(('Fig1','conceptual; estimand registry'))
 fig,ax=plt.subplots(figsize=(7.1,6.1));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 box(ax,.30,.815,.40,.15,'Declared controlled task','Fixed X; outcome law; model and folds')
 ax.plot([.50,.50],[.80,.755],color='#526168');ax.plot([.205,.795],[.755,.755],color='#526168');arrow(ax,(.205,.755),(.205,.71));arrow(ax,(.795,.755),(.795,.71))
 box(ax,.035,.55,.34,.15,'Calibration reference IDs','L=0; calculate paired scores')
 box(ax,.625,.55,.34,.15,'Separate evaluation IDs','New draws\nL=0 and added candidates')
 arrow(ax,(.205,.53),(.205,.47));arrow(ax,(.795,.53),(.795,.47))
 box(ax,.035,.30,.34,.15,'Estimate one cutoff','Rank rule from calibration only')
 box(ax,.625,.30,.34,.15,'Evaluate the fixed cutoff','Strict score > cutoff\nCounts and intervals')
 arrow(ax,(.395,.375),(.605,.375))
 box(ax,.08,.025,.84,.20,'The paired score in every reference','B: train-only selector   |   C: full-label selector\nBoth: train-only scaling and model; identical held-out folds\nScore = metric(C) − metric(B); expected status is never an input.','#EDF0F5')
 ax.text(.5,.26,'L=0 removes appended candidates; C still accesses full outcomes.',ha='center',fontsize=8.5,color='#8C3B2C')
 save(fig,out,'Fig2_scope_and_scoring');manifest.append(('Fig2','conceptual; core.py and estimand registry'))
 e1=src/'e1/main/real_failure_recovery.csv'
 if e1.exists():
  d=pd.read_csv(e1);fig,ax=plt.subplots(1,2,figsize=(7.2,3.1),gridspec_kw={'width_ratios':[1.2,1]})
  for j,model in enumerate(['random_forest','logistic_regression']):
   for i,scheme in enumerate(['record_wise','subject_wise']):
    r=d[(d.model==model)&(d.scheme==scheme)].iloc[0];x=j+(i-.5)*.25
    ax[0].errorbar(x,r.record_weighted_oof_auc_mean,yerr=[[r.record_weighted_oof_auc_mean-r.cluster_ci_lower],[r.cluster_ci_upper-r.record_weighted_oof_auc_mean]],fmt='o',color=C[i],capsize=3,label=['Record split','Subject split'][i] if j==0 else None)
  ax[0].set(xticks=[0,1],xticklabels=['RF500','Logistic'],ylim=(.45,1.02),ylabel='Recording-weighted OOF AUROC',title='A  Real disease labels');ax[0].legend(frameon=False,fontsize=8);ax[0].axhline(.5,color='#aaa',ls=':',lw=1)
  # Contrast values below are recomputed from raw summary table, not synthetic plotted points.
  for j,model in enumerate(['random_forest','logistic_regression']):
   a=d[(d.model==model)&(d.scheme=='record_wise')].iloc[0];b=d[(d.model==model)&(d.scheme=='subject_wise')].iloc[0];ax[1].bar(j,a.record_weighted_oof_auc_mean-b.record_weighted_oof_auc_mean,color=C[j],width=.5)
  ax[1].set(xticks=[0,1],xticklabels=['RF500','Logistic'],ylabel='Record minus subject AUROC',title='B  Change after split repair',ylim=(0,.20))
  fig.text(.5,-.02,'195 recordings; 32 subjects; 30 twofold repeats. Intervals: 5,000 subject bootstraps of fixed OOF predictions.',ha='center',fontsize=7);fig.tight_layout();save(fig,out,'Fig3_real_failure_recovery');manifest.append(('Fig3','results_summary/e1/main/real_failure_recovery.csv'))
 e=src/'e3e4/operating_rates.csv'
 if e.exists():
  d=pd.read_csv(e);fig,ax=plt.subplots(1,2,figsize=(7.2,3.4))
  bases=['microbiome_crc','wheat','rice'];labels=['CRC reference','Wheat','Rice']
  for j,b in enumerate(bases):
   r=d[(d.condition_id==b+'_baseline')&(d.L==0)].iloc[0]
   for i,(stat,lo,hi) in enumerate([('fixed_rate','fixed_ci_lower','fixed_ci_upper'),('rate','ci_lower','ci_upper')]):
    ax[0].errorbar(j+(i-.5)*.20,r[stat],yerr=[[r[stat]-r[lo]],[r[hi]-r[stat]]],fmt='o',capsize=3,color=C[i],label=['Strict T > 0.05 (analyst-defined)','Independent calibration'][i] if j==0 else None)
   v=d[(d.condition_id==b+'_baseline')&(d.L>0)].sort_values('L');ax[1].errorbar(v.L,v.rate,yerr=[v.rate-v.ci_lower,v.ci_upper-v.rate],marker=['o','s','^'][j],color=C[j],capsize=2,label=labels[j])
  ax[0].axhline(.05,ls=':',color='#777');ax[0].set(xticks=range(3),xticklabels=labels,ylim=(-.02,1.04),ylabel='Baseline-reference exceedance',title='A  Two explicitly defined cutoffs');ax[0].legend(fontsize=7,frameon=False)
  ax[1].axhline(.8,ls=':',color='#777');ax[1].set(xlabel='Additional independent candidates L',ylabel='Reference exceedance',ylim=(-.02,1.04),title='B  Extra-candidate response');ax[1].legend(fontsize=7,frameon=False);fig.tight_layout();save(fig,out,'Fig5_reference_operating_rates');manifest.append(('Fig5','results_summary/e3e4/operating_rates.csv'))
  tr=pd.read_csv(src/'e3e4/calibration_transfer_matrix.csv');fig,axs=plt.subplots(1,2,figsize=(9,5.4))
  for ax,task,title in zip(axs,['microbiome','regression'],['A  Classification conditions','B  Regression conditions']):
   ts=[x for x in tr.target_theta.unique() if ('microbiome' in x)==(task=='microbiome')];ss=[x for x in tr.source_theta.unique() if x in ts];p=tr.pivot(index='target_theta',columns='source_theta',values='observed_rate').reindex(index=ts,columns=ss)
   im=ax.imshow(p,aspect='equal',vmin=0,vmax=1,cmap='cividis');lab=lambda x:x.replace('microbiome_crc_','CRC ').replace('wheat_','W ').replace('rice_','R ').replace('dimension_matched','n/p matched').replace('linear_signal05','linear s0.5')
   ax.set(xticks=range(len(ss)),xticklabels=[lab(x) for x in ss],yticks=range(len(ts)),yticklabels=[lab(x) for x in ts],title=title,xlabel='Source calibration',ylabel='Target evaluation');ax.tick_params(labelsize=13);ax.xaxis.label.set_size(13);ax.yaxis.label.set_size(13);ax.title.set_size(13);plt.setp(ax.get_xticklabels(),rotation=65,ha='right')
   for j,t in enumerate(ts):
    if t in ss:ax.add_patch(plt.Rectangle((ss.index(t)-.5,j-.5),1,1,fill=False,edgecolor='white',lw=.7))
  fig.subplots_adjust(left=.1,right=.84,bottom=.35,wspace=1.15);cax=fig.add_axes([.91,.33,.018,.45]);cb=fig.colorbar(im,cax=cax,label='Held-out baseline-reference\nexceedance rate');cb.ax.tick_params(labelsize=13);cb.ax.yaxis.label.set_size(13);save(fig,out,'Fig6_calibration_transport');manifest.append(('Fig6','results_summary/e3e4/calibration_transfer_matrix.csv'))
  d=pd.read_csv(src/'e3e4/mixed_response.csv');fig,axs=plt.subplots(2,4,figsize=(7.1,5.1))
  for col,(base,model) in enumerate([(b,m) for b in ['wheat','microbiome_crc'] for m in ['linear','rf']]):
   for row,stat in enumerate(['T_scope','T_split']):
    ax=axs[row,col]
    for eta,ls in [(0,'-'),(.2,'--'),(.8,':')]:
     for enc,mk,alpha in [(False,'o',.9),(True,'s',.65)]:
      v=d[(d.base==base)&(d.model==model)&(d.statistic==stat)&(d.eta==eta)&(d.encoding==enc)].sort_values('L');ax.plot(v.L,v['mean'],ls,marker=mk,markersize=3,alpha=alpha,color=C[[0,.2,.8].index(eta)],label=f'η={eta:g}, '+('X+Z' if enc else 'X'))
    ax.axhline(0,color='#aaa',lw=.6);ax.set_title(('Wheat' if base=='wheat' else 'CRC')+' / '+model,fontsize=9)
    if col==0:ax.set_ylabel('Scope contrast' if row==0 else 'Split contrast')
    if row==1:ax.set_xlabel('Candidates L')
    ax.tick_params(labelsize=9);ax.set_xticks([0,2500])
  handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False,fontsize=9);fig.tight_layout(rect=(0,.14,1,1));save(fig,out,'Fig7_mixed_mechanisms');manifest.append(('Fig7','results_summary/e3e4/mixed_response.csv; intervals in source/supplement'))
 e2=src/'e2/per_study_contrasts.csv'
 if e2.exists():
  import json
  d=json.loads((src/'e2/e2_results.json').read_text());fig,ax=plt.subplots(1,2,figsize=(7.2,3.2));names=['Zeller','Baxter','Zackular'];keys=['crc_zeller','crc_baxter','crc_zackular'];v=[d['rf_minus_logistic_B_by_study'][k] for k in keys];ax[0].bar(names,v,color=[C[0] if x>=0 else C[1] for x in v],width=.6);ax[0].axhline(0,color='#444',lw=.8);ax[0].set(ylabel='RF120 minus logistic AUROC',title='A  Unseen-study model comparison')
  v=[x['observed'] for x in d['diagnostics']];ax[1].bar(['Logistic','RF120'],v,color=C[:2],width=.55);ax[1].axhline(.02,color='#444',ls='--',label='Analyst-defined cutoff = 0.02');ax[1].legend(fontsize=7,frameon=False);ax[1].set(ylabel='Full-label minus fold-internal AUROC',title='B  Known selector-scope contrast');fig.text(.5,-.025,'468 original-label samples; 3 studies. Both C workflows access held-out labels, regardless of the numeric response.',ha='center',fontsize=7);fig.tight_layout();save(fig,out,'Fig4_real_CRC_decision');manifest.append(('Fig4','results_summary/e2/e2_results.json; per_study_contrasts.csv'))
 pd.DataFrame(manifest,columns=['figure','source']).to_csv(out/'figure_source_manifest.csv',index=False)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);a=ap.parse_args();build(a.root)
