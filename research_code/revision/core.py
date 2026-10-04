"""Task-conditional reference computations; no expected detector status is accepted."""
import hashlib,math,warnings
import numpy as np
from scipy.special import expit
from scipy.optimize import brentq
from sklearn.linear_model import Ridge,LogisticRegression
from sklearn.ensemble import RandomForestRegressor,RandomForestClassifier
from sklearn.model_selection import KFold,StratifiedKFold,GroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score,roc_auc_score
from sklearn.exceptions import ConvergenceWarning

def seed(*parts):return int(hashlib.sha256('|'.join(map(str,parts)).encode()).hexdigest()[:8],16)
def corr_rank(X,y,k):
 yc=y-y.mean();xc=X-X.mean(0);den=np.sqrt((xc*xc).sum(0)*float(yc@yc));v=np.divide(xc.T@yc,den,out=np.zeros(X.shape[1]),where=den>0)
 return np.argsort(-np.abs(v),kind='stable')[:min(k,X.shape[1])]
def metric(y,p,task):
 if not np.isfinite(p).all():raise ValueError('nonfinite predictions')
 if task=='binary' and len(set(y))<2:raise ValueError('undefined one-class AUROC')
 return float(roc_auc_score(y,p) if task=='binary' else r2_score(y,p))
def model(task,family,rs):
 if family=='linear':return LogisticRegression(C=1,max_iter=2000,solver='lbfgs') if task=='binary' else Ridge(alpha=1)
 cls=RandomForestClassifier if task=='binary' else RandomForestRegressor
 return cls(n_estimators=120,max_features=1.0,min_samples_leaf=2,n_jobs=1,random_state=rs)
def predict_fit(X,y,tr,te,sel,task,family,rs):
 sc=StandardScaler();a=sc.fit_transform(X[tr][:,sel]);b=sc.transform(X[te][:,sel]);m=model(task,family,rs)
 with warnings.catch_warnings():
  warnings.simplefilter('error',ConvergenceWarning);m.fit(a,y[tr])
 return m.predict_proba(b)[:,1] if task=='binary' else m.predict(b)
def scope_stat(X,y,task='regression',family='linear',k=20,folds=5,split_seed=1,return_predictions=False):
 sp=StratifiedKFold(folds,shuffle=True,random_state=split_seed) if task=='binary' else KFold(folds,shuffle=True,random_state=split_seed)
 full=corr_rank(X,y,k);pb=np.full(len(y),np.nan);pc=pb.copy();fitlog=[]
 for f,(tr,te) in enumerate(sp.split(X,y)):
  sel=corr_rank(X[tr],y[tr],k);rs=seed(split_seed,'model',f)
  pb[te]=predict_fit(X,y,tr,te,sel,task,family,rs);pc[te]=predict_fit(X,y,tr,te,full,task,family,rs)
  fitlog.append({'fold':f,'n_train':len(tr),'n_test':len(te),'B_selector_fit_ids':tr.tolist(),'C_selector_fit_ids':list(range(len(y))),'transform_fit_ids':tr.tolist(),'test_ids':te.tolist()})
 out={'T_scope':metric(y,pc,task)-metric(y,pb,task),'score_B':metric(y,pb,task),'score_C':metric(y,pc,task)}
 if return_predictions:out.update(pred_B=pb,pred_C=pc,fitlog=fitlog)
 return out

def split_stat(X,y,g,task,family,k,split_seed):
 # Random folds have EXACT same test-size vector as group folds; train sizes match.
 gs=list(GroupKFold(n_splits=min(5,len(set(g)))).split(X,y,g));perm=np.random.default_rng(split_seed).permutation(len(y));pR=np.full(len(y),np.nan);pG=pR.copy();pC=pR.copy();full=corr_rank(X,y,k);start=0;logs=[]
 for f,(trG,teG) in enumerate(gs):
  teR=perm[start:start+len(teG)];start+=len(teG);trR=np.setdiff1d(np.arange(len(y)),teR)
  for tr,te,p in [(trG,teG,pG),(trR,teR,pR)]:
   if task=='binary' and len(set(y[tr]))<2:raise ValueError('single-class training group partition')
   sel=corr_rank(X[tr],y[tr],k);p[te]=predict_fit(X,y,tr,te,sel,task,family,seed(split_seed,'splitmodel',f))
  assert len(trG)==len(trR) and len(teG)==len(teR) and not set(g[trG])&set(g[teG])
  pC[teR]=predict_fit(X,y,trR,teR,full,task,family,seed(split_seed,'splitmodel',f))
  logs.append({'fold':f,'train_n':len(trG),'test_n':len(teG),'random_train_ids':trR.tolist(),'random_test_ids':teR.tolist(),'group_train_ids':trG.tolist(),'group_test_ids':teG.tolist(),'B_selected_features':corr_rank(X[trR],y[trR],k).tolist(),'C_selected_features':full.tolist(),'group_overlap_count':len(set(g[trG])&set(g[teG]))})
 return {'T_split':metric(y,pR,task)-metric(y,pG,task),'score_random_matched':metric(y,pR,task),'score_group':metric(y,pG,task),'fold_sizes':logs,'T_scope':metric(y,pC,task)-metric(y,pR,task),'score_B':metric(y,pR,task),'score_C':metric(y,pC,task),'pred_B':pR.tolist(),'pred_C':pC.tolist(),'pred_group':pG.tolist(),'observed_y':y.tolist(),'scope_folds':len(gs),'split_folds':len(gs)}

def generate(X,g,c,instance_seed,L=0,eta=0,encoding=False):
 # Streams depend on instance only; levels are paired; encoder independent of outcome offsets.
 z=(X-X.mean(0))/(X.std(0)+1e-12);n,p=z.shape
 rr=np.random.default_rng(seed(instance_seed,'architecture'));idx=rr.choice(p,min(10,p),replace=False);w=rr.normal(size=len(idx))
 s=z[:,idx]@w if c.get('architecture','linear')=='linear' else np.sin(z[:,idx[:5]]).sum(1)+z[:,idx[0]]*z[:,idx[1]]
 s=(s-s.mean())/(s.std()+1e-12)
 _,gi=np.unique(g,return_inverse=True);ng=len(set(g));mu=np.random.default_rng(seed(instance_seed,'group_outcome')).normal(size=ng);ss=mu[gi];ss=(ss-ss.mean())/(ss.std()+1e-12)
 signal=c.get('signal',0)*s;noise_rng=np.random.default_rng(seed(instance_seed,'outcome_noise'))
 if c['task']=='binary':
  logit=signal+eta*ss;prev=c.get('prevalence',.5);intercept=brentq(lambda a:expit(logit+a).mean()-prev,-40,40);prob=expit(logit+intercept);y=(noise_rng.random(n)<prob).astype(int)
 else:
  eps=noise_rng.normal(size=n)
  if c.get('noise')=='heavy':eps=noise_rng.standard_t(3,size=n)/np.sqrt(3)
  if c.get('noise')=='heteroscedastic':eps=eps*(.5+np.abs(s));eps=eps/np.sqrt(np.mean((.5+np.abs(s))**2))
  y=signal+eta*ss+eps
 features=[X]
 if L:features.append(np.random.default_rng(seed(instance_seed,'candidates')).normal(size=(2500,n)).T[:,:L])
 if encoding:
  center=np.random.default_rng(seed(instance_seed,'encoding_centers')).normal(size=(ng,10));noise=np.random.default_rng(seed(instance_seed,'encoding_noise')).normal(scale=.2,size=(n,10));features.append(center[gi]+noise)
 xa=np.column_stack(features)
 return xa,y,{'realized_y_mean':float(y.mean()),'realized_y_sd':float(y.std()),'base_signal_sd':float(signal.std()),'group_signal_sd':float((eta*ss).std()),'n':n,'p_original':p,'p_total':xa.shape[1],'L_over_p':L/p,'k_over_p_total':c.get('k',20)/xa.shape[1]}

def threshold(v,alpha=.05):
 v=np.asarray(v);j=math.ceil((len(v)+1)*(1-alpha));return float(np.sort(v)[j-1]) if j<=len(v) else float('inf')
def wilson(k,n,z=1.959963984540054):
 if n==0:return (float('nan'),float('nan'))
 p=k/n;d=1+z*z/n;m=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n));return ((p+z*z/(2*n)-m)/d,(p+z*z/(2*n)+m)/d)

def provenance_status(fitlog):
 if fitlog is None:return 'INSUFFICIENT_PIPELINE_PROVENANCE'
 violations=[set(r['fit_ids'])&set(r['test_ids']) for r in fitlog]
 return 'PROHIBITED_FIT_SCOPE_OBSERVED' if any(violations) else 'NO_PROHIBITED_OVERLAP_IN_RECORDED_SCOPE'
