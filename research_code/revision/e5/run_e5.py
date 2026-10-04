"""Frozen authentic-outcome, conditional split-randomization experiment."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import argparse,hashlib,json,platform,sys,datetime
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import scipy,sklearn,pyreadr
from sklearn.model_selection import KFold
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

X=Y=None
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def seed(role,i,stream):
    s=f'patterns-v14-e5|20260908|{role}|{i}|{stream}'
    return int.from_bytes(hashlib.sha256(s.encode()).digest()[:8],'little')%2**32
def initialize(raw):
    global X,Y
    d=pyreadr.read_r(raw);X=d['wheat.X'].to_numpy(float);Y=d['wheat.Y'].to_numpy(float)
    assert X.shape==(599,1279) and Y.shape==(599,4)
    assert list(d['wheat.Y'].columns)==['1','2','4','5']
    assert np.isfinite(X).all() and np.isfinite(Y).all()
def correlations(a,y):
    a=a-a.mean(axis=0);y=y-y.mean(axis=0)
    den=np.sqrt((a*a).sum(axis=0)[:,None]*(y*y).sum(axis=0)[None,:])
    return np.divide(np.abs(a.T@y),den,out=np.zeros_like(den),where=den>0)
def reference(arg):
    role,i,out=arg
    fs=seed(role,i,'folds');cs=seed(role,i,'candidates')
    folds=list(KFold(5,shuffle=True,random_state=fs).split(X))
    levels=[0] if role=='calibration' else [0,400,1000,2500]
    a=X if role=='calibration' else np.column_stack((X,np.random.default_rng(cs).normal(size=(2500,599)).T))
    fullcorr=correlations(a,Y)
    predictions=np.zeros((len(levels),4,2,599));selected=np.zeros((len(levels),4,2,5,20),dtype=np.int32)
    foldids=np.zeros(599,dtype=np.int16)
    for f,(tr,te) in enumerate(folds):
        foldids[te]=f
        at=a[tr];var=at.var(axis=0);bcorr=correlations(at,Y[tr])
        for li,L in enumerate(levels):
            eligible=np.flatnonzero(var[:1279+L]>1e-12)
            for ph in range(4):
                for arm,corr in enumerate((bcorr,fullcorr)):
                    ids=eligible[np.argsort(-corr[eligible,ph],kind='stable')[:20]]
                    selected[li,ph,arm,f]=ids
                    xt=a[np.ix_(tr,ids)];mu=xt.mean(0);sd=xt.std(0)
                    model=Ridge(alpha=1,solver='cholesky').fit((xt-mu)/sd,Y[tr,ph])
                    predictions[li,ph,arm,te]=model.predict((a[np.ix_(te,ids)]-mu)/sd)
    ref=f'{role}_{i:03d}'
    np.savez_compressed(Path(out)/'predictions'/f'{ref}.npz',predictions=predictions,selected=selected,fold_ids=foldids,levels=np.array(levels),y=Y)
    rows=[]
    for li,L in enumerate(levels):
        for ph,label in enumerate(['1','2','4','5']):
            b=r2_score(Y[:,ph],predictions[li,ph,0]);c=r2_score(Y[:,ph],predictions[li,ph,1])
            rows.append(dict(reference_id=ref,role=role,replicate=i,phenotype=label,L=L,R2_B=b,R2_C=c,T_scope=c-b,fold_seed=fs,candidate_seed=cs,n=599,model='ridge_alpha1',status='COMPLETE'))
    return rows
def main():
    p=argparse.ArgumentParser();p.add_argument('--raw',required=True);p.add_argument('--protocol',required=True);p.add_argument('--freeze',required=True);p.add_argument('--out',required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=False);(out/'predictions').mkdir()
    freeze=json.loads(Path(a.freeze).read_text());assert sha(a.raw)==freeze['raw_sha256'];assert sha(a.protocol)==freeze['protocol_sha256']
    env=dict(python=sys.version,platform=platform.platform(),architecture=platform.machine(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,pyreadr=pyreadr.__version__)
    manifest=dict(status='RUNNING',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),environment=env,code_sha256=sha(__file__),raw_sha256=sha(a.raw),protocol_sha256=sha(a.protocol),reference_n=200,phenotype_n=4,planned_score_rows=2000,planned_model_fits=20000,failed=[])
    (out/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    jobs=[(role,i,str(out)) for role in ['calibration','evaluation'] for i in range(100)];rows=[]
    try:
        with ProcessPoolExecutor(a.workers,initializer=initialize,initargs=(a.raw,)) as pool:
            for j,result in enumerate(pool.map(reference,jobs)):
                rows.extend(result)
                if (j+1)%20==0: print(f'Completed {j+1}/200 reference realizations',flush=True)
        df=pd.DataFrame(rows).sort_values(['role','replicate','phenotype','L']);df.to_csv(out/'scores.csv',index=False,float_format='%.17g')
        assert len(df)==2000 and len(df.reference_id.unique())==200
        manifest.update(status='COMPLETE',completed_model_fits=20000,successful_score_rows=len(df),scores_sha256=sha(out/'scores.csv'))
    except Exception as e:
        manifest['status']='FAILED';manifest['failed'].append(repr(e));raise
    finally:
        manifest['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(out/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
