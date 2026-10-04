"""Actual model fitting entry point; each result is a new reference computation."""
from pathlib import Path
import argparse,json,time,sys,os,platform,hashlib,traceback
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from core import seed,scope_stat,split_stat,generate
CACHE={}
def init(data):
 global CACHE
 for n in ['wheat','microbiome_crc','rice']:
  p=Path(data)/f'{n}.npz'
  if p.exists():
   z=np.load(p);CACHE[n]=(z['X'],z['g'])
def work(job):
 t=time.time();c=job['condition'];base=c['base'];X,g=CACHE[base]
 if c.get('sample_limit'):
  ix=np.sort(np.random.default_rng(seed(base,'sample_subset')).choice(len(X),c['sample_limit'],replace=False));X=X[ix];g=g[ix]
 if c.get('feature_limit'):
  ix=np.sort(np.random.default_rng(seed(base,'feature_subset')).choice(X.shape[1],c['feature_limit'],replace=False));X=X[:,ix]
 out={**job,'condition_id':c['id'],'base':base,'model':c['family'],'task':c['task'],'k':c.get('k',20),'folds':c.get('folds',5),'status':'RUNNING'}
 try:
  xa,y,realized=generate(X,g,c,job['instance_seed'],job['L'],job.get('eta',0),job.get('encoding',False))
  out.update(realized)
  if job['experiment']=='E4':
   out.update(split_stat(xa,y,g,c['task'],c['family'],c.get('k',20),job['split_seed']))
  else:
   out.update(scope_stat(xa,y,c['task'],c['family'],c.get('k',20),c.get('folds',5),job['split_seed']))
  out['status']='COMPLETE'
 except Exception as e:out.update(status='FAILED',failure_reason=repr(e),traceback=traceback.format_exc())
 out['wall_seconds']=time.time()-t;return out

def jobs(config,pilot=False):
 out=[];prefix='pilot_v1' if pilot else 'confirmation_v1'
 for c in config['conditions']:
  if pilot:counts={'calibration':2,'evaluation':2,'alternative':1}
  else:counts=c['counts']
  for role,n in counts.items():
   levels=c.get('alternative_L',[2500]) if role=='alternative' else [0]
   for i in range(n):
    # Alternative levels paired; calibration/evaluation/alternative IDs disjoint.
    rid=f"{prefix}_E3_{c['id']}_{role}_{i:04d}"
    for L in levels:
     out.append({'experiment':'E3','condition':c,'role':role,'reference_id':rid,'instance_id':rid+f'_L{L}','instance_seed':seed(rid),'split_seed':seed(rid,'split'),'L':L,'eta':0,'encoding':False})
 if config.get('mixed'):
  for base in ['wheat','microbiome_crc']:
   for family in ['linear','rf']:
    c={'id':f'{base}_{family}_mixed','base':base,'task':'regression' if base=='wheat' else 'binary','family':family,'k':20,'folds':5,'signal':.5,'architecture':'linear','prevalence':.5}
    for i in range(1 if pilot else config['mixed']['replicates']):
     rid=f'{prefix}_E4_{base}_{i:04d}'
     for L in config['mixed']['L']:
      for eta in config['mixed']['eta']:
       for enc in [False,True]:out.append({'experiment':'E4','condition':c,'role':'evaluation','reference_id':rid,'instance_id':rid+f'_{family}_L{L}_eta{eta}_enc{int(enc)}','instance_seed':seed(rid),'split_seed':seed(rid,'split'),'L':L,'eta':eta,'encoding':enc})
 return out

def main():
 a=argparse.ArgumentParser();a.add_argument('--config',required=True);a.add_argument('--data',required=True);a.add_argument('--out',required=True);a.add_argument('--workers',type=int,default=4);a.add_argument('--pilot',action='store_true');args=a.parse_args()
 out=Path(args.out)
 if out.exists():raise SystemExit('Refusing to overwrite run directory')
 out.mkdir(parents=True);config=json.loads(Path(args.config).read_text());todo=jobs(config,args.pilot)
 import sklearn,scipy,pandas
 meta={'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'command':sys.argv,'pilot':args.pilot,'config':config,'config_sha256':hashlib.sha256(Path(args.config).read_bytes()).hexdigest(),'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob('*.py')},'data_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(args.data).glob('*.npz')},'platform':platform.platform(),'python':sys.version,'numpy':np.__version__,'sklearn':sklearn.__version__,'scipy':scipy.__version__,'pandas':pandas.__version__,'workers':args.workers,'n_attempted':len(todo),'status':'RUNNING'}
 (out/'run_manifest.json').write_text(json.dumps(meta,indent=2));start=time.time();failed=0
 with (out/'raw_results.jsonl').open('w') as f,ProcessPoolExecutor(max_workers=args.workers,initializer=init,initargs=(args.data,)) as pool:
  futures=[pool.submit(work,j) for j in todo]
  for i,fu in enumerate(as_completed(futures),1):
   result=fu.result();failed+=result['status']!='COMPLETE';f.write(json.dumps(result,allow_nan=False)+'\n');f.flush()
   if i%50==0 or i==len(todo):print(f'{i}/{len(todo)} complete; failures={failed}; elapsed={time.time()-start:.1f}s',flush=True)
 meta.update(status='COMPLETE' if failed==0 else 'FAILED',n_failed=failed,wall_seconds=time.time()-start,finished_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),raw_sha256=hashlib.sha256((out/'raw_results.jsonl').read_bytes()).hexdigest());(out/'run_manifest.json').write_text(json.dumps(meta,indent=2))
 if failed:raise SystemExit(1)
if __name__=='__main__':main()
