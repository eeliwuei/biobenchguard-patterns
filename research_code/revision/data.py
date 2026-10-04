"""Revision loaders. Raw feature matrices returned; fold-specific fitting lives in pipelines.
No biological outcome is synthesized or substituted by this loader.
"""
from pathlib import Path
from collections import defaultdict
import hashlib,json
import numpy as np
import pandas as pd

DEFAULT_ROOT=Path(__file__).resolve().parents[2]/'data_manifest'/'raw'

def sha256(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def load_wheat(raw):
 import pyreadr
 from sklearn.decomposition import PCA
 from sklearn.cluster import KMeans
 p=Path(raw)/'crop/data/wheat.RData';d=pyreadr.read_r(str(p))
 X=d['wheat.X'].to_numpy(dtype=float);Y=d['wheat.Y'];y=Y.iloc[:,0].to_numpy(dtype=float)
 pcs=PCA(n_components=10,svd_solver='full').fit_transform(X-X.mean(0))
 g=KMeans(n_clusters=5,n_init=10,random_state=20260908).fit_predict(pcs).astype(str)
 return X,y,g,{'dataset':'wheat','task':'regression','n':len(y),'p':X.shape[1],'features':list(d['wheat.X'].columns),'sample_ids':list(map(str,Y.index)),'phenotype':'wheat.Y column 1; standardized grain yield in BGLR environment 1','group_source':'revision imposed KMeans5 on PCA10 of X; descriptive geometry, no natural group claim','raw_files':{str(p.relative_to(raw)):sha256(p)}}

def genus_of(t):
 parts=t.split(';');g=next((p[3:] for p in parts if p.startswith('g__')),'')
 if g:return 'g__'+g
 f=next((p[3:] for p in parts if p.startswith('f__')),'')
 return 'f__'+f+';g__?' if f else 'unclassified'

def load_crc(raw):
 frames=[];ys=[];gs=[];ids=[];hashes={};counts={}
 for st in ['crc_zeller','crc_baxter','crc_zackular']:
  base=Path(raw)/'microbiome/data'/f'{st}_results'
  p=base/'RDP'/f'{st}.otu_table.100.denovo.rdp_assigned';mp=base/f'{st}.metadata.txt'
  hashes[str(p.relative_to(raw))]=sha256(p);hashes[str(mp.relative_to(raw))]=sha256(mp)
  with p.open(errors='strict') as f:
   samples=f.readline().rstrip('\n').split('\t')[1:];gc=defaultdict(lambda:np.zeros(len(samples)))
   for line in f:
    v=line.rstrip('\n').split('\t')
    if len(v)!=len(samples)+1:raise ValueError('OTU row shape mismatch')
    gc[genus_of(v[0])]+=np.asarray(v[1:],float)
  G=pd.DataFrame(gc,index=samples);tot=G.sum(1)
  if (tot<=0).any():raise ValueError('Zero-read sample; explicit exclusion required')
  rel=G.div(tot,axis=0)
  md=pd.read_csv(mp,sep='\t',dtype=str,encoding='latin-1');idcol=max(md.columns,key=lambda c:md[c].isin(set(samples)).sum())
  if md[idcol].duplicated().any():raise ValueError('Duplicate metadata ID')
  md=md.set_index(idcol);dc=next(c for c in md if c.lower() in ['diseasestate','disease_state','disease'])
  lab=md[dc].reindex(rel.index);keep=lab.isin(['H','CRC']);rel=rel.loc[keep];lab=lab.loc[keep]
  y=(lab=='CRC').astype(int).to_numpy();frames.append(rel);ys.append(y);gs.extend([st]*len(y));ids.extend([st+':'+str(i) for i in rel.index])
  counts[st]={'n':len(y),'CRC':int(y.sum()),'H':int((y==0).sum()),'excluded_other_or_missing_outcome':int((~keep).sum())}
 shared=sorted(set.intersection(*[set(f.columns) for f in frames]));X=np.vstack([f[shared].to_numpy() for f in frames]);y=np.concatenate(ys)
 return X,y,np.array(gs),{'dataset':'microbiome_crc','task':'binary','n':len(y),'p':len(shared),'features':shared,'sample_ids':ids,'study_counts':counts,'phenotype':'original DiseaseState H versus CRC','group_source':'original study identity','raw_files':hashes,'processing':'genus sum; row-wise relative abundance; shared taxonomy vocabulary, no label-dependent filtering; no global learned scaling'}

def load_rice(raw):
 base=Path(raw)/'rice/data';ped=base/'RiceDiversity_44K_Genotypes_PLINK/sativas413.ped';germ=base/'germplasm.csv'
 md=pd.read_csv(germ,skiprows=1,on_bad_lines='skip');idcol=next(c for c in md if 'nsftv' in c.lower());spcol=next(c for c in md if c.lower().startswith('sub-pop'))
 md[idcol]=pd.to_numeric(md[idcol],errors='coerce');mp=dict(zip(md[idcol],md[spcol]))
 rows=[l.split() for l in ped.open() if l.strip()];ids=[int(r[1]) for r in rows];a=np.array([r[6:] for r in rows]).reshape(len(rows),-1,2)
 X=np.full(a.shape[:2],np.nan)
 for j in range(a.shape[1]):
  q=a[:,j,:];vals,c=np.unique(q[q!='0'],return_counts=True)
  if len(vals)<2:continue
  minor=vals[np.argmin(c)];v=(q==minor).sum(1).astype(float);v[(q=='0').any(1)]=np.nan;X[:,j]=v
 g=np.array([str(mp.get(i,'NA')) for i in ids]);good=np.array([s not in ['NA','nan'] for s in g]);X=X[good];g=g[good];ids=np.array(ids)[good]
 # Match historical min-group>=10 population restriction, recorded explicitly.
 vc=pd.Series(g).value_counts();good=np.isin(g,vc[vc>=10].index);X=X[good];g=g[good];ids=ids[good]
 keep=np.isfinite(X).any(0)&(np.nanstd(X,axis=0)>1e-9);X=X[:,keep]
 # Geometry-only deterministic mean imputation; revision real tasks never use this matrix.
 mean=np.nanmean(X,0);i,j=np.where(~np.isfinite(X));X[i,j]=mean[j]
 return X,np.full(len(X),np.nan),g,{'dataset':'rice','task':'regression','n':len(X),'p':X.shape[1],'sample_ids':ids.astype(str).tolist(),'group_source':'original subpopulation, restricted to groups with >=10 entries','phenotype':'not loaded; reference geometry only','raw_files':{str(ped.relative_to(raw)):sha256(ped),str(germ.relative_to(raw)):sha256(germ)}}

def load_base(name,raw):
 return {'wheat':load_wheat,'microbiome_crc':load_crc,'rice':load_rice}[name](Path(raw))

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--raw',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 for name in ['wheat','microbiome_crc','rice']:
  X,y,g,m=load_base(name,a.raw);np.savez_compressed(out/f'{name}.npz',X=X,y=y,g=g);(out/f'{name}.json').write_text(json.dumps(m,indent=2));print(name,X.shape,len(set(g)),flush=True)
