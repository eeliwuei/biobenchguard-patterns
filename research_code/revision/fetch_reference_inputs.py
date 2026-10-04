"""Acquire exact-hash wheat/rice/UCI inputs; no third-party raw data bundled."""
from pathlib import Path
import argparse,hashlib,io,json,tarfile,urllib.request,zipfile
from datetime import datetime,timezone

SOURCES=[
 ('crop/data/wheat.RData','https://raw.githubusercontent.com/gdlc/BGLR-R/master/data/wheat.RData','8710523389007dd8b7c88a5deb74ee16aea0bae158496bc11f0626fb43b5ea21','direct',None,None),
 ('rice/data/RiceDiversity_44K_Genotypes_PLINK/sativas413.ped','https://raw.githubusercontent.com/tony-travis/PIQUE/097a945f1a715078a9f877469ea60f7497f7375d/test/RiceDiversity_44K_Genotypes_PLINK.tar.xz','b25cb0225bbc61fba0a4430864d740fda0af721e224f32768bbc832b222b1b0a','tar','RiceDiversity_44K_Genotypes_PLINK/sativas413.ped','0bdd54e44a831b197fcf0bdf5d1cec7055c721cb5ab94d7c33a005e889ab62a1'),
 ('rice/data/germplasm.csv','https://web.archive.org/web/20180408144810id_/http://ricediversity.org:80/data/sets/44kgwas/RiceDiversity.44K.germplasm.csv','5cfc7bdcdc755581366eb7c66ce0f576a9013e04fa6884aed566cd18b9a30df6','direct',None,None),
 ('e1/parkinsons.data','https://archive.ics.uci.edu/static/public/174/parkinsons.zip','a36116c2deac07895b79ce97a2f99a1b3720ac4892995ff5c3b9ccb397ed2bc5','zip','parkinsons.data','9baccfbcdac15d6d9347b37c2a352fb31d65c52e367df5a47d7543a55c75b0ab')]
def digest(data):return hashlib.sha256(data).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',required=True);a=p.parse_args();root=Path(a.raw_root);root.mkdir(parents=True,exist_ok=True);records=[]
 for rel,url,expected,kind,member,archive_expected in SOURCES:
  dest=root/rel
  if dest.exists():
   if not dest.is_file() or digest(dest.read_bytes())!=expected:raise ValueError(f'Refusing inconsistent existing input: {rel}')
   records.append(dict(path=rel,status='EXISTING_HASH_VERIFIED',sha256=expected));continue
  payload=urllib.request.urlopen(url,timeout=60).read()
  if archive_expected and digest(payload)!=archive_expected:raise ValueError(f'Archive hash mismatch: {rel}')
  if kind=='tar':
   with tarfile.open(fileobj=io.BytesIO(payload),mode='r:xz') as archive:
    stream=archive.extractfile(member)
    if stream is None:raise ValueError(f'Required regular archive member absent: {member}')
    data=stream.read()
  elif kind=='zip':
   with zipfile.ZipFile(io.BytesIO(payload)) as archive:data=archive.read(member)
  else:data=payload
  if digest(data)!=expected:raise ValueError(f'Input hash mismatch: {rel}')
  dest.parent.mkdir(parents=True,exist_ok=True)
  with dest.open('xb') as f:f.write(data)
  records.append(dict(path=rel,url=url,bytes=len(data),sha256=expected,status='DOWNLOADED_HASH_VERIFIED'));print(f'Verified {rel}',flush=True)
 receipt=dict(retrieved_utc=datetime.now(timezone.utc).isoformat(),sources=records,scope='Acquisition does not relicense source data; retain source terms.')
 (root/'reference_acquisition_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':main()
