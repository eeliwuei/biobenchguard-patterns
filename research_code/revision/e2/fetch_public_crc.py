"""Fetch only the six CRC source members, verifying original Zenodo archive MD5.

Example: python fetch_public_crc.py --destination /new/revision/data_manifest/raw
Creates microbiome/data/<study>_results under the supplied destination, refuses
to overwrite a differing existing file, and retains a download provenance JSON.
Data license: CC BY-NC 4.0 (Zenodo record 840333), separate from code license.
"""
import argparse, hashlib, io, json, tarfile, urllib.request
from pathlib import Path

ap=argparse.ArgumentParser();ap.add_argument('--destination',required=True);a=ap.parse_args()
dest=Path(a.destination);dest.mkdir(parents=True,exist_ok=True)
url='https://zenodo.org/api/records/840333'
with urllib.request.urlopen(url,timeout=60)as r:record=json.load(r)
assert record['id']==840333 and record['metadata']['license']['id']=='cc-by-nc-4.0'
manifest={'record':url,'license':record['metadata']['license'],'files':[]}
for study in ['crc_zeller','crc_baxter','crc_zackular']:
    name=study+'_results.tar.gz';entry=next(x for x in record['files']if x['key']==name)
    with urllib.request.urlopen(entry['links']['self'],timeout=120)as r:payload=r.read()
    assert 'md5:'+hashlib.md5(payload).hexdigest()==entry['checksum']
    wanted=[study+'_results/'+study+'.metadata.txt',study+'_results/RDP/'+study+'.otu_table.100.denovo.rdp_assigned']
    with tarfile.open(fileobj=io.BytesIO(payload),mode='r:gz')as archive:
        for member_name in wanted:
            member=archive.getmember(member_name);assert member.isfile()
            data=archive.extractfile(member).read();path=dest/'microbiome/data'/member_name
            path.parent.mkdir(parents=True,exist_ok=True)
            digest=hashlib.sha256(data).hexdigest()
            if path.exists():assert hashlib.sha256(path.read_bytes()).hexdigest()==digest,'Existing input differs; choose a fresh destination'
            else:path.write_bytes(data)
            manifest['files'].append({'relative_path':str(path.relative_to(dest)),'sha256':digest,'archive':name,'archive_md5':entry['checksum'],'url':entry['links']['self']})
    print('Verified',name,flush=True)
(dest/'microbiome_public_fetch_manifest.json').write_text(json.dumps(manifest,indent=2))
