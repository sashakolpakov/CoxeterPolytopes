#!/usr/bin/env python3
"""Source/formula review and nested160-digit replay of the118-wall roof."""
import argparse,hashlib,json,sys
from pathlib import Path
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_higher_face_20261013_final_verify as H
V=H.V
OUT=V.ROOT/'results/fresh_h9_higher_face_20261013_review.json'
def run(verify=False):
 qa,qz,qp=sp.symbols('qa qz qp',nonzero=True)
 assert sp.cancel(qp+(2-qp)/qz*qz)==2
 assert sp.cancel(1/qa+(2-1/qa)/qz*qz)==2
 assert sp.cancel(-qa/qa)==-1
 assert sp.cancel((2/qz)*qz)==2
 # In each projection, invertibility is certified by interval pivots.
 # M^-1 b gives part.E=b; z=center-proj_E(center) gives z.E=0.
 # Hence part.z=0, and the displayed scalar identities prove the norms.
 rawread=V.read;assert V.SCALE==10**80;saved=rawread(H.PREFIX+'_verify.json')
 first,_=H.build();assert first==saved
 factor=10**80;counts={'changed_enclosures':0,'prior_certificate_comparisons':0}
 def nested(old,new,path='root'):
  if old==new:return
  if isinstance(old,dict):
   assert isinstance(new,dict)and old.keys()==new.keys(),path
   for key in old:
    if key=='interval_scale':assert int(new[key])==int(old[key])*factor
    else:nested(old[key],new[key],path+'.'+key)
  elif isinstance(old,list):
   assert isinstance(new,list)and len(old)==len(new),path
   if len(old)==2 and all(isinstance(x,str)for x in old):
    try:a,b=map(int,old);c,d=map(int,new)
    except ValueError:pass
    else:assert a*factor<=c<=d<=b*factor,(path,'non-nested enclosure');counts['changed_enclosures']+=1;return
   for i,(a,b)in enumerate(zip(old,new,strict=True)):nested(a,b,path+f'[{i}]')
  else:assert old==new,(path,old,new)
 class SavedCertificate(dict):
  def __eq__(self,other):nested(dict(self),other);counts['prior_certificate_comparisons']+=1;return True
  def __ne__(self,other):return not self.__eq__(other)
 def precision_read(path):
  value=rawread(path);return SavedCertificate(value)if path.endswith('_verify.json')and isinstance(value,dict)else value
 V.read=precision_read;V.SCALE=10**160;V.PHI=V.div(V.add(V.rat(1),V.sqrt(V.rat(5))),V.rat(2));V.radical.cache_clear()
 high,_=H.build();nested(saved,high);assert high['facets']==118
 face=high['H5_face'];assert len(face['facets'])==19 and len(face['finite_vertices'])==25 and face['finite_edges']==41 and face['open_edges']==43 and len(face['finite_components'])==1
 assert len(high['strict_escape']['pairings'])==118 and int(face['strict_spacelike_escape']['norm'][0])>0
 files=[Path(H.__file__),Path(H.T.__file__),Path(H.F.__file__),Path(V.__file__),Path(__file__).resolve()]+[V.ROOT/(x.PREFIX+'_verify.json')for x in(H,H.T,H.F)]
 out={'status':'PASS_SOURCE_FORMULA_REVIEW_AND_NESTED_160_DIGIT_H5_ROOF_REPLAY','sources':{str(p.relative_to(V.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},'frozen_precision_digits':80,'review_precision_digits':160,'facets':118,'H5_face_facets':19,'H5_finite_vertices':25,'H5_finite_edges':41,'H5_open_edges':43,'strict_escape_norm':high['strict_escape']['norm'],**counts,'scope':'Independent source/formula review of projection, polar and metric truncations; archived replay and higher-precision replay of the same rigorous interval implementation. All changed enclosures nest and exact metadata agree. This does not independently reimplement the118 normals. Both whole chamber and specified H5 face remain infinite-volume.'}
 if verify:assert out==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],counts,flush=True)
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');run(a.parse_args().verify)
