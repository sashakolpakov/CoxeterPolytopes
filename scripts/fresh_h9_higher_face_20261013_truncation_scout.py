#!/usr/bin/env python3
"""Positive two-plane truncations of the new H5 roof's hyperideal vertices."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.linalg import null_space
from scripts import fresh_h9_higher_face_20261013_linked_scout as L
S=L.S

def run():
 _,_,_,_,_,names,_,_,R,G,c=S.load();d=S.V.read(S.PREFIX+'_selected_scout.json');rr=np.array([x['normal']for x in d['recipes']]);R=np.vstack([R,rr]);names=names+[f'roof{i}'for i in range(len(rr))];fixed=[n for n in S.COMMON if n!='seed1'];pd=S.V.read(S.PREFIX+'_polar_scout.json');out=[];stats=[];seen=set()
 for vi,row in enumerate(sorted(pd['rows'],key=lambda x:-x['norm'])):
  inc=row['incident'];x=np.array(row['point']);anchors=[z['wall']for z in row['bad']]
  if not anchors:continue
  naccept=0
  for release in inc:
   if release in fixed:continue
   edge=[n for n in inc if n!=release]
   if len(edge)!=8:continue
   E=R[[names.index(n)for n in edge]];B=null_space(E@G);M=B.T@G@B
   if np.linalg.eigvalsh(M)[0]<1e-8:continue
   B=B@np.linalg.inv(np.linalg.cholesky(M).T)
   for anchor in anchors:
    a=R[names.index(anchor)];p=B.T@G@a;q=p@p
    if q<1e-10:continue
    z=np.array([-p[1],p[0]])/np.sqrt(q)
    for m in(2,3,4,5,6,8,10,12,16,24,32):
     want=-2*np.cos(np.pi/m);rad=2-want*want/q
     if rad<1e-10:continue
     for branch in(-1,1):
      r=B@(want/q*p+branch*np.sqrt(rad)*z);ps=R@G@r
      if r@G@x<1e-7 or r@G@c>=0 or not all(L.gate(t)for t in ps):continue
      key=tuple(np.round(r,8))
      if key in seen:continue
      seen.add(key);out.append({'vertex':vi,'target':row['point'],'right':edge,'release':release,'anchor':anchor,'label':m,'branch':branch,'normal':r.tolist(),'target_cut':float(r@G@x)});naccept+=1
  stats.append({'vertex':vi,'accepted':naccept,'target_norm':row['norm']})
  print(stats[-1],flush=True)
 (S.V.ROOT/(S.PREFIX+'_truncation_scout.json')).write_text(json.dumps({'status':'NUMERICAL_NEW_ROOF_TRUNCATIONS','source':S.PREFIX+'_polar_scout.json','stats':stats,'candidates':out},indent=2)+'\n');print('total',len(out),flush=True)
if __name__=='__main__':run()
