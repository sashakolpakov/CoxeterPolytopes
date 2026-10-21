#!/usr/bin/env python3
"""Finite-angle endpoint proposals on the certified 94-wall tetrahedral frontier."""
import json,sys,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_compact_next_20261013_verify as C
V=C.V;K=C.K
PREFIX='data/enumerations/fresh_h9_face_bridge_20261013'

def setup():
 d=V.read(C.PREFIX+'_caps.json');G,I,h,g,old,names,caps,center,_=C.reconstruct(d)
 ff=lambda x:sum(x)/(2*V.SCALE)
 R=np.array([[ff(x)for x in r]for r in old+caps]);names+=C.NEW
 Q=np.array([[ff(V.golden(x))for x in row]for row in G]);hh=np.array([ff(V.golden(x))for x in h]);cc=np.array([ff(x)for x in center])
 y=np.array([ff(V.golden(K.parse(x)))for x in d['escape']])
 prior=V.read(C.PREFIX+'_verify.json')
 return Q,R,names,hh,cc,y,prior

def run():
 Q,R,names,h,c,y,prior=setup();cov=R@Q;rows=[];stats=[];seen=set();labels=(2,3,4,5,6,8,10,12)
 edge_rows=[e for e in prior['new_vertex_edge_frontier']['edges']if e['status']=='open_to_light_cone']
 for ei,edge in enumerate(edge_rows):
  ids=[names.index(n)for n in edge['common_facets']];M=R[ids]@Q@R[ids].T;Mi=np.linalg.inv(M)
  profiles=[np.zeros(8)]
  for i in range(8):
   for m in(3,4,5,6,8,10,12):
    p=np.zeros(8);p[i]=-2*np.cos(np.pi/m)
    if p@Mi@p<2-1e-8:profiles.append(p)
  for i,j in itertools.combinations(range(8),2):
   for m,n in itertools.product((3,4,5),repeat=2):
    p=np.zeros(8);p[i]=-2*np.cos(np.pi/m);p[j]=-2*np.cos(np.pi/n)
    if p@Mi@p<2-1e-8:profiles.append(p)
  stat={'edge':ei,'profiles':len(profiles),'quadratics':0,'compatible':0,'target_cut':0};pp=np.array(profiles)
  for ai in range(len(R)):
   if ai in ids:continue
   A=cov[ids+[ai]];u,s,v=np.linalg.svd(A)
   if s[-1]<1e-8:continue
   inv=A.T@np.linalg.inv(A@A.T);direction=v[-1];qa=direction@Q@direction
   if abs(qa)<1e-10:continue
   for label in labels:
    rhs=np.column_stack([pp,np.full(len(pp),0 if label==2 else -2*np.cos(np.pi/label))]);bases=rhs@inv.T
    qb=2*bases@Q@direction;qc=np.einsum('ij,jk,ik->i',bases,Q,bases)-2;disc=qb*qb-4*qa*qc;stat['quadratics']+=len(profiles)
    for pi in np.flatnonzero(disc>1e-10):
     for branch in(-1,1):
      t=(-qb[pi]+branch*np.sqrt(disc[pi]))/(2*qa);r=bases[pi]+t*direction
      if r@h<=1e-8:continue
      pair=r@cov.T
      if pair.max()>1e-6:continue
      finite=np.flatnonzero(pair>-2+1e-6);labs={};valid=True
      for j in finite:
       if abs(pair[j])<1e-6:labs[int(j)]=2;continue
       z=next((m for m in labels[1:]if abs(pair[j]+2*np.cos(np.pi/m))<1e-6),None)
       if z is None:valid=False;break
       labs[int(j)]=z
      if not valid:continue
      # Parallel pairs and affine circuits are handled by the subsequent gate.
      key=tuple(np.round(r,7))
      if key in seen:continue
      seen.add(key);stat['compatible']+=1;cut=-r@y;stat['target_cut']+=int(cut>1e-8)
      rows.append({'edge_index':ei,'edge':edge,'profile_pairings':pp[pi].tolist(),'profile_labels':[2 if abs(x)<1e-8 else min(labels[1:],key=lambda m:abs(x+2*np.cos(np.pi/m)))for x in pp[pi]],'anchor':names[ai],'anchor_label':label,'normal':r.tolist(),'finite_profile':{names[j]:m for j,m in labs.items()},'parallel_walls':[names[j]for j in np.flatnonzero(abs(pair+2)<1e-6)],'target_cut':float(cut),'center':float(r@h)})
  stats.append(stat);print(stat,flush=True)
  out={'status':'NUMERICAL_FINITE_ENDPOINT_BRIDGE_SCOUT','source':C.PREFIX+'_caps.json','source_sha256':V.sha(C.PREFIX+'_caps.json'),'stats':stats,'candidates':rows,'scope':'24 certified open edges; finite spherical profiles with up to two nonright incident pairings; one additional old anchor angle; numeric proposals only.'}
  (V.ROOT/(PREFIX+'_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
 print('DONE',len(rows),'CUTTERS',sum(r['target_cut']>1e-8 for r in rows),flush=True)

if __name__=='__main__':run()
