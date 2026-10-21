#!/usr/bin/env python3
"""Joint finite-side endpoint search over the eight compact prism vertices."""
import json,sys
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.linalg import null_space
from scipy.spatial import HalfspaceIntersection
from scipy.optimize import linprog
import z3
from scripts import fresh_h9_higher_face_20261013_scout as S
T=['seed6','crosspair0','nextH4','compact0'];BASES=['wall17','facePrism']
LABELS=list(range(3,33))+[40,48,60,64,96,128]
def gate(p):
 if p>1e-7:return False
 if p<-2+1e-8 or abs(p)<1e-7:return True
 ang=np.pi/np.arccos(np.clip(-p/2,-1,1));m=int(round(ang))
 return m>=2 and abs(p+2*np.cos(np.pi/m))<1e-7
def geom(normals,Gn,cn,fixed):
 F=normals[fixed];proj=cn-(cn@Gn@F.T)@np.linalg.inv(F@Gn@F.T)@F;proj/=np.sqrt(-proj@Gn@proj)
 B=null_space(np.vstack([F@Gn,proj@Gn]));M=B.T@Gn@B;B=B@np.linalg.inv(np.linalg.cholesky(M).T)
 A=normals@Gn@B;b=normals@Gn@proj;keep=np.linalg.norm(A,axis=1)>1e-7;H=np.c_[A[keep],b[keep]]
 try:
  vv=HalfspaceIntersection(H,np.zeros(B.shape[1]),qhull_options='Qx').intersections
  q=(vv*vv).sum(axis=1)-1
  # Separate bounding LPs guard Qhull's possible output on unbounded input.
  bounded=all(linprog(s*e,A_ub=H[:,:-1],b_ub=-H[:,-1],bounds=[(None,None)]*len(e),method='highs').status==0 for e in np.eye(B.shape[1])for s in(-1,1))
  return {'bounded_affine_section':bounded,'projective_vertices':len(vv),'spacelike':int((q>1e-8).sum()),'max_norm':float(q.max()),'best_point':(proj+B@vv[np.argmax(q)]).tolist()}
 except Exception as e:return {'geometry_error':str(e)[:250]}
def run():
 _,_,_,_,_,names,_,_,R,Gn,cn=S.load();out=[]
 for free in S.COMMON:
  f=R[names.index(free)];fixed=[names.index(n)for n in S.COMMON if n!=free];corners=[]
  for k in range(4):
   for base in BASES:corners.append([n for j,n in enumerate(T)if j!=k]+[base]+[n for n in S.COMMON if n!=free])
  caps=[];lookup=set()
  for ci,edge in enumerate(corners):
   A=R[[names.index(n)for n in edge]];M=A@Gn@A.T;assert np.linalg.eigvalsh(M)[0]>0
   for anchor in [n for n in T+BASES if n not in edge]:
    a=R[names.index(anchor)];v=a-(a@Gn@A.T)@np.linalg.inv(M)@A;q=v@Gn@v
    if q>=-1e-8:continue
    for m in LABELS:
     p=-2*np.cos(np.pi/m);rr=p/q*v-np.sqrt(1-p*p/(2*q))*f;ps=R@Gn@rr
     if rr@Gn@cn>=-1e-7 or not all(gate(t)for t in ps):continue
     key=tuple(np.round(rr,9))
     if key in lookup:continue
     lookup.add(key);cover=[]
     for j,ss in enumerate(corners):
      ids=[names.index(n)for n in ss];mat=np.block([[R[ids]@Gn@R[ids].T,(R[ids]@Gn@rr)[:,None]],[(rr@Gn@R[ids].T)[None,:],np.array([[2.]])]])
      if np.linalg.eigvalsh(mat)[0]>1e-7:cover.append(j)
     caps.append({'corner':ci,'edge':edge,'anchor':anchor,'label':m,'normal':rr.tolist(),'covers_old_vertical_edges':cover})
  rr=np.array([c['normal']for c in caps]);P=rr@Gn@rr.T;conflicts=[(i,j)for i,j in combinations(range(len(caps)),2)if not gate(P[i,j])]
  solver=z3.Solver();solver.set(timeout=30000);X=[z3.Bool(f'{free}_{i}')for i in range(len(caps))]
  for i,j in conflicts:solver.add(z3.Or(z3.Not(X[i]),z3.Not(X[j])))
  for j in range(8):solver.add(z3.Or([X[i]for i,c in enumerate(caps)if j in c['covers_old_vertical_edges']]))
  status=solver.check();row={'free':free,'old_H5_geometry':geom(R,Gn,cn,fixed),'corners':corners,'candidates':caps,'conflict_pairs':conflicts,'coverage_counts':[sum(j in c['covers_old_vertical_edges']for c in caps)for j in range(8)],'cover_solver_status':str(status)}
  if status==z3.sat:
   model=solver.model();ids=[i for i,x in enumerate(X)if z3.is_true(model.eval(x))];row['selected']=ids;row['selected_H5_geometry']=geom(np.vstack([R,rr[ids]]),Gn,cn,fixed)
  print(free,len(caps),row['coverage_counts'],status,row['old_H5_geometry'],flush=True);out.append(row)
 (S.V.ROOT/(S.PREFIX+'_linked_scout.json')).write_text(json.dumps({'status':'NUMERICAL_LINKED_PRISM_ENDPOINT_SEARCH','labels':LABELS,'cases':out,'scope':'Finite label scout and numerical pair/face geometry only. No exhaustive completion obstruction.'},indent=2)+'\n')
if __name__=='__main__':run()
