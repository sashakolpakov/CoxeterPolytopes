#!/usr/bin/env python3
"""18-short-wall residual3-space with two finite/parallel angle equations."""
import itertools,json
from pathlib import Path
from collections import Counter
import numpy as np
from fresh_h20_targeted_escape_20261013_second_volume import extended_setup
ROOT=Path(__file__).resolve().parents[1]
TARGET='data/enumerations/fresh_h20_targeted_escape_20261013_bank_geometry.json'
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_edge_directions.json'
LABELS=[0,2,3,4,5,6,8,10,12] # 0 means parallel, not a finite angle.

def run():
 Q,R,h,*_=extended_setup();cov=R@Q;ns=np.einsum('ij,jk,ik->i',R,Q,R)
 g=json.loads((ROOT/TARGET).read_text());target=g['whole_bank_geometry'];w=np.array(target['point'])
 active=[i for i in target['active_indices']if i<24];assert len(active)==17
 profiles=list(itertools.product(LABELS,repeat=2));cos=np.array([[1 if m==0 else np.cos(np.pi/m)for m in row]for row in profiles])
 finite=np.cos(np.pi/np.arange(2,121));counts=Counter();caps=[];seen=set();supports=[]
 for extra in range(24):
  if extra in active:continue
  S=sorted(active+[extra]);M=R[S]@Q@R[S].T
  eig=np.linalg.eigvalsh(M)
  if eig[0]<1e-7:continue
  _,_,vh=np.linalg.svd(cov[S],full_matrices=True);N=vh[18:].T;H=N.T@Q@N
  T=cov@N;targetcov=w@Q@N;hc=h@Q@N;other=[i for i in range(87)if i not in S]
  supports.append(S);start=len(caps)
  for a,b in itertools.combinations(other,2):
   A=T[[a,b]];d=np.cross(A[0],A[1]);length=np.linalg.norm(d)
   if length<1e-9:continue
   d/=length
   aa=d@H@d;pinv=A.T@np.linalg.inv(A@A.T)
   rhs=-cos*np.sqrt(2*ns[[a,b]]);initial=rhs@pinv.T
   bb=2*initial@H@d;cc=np.einsum('ij,jk,ik->i',initial,H,initial)-2
   counts['quadratics']+=len(profiles)
   disc=bb*bb-4*aa*cc
   branches=(0,)if abs(aa)<1e-10 else(-1,1)
   for branch in branches:
    rows=np.flatnonzero(abs(bb)>1e-10)if branch==0 else np.flatnonzero(disc>=-1e-8)
    t=-cc[rows]/bb[rows]if branch==0 else(-bb[rows]+branch*np.sqrt(np.maximum(0,disc[rows])))/(2*aa)
    u=initial[rows]+t[:,None]*d
    good=(u@hc<-1e-7)&(u@targetcov>1e-7);rows=rows[good];u=u[good]
    if not len(rows):continue
    counts['oriented_cutters']+=len(rows)
    ratios=-(u@T.T)/np.sqrt(2*ns)
    ix=np.clip(np.searchsorted(finite,ratios),1,len(finite)-1)
    error=np.minimum(abs(ratios-finite[ix]),abs(ratios-finite[ix-1]));error=np.where(ratios>=1-1e-7,0,error)
    for j in np.flatnonzero(np.max(error,axis=1)<1e-7):
     r=N@u[j];assert abs(r@Q@r-2)<1e-5
     key=tuple(np.round(r/np.linalg.norm(r),9))
     if key in seen:continue
     seen.add(key);caps.append({'short_orthogonal_support':S,'anchors':[a,b],'labels':list(profiles[rows[j]]),
      'normal':r.tolist(),'target_pairing':float(r@Q@w),'ratios':ratios[j].tolist()})
  print({'extra':extra,'new_caps':len(caps)-start,'caps':len(caps),'counts':dict(counts)},flush=True)
 out={'status':'NUMERICAL_18_SHORT_TWO_ANCHOR_TARGETED_SCOUT','target_source':TARGET,'labels':LABELS,
  'parallel_label':0,'supports':supports,'counts':dict(counts),'caps':caps,
  'scope':'Extend17activeoriginalshortwalls ofthefullvertex-bankescape to18PDwalls. Prescribe twofiniteorparallelangles againstallremaining87walls, solvebothnormquadraticbranchesincludinglinearcases, andtestall87relations. Allcapsmustcuttheonecommonnumericaltarget. Numericproposalsonly.'}
 (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print({'done':True,'caps':len(caps)})

if __name__=='__main__':run()
