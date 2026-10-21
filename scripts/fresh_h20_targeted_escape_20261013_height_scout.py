#!/usr/bin/env python3
"""Target-directed single-height caps from the next86-wall polar projection."""
import itertools,json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fresh_h20_targeted_escape_20261013 import load,vector,ev

ROOT=Path(__file__).resolve().parents[1]
CAP='data/enumerations/fresh_h20_targeted_escape_20261013_exact.json'
VOLUME='data/enumerations/fresh_h20_targeted_escape_20261013_volume.json'
OUT=ROOT/'data/enumerations/fresh_h20_targeted_escape_20261013_height_scout.json'
LABELS=[2,3,4,5,6,8,10,12]

def setup():
 b,p,e,Q,R,h,_=load();cap=json.loads((ROOT/CAP).read_text());vol=json.loads((ROOT/VOLUME).read_text())
 r=np.array([ev(x[0])+np.sqrt(10)*ev(x[1])for x in cap['normal']]);R=np.vstack([R,r])
 w=vector(vol['certificate']['vector']);return b,p,Q,R,h,w,vol

def height_candidates(v,Q,R,h,w,labels=LABELS):
 cov=R@Q;norms=np.einsum('ij,jk,ik->i',R,Q,R)
 spatial=v.copy();spatial[-1]=0;q0=spatial@Q@spatial
 a=cov@spatial;z=cov[:,-1]
 c2=np.cos(np.pi/np.array(labels))**2
 qa=z[:,None]**2-3*norms[:,None]*c2
 qb=2*a[:,None]*z[:,None]*np.ones_like(c2)
 qc=a[:,None]**2-q0*norms[:,None]*c2
 disc=qb*qb-4*qa*qc
 proposals=[]
 for anchor in range(len(R)):
  for j,m in enumerate(labels):
   aa,bb,cc=qa[anchor,j],qb[anchor,j],qc[anchor,j]
   if abs(aa)<1e-9:
    ts=[]if abs(bb)<1e-9 else[-cc/bb]
   elif disc[anchor,j]<-1e-8:continue
   else:ts=[(-bb+s*np.sqrt(max(0,disc[anchor,j])))/(2*aa)for s in(-1,1)]
   for t in ts:
    norm=q0+3*t*t
    if norm<=1e-9:continue
    r=spatial.copy();r[-1]=t;r*=np.sqrt(2/norm)
    if r@Q@h>=-1e-7 or r@Q@w<=1e-7:continue
    if (a[anchor]+z[anchor]*t)>1e-7:continue
    proposals.append((anchor,m,t,r))
 if not proposals:return [],None
 rr=np.array([t[3]for t in proposals]);ratios=-(rr@cov.T)/np.sqrt(2*norms)
 finite=np.cos(np.pi/np.arange(2,121));ix=np.clip(np.searchsorted(finite,ratios),1,len(finite)-1)
 errors=np.minimum(abs(ratios-finite[ix]),abs(ratios-finite[ix-1]));errors=np.where(ratios>=1-1e-7,0,errors)
 maximum=np.max(errors,axis=1);out=[]
 for i in np.flatnonzero(maximum<1e-7):
  anchor,m,t,r=proposals[i];out.append({'anchor':anchor,'label':m,'height':float(t),'normal':r.tolist(),
   'escape_pairing':float(r@Q@w),'ratios':ratios[i].tolist()})
 k=int(np.argmin(maximum));return out,{'error':float(maximum[k]),'bad_indices':np.flatnonzero(errors[k]>1e-7).tolist(),
   'anchor':proposals[k][0],'label':proposals[k][1]}

def run():
 b,p,Q,R,h,w,vol=setup();cov=R@Q;norms=np.einsum('ij,jk,ik->i',R,Q,R)
 rng=np.random.default_rng(86107);records=[];seen=set();found=[]
 for i in range(33):
  if i==0:ids=vol['certificate']['numerical_active_indices']
  else:
   lp=linprog(rng.normal(size=21),A_ub=cov/np.sqrt(norms[:,None]),b_ub=np.zeros(86),
    A_eq=(h@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
   if not lp.success or lp.x@Q@lp.x<1e-10:continue
   ids=np.argsort(abs(cov@lp.x)/np.sqrt(norms))[:20].tolist()
  ids=tuple(sorted(ids))
  if ids in seen:continue
  seen.add(ids);_,s,vh=np.linalg.svd(cov[list(ids)],full_matrices=True)
  if s[-1]<1e-8:continue
  v=vh[-1];v*=(-1 if v@Q@h>0 else 1)
  if abs(v[-1])<1e-8:continue
  v/=abs(v[-1]);caps,best=height_candidates(v,Q,R,h,w)
  records.append({'active':list(ids),'polar':v.tolist(),'candidates':len(caps),'best_failed_profile':best})
  for cap in caps:cap['projection_index']=len(records)-1;found.append(cap)
  print({'tested':len(records),'candidates':len(found),'best':best},flush=True)
 OUT.write_text(json.dumps({'status':'NUMERICAL_TARGETED_HEIGHT_REPLACEMENT_SCOUT','labels':LABELS,
   'projections':records,'candidates':found,'scope':'Bounded33LP/known-gap polar projections; vary only lastcoordinate. Every accepted row cuts the one frozen86-wall strictescape and numerically passes all86pairgates. No exactcap or exhaustive claim.'},indent=2)+'\n')

if __name__=='__main__':run()
