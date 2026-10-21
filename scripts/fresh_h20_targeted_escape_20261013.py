#!/usr/bin/env python3
"""Three simultaneous finite-angle repairs targeting the surviving H20 escape.

Numerical proposals only. Keep the17 active mirrors perpendicular to the sole
offending wall53 and release all three of its active neighbors22,34,44.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
POLAR='data/enumerations/fresh_dimension_bridge_polar_20261013.json'
EXT='data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json'
OUT=ROOT/'data/enumerations/fresh_h20_targeted_escape_20261013_scout.json'

def read(p):return json.loads((ROOT/p).read_text())
def ev(x):return float(Fraction(x[0]))+np.sqrt(6)*float(Fraction(x[1]))
def vector(x):return np.array([ev(t)for t in x])
def load():
 b=read(BASE);p=read(POLAR);e=read(EXT)
 Q=np.array([vector(r)for r in b['form']])
 R=[vector(r)for r in b['normals']]+[vector(r['vector'])for r in p['cuts']]
 for c in e['caps']:
  assert c['field']['kind']=='quadratic'and c['field']['base']=={'kind':'sqrt6'}
  d=np.sqrt(ev(c['field']['radicand']))
  R.append(np.array([ev(x[0])+d*ev(x[1])for x in c['vector']]))
 return b,p,e,Q,np.array(R),vector(b['strict_center']),vector(p['strict_escape']['vector'])

def run(labels):
 b,p,e,Q,R,h,w=load();cov=R@Q;ns=np.einsum('ij,jk,ik->i',R,Q,R)
 active=p['strict_escape']['source_active_indices'];bad=53
 released=[i for i in active if abs(cov[i]@R[bad])>1e-8]
 assert released==[22,34,44]
 kept=[i for i in active if i not in released];assert len(kept)==17
 local=[bad]+released
 profiles=list(itertools.product(labels,repeat=3));cos=np.cos(np.pi/np.array(profiles))
 finite=np.cos(np.pi/np.arange(2,121));counts=Counter();found=[];nearest=[]
 for anchors in itertools.combinations(local,3):
  ids=kept+list(anchors);A=cov[ids];u,s,vh=np.linalg.svd(A,full_matrices=True)
  assert s[-1]>1e-8
  direction=vh[-1];pinv=np.linalg.pinv(A)[:,-3:]
  initial=(-cos*np.sqrt(2*ns[list(anchors)]))@pinv.T
  aa=direction@Q@direction;bb=2*initial@Q@direction
  cc=np.einsum('ij,jk,ik->i',initial,Q,initial)-2;disc=bb*bb-4*aa*cc
  assert abs(aa)>1e-10
  counts['quadratics']+=len(profiles)
  for sign in(-1,1):
   rows=np.flatnonzero(disc>=-1e-9)
   ts=(-bb[rows]+sign*np.sqrt(np.maximum(0,disc[rows])))/(2*aa)
   rr=initial[rows]+ts[:,None]*direction
   take=(rr@Q@h<-1e-7)&(rr@Q@w>1e-7)
   rows=rows[take];rr=rr[take];counts['oriented_escape_cutters']+=len(rows)
   if not len(rows):continue
   ratios=-(rr@cov.T)/np.sqrt(2*ns)
   ix=np.clip(np.searchsorted(finite,ratios),1,len(finite)-1)
   errors=np.minimum(abs(ratios-finite[ix]),abs(ratios-finite[ix-1]))
   errors=np.where(ratios>=1-1e-7,0,errors)
   allok=np.max(errors,axis=1)<1e-7
   for j in np.argsort(np.max(errors,axis=1))[:5]:
    nearest.append({'anchors':list(anchors),'labels':list(profiles[rows[j]]),
      'maximum_pair_error':float(max(errors[j])),
      'bad_indices':np.flatnonzero(errors[j]>1e-7).tolist(),
      'escape_pairing':float(rr[j]@Q@w),'normal':rr[j].tolist(),
      'ratios':ratios[j].tolist()})
   for j in np.flatnonzero(allok):
    r=rr[j];assert abs(r@Q@r-2)<1e-6
    found.append({'anchors':list(anchors),'labels':list(profiles[rows[j]]),
       'normal':r.tolist(),'escape_pairing':float(r@Q@w),'ratios':ratios[j].tolist()})
  print({'anchors':anchors,'counts':dict(counts),'accepted':len(found)},flush=True)
 out={'status':'NUMERICAL_THREE_ANGLE_TARGETED_ESCAPE_SCOUT','source_hashes':{
  path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()for path in(BASE,POLAR,EXT)},
  'target':{'strict_escape_source':POLAR,'polar_index':27,'bad_wall':bad,'released':released,'kept':kept},
  'labels':labels,'counts':dict(counts),'candidates':found,
  'nearest':sorted(nearest,key=lambda r:r['maximum_pair_error'])[:12],
  'scope':'Keep17 specified orthogonalities; prescribe three angles at each three-of-four subset of walls53,22,34,44; solve both norm2 quadratic branches. Require cutting the saved strict escape and test all85 old pairs numerically against labels2..120 or dotted. Only these finite profiles are enumerated; no exact cap or universal obstruction claim.'}
 OUT.write_text(json.dumps(out,indent=2)+'\n')
 print({'status':out['status'],'candidates':len(found),'best_max_error':out['nearest'][0]['maximum_pair_error']if nearest else None})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--labels',type=int,nargs='+',default=[2,3,4,5,6,8,10,12]);args=p.parse_args();run(args.labels)
