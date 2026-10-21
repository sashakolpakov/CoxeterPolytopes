#!/usr/bin/env python3
"""Bounded numerical mixed-E7/source-edge event scout; no exact acceptance."""
import argparse,json,math,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from hashlib import sha256
import numpy as np
import mpmath as mp
import fresh_h20_e7_intrinsic_caps_20261014_cusps_review as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_causal_residuals_20261014'
BASE='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
FACE='results/fresh_h20_e7_intrinsic_caps_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def data():
 d,f=map(K.read,[BASE,FACE]);Q=[list(map(F,row))for row in d['form']];N=list(map(K.vector,d['normals']));H=K.vector(d['strict_center']);E=K.vector(d['strict_escape']['point'])
 cache={};sp=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
 def dot(v,w):
  key=tuple(w)
  if key not in cache:cache[key]=[tuple(sum(q*w[j][k]for j,q in row)for k in range(8))for row in sp]
  return K.total(K.mul(a,b)for a,b in zip(v,cache[key]))
 assert len(N)==87 and all(K.sign(dot(n,E))<0 for n in N)and K.sign(dot(E,E))>0
 assert all(K.sign(dot(n,H))<0 for n in N)and K.sign(dot(H,H))<0
 S=d['E7_support'];G=[[K.rational(dot(N[i],N[j]))for j in S]for i in S];K.pd(G);I=K.inv(G)
 rows=f['face_finite_vertices']+f['face_ideal_vertices'];assert len(rows)==46
 V=[K.vector(row['ambient_point'])for row in rows]
 C=[list(map(F,row))for row in f['face_form']];RV=[list(map(F,row['point']))for row in rows]
 VG=[[sum(v[i]*C[i][j]*w[j]for i in range(14)for j in range(14))for w in RV]for v in RV]
 assert all(VG[i][i]<0 for i in range(44))and VG[44][44]==VG[45][45]==0
 assert all(VG[i][j]<0 for i in range(46)for j in range(i))
 FG=[list(map(F,row))for row in f['face_gram']];edges=[]
 for i,j in combinations(range(46),2):
  ids=sorted(set(rows[i]['incident'])&set(rows[j]['incident']))
  if len(ids)!=12:continue
  try:K.pd([[FG[a][b]for b in ids]for a in ids])
  except AssertionError:continue
  edges.append([i,j,ids])
 mp.mp.dps=70;rad=[mp.sqrt(q)for q in K.RAD]
 def val(x):return float(sum(mp.mpf(a.numerator)/a.denominator*b for a,b in zip(x,rad)))
 VP=np.array([[val(dot(v,n))for n in N]for v in V]);VE=np.array([val(dot(v,E))for v in V])
 SP=np.array([[val(dot(N[i],n))for n in N]for i in S]);SE=np.array([val(dot(N[i],E))for i in S])
 return d,f,S,np.array(I,float),np.array([[val(x)for x in n]for n in N]),np.array([[val(x)for x in v]for v in V]),np.array(Q,float),np.array(VG,float),VP,VE,SP,SE,edges
def run(sizes,single_labels=None):
 start=time.time();d,f,S,J,N,V,Q,VG,VP,VE,SP,SE,edges=data();norms=np.array(d['norms'],float)
 profiles=[]
 for size in sizes:
  for support in combinations(range(7),size):
   c=np.zeros(7);c[list(support)]=1;wc=c@J
   profiles.append(dict(support=list(support),contact_labels=[3]*size,weights=wc,mass=float(wc@c-2),pairs=wc@SP,escape=float(wc@SE)))
 if single_labels:
  for label in single_labels:
   c=np.zeros(7);c[6]=2*np.cos(np.pi/label);wc=c@J
   profiles.append(dict(support=[6],contact_labels=[label],weights=wc,mass=float(wc@c-2),pairs=wc@SP,escape=float(wc@SE)))
 counts=Counter();candidates=[];near=[];labels=[2,3,4,6]
 contacts=np.sqrt(2*norms)[None,:]*np.cos(np.pi/np.array(labels))[:,None];contacts[0]=0
 cosines=np.cos(np.pi/np.arange(2,65));seen=set();best=-math.inf
 for pi,p in enumerate(profiles):
  wp=p['pairs'];we=p['escape'];w=p['weights']@N[S]
  for ei,(iv,jv,support)in enumerate(edges):
   aa,bb=VP[iv].copy(),VP[jv].copy();swap=np.abs(aa)<np.abs(bb)
   alpha=np.where(swap,bb,aa);beta=np.where(swap,aa,bb)
   valid=np.abs(alpha)>1e-12
   qa=np.where(swap,VG[jv,jv],VG[iv,iv]);qb=np.where(swap,VG[iv,iv],VG[jv,jv]);qab=VG[iv,jv]
   for li,label in enumerate(labels):
    u=np.zeros(87);v=np.zeros(87);u[valid]=(wp[valid]-contacts[li,valid])/alpha[valid];v[valid]=-beta[valid]/alpha[valid]
    A=qa*v*v+2*qab*v+qb;B=2*qa*u*v+2*qab*u;C=qa*u*u+p['mass'];disc=B*B-4*A*C
    for branch in [-1,1]:
     mask=valid&(np.abs(A)>1e-10)&(disc>=-1e-8)
     x=np.zeros(87);x[mask]=(-B[mask]+branch*np.sqrt(np.maximum(disc[mask],0)))/(2*A[mask])
     if branch==1:
      linear=valid&(np.abs(A)<=1e-10)&(np.abs(B)>1e-10);x[linear]=-C[linear]/B[linear];mask|=linear
     y=u+v*x;a=np.where(swap,x,y);b=np.where(swap,y,x)
     counts['quadratic_branches']+=int(sum(mask));mask&=(a>=-1e-8)&(b>=-1e-8)
     counts['positive_cone_branches']+=int(sum(mask));cut=-we+a*VE[iv]+b*VE[jv]
     if mask.any():best=max(best,float(np.max(cut[mask])))
     mask&=cut>1e-7;counts['escape_cutters']+=int(sum(mask))
     for anchor in np.flatnonzero(mask):
      pair=-wp+a[anchor]*VP[iv]+b[anchor]*VP[jv]
      if np.max(pair)>1e-7:continue
      counts['acute_cutters']+=1;cos=-pair/np.sqrt(2*norms)
      finite=(cos>1e-7)&(cos<1-1e-7);errors=np.min(np.abs(cos[finite,None]-cosines[None,:]),axis=1)
      bad=np.flatnonzero(finite)[errors>1e-7]
      rec=dict(profile=pi,E7_positions=p['support'],E7_labels=p['contact_labels'],source_edge=[iv,jv],source_edge_support=support,anchor=int(anchor),label=label,branch=branch,a=float(a[anchor]),b=float(b[anchor]),escape_pairing=float(cut[anchor]))
      if len(bad):
       if len(near)<30:near.append(dict(rec,bad_walls=bad.tolist(),bad_cosines=cos[bad].tolist()))
       continue
      r=-w+a[anchor]*V[iv]+b[anchor]*V[jv];key=tuple(np.round(r,8))
      if key in seen:continue
      seen.add(key);assert abs(r@Q@r-2)<1e-4
      candidates.append(dict(rec,normal=r.tolist(),pairings=pair.tolist()))
  print('profile',pi,p['support'],'candidates',len(candidates),'counts',dict(counts),flush=True)
 paths=[BASE,FACE,f'scripts/{TAG}_scout.py','scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps_review.py']
 out=dict(status='BOUNDED_NUMERICAL_MIXED_E7_SOURCE_EDGE_SCOUT',sources={x:sha(x)for x in paths},source_edges=edges,profile_sizes=sizes,single_wall65_labels=single_labels,E7_support=S,profile_count=len(profiles),anchor_labels=labels,pair_test_labels=list(range(2,65)),counts=dict(counts),best_proposed_escape_pairing=best,candidates=candidates,rejected_examples=near,elapsed_seconds=time.time()-start,scope='Numerical candidates only. Listed E7 profiles; residual in positive spans of source-edge endpoints; anchors2/3/4/6. Not a continuous all-profile or all-label exclusion.')
 suffix='single_scout'if single_labels else'scout'
 (ROOT/f'results/{TAG}_{suffix}.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--sizes',default='2,3');ap.add_argument('--single-labels',default='');a=ap.parse_args();d=run([int(x)for x in a.sizes.split(',')if x],[int(x)for x in a.single_labels.split(',')if x]or None);print(d['status'],len(d['candidates']))
