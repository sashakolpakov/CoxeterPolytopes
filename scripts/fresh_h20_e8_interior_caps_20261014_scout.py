#!/usr/bin/env python3
"""Bounded numerical source-edge/interior E8 cap event scout.
No exact acceptance or completeness claim for all intrinsic caps.
"""
import json,math,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from hashlib import sha256
from collections import Counter
import numpy as np
import mpmath as mp
import fresh_h20_contact_closure_20261014_replacement as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_interior_caps_20261014_scout'
BASE='results/fresh_h20_e8_cusp_completion_20261014.json';FACE='results/fresh_h20_contact_closure_20261014_h12_replacement.json';SOURCE='results/fresh_h12_a12_source_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def data():
 d,f,s=[json.loads((ROOT/p).read_text())for p in [BASE,FACE,SOURCE]]
 for row in [d,f,s]:
  for p,h in row['sources'].items():assert sha(p)==h
 mp.mp.dps=70;rad=list(map(mp.sqrt,[1,2,79,158]))
 def real(x):return float(sum(mp.mpf(a.numerator)/a.denominator*b for a,b in zip(x,rad)))
 Q=[[F(x)for x in row]for row in d['form']];N=[[tuple(map(F,x))for x in row]for row in d['normals']];H=[tuple(map(F,x))for x in d['strict_center']];escape=[tuple(map(F,x))for x in d['strict_escape']['point']]
 cache={}
 def dot(x,y):
  key=tuple(y)
  if key not in cache:cache[key]=[tuple(sum(Q[i][j]*y[j][k]for j in range(21))for k in range(4))for i in range(21)]
  out=K.ZZ
  for a,b in zip(x,cache[key]):out=K.add(out,K.mul(a,b))
  return out
 E=f['E8_support'];EI=np.array(f['E8_inverse'],float)
 vertices=f['H12_face_finite_vertices']+f['H12_face_ideal_vertices'];by_inc={tuple(v['source_incident']):v for v in vertices}
 V=[[tuple(map(F,x))for x in by_inc[tuple(v['incident'])]['point']]for v in s['projective_vertices']]
 VS=np.array([v['point']for v in s['projective_vertices']],dtype=object)
 VS=np.array([[float(F(x))for x in row]for row in VS]);Vgram=VS@np.diag([-1]+[1]*12)@VS.T
 vp=np.array([[real(dot(v,n))for n in N]for v in V]);ve=np.array([real(dot(v,escape))for v in V])
 ep=np.array([[real(dot(N[i],n))for n in N]for i in E]);ee=np.array([real(dot(N[i],escape))for i in E])
 G=np.array(s['gram'],float);edges=[]
 for i,j in combinations(range(37),2):
  support=sorted(set(s['projective_vertices'][i]['incident'])&set(s['projective_vertices'][j]['incident']))
  if len(support)==11 and np.linalg.eigvalsh(G[np.ix_(support,support)])[0]>1e-8:edges.append((i,j,support))
 assert len(edges)==234
 return d,f,s,np.array([[real(x)for x in n]for n in N]),np.array([[real(x)for x in v]for v in V]),np.array(Q,float),EI,ep,ee,vp,ve,Vgram,edges
def run():
 start=time.time();d,f,s,N,V,Q,EI,ep,ee,vp,ve,VG,edges=data();norms=np.array(d['norms'],float);E=f['E8_support'];profiles=[]
 for k in [1,2]:
  for support in combinations(range(8),k):
   c=np.zeros(8);c[list(support)]=1;wc=c@EI;mass=float(wc@c-2)
   profiles.append(dict(E8_positions=list(support),E8_walls=[E[i]for i in support],coefficients=wc,mass=mass,old_pairs=wc@ep,escape_pair=float(wc@ee)))
 labels=[2,3,4,6];contacts=np.array([0,.5,math.sqrt(.5),math.sqrt(.75)])[:,None]*np.sqrt(2*norms)[None,:]
 candidates=[];counts=Counter();excluded=[];classification_cosines=np.cos(np.pi/np.arange(2,25))
 for pi,p in enumerate(profiles):
  if p['mass']<1e-9:counts['null_profile_skipped']+=1;continue
  wp=p['old_pairs'];wE=p['escape_pair'];w=p['coefficients']@N[E]
  for ei,(iv,jv,support)in enumerate(edges):
   av,bv=vp[iv],vp[jv];q1,q2,q12=VG[iv,iv],VG[jv,jv],VG[iv,jv]
   for anchor in range(79):
    alpha,beta=av[anchor],bv[anchor]
    if abs(alpha)<1e-10 and abs(beta)<1e-10:continue
    swapped=abs(alpha)<abs(beta)
    if swapped:alpha,beta=beta,alpha;qa,qb=q2,q1
    else:qa,qb=q1,q2
    for label,contact in zip(labels,contacts[:,anchor]):
     D=wp[anchor]-contact
     u=D/alpha;v=-beta/alpha
     aa=qa*v*v+2*q12*v+qb;bb=2*qa*u*v+2*q12*u;cc=qa*u*u+p['mass']
     if abs(aa)<1e-10:
      roots=[]if abs(bb)<1e-10 else[-cc/bb]
     else:
      disc=bb*bb-4*aa*cc
      if disc<-1e-8:continue
      sq=math.sqrt(max(0,disc));roots=[(-bb+sq)/(2*aa),(-bb-sq)/(2*aa)]
     for v2 in roots:
      v1=u+v*v2;a,b=(v2,v1)if swapped else(v1,v2)
      counts['quadratic_branches']+=1
      if a<-1e-8 or b<-1e-8:continue
      counts['positive_cone_branches']+=1;cut=-wE+a*ve[iv]+b*ve[jv]
      if cut<=1e-7:continue
      counts['escape_cutters']+=1;pair=-wp+a*av+b*bv
      if np.max(pair)>1e-7:continue
      counts['acute_cutters']+=1
      cosine=-pair/np.sqrt(2*norms);finite=(cosine<1-1e-7)&(cosine>1e-7)
      nearest=np.min(np.abs(cosine[finite,None]-classification_cosines[None,:]),axis=1)
      bad=np.flatnonzero(finite)[nearest>1e-7]
      if len(bad):
       if len(excluded)<50:excluded.append(dict(profile=pi,edge=ei,anchor=anchor,label=label,bad_walls=bad.tolist(),bad_cosines=cosine[bad].tolist(),escape_pairing=cut))
       continue
      r=-w+a*V[iv]+b*V[jv]
      assert abs(r@Q@r-2)<1e-5
      candidates.append(dict(profile=pi,E8_positions=p['E8_positions'],source_edge=[iv,jv],source_edge_support=support,anchor=anchor,label=label,a=a,b=b,normal=r.tolist(),pairings=pair.tolist(),escape_pairing=cut))
 paths=[BASE,FACE,SOURCE,f'scripts/{TAG}.py','scripts/fresh_h20_contact_closure_20261014_replacement.py']
 return dict(status='NUMERICAL_BOUNDED_SOURCE_EDGE_EVENT_SCOUT',sources={p:sha(p)for p in paths},counts=dict(counts),source_edges=len(edges),E8_profiles=len(profiles),anchor_labels=labels,candidate_pair_label_test=list(range(2,25)),candidates=candidates,rejected_examples=excluded,elapsed_seconds=time.time()-start,scope='One/two pi3 E8 contacts, y in positive span of the endpoints of234 source edges, one retained79 anchor at2/3/4/6. Numerical proposals only; no all-profile or all-label nonexistence theorem. Null-residual single profile skipped as cusp-only.')
if __name__=='__main__':
 d=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['counts'],'candidates',len(d['candidates']))
