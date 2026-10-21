#!/usr/bin/env python3
"""Two dotted distances coupled to an existing cap angle on the label-4 seed.

Floating arithmetic proposes a bounded scored list only. Exact Q(phi)(sqrt d)
identities certify all accepted roots and strict spacelike common escapes.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations,product
import json
from pathlib import Path
import heapq
import numpy as np
from scipy.optimize import linprog
from scripts.fresh_h8_alternate_caps_20261013 import *
from scripts.fresh_h8_release_h9_free_caps_20261013 import EZ,EO,eadd,esub,escale,emul,esign,ereal,eencode
BASE=ROOT/'data/enumerations/fresh_h8_alternate_4_4_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h8_alternate_coupled_20261013.json'
# These anchors pass a numerical necessary no-affine-subdiagram screen.
# This selection is a proposal heuristic, not an exact affine-exclusion claim.
ANCHORS=(7,10,13,16,19,26,27,28,35,39,40,41,44,46,59,60,61,62,67,70,71,77,82,92,97,98,100,101,103,104,105,106)

def proposals(I,old,center,y,max_support=3,per_anchor=3):
 In=np.array([[real(x)for x in r]for r in I]);An=np.array([[real(x)for x in old[i]]for i in ANCHORS]);yn=np.array([real(x)for x in y]);cn=np.array([real(x)for x in center]);Iy=In@yn;Ic=In@cn
 Ms=np.tile(np.array([0.,1.,real(P)]),len(ANCHORS));refs=np.repeat(An,3,axis=0);anchorids=np.repeat(np.array(ANCHORS),3);targets=np.tile(np.arange(3),len(ANCHORS));heaps={i:[]for i in ANCHORS};serial=0;admissible=0
 for u,v in combinations(range(10),2):
  viable=abs(refs[:,v])>1e-8
  if not viable.any():continue
  ar=refs[viable];ai=anchorids[viable];ti=targets[viable];mm=Ms[viable]
  t=-ar[:,u]/ar[:,v];w1=np.zeros((len(t),10));w1[:,u]=1;w1[:,v]=t;A=np.einsum('bi,ij,bj->b',w1,In,w1)
  rest=[j for j in range(10)if j not in (u,v)]
  for k in range(max_support+1):
   for supp in combinations(rest,k):
    for labels in product((0,1),repeat=k):
     f=np.zeros(10)
     for j,l in zip(supp,labels):f[j]=(1.,real(P))[l]
     L=(mm-ar@f)/ar[:,v];w0=np.broadcast_to(f,(len(t),10)).copy();w0[:,v]=L
     B=2*np.einsum('bi,ij,bj->b',w1,In,w0);C=np.einsum('bi,ij,bj->b',w0,In,w0)-2;D=B*B-4*A*C
     good=(D>1e-10)&(abs(A)>1e-10)
     if not good.any():continue
     for branch in (-1,1):
      x=np.full(len(t),np.nan);x[good]=(-B[good]+branch*np.sqrt(D[good]))/(2*A[good]);z=L+t*x
      ws=w0+x[:,None]*w1;cut=ws@Iy;orientation=ws@Ic
      ok=good&(x>2+1e-7)&(z>2+1e-7)&(cut>1e-7)&(orientation< -1e-7)
      admissible+=int(ok.sum())
      for j in np.flatnonzero(ok):
       # Normalize the violation to avoid rewarding mere huge coefficients.
       score=float(cut[j]/np.linalg.norm(ws[j]@In));anchor=int(ai[j]);serial+=1
       item=(score,serial,(u,v,supp,labels,anchor,int(ti[j]),branch))
       h=heaps[anchor]
       if len(h)<per_anchor:heapq.heappush(h,item)
       elif item[0]>h[0][0]:heapq.heapreplace(h,item)
 selected=sorted([z for h in heaps.values()for z in h],reverse=True);out=[];seen=set()
 for score,_,(u,v,supp,labels,anchor,target,branch)in selected:
  f=[Z]*10
  for j,l in zip(supp,labels):f[j]=(O,P)[l]
  ref=old[anchor];M=(Z,O,P)[target];L=div(sub(M,dot(ref,f)),ref[v]);T=div(neg(ref[u]),ref[v]);w0=list(f);w0[v]=L;w0=tuple(w0);w1=tuple(O if j==u else T if j==v else Z for j in range(10))
  A=pair(w1,I,w1);B=scale(pair(w0,I,w1),2);C=sub(pair(w0,I,w0),scale(O,2));D=sub(mul(B,B),scale(mul(A,C),4));assert sign(D)>0
  x=(div(neg(B),scale(A,2)),div(scale(O,branch),scale(A,2)))
  wa=tuple(add(a,mul(b,x[0]))for a,b in zip(w0,w1));wb=tuple(mul(b,x[1])for b in w1)
  pa=tuple(neg(a)for a in wa);pb=tuple(neg(b)for b in wb);ca,cb=mv(I,pa),mv(I,pb)
  key=(D,ca,cb)
  if key in seen:continue
  seen.add(key)
  assert all(esign((a,b),D)<=0 for a,b in zip(pa,pb))
  assert esign(esub((wa[u],wb[u]),(scale(O,2),Z)),D)>0 and esign(esub((wa[v],wb[v]),(scale(O,2),Z)),D)>0
  norm=(add(dot(pa,ca),mul(D,dot(pb,cb))),scale(dot(pa,cb),2));assert esign(esub(norm,(scale(O,2),Z)),D)==0
  assert esign(eadd((dot(pa,ref),dot(pb,ref)),(M,Z)),D)==0
  assert esign((dot(ca,center),dot(cb,center)),D)>0
  assert esign((dot(ca,y),dot(cb,y)),D)<0
  # Every finite seed angle is checked exactly; free pairs are strictly dotted.
  for j,(a,b)in enumerate(zip(pa,pb)):
   sq=emul((a,b),(a,b),D);ratio=emul(sq,(div(scale(O,2),I_G[j][j]),Z),D)
   assert esign(esub(ratio,(scale(O,4),Z)),D)>0 or any(esign(esub(ratio,(q,Z)),D)==0 for q in (Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P)))
  out.append({'free_seed_indices':[u,v],'fixed_pairing_magnitudes':vec(f),'anchor_base_index':anchor,'anchor_pairing_magnitude':enc(M),'norm_polynomial_descending':vec([A,B,C]),'discriminant':enc(D),'free_x':eencode(x),'pairings_rational':vec(pa),'pairings_radical':vec(pb),'coefficients_rational':vec(ca),'coefficients_radical':vec(cb),'target_escape_pairing':eencode((dot(ca,y),dot(cb,y))),'numeric_score':score})
 return out,admissible

def escape(I,old,caps,center,tries=400):
 records=[(parse(r['discriminant']),tuple(map(parse,r['coefficients_rational'])),tuple(map(parse,r['coefficients_radical'])))for r in caps]
 A=np.array([[real(x)for x in r]for r in old]+[[real(x)+real(d)**.5*real(z)for x,z in zip(a,b)]for d,a,b in records]);A/=np.linalg.norm(A,axis=1)[:,None]
 qn=np.array([[real(x)for x in r]for r in I]);cn=np.array([real(x)for x in center]);cn/=sum(cn);rng=np.random.default_rng(20261013)
 for trial in range(tries):
  p=linprog(rng.normal(size=10),A_ub=-A,b_ub=np.zeros(len(A)),A_eq=np.ones((1,10)),b_eq=[1],bounds=(0,None),method='highs');assert p.success,p.message
  if p.x@qn@p.x<=1e-8:continue
  ep=.01
  while(p.x+ep*cn)@qn@(p.x+ep*cn)<=0:ep/=2
  target=p.x+ep*cn
  for den in (10**8,10**10,10**12,10**14):
   y=tuple(scale(O,F(round(float(t)*den),den))for t in target);q=pair(y,I,y)
   if sign(q)<=0 or any(sign(x)<=0 for x in y):continue
   if any(sign(dot(a,y))<=0 for a in old):continue
   if any(esign((dot(a,y),dot(b,y)),d)<=0 for d,a,b in records):continue
   cp=pair(y,I,center)
   if sign(cp)>0:
    # Add a positive multiple of the strict timelike center, reversing
    # the center pairing while preserving the positive Lorentz norm.
    lam=div(scale(cp,-2),pair(center,I,center))
    y=tuple(add(x,mul(lam,c))for x,c in zip(y,center));q=pair(y,I,y);cp=pair(y,I,center)
   elif sign(cp)==0:
    lam=F(1)
    while sign(pair(tuple(add(x,scale(c,lam))for x,c in zip(y,center)),I,tuple(add(x,scale(c,lam))for x,c in zip(y,center))))<=0:lam/=2
    y=tuple(add(x,scale(c,lam))for x,c in zip(y,center));q=pair(y,I,y);cp=pair(y,I,center)
   assert sign(q)>0 and sign(cp)<0
   return {'status':'EXACT_STRICT_SPACELIKE_ESCAPE','y':vec(y),'norm':enc(q),'center_pairing':enc(cp),'trial':trial,'all_caps_strict':True}
 return {'status':'UNDECIDED_NO_NUMERIC_ESCAPE','trials':tries}

def run(rounds=3):
 global I_G
 data=json.loads(BASE.read_text());I=matrix(data['seed_gram_inverse']);I_G=matrix(data['seed_gram']);old=[tuple(map(parse,r['seed_coefficients']))for r in data['caps']];center=tuple(map(parse,data['strict_center_y']));y=tuple(map(parse,data['volume_gate']['y']));caps=[];history=[]
 for iteration in range(rounds):
  new,count=proposals(I,old,center,y);caps.extend(new);gate=escape(I,old,caps,center);row={'round':iteration,'numerically_admissible_coupled_solutions':count,'exact_new_caps':len(new),'total_new_caps':len(caps),'target_escape_y':vec(y),'volume_gate':gate};history.append(row)
  out={'schema':'alternate-H8-two-dotted-coupled-caps-v1','base_source':str(BASE.relative_to(ROOT)),'fixed_pairing_maximum_support':3,'kept_proposals_per_anchor':3,'anchor_base_indices':list(ANCHORS),'caps':caps,'rounds':history,'scope':'Bounded, numerically proposed cap construction on unused label-4/omit-4 seed. All accepted normals and common escapes exactly certified; no exhaustive cap or completion claim.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in row.items()if k!='target_escape_y'},flush=True)
  if gate['status']!='EXACT_STRICT_SPACELIKE_ESCAPE':break
  y=tuple(map(parse,gate['y']))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--rounds',type=int,default=3);a=p.parse_args();run(a.rounds)
