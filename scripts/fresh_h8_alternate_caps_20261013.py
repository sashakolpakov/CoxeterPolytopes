#!/usr/bin/env python3
"""Bounded exact cap proposals for three previously unused released H8 lifts."""
import argparse
from fractions import Fraction as F
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scripts.fresh_h4_blocks_independent_verify import Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,pair,matrix,inverse,rref
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h8_lift_h9_20261013.json'
real=lambda x:float(x[0])+(1+5**.5)/2*float(x[1])
vec=lambda r:[enc(x)for x in r]
mat=lambda a:[vec(r)for r in a]

def seed(label,deleted):
 raw=next(s for s in json.loads(SOURCE.read_text())['seeds']if s['new_label']==label)
 Q=matrix(raw['host_gram']);allroots=matrix(raw['wall_normals']);ids=[i for i in range(11)if i!=deleted]
 R=tuple(allroots[i]for i in ids);G=tuple(tuple(pair(r,Q,s)for s in R)for r in R);I=inverse(G)
 h=tuple(map(parse,raw['strict_center']));y=tuple(neg(pair(h,Q,r))for r in R)
 assert all(sign(t)>0 for t in y) and sign(pair(y,I,y))<0
 return Q,R,G,I,y,ids

def strict_escape(I,A,center,tries=200):
 # A contains coefficients c=G^-1 p; an outward cap inequality is c.y>=0.
 assert all(sign(x)>0 for x in center) and sign(pair(center,I,center))<0
 assert all(sign(dot(a,center))>0 for a in A)
 qn=np.array([[real(x)for x in r]for r in I]);an=np.array([[real(x)for x in r]for r in A]);cn=np.array([real(x)for x in center]);cn/=sum(cn)
 if len(A):an/=np.linalg.norm(an,axis=1)[:,None]
 rng=np.random.default_rng(20261013)
 for trial in range(tries):
  res=linprog(rng.normal(size=10),A_ub=-an if len(A)else None,b_ub=np.zeros(len(A))if len(A)else None,A_eq=np.ones((1,10)),b_eq=[1],bounds=(0,None),method='highs')
  assert res.success,res.message
  if res.x@qn@res.x<=1e-8:continue
  eps=.01
  while (res.x+eps*cn)@qn@(res.x+eps*cn)<=0:eps/=2
  target=res.x+eps*cn
  for den in (10**8,10**10,10**12):
   y=tuple(scale(O,F(round(float(x)*den),den))for x in target);q=pair(y,I,y)
   if any(sign(t)<=0 for t in y)or sign(q)<=0:continue
   values=[dot(a,y)for a in A]
   if any(sign(x)<=0 for x in values):continue
   cp=pair(y,I,center)
   if sign(cp)>=0:continue
   return {'status':'EXACT_STRICT_SPACELIKE_ESCAPE','y':vec(y),'norm':enc(q),'center_pairing':enc(cp),'trial':trial,'all_cap_values_strict':True}
 return {'status':'UNDECIDED_NO_ESCAPE','trials':tries}

def run(label,deleted,max_support,output):
 Q,R,G,I,center,ids=seed(label,deleted);In=np.array([[real(x)for x in r]for r in I])
 magnitudes=(O,P,scale(O,2),add(O,P),scale(P,2));vals=np.array([real(x)for x in magnitudes]);caps=[];counts=[];seen=set()
 for k in range(1,max_support+1):
  choices=np.array(list(product(range(len(vals)),repeat=k)),dtype=int);pp=vals[choices];prop=0;orientation=0;invalid=0
  for support in combinations(range(10),k):
   squares=np.einsum('bi,ij,bj->b',pp,In[np.ix_(support,support)],pp)
   for ci in np.flatnonzero((abs(squares-1)<1e-7)|(abs(squares-2)<1e-7)):
    prop+=1;p=[Z]*10
    for j,m in zip(support,choices[ci]):p[j]=neg(magnitudes[m])
    p=tuple(p);c=mv(I,p);q=dot(p,c)
    if q not in (O,scale(O,2)):continue
    if any(sign(sub(div(scale(mul(x,x),4),mul(q,G[j][j])),scale(O,4)))<0 and div(scale(mul(x,x),4),mul(q,G[j][j]))not in (Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P))for j,x in enumerate(p)):
     invalid+=1;continue
    if sign(dot(c,center))<=0:orientation+=1;continue
    if c in seen:continue
    seen.add(c);r=tuple(total(mul(a,w[j])for a,w in zip(c,R))for j in range(10))
    assert pair(r,Q,r)==q and tuple(pair(r,Q,w)for w in R)==p
    caps.append({'normal':vec(r),'norm':enc(q),'seed_pairings':vec(p),'seed_coefficients':vec(c),'support':k})
  row={'support':k,'numeric_norm_proposals':prop,'incompatible_angles':invalid,'excluded_by_fixed_center':orientation,'caps_so_far':len(caps)};counts.append(row);print(label,deleted,row,flush=True)
 A=[tuple(map(parse,c['seed_coefficients']))for c in caps]
 out={'schema':'alternate-released-H8-cap-bank-v1','source':str(SOURCE.relative_to(ROOT)),'seed_new_label':label,'deleted_original_wall':deleted,'retained_original_walls':ids,'host_gram':mat(Q),'seed_normals':mat(R),'seed_gram':mat(G),'seed_gram_inverse':mat(I),'strict_center_y':vec(center),'seed_rank':len(rref(R)[1]),'dual_vertex_norms':vec([I[j][j]for j in range(10)]),'maximum_support':max_support,'pairing_magnitudes':vec(magnitudes),'counts':counts,'caps':caps,'volume_gate':strict_escape(I,A,center),'scope':'Finite numeric proposal dictionary; every archived cap and escape exact. Fixed center restricts accepted caps; no exhaustion or general completion claim.'}
 output.write_text(json.dumps(out,indent=2)+'\n');print({'seed':[label,deleted],'caps':len(caps),'gate':out['volume_gate']['status']},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--label',type=int,required=True);p.add_argument('--delete',type=int,required=True);p.add_argument('--support',type=int,default=5);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.label,a.delete,a.support,a.output)
