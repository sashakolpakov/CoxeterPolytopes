#!/usr/bin/env python3
"""Bounded numerical proposal of PD18 coupled polar repairs; no acceptance theorem."""
import json
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_dimension_construction_routes_20261014.json'
def val(t):return float(Fraction(t['a']))+float(Fraction(t['b']))*np.sqrt(2)
def run():
 d=json.loads(SOURCE.read_text());Q=np.array([[val(x)for x in r]for r in d['form']]);N=np.array([[val(x)for x in r]for r in d['normals']]);h=np.array([float(Fraction(x['a'])+Fraction(x['b']))for x in d['strict_center']]);v=np.array([val(x)for x in d['strict_escape']['boundary_point']]);E=np.array([val(x)for x in d['strict_escape']['point']]);G=N@Q@N.T
 supports=[[0,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,27,33],[0,2,3,4,5,6,7,8,9,10,11,13,14,15,16,19,27,35],[0,2,3,4,5,6,7,8,10,11,12,13,14,15,17,19,27,43],[0,2,3,5,6,7,8,9,10,11,12,13,14,16,17,19,27,37]]
 caps=[];recipes=[];vertices=[]
 for si,S in enumerate(supports):
  B=G[np.ix_(S,S)];h0=h-N[S].T@np.linalg.solve(B,N[S]@Q@h)
  # Independent nullspace direction perpendicular to support, v and h0.
  _,_,V=np.linalg.svd(np.vstack([N[S]@Q,v@Q,h0@Q]));w=V[-1];w/=np.sqrt(w@Q@w)
  a=h0@Q@h0;b=2*v@Q@h0;c=v@Q@v;t0=(-b-np.sqrt(b*b-4*a*c))/(2*a)
  assert t0>0
  for bits in range(4,35):
   t=np.floor(t0*2**bits)/2**bits;m=v+t*h0;q=m@Q@m
   if q<=0:continue
   r1=m+np.sqrt(q)*w;r2=m-np.sqrt(q)*w;arr=np.array([r1,r2]);oldpairs=arr@Q@N.T
   others=[i for i in range(60)if i not in S]
   if np.max(oldpairs[:,others]/np.sqrt(2*q*np.diag(G)[others]))>=-1-1e-7:continue
   if np.min(arr@Q@E)<=1e-7:continue
   break
  else:raise AssertionError('no dyadic cap')
  p=h0-(h0@Q@m)/q*m
  assert p@Q@p<0
  vertices.append(p.tolist());recipes.append(dict(support=S,dyadic_t=str(Fraction(int(round(t*2**bits)),2**bits)),bits=bits,limiting_t=float(t0),q=float(q),target_min=float(np.min(arr@Q@E)),old_dotted_max=float(np.max(oldpairs[:,others]/np.sqrt(2*q*np.diag(G)[others])))))
  caps.extend([r1,r2])
 base_caps=np.array(caps);caps=[]
 for aa,bb in [(1,1),(1,-1),(-1,1),(-1,-1)]:
  cc=base_caps.copy();cc[:,19]*=aa;cc[:,20]*=bb;caps.extend(cc)
 A=np.vstack([N,caps]);n=np.einsum('ij,jk,ik->i',A,Q,A);pairs=A@Q@A.T/np.sqrt(n[:,None]*n[None,:]);newbad=[]
 for i in range(60,len(A)):
  for j in range(i):
   p=pairs[i,j]
   if abs(p)<1e-5 or p<=-1+1e-6:continue
   newbad.append([i,j,float(p)])
 # Projective slice by original strict timelike center; all vertices bounded.
 H=h/np.sqrt(-h@Q@h);_,_,V=np.linalg.svd((H@Q)[None,:]);T=V[1:].T
 # y=H+Tz gives (y,H)=-1. Euclidean transform makes spatial norm Euclidean.
 L=np.linalg.cholesky(T.T@Q@T);T=T@np.linalg.inv(L.T)
 rng=np.random.default_rng(20261014);esc=[]
 for trial in range(64):
  direction=rng.normal(size=20);direction/=np.linalg.norm(direction)
  lp=linprog(-direction,A_ub=A@Q@T,b_ub=-A@Q@H,bounds=[(None,None)]*20,method='highs')
  if not lp.success:continue
  for _ in range(3):
   z=lp.x;direction=z/np.linalg.norm(z)
   lp=linprog(-direction,A_ub=A@Q@T,b_ub=-A@Q@H,bounds=[(None,None)]*20,method='highs')
  p=H+T@lp.x;q=float(p@Q@p)
  if q>1e-7:esc.append(dict(point=p.tolist(),q=q,max_pair=float(np.max(A@Q@p))))
 out=dict(status='NUMERICAL_COUPLED_PD18_PROPOSALS_ONLY',recipes=recipes,base_caps=base_caps.tolist(),caps=np.array(caps).tolist(),rational_center=h.tolist(),finite_vertex_proposals=vertices,new_pair_failures=newbad,relaxed_full_bank_escapes=sorted(esc,key=lambda x:-x['q'])[:5],scope='Four specified support faces, dyadic heights; numerical scout only. Thirty-two caps retained jointly only if new_pair_failures empty. LP is proposal search, not volume certificate.')
 (ROOT/'results/fresh_h20_coupled_lift_20261014_scout.json').write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'caps',len(caps),'bad',newbad,'escapes',[(e['q'])for e in out['relaxed_full_bank_escapes']]);print(recipes)
if __name__=='__main__':run()
