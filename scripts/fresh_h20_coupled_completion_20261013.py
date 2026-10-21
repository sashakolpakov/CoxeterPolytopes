#!/usr/bin/env python3
"""Construct a full affine link at the new18/96 cusp (numerical height scout)."""
import argparse,json,sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fresh_h20_gap_next_20261013_exact import setup as setup96, read, sha, ROOT, kernel
from fresh_dimension_bridge_20261013_verify import E,Z,O,vec,parse,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv
from fresh_h20_edge_family_completion_20261013_cusp_repairs import base_root,sumvec
SOURCE='data/enumerations/fresh_h20_gap_next_20261013_exact.json'
DIAG='data/enumerations/fresh_h20_gap_next_20261013_cusp_diagnostic.json'
OUT='data/enumerations/fresh_h20_coupled_completion_20261013_scout.json'
F10=Extension(B,E(F(10)))
D8IDS=[[0,8,7,16,4,19,36,10],[2,14,13,26,11,20,5,12]]
COEFF=[2,1,2,2,2,2,1,1]

def num(f,x):
 if f is B:return float(x.a)+float(x.b)*np.sqrt(6)
 return num(f.base,x[0])+np.sqrt(num(f.base,f.d))*num(f.base,x[1])
def setup():
 Q,h,hc,roots=setup96();d=read(SOURCE)
 for c in d['new_caps']:
  if c['field']['radicand']==['1','0']:
   assert all(parse(x[1])==Z for x in c['normal']);roots.append((B,tuple(parse(x[0])for x in c['normal'])))
  else:
   f=Extension(B,parse(c['field']['radicand']));roots.append((f,tuple(map(f.parse,c['normal']))))
 assert len(roots)==98
 return Q,h,hc,roots

def frames(Q,h,hc,roots):
 d=read(SOURCE);r96=roots[96][1];delta=tuple(F10.add(x,F10.lift(y))for x,y in zip(r96,base_root(*roots[18]),strict=True))
 theta=[]
 for ids in D8IDS:
  bs=[tuple(x/(E(F(3))if i>=24 else O)for x in base_root(*roots[i]))for i in ids]
  z=sumvec([(E(F(c)),r)for c,r in zip(COEFF,bs,strict=True)])
  assert dot(z,mv(Q,z))==E(F(2))
  ps=[dot(z,mv(Q,r))for r in bs];assert ps.count(O)==1 and ps.count(Z)==7
  theta.append(z)
 diag=read(DIAG);p=tuple(E(F(107,327))*x for x in vec(diag['recession']['vector']))
 ids=[0,2,4,5,7,8,10,11,12,13,14,16,19,20,26,36]
 p1=vec(d['new_caps'][0]['polar_certificate']['vector']);a=base_root(*roots[3]);b=tuple(x/E(F(3))for x in base_root(*roots[42]))
 rows=[mv(Q,r)for r in [base_root(*roots[i])for i in ids]+[base_root(*roots[18]),p1,tuple(x+y for x,y in zip(a,b,strict=True)),p]]
 q,piv=kernel(rows);q=tuple(E(F(355,327))*x for x in q)
 assert dot(p,mv(Q,p))==E(F(5))and dot(q,mv(Q,q))==O and dot(p,mv(Q,q))==Z
 return delta,theta,p,q

def run():
 Qe,he,hc,roots=setup();delta,theta,p,q=frames(Qe,he,hc,roots)
 Q=np.array([[num(B,x)for x in row]for row in Qe]);h=np.array([num(B,x)for x in he]);R=np.array([[num(f,x)for x in r]for f,r in roots]);D=np.array([num(F10,x)for x in delta]);wn=np.einsum('ij,jk,ik->i',R,Q,R);dp=R@Q@D
 target=np.array([num(B,x)for x in vec(read(SOURCE)['strict_escape']['vector'])]);targetc=Q@target
 dirs=[tuple(-x for x in t)for t in theta]+[p,tuple(-x for x in p),q,tuple(-x for x in q)]
 rows=[];new=[]
 for v in dirs:
  z=np.array([num(B,x)for x in v]);n=z@Q@z;zp=R@Q@z;ix=np.flatnonzero(dp<-1e-7)
  hs=(zp[ix]+np.sqrt(n*wn[ix]))/-dp[ix];H=max(0,float(hs.max()),float((z@Q@h)/(-D@Q@h)))
  # Rational outward height with a small margin avoids accidental extra cusps.
  rat=F(int(np.ceil(H*10000))+1,10000);r=z+float(rat)*D;new.append(r)
  rows.append({'base_normal':[x.enc()for x in v],'base_norm':n,'minimum_all_dotted_height':H,'height':str(rat),'limiting_wall':int(ix[np.argmax(hs)]),'target_pairing':float(r@targetc),'normal':r.tolist()})
 U=np.vstack([R,new]);norms=np.einsum('ij,jk,ik->i',U,Q,U);P=U@Q@U.T;cs=-P/np.sqrt(norms[:,None]*norms[None,:]);bad=(cs<-1e-7)|((cs<1-1e-7)&(np.min(abs(cs[...,None]-np.array([0,.5,np.sqrt(.5),np.sqrt(.75)])),axis=-1)>1e-6));np.fill_diagonal(bad,False)
 print('height rows',[(r['height'],r['limiting_wall'],r['target_pairing'])for r in rows],flush=True);print('pair conflicts',np.argwhere(bad).tolist(),flush=True)
 rng=np.random.default_rng(202610131);A=U@Q/np.sqrt(norms[:,None]);hn=h/np.sqrt(-h@Q@h);trials=[];escape=None
 for k in range(300):
  ans=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
  if not ans.success:continue
  qv=float(ans.x@Q@ans.x);trials.append({'attempt':k,'norm':qv})
  if qv>1e-7:
   escape={'vector':ans.x.tolist(),'norm':qv,'active_indices':np.flatnonzero(abs(A@ans.x)<1e-7).tolist()};break
 out={'status':'NUMERICAL_104_WALL_FULL_CUSP_COMPLETION','source':SOURCE,'source_sha256':sha(SOURCE),'new_caps':rows,'pair_conflicts':np.argwhere(bad).tolist(),'escape':escape,'trials':trials}
 (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print('escape',escape,flush=True)
if __name__=='__main__':run()
