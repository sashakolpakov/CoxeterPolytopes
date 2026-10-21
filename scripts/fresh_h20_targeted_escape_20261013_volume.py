#!/usr/bin/env python3
"""Bounded LP proposals and exact strict-escape gate for the new86-wall cone."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fresh_h20_targeted_escape_20261013 import load,vector,ev,BASE,POLAR,EXT
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot

ROOT=Path(__file__).resolve().parents[1]
CAP='data/enumerations/fresh_h20_targeted_escape_20261013_exact.json'
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_volume.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def setup():
 b,p,e,Qn,Rn,hn,wn=load();new=read(CAP)
 r=np.array([ev(x[0])+np.sqrt(10)*ev(x[1])for x in new['normal']]);Rn=np.vstack([Rn,r])
 Q=tuple(map(vec,b['form']));h=vec(b['strict_center']);hc=mv(Q,h)
 roots=[(B,r)for r in list(map(vec,b['normals']))+[vec(x['vector'])for x in p['cuts']]]
 for c in e['caps']+[{'field':new['field'],'vector':new['normal']}]:
  f=Extension(B,parse(c['field']['radicand']));roots.append((f,tuple(map(f.parse,c['vector']))))
 return Qn,Rn,hn,Q,h,hc,roots

def check(v,Q,hc,roots):
 cv=mv(Q,v);q=dot(v,cv);time=dot(v,hc)
 assert q.sign()>0 and time.sign()<0
 pairs=[]
 for f,r in roots:
  p=fdot(f,r,cv if f is B else tuple(map(f.lift,cv)))
  assert f.sign(p)<0;pairs.append(f.encode(p))
 return {'vector':[x.enc()for x in v],'norm':q.enc(),'center_pairing':time.enc(),'strict_wall_pairings':pairs}

def run(attempts,seed):
 Qn,Rn,hn,Q,h,hc,roots=setup();n=np.einsum('ij,jk,ik->i',Rn,Qn,Rn)
 cov=Rn@Qn/np.sqrt(n[:,None]);slice_h=hn/np.sqrt(-hn@Qn@hn)
 rng=np.random.default_rng(seed);tries=[];cert=None
 for i in range(attempts):
  result=linprog(rng.normal(size=21),A_ub=cov,b_ub=np.zeros(86),
   A_eq=(slice_h@Qn)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
  if not result.success:tries.append({'attempt':i,'lp_status':int(result.status)});continue
  q=result.x@Qn@result.x
  ids=np.flatnonzero(abs(cov@result.x)<1e-7).tolist()
  tries.append({'attempt':i,'numeric_norm':float(q),'active_indices':ids})
  if q<1e-7:continue
  scale=max(abs(result.x));base=tuple(E(F(round(x/scale*10**9),10**9))for x in result.x)
  for power in range(12,61):
   v=tuple(x+E(F(1,2**power))*z for x,z in zip(base,h,strict=True))
   try:cert=check(v,Q,hc,roots)
   except AssertionError:continue
   cert.update({'source_attempt':i,'perturbation_power':power,'numerical_active_indices':ids});break
  if cert:break
 return {'status':'EXACT_STRICT_ESCAPE_FOR_86_WALL_H20'if cert else'BOUNDED_NUMERICAL_86_WALL_VOLUME_SCOUT_UNDECIDED',
  'source_hashes':{p:sha(p)for p in(BASE,POLAR,EXT,CAP)},'seed':seed,'attempt_limit':attempts,
  'attempts':tries,'certificate':cert,
  'scope':'LP proposes points only. Every saved escape is checked exactly in separate radical extensions against all86walls, with positive norm and future orientation. Failure to find an escape is not a finite-volume certificate.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--attempts',type=int,default=200);p.add_argument('--seed',type=int,default=20262086);a=p.parse_args()
 if a.verify:
  out=read(OUT);assert all(sha(p)==h for p,h in out['source_hashes'].items())
  _,_,_,Q,h,hc,roots=setup();c=out['certificate'];assert c
  actual=check(vec(c['vector']),Q,hc,roots)
  assert all(actual[k]==c[k]for k in actual)
 else:
  out=run(a.attempts,a.seed);(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':out['status'],'attempts':len(out['attempts']),'norm':out['certificate']['norm']if out['certificate']else None}))
