#!/usr/bin/env python3
"""Numerical cone-projection screen of the bounded noncutting endpoint bank."""
import argparse,json
from collections import Counter
from hashlib import sha256
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_ridge_pair_20261016_polar as P
ROOT=P.ROOT;TAG='fresh_nonarithmetic_ridge_pair_20261016_compatibility'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 Q,N,H,E=P.geometry();dom=P.P.B.read(P.DOM);ends=P.P.B.read(P.END);cov=[Q*n for n in N]+[-Q*E];hcov=Q*H
 rays=[]
 for row in dom['records']:
  v=mp.matrix(P.P.null([list(cov[j])for j in row['active']]));v/=-(v.T*hcov)[0]
  rays.append(v)
 endpoints=[]
 for i,x in enumerate(ends['candidates']):
  u=mp.matrix([mp.mpf(z)for z in x['point']])
  if not any(max(abs(a-b)for a,b in zip(u,v))<mp.mpf('1e-25')for _,v in endpoints):endpoints.append((i,u))
 results=[]
 for i,u in endpoints:
  vals=[];orientation=[]
  hu=(H.T*Q*u)[0]
  for v in rays:
   q=(v.T*Q*v)[0];g=(v.T*Q*u)[0];vals.append(q-g*g/2);orientation.append(-1-g*hu/2)
  positive=[j for j,v in enumerate(vals)if v>mp.mpf('1e-25')];mixed=[j for j,v in enumerate(orientation)if v>=-mp.mpf('1e-25')]
  results.append(dict(endpoint_index=i,recipe={k:v for k,v in ends['candidates'][i].items()if k not in['point','bad']},positive_projection_ray_indices=positive,nonfuture_projection_ray_indices=mixed,maximum_projection_norm=mp.nstr(max(vals),40),minimum_projection_norm=mp.nstr(min(vals),40)))
  print('endpoint',i,'positive projections',len(positive),'max',mp.nstr(max(vals),12),flush=True)
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',P.DOM,P.END,P.CL,P.CO,P.P.BATCH]
 return dict(status='NUMERICAL_CUTTER_CONE_PROJECTION_SCREEN_OF_ENDPOINT_BANK',sources={p:sha(p)for p in paths},cone_ray_count=len(rays),distinct_endpoints=len(endpoints),endpoints_with_positive_projection=sum(bool(r['positive_projection_ray_indices'])for r in results),records=results,scope='Numerical high-precision diagnostic only. A positive projected norm at a cutter-cone generator is a necessary opening for a finite contact, not a Coxeter pair or geometric completion. Negative numerical results require exact interval replay before any obstruction claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved projection screen')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['distinct_endpoints'],d['endpoints_with_positive_projection'])
