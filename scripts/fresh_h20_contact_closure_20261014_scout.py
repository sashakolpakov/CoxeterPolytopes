#!/usr/bin/env python3
"""Small LP proposals for exact dual certificates; no numerical acceptance."""
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_contact_closure_20261014'
def run():
 D=json.loads((ROOT/'results/fresh_h20_coupled_lift_20261014.json').read_text());src=json.loads((ROOT/'results/fresh_dimension_construction_routes_20261014.json').read_text());Q=np.array([[float(F(x))for x in row]for row in D['form']]);N=np.array([[float(F(x['a']))+float(F(x['b']))*np.sqrt(2)for x in row]for row in src['normals']]);E=np.array([float(F(x))for x in D['strict_escape']['point']]);out=[]
 for sg in product([1,-1],repeat=2):
  flip=np.ones(21);flip[19:]=sg
  for si,rec in enumerate(D['base_recipes']):
   T=N[rec['support']]*flip
   _,_,V=np.linalg.svd(T@Q);X=V[18:].T;A=N@Q@X;b=E@Q@X
   lp=linprog(np.ones(60),A_eq=A.T,b_eq=b,bounds=(0,None),method='highs')
   ids=[i for i,x in enumerate(lp.x)if x>1e-7]if lp.success else []
   if lp.success:assert len(ids)==3
   out.append(dict(template=si,signs=sg,dual_walls=ids,dual_feasible=bool(lp.success)))
 # Release one support wall and ask only whether the same dual-cone obstruction survives.
 released=[]
 for row in out:
  sg=np.ones(21);sg[19:]=row['signs'];S=D['base_recipes'][row['template']]['support']
  for j in S:
   T=N[[i for i in S if i!=j]]*sg;_,_,V=np.linalg.svd(T@Q);X=V[17:].T
   lp=linprog(np.ones(60),A_eq=(N@Q@X).T,b_eq=E@Q@X,bounds=(0,None),method='highs')
   released.append(dict(template=row['template'],signs=row['signs'],released=j,dual_feasible=bool(lp.success)))
 return dict(status='NUMERICAL_DUAL_PROPOSALS_ONLY',dual_proposals=out,released_support_relaxation=released,scope='Only dual supports become exact in the separate primary checker; the PD17 release feasibility screen is numerical, not a complete cap search.')
if __name__=='__main__':
 d=run();(ROOT/f'results/{TAG}_scout.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],len(d['dual_proposals']),'PD17 not dual obstructed',sum(not x['dual_feasible']for x in d['released_support_relaxation']))
