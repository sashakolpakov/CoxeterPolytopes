#!/usr/bin/env python3
"""Scout a stronger linear relaxation for the 51 unbounded I2 label tails."""
import json
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];P='data/enumerations/fresh_h9_component_bridge_20261013'

def run():
 d=json.loads((ROOT/(P+'_faces_scout.json')).read_text());cert=json.loads((ROOT/(P+'_verify.json')).read_text());cap=json.loads((ROOT/(P+'_caps.json')).read_text());theorem=json.loads((ROOT/(P+'_continuous_verify.json')).read_text())
 G=np.array(d['gram']);R=np.vstack([d['normals'],[[sum(map(int,x))/(2*10**80)for x in r]for r in cert['normal_enclosures']]])
 names=d['names']+[r['name']for r in cert['recipes']];center=np.array(d['center']);phi=(1+np.sqrt(5))/2;y=np.array([float(Fraction(x['a']))+phi*float(Fraction(x['b']))for x in cap['escape']]);w=-np.linalg.solve(G,y)
 edges=[r['edge']for r in cert['new_vertex_frontier']if len(r['end_walls'])==1];cov=R@G;rows=[]
 for tail in theorem['uncovered_I2_tails']:
  ss=edges[tail['edge']];ids=[names.index(n)for n in ss];A=R[ids];Mi=np.linalg.inv(A@G@A.T);proj=lambda v:v-(A@G@v)@Mi@A
  tv=proj(center);tv/=np.sqrt(-tv@G@tv);ev=proj(w);ev+=(ev@G@tv)*tv;ev/=np.sqrt(ev@G@ev)
  # Scale the candidate so its sole nonzero edge pairing is -1. This
  # removes the angle parameter entirely; the LP also drops its norm bound.
  base=-R[names.index(tail['wall'])]/2;outside=[j for j in range(104)if j not in ids]
  ac=np.c_[cov[outside]@tv,cov[outside]@ev];bc=-(cov[outside]@base);objective=np.array([tv@y,ev@y]);constant=float(base@y)
  res=linprog(objective,A_ub=ac,b_ub=bc,bounds=[(None,None)]*2,method='highs')
  row={**tail,'lp_status':int(res.status)}
  if res.success:
   row.update({'target_minimum':float(res.fun+constant),'solution':res.x.tolist(),'residual_norm':float(res.x[1]**2-res.x[0]**2),'active_walls':[names[outside[i]]for i in np.flatnonzero(abs(res.ineqlin.residual)<1e-7)],'dual_support':[(names[outside[i]],float(-v))for i,v in enumerate(res.ineqlin.marginals)if abs(v)>1e-8]})
  rows.append(row)
 out={'status':'NUMERICAL_LINEAR_RELAXATION_OF_ALL_I2_HEIGHTS','rows':rows,'scope':'Only numerical LP proposals. Scaling pair(a,r)=-1 removes mu=2cos(pi/m); dropping the residual norm constraint q in(0,1/6] enlarges every m>=6 family. A positive exact dual lower bound would exclude every real coupling here; no exact proof yet.'}
 (ROOT/(P+'_dihedral_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
 print([(r['edge'],r['wall'],r.get('target_minimum'),r.get('dual_support'))for r in rows],flush=True)

if __name__=='__main__':run()
