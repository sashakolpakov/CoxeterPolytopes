#!/usr/bin/env python3
"""Numerical full-H12 residual optimization; no exact exclusion by itself."""
import json
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from hashlib import sha256
import numpy as np
import cvxpy as cp

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_completion_gate_20261014'
PATH='results/fresh_h20_e8_cusp_completion_20261014.json'
def value(a):return sum(float(F(x))*np.sqrt(d)for x,d in zip(a,[1,2,79,158]))
def run():
 d=json.loads((ROOT/PATH).read_text());f=json.loads((ROOT/'results/fresh_h20_contact_closure_20261014_e8_face.json').read_text())
 Q=np.array([[float(F(x))for x in row]for row in d['form']]);N=np.array([[value(x)for x in r]for r in d['normals']]);E=np.array([value(x)for x in d['strict_escape']['point']]);J=d['E8_support'];B=N[J]@Q@N[J].T;I=np.linalg.inv(B)
 A=[i-1 for i in f['orthogonal_face_walls']];T=np.array([value(x)for x in f['normalized_timelike_basis']]);C=N[A]@Q@N[A].T;L=np.linalg.cholesky(C)
 spatial=np.linalg.solve(L,N[A]@Q@N.T).T;time=np.sqrt(13)*(N@Q@T);M=np.column_stack([spatial,time]);M[np.max(np.abs(M),axis=1)<1e-8]=0
 ep=np.r_[np.linalg.solve(L,N[A]@Q@E),np.sqrt(13)*T@Q@E]
 contributions=-I@(N[J]@Q@E);wallweights=I@(N[J]@Q@N.T)
 scales=np.maximum(np.linalg.norm(M,axis=1),1);x=cp.Variable(13);rhs=cp.Parameter(79);radius=cp.Parameter(nonneg=True)
 linear=M/scales[:,None]@x<=rhs;cone=cp.SOC(x[12],cp.hstack([x[:12],radius]));problem=cp.Problem(cp.Maximize(ep@x),[linear,cone])
 rows=[]
 profiles=[]
 for size in [1,2]:
  for support in combinations(range(8),size):
   c=np.zeros(8);c[list(support)]=1;profiles.append((support,c))
 for i in range(8):
  c=np.zeros(8);c[i]=np.sqrt(2);profiles.append(((i,),c))
 for support,c in profiles:
  D=float(c@I@c);rhs.value=(c@wallweights)/scales;radius.value=np.sqrt(max(0,D-2))
  problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-10,tol_gap_rel=1e-10,max_iter=250)
  assert problem.status in ['optimal','optimal_inaccurate']
  y=x.value;objective=float(contributions@c+problem.value);norm=float(y[:12]@y[:12]-y[12]**2);dual=linear.dual_value/scales
  rows.append(dict(support=list(support),contact_vector=c.tolist(),w_norm=D,relaxed_max_target=objective,residual_norm=norm,target_residual_norm=2-D,orthonormal_face_coordinates=y.tolist(),rational_face_coordinates=np.r_[np.linalg.solve(L.T,y[:12]),np.sqrt(13)*y[12]].tolist(),active_walls=np.flatnonzero(dual>1e-5).tolist(),positive_duals=[[i,float(t)]for i,t in enumerate(dual)if t>1e-7],max_acute_residual=float(np.max(M@y-c@wallweights)),cone_dual=[np.asarray(z).tolist()for z in cone.dual_value],status=problem.status))
  print(support,round(objective,8),'norm',round(norm,8),'active',rows[-1]['active_walls'],flush=True)
 out=dict(status='NUMERICAL_ALL_Y_ACUTE_SOCP_SUPERSET_PROFILE_SCREEN',source_sha256=sha256((ROOT/PATH).read_bytes()).hexdigest(),profiles=rows,scope='36one/two-node pi3 profiles and8one-node pi4 profiles. Residual norm inequality is a convex relaxation; numerical duals are proposals for exact bounds, not certificates.')
 (ROOT/f'results/{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n')
 return out
if __name__=='__main__':run()
