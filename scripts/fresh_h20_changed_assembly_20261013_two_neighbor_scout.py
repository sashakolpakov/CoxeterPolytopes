#!/usr/bin/env python3
"""Locate target-cutting residual spheres after freeing two neighbor rows."""
import json,numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT
from fresh_h20_coupled_completion_20261013 import read,num
from fresh_h20_coupled_completion_20261013_second_exact import OUT as SOURCE
from fresh_dimension_bridge_two_angle_20261013_independent import B
OUT='data/enumerations/fresh_h20_changed_assembly_20261013_two_neighbor_scout.json'
Q,R,h,nei,*_=setup();cov=R@Q;qn=np.einsum('ij,jk,ik->i',R,Q,R);x=np.array([num(B,B.parse(z))for z in read(SOURCE)['strict_escape']['vector']]);out=[]
for j in[9,3,11,40]:
 for k in nei:
  if k==j:continue
  I=[i for i in nei if i not in(j,k)]+[175];J=[i for i in range(176)if i not in I]
  N=null_space(cov[I]);ev,E=np.linalg.eigh(N.T@Q@N)
  if min(ev)<1e-8:continue
  N=N@E@np.diag(1/np.sqrt(ev));c=R[j]-N@(N.T@Q@R[j]);rho=qn[j]-c@Q@c
  if rho<1e-8:continue
  L=N*np.sqrt(rho);A=cov[J]@L;b=-cov[J]@c;scale=np.maximum(np.linalg.norm(A,axis=1),1);A=A/scale[:,None];b=b/scale
  tv=x@Q@L;tc=x@Q@c
  ans=minimize(lambda u:-tv@u,np.zeros(len(N.T)),jac=lambda u:-tv,constraints=[{'type':'ineq','fun':lambda u:b-A@u,'jac':lambda u:-A},{'type':'ineq','fun':lambda u:1-u@u,'jac':lambda u:-2*u}],method='SLSQP',options={'maxiter':500,'ftol':1e-10})
  value=float(tc+tv@ans.x);r=c+L@ans.x
  row={'copied_neighbor':j,'released_second':k,'zero_fixed_indices':I,'residual_dimension':N.shape[1],'radius_squared':float(rho),'solver_success':bool(ans.success),'target_upper_proposal':value,'sphere_parameter_norm':float(ans.x@ans.x),'max_acute_residual':float(max(A@ans.x-b)),'normal_proposal':r.tolist()}
  out.append(row)
  print(j,k,'success',ans.success,'target',value,'u2',row['sphere_parameter_norm'],flush=True)
(ROOT/OUT).write_text(json.dumps({'status':'NUMERICAL_TWO_NEIGHBOR_TARGET_FEASIBILITY','source':SOURCE,'families':out},indent=2)+'\n')
