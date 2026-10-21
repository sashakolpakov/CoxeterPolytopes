#!/usr/bin/env python3
"""Numerical separation of three remaining cusp-normal support patterns.

This convex relaxation uses all current caps through constraint generation.
It is a proposal/upper-bound scout, never an exact exclusion certificate.
"""
import json
from fractions import Fraction as F
import numpy as np
import cvxpy as cp
from fresh_h22_coupled_cusp_20261013_scout import setup,P,PRE
from fresh_h22_new_cusp_completion_20261013_verify import setup as cap_setup
from fresh_h22_new_cusp_completion_20261013_tight import shifted,DEN

def run():
 Q,W,G,h,X,loc,cs,fin,gi,U,ds,dt,ws,wt,dw,mask,ps,pt,dp=setup()
 supports=[[],[28],[1]];wbar=ws-wt[:,None]*ds/dt
 bases=[];rhos=[]
 for support in supports:
  p=np.array([-1. if k in support else 0. for k in fin])
  bases.append((gi@p)@wbar[fin]);rhos.append(2-p@gi@p)
 bases=np.array(bases).T
 delta=np.array([float(F(x))for x in loc['null_ray']]);delta[-1]*=np.sqrt(11)
 oldframe=np.column_stack([X@U[:21],delta[:-1],X@bases[:21]])
 heightframe=np.r_[U[-1],delta[-1],bases[-1]]
 blocks=[]
 def add(name,old,z):
  M=old@(Q@oldframe);z=np.broadcast_to(z,len(old)).copy()
  blocks.append((name,M,z));print('prepared',name,len(old),flush=True)
 _,_,_,_,_,_,D,V,B,A,ids,sizes=cap_setup()
 zz=np.load(P/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz')
 add('old_null',zz['old_null_projections'],zz['four_heights']/4)
 zz=np.load(P/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')
 add('old_timelike',zz['twice_old_projections']/2,np.sqrt(2.5))
 zz=np.load(P/'data/enumerations/fresh_h22_gap_face_orbit_20261013.npz')
 add('old_k6',zz['six_old_projections']/6,np.sqrt(11))
 for name,old,t,ix,index in shifted(D,V,B,A,ids,sizes,Q,True):add(name,old/DEN,t*np.sqrt(11)/10)
 add('center',h[None,:],4.)
 del D,V,B,A,ids,sizes,old,zz
 lhs0=np.column_stack([ws[mask]@U,dw[mask]])
 obj=np.r_[U.T@ps,dp];rows=[]
 for bi,(support,rho)in enumerate(zip(supports,rhos,strict=True)):
  base=bases[:,bi];rhs0=-ws[mask]@base
  active=set();iterations=[];solution=None
  for step in range(50):
   chosen=sorted(active);extra=[];erhs=[]
   for g,i,s in chosen:
    _,M,z=blocks[g];row=M[i]+s*z[i]*heightframe
    scale=max(1.,np.linalg.norm(row[:5]));extra.append(row[:5]/scale);erhs.append(-row[5+bi]/scale)
   lhs=np.vstack([lhs0,np.array(extra).reshape((-1,5))]);rhs=np.r_[rhs0,erhs]
   v=cp.Variable(5);cons=[lhs@v<=rhs,cp.norm(v[:4])<=np.sqrt(rho)]
   prob=cp.Problem(cp.Maximize(obj@v+base@ps),cons)
   prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=500)
   assert prob.status in ('optimal','optimal_inaccurate'),prob.status
   vv=v.value;value=float(prob.value);print('support',support,'step',step,'upper',value,'active',len(active),flush=True)
   rec={'iteration':step,'upper':value,'status':prob.status,'active':len(active)};iterations.append(rec)
   solution={'variables':vv.tolist(),'upper':value,'old_duals':cons[0].dual_value[:len(lhs0)].tolist(),
             'extra_duals':cons[0].dual_value[len(lhs0):].tolist(),'sphere_dual':float(cons[1].dual_value),
             'active_caps':[[blocks[g][0],i,s]for g,i,s in chosen]}
   if value< -1e-7:break
   new=[];worst=0.;co=np.r_[vv,np.eye(3)[bi]];height=float(heightframe@co)
   for g,(name,M,z)in enumerate(blocks):
    vals=M@co+z*abs(height);worst=max(worst,float(max(vals)))
    wrong=np.flatnonzero(vals>2e-6)
    if len(wrong):
     wrong=wrong[np.argsort(vals[wrong])[-50:]];new.extend((g,int(i),1 if height>=0 else -1)for i in wrong)
   rec['worst_cap_acute_violation']=worst
   if not new:break
   before=len(active);active.update(new);assert len(active)>before
  else:raise RuntimeError('constraint-generation iteration bound')
  rows.append({'support':support,'residual_norm':float(rho),'iterations':iterations,'solution':solution})
  out={'status':'NUMERICAL_NEW_CUSP_OTHER_SUPPORTS_FULL_CAP_ACUTE_RELAXATION','rows':rows,
       'scope':'Only support empty, old28, or old1, with nonzero prescribed pairing -1. Norm sphere relaxed to ball. Selected whole-current-chamber acute constraints give numerical upper bounds; no exact exclusion or finite-angle compatibility is claimed.'}
  (P/'data/enumerations'/f'{PRE}_other_scout.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':run()
