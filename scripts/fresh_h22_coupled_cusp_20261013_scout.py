#!/usr/bin/env python3
"""Numerical coupled affine-A20 cusp completion; exact certification separate.

Join the five finite chains at the newly exposed cusp by five new norm-two
bridges, leaving the existing affine A1 perpendicular. Search each possible
individual bridge over the complete four-dimensional tangent sphere using a
convex-ball relaxation and retained-wall dotted constraints.
"""
import argparse
import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import numpy as np
import scipy.linalg as la
from scipy.optimize import minimize
from fresh_h22_borcherds_extension_20261013 import prepare

P=Path(__file__).resolve().parents[1]
PRE='fresh_h22_coupled_cusp_20261013'

def setup():
 Q,W,G,h,X,A,b,Qn,_=prepare()
 local=json.loads((P/'results/fresh_h22_lower_new_cusp_20261013.json').read_text())['representatives'][0]
 delta=np.array([float(F(x))for x in local['null_ray']]);delta[-1]*=np.sqrt(11)
 ds=np.r_[delta[:-1]@Q@X,delta[-1]];dt=-delta[:-1]@Q@h/np.sqrt(14)
 ws=np.column_stack([A,np.zeros(210)]);wt=b
 assert abs(ds@ds-dt*dt)<1e-5 and abs(dt-72/np.sqrt(14))<1e-6
 # Horospherical projections have zero time coordinate and unchanged Gram.
 wbar=ws-wt[:,None]*ds/dt
 ids=local['incident_old_wall_indices'];cs=[]
 for comp in local['components']:
  if comp['type'].startswith('finite'):
   js=[int(x[3:])for x in comp['walls']];sub=G[np.ix_(js,js)]
   if len(js)==1:ordered=js
   else:
    ends=np.flatnonzero(np.sum(sub==-1,axis=1)==1);ordered=[js[int(ends[0])]]
    while len(ordered)<len(js):
     nxt=[j for j in js if G[ordered[-1],j]==-1 and j not in ordered];assert len(nxt)==1;ordered+=nxt
   cs.append(ordered)
 fin=sum(cs,[]);gi=np.linalg.inv(G[np.ix_(fin,fin)])
 # An orthonormal frame of the four unused Euclidean tangent directions.
 U=la.null_space(np.vstack([ds,wbar[ids]]));assert U.shape==(22,4)
 assert np.max(abs(U.T@U-np.eye(4)))<1e-8
 pjson=json.loads((P/'results/fresh_h22_new_cusp_completion_20261013_lower_escape.json').read_text())
 pp=np.array(pjson['numerator'])/pjson['denominator'];ps=np.r_[pp[:-1]@Q@X,pp[-1]];pt=-pp[:-1]@Q@h/np.sqrt(14)
 dp=ds@ps-dt*pt
 dw=ws@ds-wt*dt
 assert max(dw)<1e-5
 mask=dw<-1e-5
 return Q,W,G,h,X,local,cs,fin,gi,U,ds,dt,ws,wt,dw,mask,ps,pt,dp

def run(acute=False):
 Q,W,G,h,X,local,cs,fin,gi,U,ds,dt,ws,wt,dw,mask,ps,pt,dp=setup()
 print('chains',cs,'target.delta',dp,flush=True)
 rows=[];target_u=U.T@ps
 for ci,cj in combinations(range(5),2):
  for wi in sorted(set((cs[ci][0],cs[ci][-1]))):
   for wj in sorted(set((cs[cj][0],cs[cj][-1]))):
    p=np.array([-1. if k in(wi,wj)else 0. for k in fin]);coef=gi@p
    wbar=ws-wt[:,None]*ds/dt
    aa=coef@wbar[fin];rho=2-p@gi@p
    assert rho>0 and abs(aa@aa+rho-2)<1e-6
    lhs=np.column_stack([ws[mask]@U,dw[mask]])
    rhs=(np.zeros(sum(mask)) if acute else -np.sqrt(2*np.diag(G)[mask]))-ws[mask]@aa
    objective=np.r_[target_u,dp]
    init=np.r_[np.zeros(4),max(1,float(np.max(rhs/dw[mask])))+1]
    constraints=[{'type':'ineq','fun':lambda v,rhs=rhs,lhs=lhs:rhs-lhs@v,
                  'jac':lambda v,lhs=lhs:-lhs},
                 {'type':'ineq','fun':lambda v,rho=rho:rho-v[:4]@v[:4],
                  'jac':lambda v:np.r_[-2*v[:4],0]}]
    sol=minimize(lambda v:-objective@v,init,jac=lambda v:-objective,
                 bounds=[(None,None)]*5,constraints=constraints,
                 method='SLSQP',options={'ftol':1e-10,'maxiter':1000})
    maxviol=float(np.max(lhs@sol.x-rhs));normerr=float(sol.x[:4]@sol.x[:4]-rho)
    basepair=float(aa@ps);value=float(basepair+objective@sol.x)
    row={'components':[ci,cj],'endpoints':[wi,wj],'rho_squared':float(rho),
         'success':bool(sol.success),'message':sol.message,'value':value,
         'u_and_shift':sol.x.tolist(),'max_violation':maxviol,'sphere_excess':normerr,
         'base_spatial':aa.tolist(),'active_old':[int(i)for i in np.flatnonzero(mask)[np.flatnonzero(abs(lhs@sol.x-rhs)<1e-6)]]}
    rows.append(row);print('bridge',wi,wj,'cut',value,'sphere',normerr,'ok',sol.success,flush=True)
 rows.sort(key=lambda r:r['value'],reverse=True)
 out={'status':'NUMERICAL_AFFINE_A20_COUPLED_CUSP_BRIDGE_RELAXATION','chains':cs,
      'fixed_cusp':'lower delta4','target_delta_pairing':float(dp),'tangent_frame':U.tolist(),
      'delta_spatial':ds.tolist(),'delta_time':float(dt),'candidates':rows,
      'constraint_mode':'acute' if acute else 'dotted',
      'scope':'All32 chain-endpoint pairs, each with a numerical convex-ball relaxation of its tangent sphere and all real shifts. Only original210 nonincident constraints in the stated mode are included. Positive values are proposals, negative numerical values are not exact exclusions. Whole current cap family and coherent five-wall completion remain to be checked.'}
 suffix='_acute' if acute else ''
 (P/'data/enumerations'/f'{PRE}{suffix}_scout.json').write_text(json.dumps(out,indent=2)+'\n')
 print('best',rows[0]['value'],rows[0]['endpoints'],flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--acute',action='store_true')
 run(parser.parse_args().acute)
