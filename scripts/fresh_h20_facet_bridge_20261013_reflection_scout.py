#!/usr/bin/env python3
"""Numerical coherent reflected-subchamber replacement of facet175."""
import json,numpy as np
from scipy.optimize import linprog
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT
from fresh_h20_coupled_completion_20261013_orbit_exact import OUT as SOURCE
from fresh_h20_coupled_completion_20261013 import read,num
OUT='data/enumerations/fresh_h20_facet_bridge_20261013_reflection_scout.json'
Q,R,h,ids,Qe,he,hc,roots,f,g=setup();R=R[:175];P=np.array([num(f,f.parse(x))for x in read(SOURCE)['final_hyperideal_point']['vector']]);cs=np.array([0]+[np.cos(np.pi/m)for m in range(3,121)]);ns=np.einsum('ij,jk,ik->i',R,Q,R);caps=[];rejected=[];fixed=[]
for i in range(175):
 pp=P@Q@R[i]
 if abs(pp)<1e-7:fixed.append(i);continue
 r=R[i]-2*pp/(P@Q@P)*P
 if r@Q@P<1e-7 or r@Q@h>=0:continue
 z=-R@Q@r/np.sqrt(ns[i]*ns);bad=np.flatnonzero(~((z>=1-1e-7)|(np.min(abs(z[:,None]-cs),axis=1)<1e-7))).tolist();row={'original':i,'normal':r.tolist(),'finite_neighbors':np.flatnonzero(z<1-1e-7).tolist()}
 if not bad:caps.append(row)
 else:rejected.append(dict(row,bad=bad,bad_cosines=[float(z[k])for k in bad]))
print('new',len(caps),'rejected',len(rejected),'fixed',len(fixed),flush=True)
C=np.array([x['normal']for x in caps]);U=np.vstack([R,C]);nu=np.einsum('ij,jk,ik->i',U,Q,U);A=U@Q/np.sqrt(nu[:,None]);hn=h/np.sqrt(-h@Q@h);rng=np.random.default_rng(2026101322);trials=[];escape=None
for k in range(200):
 x=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
 if x.success:
  q=float(x.x@Q@x.x);trials.append(q)
  if q>1e-7:escape={'vector':x.x.tolist(),'norm':q,'active_indices':np.flatnonzero(abs(A@x.x)<1e-7).tolist(),'attempt':k};break
print('escape',escape,'trials',len(trials),flush=True)
(ROOT/OUT).write_text(json.dumps({'status':'NUMERICAL_COHERENT_POLAR_REFLECTION_REPLACEMENT','source':SOURCE,'retained':175,'reflected':caps,'rejected':rejected,'fixed':fixed,'wall_count':len(U),'trials':trials,'escape':escape},indent=2)+'\n')
