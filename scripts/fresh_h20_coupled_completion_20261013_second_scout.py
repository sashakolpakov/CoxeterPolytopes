#!/usr/bin/env python3
"""Numerical coupled repair at the complete first-cap orbit's remaining gap."""
import json,numpy as np
from fresh_h20_coupled_completion_20261013_dotted_exact import setup,OUT as FIRST
from fresh_h20_coupled_completion_20261013 import num,read,ROOT
from fresh_dimension_bridge_20261013_verify import parse
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension
SRC='data/enumerations/fresh_h20_coupled_completion_20261013_orbit_scout.json'
OUT='data/enumerations/fresh_h20_coupled_completion_20261013_second_scout.json'
Qe,he,hc,roots=setup();j=read(FIRST);c=j['new_cap'];f=Extension(B,parse(c['field']['radicand']));roots.append((f,tuple(map(f.parse,c['normal']))));Q=np.array([[num(B,x)for x in r]for r in Qe]);h=np.array([num(B,x)for x in he]);R=np.array([[num(g,x)for x in r]for g,r in roots]);j=read(SRC);R=np.vstack([R,[x['normal']for x in j['accepted']]]);target=np.array(j['all_accepted_escape']['vector']);p=-target/(target@Q@R[18]);A=R[18];Bv=R[30]/3;C=R[46]/3;ns=np.einsum('ij,jk,ik->i',R,Q,R);angles=np.array([0]+[2*np.cos(np.pi/m)for m in(3,4,5,6,8,10,12)]+[2]);best=[];caps=[]
def test(r,meta):
 if r@Q@h>=0 or r@Q@target<1e-8:return
 t=-R@Q@r/np.sqrt(2*ns);good=(t>=1-1e-7)|(np.min(abs(t[:,None]-angles/2),axis=1)<1e-7);bad=np.flatnonzero(~good).tolist();rec=dict(meta,normal=r.tolist(),bad=bad,ratios=[float(t[i])for i in bad],cut=float(r@Q@target));best.append(rec)
 if not bad:caps.append(rec)
# finite leaf coefficients first
for ia,a in enumerate(angles):
 for ib,b in enumerate(angles):
  for ic,c in enumerate(angles):
   u=a*A+b*Bv+c*C;aa=p@Q@p;bb=2*p@Q@u;cc=u@Q@u-2;D=bb*bb-4*aa*cc
   if D<-1e-8:continue
   for k in((-bb+np.sqrt(max(0,D)))/(2*aa),(-bb-np.sqrt(max(0,D)))/(2*aa)):test(k*p+u,{'type':'finite_leaves','labels':[ia,ib,ic],'coefficients':[a,b,c],'k':k})
# free dotted leaf coefficients: a=c, two finite old-leaf anchors
for il,L in enumerate(angles):
 for im,M in enumerate(angles):
  u=L*p+(.75*L-.5*M)*Bv;v=2*p+A+1.5*Bv+C;aa=v@Q@v;bb=2*v@Q@u;cc=u@Q@u-2;D=bb*bb-4*aa*cc
  if D<-1e-8:continue
  for a in((-bb+np.sqrt(max(0,D)))/(2*aa),(-bb-np.sqrt(max(0,D)))/(2*aa)):
   b=1.5*a+.75*L-.5*M
   if a<2-1e-8 or b<2-1e-8:continue
   test(u+a*v,{'type':'free_dotted_leaves','labels':[il,im],'coefficients':[a,b,a],'k':2*a+L})
best.sort(key=lambda x:len(x['bad']));print('qP',p@Q@p,'success',len(caps),'best',[(x['type'],x['labels'],x['coefficients'],x['k'],x['bad'],x['ratios'])for x in best[:12]],flush=True)
(ROOT/OUT).write_text(json.dumps({'status':'NUMERICAL_SECOND_COUPLED_REPAIR','source':SRC,'polar':p.tolist(),'caps':caps,'best':best[:30]},indent=2)+'\n')
