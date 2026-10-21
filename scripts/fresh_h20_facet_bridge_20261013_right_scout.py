#!/usr/bin/env python3
"""All-right19-support facet bridges: fourteen unit-circle angle candidates."""
import numpy as np,json
from scipy.linalg import null_space
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT
Q,R,h,ids,*_=setup();N=null_space(R[ids+[175]]@Q);v,V=np.linalg.eigh(N.T@Q@N);N=N@V@np.diag(1/np.sqrt(v));P=R@Q@N;ns=np.einsum('ij,jk,ik->i',R,Q,R);a=P[18]/np.sqrt(2);norm=np.linalg.norm(a);w=np.array([-a[1],a[0]])/norm;labels=np.arange(2,121);allowed=np.cos(np.pi/labels);out=[]
for m in range(2,9):
 c=np.cos(np.pi/m)
 for sg in[-1,1]:
  coeff=-c*a/(norm*norm)+sg*np.sqrt(1-c*c/(norm*norm))*w;r=N@coeff;t=-R@Q@r/np.sqrt(ns);good=(t>=1-1e-7)|(np.min(abs(t[:,None]-allowed),axis=1)<1e-7);bad=np.flatnonzero(~good).tolist();center=float(r@Q@h)
  row={'label':m,'branch':sg,'normal':r.tolist(),'center_pairing':center,'bad':bad,'bad_ratios':[float(t[i])for i in bad]};out.append(row);print(m,sg,'center',center,'bad',[(i,round(float(t[i]),8))for i in bad],flush=True)
(ROOT/'data/enumerations/fresh_h20_facet_bridge_20261013_right_scout.json').write_text(json.dumps({'status':'NUMERICAL_ALL_RIGHT19_SUPPORT_SCOUT','support':ids+[175],'anchor':18,'anchor_projection_ratio':norm,'candidates':out},indent=2)+'\n')
