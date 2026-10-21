#!/usr/bin/env python3
"""Release one finite neighbor, retain perpendicularity to facet175, solve two angles."""
import json,numpy as np
from scipy.linalg import null_space
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT
Q,R,h,neighbors,*_=setup();ns=np.einsum('ij,jk,ik->i',R,Q,R);cov=R@Q;N=null_space(cov[neighbors+[175]]);v,V=np.linalg.eigh(N.T@Q@N);N=N@V@np.diag(1/np.sqrt(v));hp=h@Q@N;w=N@np.array([-hp[1],hp[0]])/np.linalg.norm(hp);c=h-(h@Q@R[175])*R[175];t=np.sqrt(-c@Q@c);targets=[]
for sg in[-1,1]:
 for eps in[1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]:
  x=c+sg*t*(1+eps)*w
  if np.max(cov[np.arange(176)!=175]@x)<0:break
 targets.append(x)
print('targets',[(float(x@Q@x),float(np.max(cov@x)))for x in targets],flush=True)
labels=[2,3,4,5,6,8,10,12];allowed=np.array([0]+[np.cos(np.pi/m)for m in range(3,121)]);caps=[];profiles=0;seen=set();summary=[];best=[]
for released in neighbors:
 Z=[i for i in neighbors if i!=released]+[175]
 for anchor in[18,50,158]:
  J=Z+[released,anchor];M=cov[J];u,s,vh=np.linalg.svd(M,full_matrices=True)
  if s[-1]<1e-8:continue
  n=vh[-1];inverse=np.linalg.pinv(M,rcond=1e-11)
  for m in labels:
   for k in labels:
    profiles+=1;b=np.zeros(20);b[-2]=-np.sqrt(2*ns[released])*np.cos(np.pi/m);b[-1]=-np.sqrt(2*ns[anchor])*np.cos(np.pi/k);r0=inverse@b;A=n@Q@n;BB=2*r0@Q@n;CC=r0@Q@r0-2;D=BB*BB-4*A*CC
    if D<-1e-8:continue
    if abs(A)<1e-10:ts=[-CC/BB]if abs(BB)>1e-10 else[]
    else:ts=[(-BB+sg*np.sqrt(max(0,D)))/(2*A)for sg in[-1,1]]
    for branch,tt in enumerate(ts):
     r=r0+tt*n
     if abs(r@Q@r-2)>1e-6 or r@Q@h>=-1e-6 or np.max(abs(M@r-b))>1e-7:continue
     cuts=[float(r@Q@x)for x in targets]
     z=-cov@r/np.sqrt(2*ns);good=(z>=1-1e-7)|(np.min(abs(z[:,None]-allowed),axis=1)<1e-7);bad=np.flatnonzero(~good).tolist();row={'released':released,'zero_support':Z,'anchor':anchor,'labels':[m,k],'branch':branch,'normal':r.tolist(),'cuts':cuts,'bad':bad,'bad_ratios':[float(z[i])for i in bad]}
     if len(bad)<4:best.append(row)
     if not bad:
      key=tuple(np.round(r,8))
      if key not in seen:seen.add(key);caps.append(row);print('candidate',len(caps),released,anchor,m,k,branch,cuts,flush=True)
 summary.append({'released':released,'profiles':profiles,'caps':len(caps)});print(summary[-1],flush=True)
(ROOT/'data/enumerations/fresh_h20_facet_bridge_20261013_release_scout.json').write_text(json.dumps({'status':'NUMERICAL_RELEASED_NEIGHBOR_FACET_BRIDGES','source':'data/enumerations/fresh_h20_coupled_completion_20261013_second_exact.json','labels':labels,'anchors':[18,50,158],'profiles':profiles,'targets':[x.tolist()for x in targets],'caps':caps,'near_misses':best,'summary':summary},indent=2)+'\n')
