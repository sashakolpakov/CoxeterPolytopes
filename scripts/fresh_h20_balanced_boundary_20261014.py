#!/usr/bin/env python3
"""Numerical boundary-center cap scout with finite omitted seed contacts.

Prescribes at most three mixed finite contacts among19 right-contact slots.
The remaining three seed contacts can be finite, parallel, or dotted. A fixed
centralR2 contact determines up to two real radii. No positive geometric claim
is made without a separate exact reconstruction.
"""
from itertools import combinations,product
import json,math,time,argparse
import numpy as np
from fresh_h20_balanced_20261014 import ROOT,data,sha
TAG='fresh_h20_balanced_boundary_20261014'
def run():
 st=time.time();D=data();VV=np.array([v['point']for v in D['vertices']]);v=D['vertices'][0];z=v['point'];I=v['retained'];W=D['vectors'][I];dual=np.linalg.inv(W@W.T)@W;J=dual@dual.T;zy=dual@z
 target=.4*VV;targetq=np.sum(target*target,axis=1);cs=np.array([math.cos(math.pi/m)for m in range(2,61)]);labels=[3,4,5,6];values={m:math.cos(math.pi/m)for m in labels};central=list(range(2,17))+[0]
 profiles=trials=accepted=nonendpoint=0;hits=[];best=None;accepted_rows=[]
 for k in range(4):
  for inds0 in combinations(range(19),k):
   inds=list(inds0)
   for labs in product(labels,repeat=k):
    profiles+=1;co=np.array([values[m]for m in labs]);p=float(co@J[np.ix_(inds,inds)]@co)if k else 0.;zyy=float(co@zy[inds])if k else 0.;y=co@dual[inds]if k else np.zeros(19)
    sums=np.zeros(3)
    for ii,cc in zip(inds,co):sums[0 if ii<6 else 1 if ii<12 else 2]+=cc
    for cm in central:
     trials+=1;cc=1 if cm==0 else math.cos(math.pi/cm);A=p-1;bb=-2*zyy-4*cc;C=192;dis=bb*bb-4*A*C
     if abs(A)<1e-10:roots=[-C/bb]
     elif dis>=0:roots=[(-bb+math.sqrt(dis))/(2*A),(-bb-math.sqrt(dis))/(2*A)]
     else:continue
     for R in roots:
      if R<=1.8-1e-9:continue # A central-contact ball of R<=1.8 cannot cut |x|=5.6.
      omitted=np.array([7,7,8])/R-sums;pos=np.searchsorted(cs,omitted);a=np.clip(pos,0,len(cs)-1);b=np.clip(pos-1,0,len(cs)-1);err=np.minimum(abs(omitted-cs[a]),abs(omitted-cs[b]))
      if not np.all((omitted>=1-1e-8)|(err<1e-8)):continue
      accepted+=1;nonendpoint+=int(p>1+1e-8);c=z-R*y;power=targetq+float(c@c)-2*target@c-R*R;j=int(np.argmin(power));row=dict(contacts=[[I[i],m]for i,m in zip(inds,labs)],central_label=cm,radius=R,spatial_schur_cost=p,omitted_cosines=omitted.tolist(),target_index=j,minimum_power=float(power[j]))
      accepted_rows.append(row)
      if best is None or row['minimum_power']<best['minimum_power']:best=row
      if power[j]<-1e-8:hits.append(row)
 out=dict(status='BOUNDED_NUMERICAL_BOUNDARY_CENTRAL_CONTACT_SCOUT',sources={f:sha(f)for f in['scripts/fresh_h20_balanced_20261014.py',f'scripts/{TAG}.py']},retained_support=I,prescribed_contact_limit=3,retained_finite_labels=labels,omitted_finite_labels_tested=list(range(2,61)),central_labels=central,profile_count=profiles,profile_label_trials=trials,accepted_count=accepted,nonendpoint_count=nonendpoint,strict_target_cutters=len(hits),best=best,cutters=hits,accepted_recipes=accepted_rows,elapsed_seconds=time.time()-st,scope='Bounded numerical all22seed and centralR2 pair proposals. All392 transported .4vertex targets tested. Omitted finite contacts are allowed, unlike the prior dotted-only scout. No exact cap, complete label census, joint compatibility, or finite-volume claim.')
 (ROOT/'results'/f'{TAG}.json').write_text(json.dumps(out,indent=2)+'\n');print('profiles',profiles,'accepted',accepted,'nonendpoint',nonendpoint,'cutters',len(hits),'best',best);print(hits[:3])
if __name__=='__main__':run()
