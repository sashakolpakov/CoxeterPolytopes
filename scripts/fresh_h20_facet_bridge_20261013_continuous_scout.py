#!/usr/bin/env python3
"""Positive2-plane norm circles with all variable remaining anchors strictly dotted."""
import numpy as np,json,math
from scipy.linalg import null_space
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT
Q,R,h,ids,*_=setup();ns=np.einsum('ij,jk,ik->i',R,Q,R);J=ids+[175];G=R[J]@Q@R[J].T;inv=np.linalg.inv(G);N=null_space(R[J]@Q);va,V=np.linalg.eigh(N.T@Q@N);N=N@V@np.diag(1/np.sqrt(va));P=R@Q@N;fixed=np.linalg.norm(P,axis=1)<1e-7;labels=[2,3,4,5,6,8,10,12,0];allowed=np.array([0]+[np.cos(np.pi/m)for m in range(3,121)]);caps=[];families=[];seen=set();count=0
hp=h@Q@N;w=N@np.array([-hp[1],hp[0]])/np.linalg.norm(hp);c=h-(h@Q@R[175])*R[175];tt=np.sqrt(-c@Q@c);targets=[c+sg*tt*(1+1e-4)*w for sg in[-1,1]]
def intervals(A,B,C):
 amp=math.hypot(A,B)
 if C>=amp-1e-10:return[]
 if C<=-amp+1e-10:return[(0.,2*math.pi)]
 ph=math.atan2(B,A);a=math.acos(-C/amp);lo=(ph+a)%(2*math.pi);length=2*math.pi-2*a;hi=lo+length
 return[(lo,hi)]if hi<=2*math.pi else[(lo,2*math.pi),(0.,hi-2*math.pi)]
def intersect(a,b):return[(max(lo,x),min(hi,y))for lo,hi in a for x,y in b if min(hi,y)-max(lo,x)>1e-9]
for release in ids:
 for ml in labels:
  for mr in labels:
   count+=1;rhs=np.zeros(19);rhs[-1]=-np.sqrt(2)*(np.cos(np.pi/ml)if ml else 1);rhs[ids.index(release)]=-np.sqrt(2*ns[release])*(np.cos(np.pi/mr)if mr else 1);r0=(inv@rhs)@R[J];q=2-r0@Q@r0
   if q<=1e-9:continue
   rat=-R@Q@r0/np.sqrt(2*ns)
   if not np.all((rat[fixed]>=1-1e-7)|(np.min(abs(rat[fixed,None]-allowed),axis=1)<1e-7)):continue
   UV=np.sqrt(q)*N;Puv=R@Q@UV;pr=R@Q@r0;conds=[(Puv[i,0],Puv[i,1],pr[i]+np.sqrt(2*ns[i]))for i in range(176)if not fixed[i]]
   hp=h@Q@UV;conds.append((hp[0],hp[1],h@Q@r0));I=[(0.,2*math.pi)]
   for A,B,C in conds:
    I=intersect(I,intervals(A,B,C))
    if not I:break
   if not I:continue
   for targetid,x in enumerate(targets):
    xp=x@Q@UV;TI=intersect(I,intervals(-xp[0],-xp[1],-x@Q@r0))
    for lo,hi in TI:
     th=(lo+hi)/2;r=r0+UV@np.array([np.cos(th),np.sin(th)]);key=tuple(np.round(r,8))
     if key in seen:continue
     seen.add(key);entry={'released':release,'facet_label':ml,'released_label':mr,'parameter_interval':[lo,hi],'parameter':th,'normal':r.tolist(),'target':targetid,'cut':float(r@Q@x),'norm':float(r@Q@r)};caps.append(entry);print('cap',len(caps),release,ml,mr,targetid,lo,hi,flush=True)
   families.append({'released':release,'facet_label':ml,'released_label':mr,'intervals':I})
print('profiles',count,'dotted families',len(families),'cutters',len(caps),flush=True)
(ROOT/'data/enumerations/fresh_h20_facet_bridge_20261013_continuous_scout.json').write_text(json.dumps({'status':'NUMERICAL_CONTINUOUS_DOTTED_ANCHOR_CIRCLES','profiles':count,'labels':labels,'families':families,'caps':caps,'targets':[x.tolist()for x in targets]},indent=2)+'\n')
