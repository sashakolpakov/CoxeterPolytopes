#!/usr/bin/env python3
"""Duplicate a finite neighbor through a reflection fixing the last H20 facet.

The other17 finite-neighbor pairings are inherited verbatim, including dotted
and nonright relations. This is different from the old all-right bridge plane.
One self-image Coxeter label leaves an affine circle; a bounded-angle anchor
then yields explicit proposals. Numerical only.
"""
import json,numpy as np
from scipy.linalg import null_space
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT

OUT='data/enumerations/fresh_h20_changed_assembly_20261013_neighbor_scout.json'
Q,R,h,neighbors,*_=setup();qn=np.einsum('ij,jk,ik->i',R,Q,R);cov=R@Q
I=neighbors+[175];M=cov[I];G=R[I]@Q@R[I].T;Gi=np.linalg.inv(G)
N=null_space(M);ev,E=np.linalg.eigh(N.T@Q@N);assert min(ev)>0;N=N@E@np.diag(1/np.sqrt(ev))
projection=cov@N;hproj=h@Q@N
targets=[np.array(x)for x in json.load(open(ROOT/'data/enumerations/fresh_h20_facet_bridge_20261013_release_scout.json'))['targets']]
labels=[2,3,4,5,6,8,10,12,0];allowed=np.array([0.]+[np.cos(np.pi/m)for m in range(3,121)])
circles=[];candidates=[];near=[];seen=set()
for j in neighbors:
 for label in labels:
  cosine=1 if label==0 else np.cos(np.pi/label)
  rhs=R[I]@Q@R[j];rhs[I.index(j)]=-qn[j]*cosine
  center=R[I].T@Gi@rhs;radius2=qn[j]-center@Q@center
  row={'neighbor':j,'self_image_label':'parallel'if label==0 else label,'radius_squared':float(radius2)}
  if radius2<-1e-7:row['status']='empty real circle';circles.append(row);continue
  if radius2<1e-9:row['status']='degenerate circle';circles.append(row);continue
  radius=np.sqrt(radius2);base=cov@center;amp=radius*np.linalg.norm(projection,axis=1);den=np.sqrt(qn*qn[j])
  zlo=-(base+amp)/den;zhi=-(base-amp)/den
  choices=[]
  for k in range(176):
   if k in I or amp[k]<1e-8:continue
   if zhi[k]<1-1e-6:
    mmax=int(np.floor(np.pi/np.arccos(max(0,min(1,zhi[k])))))+1
    possible=[m for m in range(2,mmax+1)if zlo[k]-1e-8<=np.cos(np.pi/m)<=zhi[k]+1e-8]
    choices.append((len(possible),k,possible))
  if not choices:row['status']='no discrete bounded anchor';circles.append(row);continue
  _,k,ks=min(choices);row.update(status='discrete anchor',anchor=k,anchor_labels=ks,anchor_cosine_interval=[float(zlo[k]),float(zhi[k])])
  circles.append(row)
  a=projection[k];an=np.linalg.norm(a);unit=a/an;perp=np.array([-unit[1],unit[0]])
  for m in ks:
   wanted=-np.sqrt(qn[j]*qn[k])*np.cos(np.pi/m);v=(wanted-base[k])/(radius*an)
   if abs(v)>1+1e-7:continue
   for sg in(-1,1):
    theta=v*unit+sg*np.sqrt(max(0,1-v*v))*perp;r=center+radius*N@theta
    if abs(r@Q@r-qn[j])>1e-5 or r@Q@h>=-1e-7:continue
    pairs=cov@r;z=-pairs/den
    ok=(z>=1-1e-7)|(np.min(abs(z[:,None]-allowed),axis=1)<1e-7)
    bad=np.flatnonzero(~ok).tolist();cuts=[float(r@Q@x)for x in targets]
    rec={'neighbor':j,'self_image_label':row['self_image_label'],'anchor':k,'anchor_label':m,'branch':sg,
         'normal':r.tolist(),'norm':float(qn[j]),'target_pairings':cuts,'bad':bad,'bad_cosines':[float(z[x])for x in bad]}
    if len(bad)<=3:near.append(rec)
    if not bad:
     key=tuple(np.round(r/np.sqrt(qn[j]),9))
     if key not in seen:seen.add(key);candidates.append(rec);print('candidate',j,label,k,m,sg,'cuts',cuts,flush=True)
 print('neighbor',j,'candidates',len(candidates),'near',len(near),flush=True)
out={'status':'NUMERICAL_REFLECTED_FINITE_NEIGHBOR_CIRCLE_SCOUT','source':'data/enumerations/fresh_h20_coupled_completion_20261013_second_exact.json',
     'neighbor_indices':neighbors,'self_image_labels':labels,'circles':circles,'compatible_candidates':candidates,'near_misses':near,
     'scope':'Preserve the original chosen neighbor row against the other17 finite neighbors and facet175; solve one self-image label and one finite-anchor equation. Floating proposals only.'}
(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
