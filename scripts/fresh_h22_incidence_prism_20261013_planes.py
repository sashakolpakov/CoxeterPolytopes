#!/usr/bin/env python3
"""Construct cap proposals from21 right angles and one prescribed angle.

Works in a numerically orthonormal Lorentz frame. Scope is selected remote
vertices/anchors; no completeness or exact-Coxeter assertion is made here.
"""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import null_space,qr
from fresh_h22_borcherds_extension_20261013 import prepare
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_incidence_prism_20261013'
def run():
 Q,W,G,h,X,A,b,Qn,_=prepare();J=np.diag([-1.]+[1.]*22)
 data=json.loads((P/'data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json').read_text())['candidates'];chosen=[c for c in data if c['twice_height_square']==4];C=np.array([c['twice_old_projection']for c in chosen],dtype=float)/2
 N=np.vstack([np.column_stack([W,np.zeros(len(W))]),np.column_stack([C,np.full(len(C),np.sqrt(2))]),np.column_stack([C,np.full(len(C),-np.sqrt(2))])]);M=np.zeros((23,23));M[:22,:22]=Q;M[-1,-1]=1
 B=np.zeros((23,23));B[:22,0]=h/np.sqrt(14);B[:22,1:22]=X;B[-1,-1]=1
 R=N@M@B@J;wn=np.einsum('ij,ij->i',R,R@J);assert max(abs(wn-np.r_[np.einsum('ij,ij->i',W,W@Q),np.full(602,2)]))<1e-7
 old=json.loads((P/('data/enumerations/'+PRE+'_scout.json')).read_text());runs=sorted([r for r in old['runs']if r.get('polar_invalid_count',0)>0],key=lambda r:r['polar_invalid_count'])[:40]
 candidates=[];seen=set();cosines=np.cos(np.pi/np.arange(2,121));labels=[2,3,4,5,6,7,8,10,12,16,20,24]
 for rr in runs:
  active=rr['active'];_,_,piv=qr(R[active].T,mode='economic',pivoting=True);base=[active[int(i)]for i in piv[:22]];point=J@B.T@M@np.array(rr['point']);new=0
  for released in base:
   fixed=[i for i in base if i!=released];E=null_space(R[fixed]@J);assert E.shape==(23,2);F=E.T@J@E
   for anchor in [released]+rr['polar_invalid_contacts']:
    a=E.T@J@R[anchor];an=float(a@a)
    if an<1e-16:continue
    v=np.array([-a[1],a[0]])/np.sqrt(an)
    for m in labels:
     p=-np.sqrt(2*wn[anchor])*np.cos(np.pi/m);y0=p*a/an;coeff=[float(v@F@v),float(2*v@F@y0),float(y0@F@y0-2)]
     if abs(coeff[0])<1e-10:
      roots=[]if abs(coeff[1])<1e-10 else[-coeff[2]/coeff[1]]
     else:
      dis=coeff[1]**2-4*coeff[0]*coeff[2]
      if dis<-1e-9:continue
      roots=[(-coeff[1]+s*np.sqrt(max(dis,0)))/(2*coeff[0])for s in(-1,1)]
     for t in roots:
      r=E@(y0+t*v)
      if r[0]<1e-8 or r@J@point<1e-7:continue
      c=-(R@J@r)/np.sqrt(2*wn)
      if min(c)<-1e-6:continue
      finite=c<1-1e-6
      if np.any(np.min(abs(c[finite,None]-cosines),axis=1)>1e-6):continue
      key=tuple(np.round(r,7))
      if key in seen:continue
      seen.add(key);new+=1;candidates.append({'vertex_trial':rr['trial'],'fixed_right':fixed,'anchor':anchor,'label':m,'canonical_normal':r.tolist(),'physical_normal':(B@r).tolist(),'norm':float(r@J@r),'escape_pair':float(r@J@point),'finite_contacts':[{'wall':int(i),'label':int(np.argmin(abs(c[i]-cosines)))+2}for i in np.flatnonzero(finite)],'minimum_dotted_excess':float(min(c[~finite]-1))if np.any(~finite)else None})
  print('vertex',rr['trial'],'bad',rr['polar_invalid_count'],'new',new,'total',len(candidates),flush=True)
 out={'status':'NUMERICAL_21_RIGHT_ONE_ANGLE_CAP_PROPOSALS','vertices_checked':len(runs),'angle_labels':labels,'candidates':candidates,'scope':'All812 numerical pair gates pass, and each proposal cuts its selected remote vertex. Exact reconstruction, mutual compatibility and volume remain required.'};(P/('data/enumerations/'+PRE+'_planes.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(candidates),flush=True)
if __name__=='__main__':run()
