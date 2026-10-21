#!/usr/bin/env python3
"""Adaptive finite-profile cuts of the current remote H5 vertex."""
import argparse,json,sys
from itertools import product,combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.linalg import null_space
from scipy.spatial import HalfspaceIntersection
from scripts import fresh_h9_higher_face_20261013_linked_scout as L
S=L.S

def vertices(R,G,c,fixed):
 F=R[fixed];p=c-(c@G@F.T)@np.linalg.inv(F@G@F.T)@F;p/=np.sqrt(-p@G@p);B=null_space(np.vstack([F@G,p@G]));M=B.T@G@B;B=B@np.linalg.inv(np.linalg.cholesky(M).T);A=R@G@B;b=R@G@p;keep=np.linalg.norm(A,axis=1)>1e-7;H=np.c_[A[keep],b[keep]];vv=HalfspaceIntersection(H,np.zeros(5),qhull_options='Qx').intersections;xs=p+vv@B.T;qs=(vv*vv).sum(axis=1)-1;order=np.argsort(-qs);return[(float(qs[i]),xs[i],[j for j,x in enumerate(R@G@xs[i])if abs(x)<1e-6])for i in order]
def noaffine(R,G,r):
 ps=R@G@r;ii=np.where(ps>-2+1e-6)[0];A=np.vstack([r,R[ii]]);M=A@G@A.T
 for mask in range(1<<len(ii)):
  ids=[0]+[i+1 for i in range(len(ii))if mask>>i&1];m=M[np.ix_(ids,ids)]
  if np.any(m<-2-1e-7):continue
  eig=np.linalg.eigvalsh(m)
  if eig[0]>-1e-7 and eig[0]<1e-7:return False
 return True

def run(steps):
 _,_,_,_,_,names,_,_,R,G,c=S.load();d=S.V.read(S.PREFIX+'_selected_scout.json');R=np.vstack([R,np.array([x['normal']for x in d['recipes']])]);names+=['roof'+str(i)for i in range(7)];aa=S.V.read(S.PREFIX+'_adaptive_scout.json')['steps'];tt=S.V.read(S.PREFIX+'_truncation_select.json')['bank'];extras=[x['normal']for x in aa]+[tt[i]['normal']for i in(0,11)];R=np.vstack([R,np.array(extras)]);names+=['adaptive0','adaptive1','roofPolar','roofMetric'];fixed=[names.index(n)for n in S.COMMON if n!='seed1'];labels=np.array(list(product((2,3,4,5,6,8),repeat=5)));bb=-2*np.cos(np.pi/labels);records=[];vv=vertices(R,G,c,fixed)
 for step in range(steps):
  proposals=[]
  for vi,(q,x,inc)in enumerate(vv[:8]):
   if q<1e-7:break
   other=[i for i in inc if i not in fixed]
   if len(other)!=5:continue
   E=R[fixed+other];M=E@G@E.T
   if abs(np.linalg.det(M))<1e-9:continue
   Mi=np.linalg.inv(M);bs=np.c_[np.zeros((len(bb),4)),bb];co=bs@Mi;part=co@E;qp=np.einsum('bi,bi->b',bs,co);idx=np.where(qp<2-1e-8)[0];rr=part[idx]+np.sqrt(2-qp[idx])[:,None]*(x/np.sqrt(q));ps=rr@G@R.T;ms=np.maximum(2,np.rint(np.pi/np.arccos(np.clip(-ps/2,-1+1e-15,1-1e-15))));finite=(ps<=1e-7)&(ps>=-2)&(np.abs(ps+2*np.cos(np.pi/ms))<1e-7);valid=np.all((ps<-2-1e-7)|(np.abs(ps)<1e-7)|finite,axis=1)&(rr@G@c<0)&np.all(np.abs(ps+2)>1e-7,axis=1)
   order=sorted(np.where(valid)[0],key=lambda k:(sum(labels[idx[k]]!=2),max(labels[idx[k]]),sum(labels[idx[k]])))[:50]
   for k in order:
    r=rr[k]
    if not noaffine(R,G,r):continue
    cuts=np.array([r@G@z[1]for z in vv]);score=sum(max(0,z[0])for z,t in zip(vv,cuts)if t>1e-6)
    proposals.append({'right_fixed':[names[j]for j in fixed],'incident':[names[j]for j in other],'profile':labels[idx[k]].tolist(),'normal':r.tolist(),'target':x.tolist(),'target_norm':q,'source_vertex':vi,'cut_score':float(score)})
  proposals.sort(key=lambda x:-x['cut_score']);cases=[]
  for p in proposals[:30]:
   Rn=np.vstack([R,np.array(p['normal'])]);newvv=vertices(Rn,G,c,fixed);cases.append((newvv[0][0],p,newvv))
  if not cases:print('STOP no capped profile',step,flush=True);break
  cases.sort(key=lambda x:x[0]);maxq,best,newvv=cases[0];best['new_name']='secondAdaptive'+str(step);best['remaining_max_norm']=maxq;best['projective_vertices']=len(newvv);best['finite_vertices']=sum(q< -1e-7 for q,_,_ in newvv);best['candidate_count']=len(proposals);best['geometry_trials']=len(cases);records.append(best);R=np.vstack([R,np.array(best['normal'])]);names.append(best['new_name']);vv=newvv;print(step,best['incident'],best['profile'],maxq,'finite',best['finite_vertices'],flush=True)
  out={'status':'NUMERICAL_ADAPTIVE_H5_PROFILE_TRUNCATIONS','source':S.PREFIX+'_selected_scout.json','steps':records,'remaining_max_norm':maxq,'remaining_point':newvv[0][1].tolist(),'remaining_incident':[names[i]for i in newvv[0][2]],'scope':'Numerical finite-profile and affine gates. No exact or exhaustive completion claim.'};(S.V.ROOT/(S.PREFIX+'_adaptive_next_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
  if maxq< -1e-7:print('NUMERICALLY COMPACT H5',flush=True);break
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--steps',type=int,default=8);a=p.parse_args();run(a.steps)
