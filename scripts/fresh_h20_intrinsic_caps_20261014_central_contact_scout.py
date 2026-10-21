#!/usr/bin/env python3
"""Finite central-ball contacts along mark-weighted face-center rays.

Floating proposal scout.  Retained seed contacts stay3/5 via R=2(1-gamma).
Each prescribed central angle solves one explicit quadratic, not a gamma grid.
"""
import argparse,json,math,time
from itertools import product,combinations
from collections import Counter
import numpy as np
from fresh_h20_intrinsic_caps_20261014 import ROOT,TAG,data,sha
from fresh_h20_intrinsic_caps_20261014_quarter_face_scout import build

def allowed(x,cs,tol=1e-8):
 out=x>=1-tol
 for c in cs:out|=abs(x-c)<tol
 return out

def run(mode,maximum):
 t=time.time();phi=(1+5**.5)/2;D=data(phi);V=np.array([v['point']for v in D['vertices']])
 starts=[0,12,15];dims=[slice(0,11),slice(11,13),slice(13,19)];sizes=[12,3,7]
 FV=[np.array([next(v['point'][sl]for v in D['vertices']if v['omitted'][f]==st+i)for i in range(n)])for f,(st,n,sl)in enumerate(zip(starts,sizes,dims))]
 adj=[set(map(int,np.flatnonzero(D['unit_Gram'][i,:12]<-.1)))for i in range(12)];cy=[0,min(adj[0])]
 while len(cy)<12:cy.append(next(i for i in adj[cy[-1]]if i!=cy[-2]))
 AA=sorted({tuple(sorted(cy[(i+j)%12]for j in range(k)))for k in range(1,5)for i in range(12)})if mode=='contiguous'else[a for k in range(1,5)for a in combinations(range(12),k)]
 BB=[a for k in range(1,4)for a in combinations(range(3),k)];MM=np.array(D['marks'][2],dtype=int);EE=[a for k in range(1,8)for a in combinations(range(7),k)if sum(MM[list(a)])<=4]
 FC=[np.array([np.average(FV[f][list(a)],axis=0,weights=np.array(D['marks'][f])[list(a)])for a in choices])for f,choices in enumerate([AA,BB,EE])]
 choices=list(product(range(len(AA)),range(len(BB)),range(len(EE))));P=np.array([np.concatenate([FC[f][a]for f,a in enumerate(key)])for key in choices]);N=np.sum(P*P,axis=1)
 S=np.array([[len(AA[a]),len(BB[b]),sum(MM[list(EE[e])])]for a,b,e in choices]);assert min(N)>4
 labels=list(range(2,maximum+1));cos=np.cos(np.pi/np.array(labels));cs=list(cos)+[1.]
 quarter=json.loads((ROOT/f'results/{TAG}_quarter_face_scout.json').read_text());_,_,_,qC,_,_,_,_=build('contiguous');qC=qC[quarter['selected']]
 oldC=np.concatenate([V,V/2,qC]);oldR=np.concatenate([np.ones(504),np.full(len(qC),1.5)]);old2=np.sum(oldC*oldC,axis=1)
 counts=[];survivors=[];seed_total=0;best=None;all_retained=[]
 for label,c in list(zip(labels,cos))+[('parallel',1.)]:
  b=8*(1+c);a=N-4;gamma=2*b/(np.sqrt(b*b+4*a*b)+b);R=2*(1-gamma);T=gamma/(1-gamma)
  ratios=np.column_stack([.5+6*T/S[:,0],phi/2*(1+3*T/S[:,1]),.5+6*T/S[:,2]])
  ok=allowed(ratios,cs).all(axis=1);indices=np.flatnonzero(ok);seed_total+=len(indices)
  row=dict(central_label=label,seed_compatible=int(len(indices)),retains_all_vertex=0,retains_all_half=0,retains_all_quarter=0,retains_all=0)
  for start in range(0,len(indices),1000):
   ii=indices[start:start+1000];C=P[ii]*gamma[ii,None];RR=R[ii];dd=np.maximum(np.sum(C*C,axis=1)[:,None]+old2[None,:]-2*C@oldC.T,0);co=(dd-RR[:,None]**2-oldR[None,:]**2)/(2*RR[:,None]*oldR[None,:]);gate=allowed(co,cs)
   vpass=gate[:,:252].all(axis=1);hpass=gate[:,252:504].all(axis=1);qpass=gate[:,504:].all(axis=1);passall=gate.all(axis=1)
   row['retains_all_vertex']+=int(vpass.sum());row['retains_all_half']+=int(hpass.sum());row['retains_all_quarter']+=int(qpass.sum());row['retains_all']+=int(passall.sum())
   conflicts=(~gate).sum(axis=1)
   for j in np.flatnonzero(passall):
    k=int(ii[j]);a0,b0,e0=choices[k];r=dict(A11=list(AA[a0]),triangle=list(BB[b0]),E6=list(EE[e0]),central_label=label,gamma=float(gamma[k]),radius=float(R[k]),center=C[j].tolist(),seed_omitted_ratios=ratios[k].tolist());all_retained.append(r)
   j=int(np.argmin(conflicts))if len(ii)else None
   if j is not None and(best is None or int(conflicts[j])<best['conflicts']):
    k=int(ii[j]);a0,b0,e0=choices[k];best=dict(A11=list(AA[a0]),triangle=list(BB[b0]),E6=list(EE[e0]),central_label=label,gamma=float(gamma[k]),radius=float(R[k]),center=C[j].tolist(),seed_omitted_ratios=ratios[k].tolist(),conflicts=int(conflicts[j]),conflict_ids=np.flatnonzero(~gate[j]).tolist(),conflict_cosines=co[j,~gate[j]].tolist())
  counts.append(row)
  if label in[2,3,4,5,6,10,maximum,'parallel']:print(row,flush=True)
 out=dict(status='NUMERICAL_CENTRAL_CONTACT_FACE_RAY_SCOUT',sources={p:sha(p)for p in[f'scripts/{TAG}.py',f'scripts/{TAG}_quarter_face_scout.py',f'scripts/{TAG}_central_contact_scout.py',f'results/{TAG}_quarter_face_scout.json']},A11_mode=mode,face_rays=len(choices),central_labels=labels+['parallel'],cross_labels_tested=labels,continuous_radius_relation='R=2(1-gamma)',central_equation='(N-4)gamma^2+8(1+c)gamma-8(1+c)=0',seed_compatible=seed_total,counts=counts,retains_all_587=len(all_retained),all_retained_candidates=all_retained,best_replacement_candidate=best,elapsed_seconds=time.time()-t,scope='Floating finite-label central-contact equations on the displayed face rays. Each positive gamma is solved analytically. No full continuous/all-label exclusion or exact Coxeter/volume claim.')
 path=ROOT/f'results/{TAG}_central_contact_{mode}_scout.json';path.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],'seed',seed_total,'all',len(all_retained),'bestconflicts',best['conflicts']if best else None,'seconds',time.time()-t,flush=True)
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['contiguous','all'],default='contiguous');p.add_argument('--max-label',type=int,default=20);a=p.parse_args();run(a.mode,a.max_label)
