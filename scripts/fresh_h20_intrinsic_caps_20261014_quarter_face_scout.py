#!/usr/bin/env python3
"""Numerical quarter-face cap covering search; no volume or exact claims."""
import argparse
from itertools import combinations,product
import json,math,time,os
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix,vstack
from fresh_h20_intrinsic_caps_20261014 import ROOT,TAG,data,sha

OUT=ROOT/f'results/{TAG}_quarter_face_scout.json'
R=1.5
def allowed(co,phi,tol=1e-8):
 return (co>=1-tol)|(abs(co)<tol)|(abs(co-.5)<tol)|(abs(co-phi/2)<tol)

def build(mode):
 phi=(1+math.sqrt(5))/2;D=data(phi);V=np.array([v['point']for v in D['vertices']]);X=V/4
 starts=[0,12,15];dims=[slice(0,11),slice(11,13),slice(13,19)]
 factorV=[]
 for f,(st,n,sl)in enumerate(zip(starts,[12,3,7],dims)):
  factorV.append(np.array([next(v['point'][sl]for v in D['vertices']if v['omitted'][f]==st+i)for i in range(n)]))
 U=D['unit_Gram'][:12,:12];adj=[set(map(int,np.flatnonzero(U[i]<-.1)))for i in range(12)]
 assert all(len(a)==2 for a in adj)
 cycle=[0,min(adj[0])]
 while len(cycle)<12:cycle.append(next(i for i in adj[cycle[-1]]if i!=cycle[-2]))
 assert len(set(cycle))==12 and cycle[0]in adj[cycle[-1]]
 if mode=='contiguous':AA=sorted({tuple(sorted(cycle[(i+j)%12]for j in range(k)))for k in range(1,5)for i in range(12)})
 else:AA=[a for k in range(1,5)for a in combinations(range(12),k)]
 BB=[a for k in range(1,4)for a in combinations(range(3),k)]
 MM=np.array(D['marks'][2],dtype=int);EE=[a for k in range(1,8)for a in combinations(range(7),k)if sum(MM[list(a)])<=4]
 FC=[np.array([np.average(factorV[f][list(a)],axis=0,weights=np.array(D['marks'][f])[list(a)])/4 for a in choices])for f,choices in enumerate([AA,BB,EE])]
 choices=list(product(range(len(AA)),range(len(BB)),range(len(EE))))
 C=np.array([np.concatenate([FC[f][a]for f,a in enumerate(key)])for key in choices])
 center2=np.sum(C*C,axis=1);x2=np.sum(X*X,axis=1)
 dist=np.maximum(center2[:,None]+x2[None,:]-2*C@X.T,0)
 cover=dist<=R*R+1e-8;use=cover.any(axis=1);before=len(C);covering_count=int(use.sum())
 C=C[use];dist=dist[use];cover=cover[use];choices=[a for a,keep in zip(choices,use)if keep]
 # Check every retained seed and all505 existing spherical caps.  These are
 # numerical gates; an eventual positive selection requires exact replay.
 seed=(D['offsets'][None,:]-C@D['vectors'].T)/R
 seed_ok=allowed(seed,phi).all(axis=1)
 oldC=np.concatenate([np.zeros((1,19)),V,V/2]);oldR=np.concatenate([[2.],np.ones(504)])
 compat=np.ones(len(C),dtype=bool);minimum=math.inf
 for a in range(0,len(C),2000):
  cc=C[a:a+2000];dd=np.maximum(np.sum(cc*cc,axis=1)[:,None]+np.sum(oldC*oldC,axis=1)[None,:]-2*cc@oldC.T,0)
  co=(dd-R*R-oldR[None,:]**2)/(2*R*oldR[None,:]);compat[a:a+len(cc)]=allowed(co,phi).all(axis=1);minimum=min(minimum,float(co.min()))
 use=seed_ok&compat;C=C[use];dist=dist[use];cover=cover[use];choices=[a for a,keep in zip(choices,use)if keep]
 recipes=[dict(A11=list(AA[a]),triangle=list(BB[b]),E6=list(EE[e]))for a,b,e in choices]
 return D,V,X,C,cover,dist,recipes,dict(A11_mode=mode,A11_cycle=cycle,A11_subsets=len(AA),triangle_subsets=len(BB),E6_subsets=len(EE),raw_caps=before,caps_covering_probe=covering_count,seed_rejected=int((~seed_ok).sum()),existing527_rejected=int((~compat).sum()),minimum_old_spherical_cosine=minimum)

def run(mode,seconds,rounds):
 t=time.time();phi=(1+math.sqrt(5))/2;D,V,X,C,cover,dist,recipes,counts=build(mode);n=len(C)
 print('PID',os.getpid(),'bank',n,counts,'coverage',np.unique(cover.sum(axis=0),return_counts=True),flush=True)
 assert n and cover.any(axis=0).all()
 ri,ci=np.where(cover.T);base=coo_matrix((np.ones(len(ri)),(ri,ci)),shape=(252,n)).tocsr()
 conflicts=set();history=[];accepted=None;last=[]
 # Integer coverage optimization with lazily generated *necessary* pair
 # conflicts; no conflict is inferred from geometry of an unselected pair.
 for it in range(rounds):
  if conflicts:
   rows=[];cols=[]
   for k,(a,b)in enumerate(sorted(conflicts)):rows.extend([k,k]);cols.extend([a,b])
   neg=coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(conflicts),n)).tocsr();A=vstack([base,neg]).tocsr()
  else:A=base
  lo=np.concatenate([np.ones(252),np.full(len(conflicts),-np.inf)]);hi=np.concatenate([np.full(252,np.inf),np.ones(len(conflicts))])
  res=milp(np.ones(n),integrality=np.ones(n),bounds=Bounds(0,1),constraints=LinearConstraint(A,lo,hi),options={'time_limit':seconds,'mip_rel_gap':0.01})
  row=dict(round=it,status=int(res.status),message=res.message,constraints=len(conflicts),elapsed_seconds=time.time()-t)
  if res.x is None:history.append(row);print(row,flush=True);break
  chosen=np.flatnonzero(res.x>.5);last=chosen.tolist();assert cover[chosen].any(axis=0).all()
  cc=C[chosen];dd=np.maximum(np.sum(cc*cc,axis=1)[:,None]+np.sum(cc*cc,axis=1)[None,:]-2*cc@cc.T,0);co=dd/(2*R*R)-1
  bad=np.argwhere(np.triu(~allowed(co,phi),1));new={(min(int(chosen[a]),int(chosen[b])),max(int(chosen[a]),int(chosen[b])))for a,b in bad}
  assert not new&conflicts
  row.update(selected=len(chosen),new_conflicts=len(new));history.append(row);print(row,flush=True)
  if not new:accepted=last;break
  conflicts|=new
 status='NUMERICAL_COMPATIBLE_PROBE_COVER_FOUND'if accepted is not None else('NUMERICAL_FINITE_BANK_INFEASIBLE'if history[-1]['status']==2 else'NUMERICAL_LAZY_SEARCH_INCOMPLETE')
 paths=[f'scripts/{TAG}.py',f'scripts/{TAG}_quarter_face_scout.py',f'results/{TAG}_golden_seed.json',f'results/{TAG}_golden_extension.json']
 out=dict(status=status,sources={p:sha(p)for p in paths},counts=counts,candidate_count=n,radius='3/2',center_scale='1/4',probe_count=252,coverage_per_probe=cover.sum(axis=0).tolist(),history=history,conflict_count=len(conflicts),selected=accepted,last_relaxed_selection=last,selected_recipes=[recipes[i]for i in accepted]if accepted is not None else[],candidate_recipes=recipes,scope='Numerical face-center bank, floating Coxeter gates and finite252-probe cover. A selected set needs exact arithmetic and full volume checks. Infeasibility concerns only this bank; time/round limits are explicitly incomplete.')
 target=OUT if mode=='contiguous'else OUT.with_name(OUT.stem+'_all_subsets.json');target.write_text(json.dumps(out,indent=2)+'\n')
 np.savez_compressed(ROOT/f'data/enumerations/{TAG}_quarter_face_{mode}.npz',centers=C,cover=cover,probe_distances=dist)
 print(status,'output',str(target),flush=True)
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['contiguous','all'],default='contiguous');p.add_argument('--seconds',type=float,default=30);p.add_argument('--rounds',type=int,default=30);a=p.parse_args();run(a.mode,a.seconds,a.rounds)
