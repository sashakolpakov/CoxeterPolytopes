#!/usr/bin/env python3
"""Numerical joint replacement cover: quarter and central-contact face caps."""
import argparse,json,math,time,os
from itertools import product,combinations
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix,vstack
from fresh_h20_intrinsic_caps_20261014 import ROOT,TAG,data,sha
from fresh_h20_intrinsic_caps_20261014_quarter_face_scout import build
from fresh_h20_intrinsic_caps_20261014_central_contact_scout import allowed

def run(seconds,rounds):
 t=time.time();phi=(1+5**.5)/2;D,V,X,Cq,_,_,Rq,counts=build('contiguous')
 dims=[slice(0,11),slice(11,13),slice(13,19)];starts=[0,12,15];sizes=[12,3,7]
 FV=[np.array([next(v['point'][sl]for v in D['vertices']if v['omitted'][f]==st+i)for i in range(n)])for f,(st,n,sl)in enumerate(zip(starts,sizes,dims))]
 cycle=counts['A11_cycle'];AA=sorted({tuple(sorted(cycle[(i+j)%12]for j in range(k)))for k in range(1,5)for i in range(12)})
 BB=[a for k in range(1,4)for a in combinations(range(3),k)];MM=np.array(D['marks'][2],dtype=int);EE=[a for k in range(1,8)for a in combinations(range(7),k)if sum(MM[list(a)])<=4]
 FC=[np.array([np.average(FV[f][list(a)],axis=0,weights=np.array(D['marks'][f])[list(a)])for a in choices])for f,choices in enumerate([AA,BB,EE])]
 choices=list(product(range(len(AA)),range(len(BB)),range(len(EE))));P=np.array([np.concatenate([FC[f][a]for f,a in enumerate(key)])for key in choices]);N=np.sum(P*P,axis=1)
 S=np.array([[len(AA[a]),len(BB[b]),sum(MM[list(EE[e])])]for a,b,e in choices]);labels=list(range(2,65));cs=list(np.cos(np.pi/np.array(labels)))+[1.]
 centers=[*Cq];radii=[1.5]*len(Cq);recipes=[dict(r,kind='quarter',gamma=.25,radius=1.5)for r in Rq]
 seed_counts=[]
 for label in[2,3,5]:
  c=math.cos(math.pi/label);b=8*(1+c);a=N-4;ga=2*b/(np.sqrt(b*b+4*a*b)+b);RR=2*(1-ga);T=ga/(1-ga)
  ratios=np.column_stack([.5+6*T/S[:,0],phi/2*(1+3*T/S[:,1]),.5+6*T/S[:,2]])
  ii=np.flatnonzero(allowed(ratios,cs).all(axis=1));seed_counts.append([label,len(ii)])
  for k in ii:
   a0,b0,e0=choices[k];centers.append(P[k]*ga[k]);radii.append(RR[k]);recipes.append(dict(kind='central',central_label=label,A11=list(AA[a0]),triangle=list(BB[b0]),E6=list(EE[e0]),gamma=float(ga[k]),radius=float(RR[k])))
 C=np.array(centers);rr=np.array(radii);probes=np.concatenate([V/4,V/8]);dist=np.maximum(np.sum(C*C,axis=1)[:,None]+np.sum(probes*probes,axis=1)[None,:]-2*C@probes.T,0);cover=dist<=rr[:,None]**2+1e-8
 use=cover.any(axis=1);C=C[use];rr=rr[use];cover=cover[use];recipes=[r for r,k in zip(recipes,use)if k]
 oldC=np.concatenate([np.zeros((1,19)),V,V/2]);oldR=np.concatenate([[2.],np.ones(504)]);passold=np.ones(len(C),dtype=bool)
 for start in range(0,len(C),1000):
  cc=C[start:start+1000];rs=rr[start:start+1000];dd=np.maximum(np.sum(cc*cc,axis=1)[:,None]+np.sum(oldC*oldC,axis=1)[None,:]-2*cc@oldC.T,0)
  co=(dd-rs[:,None]**2-oldR[None,:]**2)/(2*rs[:,None]*oldR[None,:]);passold[start:start+len(cc)]=allowed(co,cs).all(axis=1)
 C=C[passold];rr=rr[passold];cover=cover[passold];recipes=[r for r,k in zip(recipes,passold)if k];n=len(C)
 print('PID',os.getpid(),'bank',n,'quarter',sum(r['kind']=='quarter'for r in recipes),'seedcentral',seed_counts,'probe coverage',np.unique(cover.sum(axis=0),return_counts=True),flush=True)
 assert cover.any(axis=0).all()
 row,col=np.where(cover.T);base=coo_matrix((np.ones(len(row)),(row,col)),shape=(504,n)).tocsr();conflicts=set();history=[];selected=None;last=[]
 for it in range(rounds):
  if conflicts:
   ij=sorted(conflicts);rows=np.repeat(np.arange(len(ij)),2);cols=np.array(ij).flatten();extra=coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(len(ij),n)).tocsr();A=vstack([base,extra]).tocsr()
  else:A=base
  lo=np.concatenate([np.ones(504),np.full(len(conflicts),-np.inf)]);hi=np.concatenate([np.full(504,np.inf),np.ones(len(conflicts))])
  res=milp(np.ones(n),integrality=np.ones(n),bounds=Bounds(0,1),constraints=LinearConstraint(A,lo,hi),options={'time_limit':seconds,'mip_rel_gap':.01})
  rec=dict(round=it,status=int(res.status),message=res.message,existing_conflicts=len(conflicts),elapsed_seconds=time.time()-t)
  if res.x is None:history.append(rec);print(rec,flush=True);break
  ids=np.flatnonzero(res.x>.5);last=ids.tolist();assert cover[ids].any(axis=0).all();cc=C[ids];rs=rr[ids]
  dd=np.maximum(np.sum(cc*cc,axis=1)[:,None]+np.sum(cc*cc,axis=1)[None,:]-2*cc@cc.T,0);co=(dd-rs[:,None]**2-rs[None,:]**2)/(2*rs[:,None]*rs[None,:]);bad=np.argwhere(np.triu(~allowed(co,cs),1));new={(min(int(ids[a]),int(ids[b])),max(int(ids[a]),int(ids[b])))for a,b in bad};assert not new&conflicts
  rec.update(selected=len(ids),central=sum(recipes[i]['kind']=='central'for i in ids),new_conflicts=len(new));history.append(rec);print(rec,flush=True)
  if not new:selected=last;break
  conflicts|=new
 out=dict(status='NUMERICAL_JOINT_PROBE_COVER_FOUND'if selected is not None else('NUMERICAL_JOINT_BANK_INFEASIBLE'if history[-1]['status']==2 else'NUMERICAL_JOINT_LAZY_SEARCH_INCOMPLETE'),sources={p:sha(p)for p in[f'scripts/{TAG}_joint_contact_scout.py',f'scripts/{TAG}_central_contact_scout.py',f'scripts/{TAG}_quarter_face_scout.py',f'results/{TAG}_golden_extension.json']},fixed_wall_count=527,central_labels=[2,3,5],cross_labels_tested=labels,quarter_bank_count=sum(r['kind']=='quarter'for r in recipes),candidate_count=n,probe_count=504,probe_types=['quarter_vertex','eighth_vertex'],candidate_recipes=recipes,history=history,conflict_count=len(conflicts),selected=selected,selected_recipes=[recipes[i]for i in selected]if selected is not None else[],last_relaxed_selection=last,scope='Numerical finite-bank search allowing quarter-layer replacement. Fixed527 retained. Covers only504 probes; exact pair/geometry and full volume remain necessary.')
 (ROOT/f'results/{TAG}_joint_contact_scout.json').write_text(json.dumps(out,indent=2)+'\n');np.savez_compressed(ROOT/f'data/enumerations/{TAG}_joint_contact_scout.npz',centers=C,radii=rr,cover=cover,probes=probes)
 print(out['status'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=20);p.add_argument('--rounds',type=int,default=50);a=p.parse_args();run(a.seconds,a.rounds)
