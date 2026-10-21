#!/usr/bin/env python3
"""Bounded numerical Lorentzian Gram fitting after exact necessary labels."""
import argparse,json,time
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from fresh_h9_highbase_construction_20261014_cell import finite

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_surgery_20261014'
PAIRS=list(combinations(range(14),2));PI={p:i for i,p in enumerate(PAIRS)}
def read(s):return json.loads((ROOT/f'results/{TAG}_{s}.json').read_text())
def extensions(row,proposal):
 V=[tuple(i for i in range(14)if v>>i&1)for v in row['canonical_vertices']]
 vp=[tuple(PI[p]for p in combinations(v,2))for v in V]
 labels=proposal['label_vector'][:];options=[];forced=[]
 for p in proposal['unconstrained_meeting_pairs']:
  k=PI[tuple(p)];allowed=[]
  for m in range(2,8):
   labels[k]=m
   if all(finite(9,tuple(labels[j]for j in vv))for vv in vp):allowed.append(m)
  labels[k]=2
  if len(allowed)>1:options.append((p,allowed))
  else:assert allowed==[2];forced.append(p)
 out=[]
 for choice in product(*(a for p,a in options)):
  for (p,a),m in zip(options,choice):labels[PI[tuple(p)]]=m
  if all(finite(9,tuple(labels[j]for j in vv))for vv in vp):out.append(labels[:])
 return V,options,forced,out

def fit(labels,V,dotted,maxeval):
 variable=[(tuple(p),'dotted')for p in dotted]+[(p,'finite_ge7')for p,m in zip(PAIRS,labels)if m==7]
 G=np.eye(14)
 for (i,j),m in zip(PAIRS,labels):G[i,j]=G[j,i]=0 if m==2 else -np.cos(np.pi/m)
 for (i,j),kind in variable:G[i,j]=G[j,i]=-1.2 if kind=='dotted'else -np.cos(np.pi/7)
 vertex=max(V,key=lambda v:np.linalg.det(G[np.ix_(v,v)]));remaining=[i for i in range(14)if i not in vertex]
 lo=np.array([1+1e-9 if kind=='dotted'else np.cos(np.pi/7)for p,kind in variable]);hi=np.array([30 if kind=='dotted'else 1-1e-9 for p,kind in variable])
 def schur(x):
  for ((i,j),kind),a in zip(variable,x):G[i,j]=G[j,i]=-a
  A=G[np.ix_(vertex,vertex)];B=G[np.ix_(vertex,remaining)];C=G[np.ix_(remaining,remaining)]
  return C-B.T@np.linalg.solve(A,B)
 def fun(x):
  S=schur(x);z=np.sqrt(np.maximum(-S.diagonal(),0))
  return np.r_[(S+np.outer(z,z))[np.triu_indices(5,1)],np.maximum(S.diagonal(),0)]
 best=None
 for start in[1.15,2,5]:
  x=np.array([start if kind=='dotted'else .95 for p,kind in variable]);r=least_squares(fun,x,bounds=(lo,hi),max_nfev=maxeval,ftol=1e-11,xtol=1e-11,gtol=1e-11)
  score=float(np.max(np.abs(fun(r.x))))
  if best is None or score<best['residual']:
   S=schur(r.x);best=dict(residual=score,parameters=r.x.tolist(),evaluations=int(r.nfev),status=int(r.status),schur_diagonal=S.diagonal().tolist(),eigenvalues=np.linalg.eigvalsh(G).tolist())
  if score<1e-8:break
 best.update(vertex=list(vertex),variables=[dict(pair=list(p),kind=k)for p,k in variable],labels=labels)
 return best

def run(maxeval):
 inc=read('incidence');dfs=read('dfs_0_5');rows=[];t=time.time()
 for p in dfs['cases']:
  if 'label_vector'not in p:continue
  V,options,forced,ext=extensions(inc['cases'][p['case']],p)
  print('case',p['case'],'metric families',len(ext),flush=True)
  fits=[]
  for i,l in enumerate(ext):
   f=fit(l,V,p['dotted_pairs'],maxeval);fits.append(f)
   if i%20==0 or f['residual']<1e-7:print(' fit',i,'residual',f['residual'],'seconds',round(time.time()-t,2),flush=True)
  rows.append(dict(case=p['case'],free_pair_options=options,forced_right_pairs=forced,exact_finite_extensions=len(ext),fits=fits))
 paths=[f'scripts/{TAG}_metric_scout.py',f'results/{TAG}_incidence.json',f'results/{TAG}_dfs_0_5.json','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='BOUNDED_NUMERICAL_GRAM_SCOUT',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},max_evaluations_per_start=maxeval,dotted_upper_bound=30,starts=[1.15,2,5],cases=rows,scope='Only the three saved necessary labelings. Exact finite-vertex extension categories, followed by numerical rank-ten Lorentzian Schur fitting. Category>=7 is relaxed continuously. Failure is not an exact metric exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--evaluations',type=int,default=400);a=p.parse_args();o=run(a.evaluations);(ROOT/f'results/{TAG}_metric_scout.json').write_text(json.dumps(o,indent=2)+'\n');print('DONE',[(r['case'],min(f['residual']for f in r['fits']))for r in o['cases']])
