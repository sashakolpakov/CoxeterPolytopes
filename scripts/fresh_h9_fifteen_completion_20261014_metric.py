#!/usr/bin/env python3
"""Numerical rank-ten Lorentzian fit, only after a full-category survivor.

The six-dimensional Schur complement must be negative rank one. Continuous
category>=7 parameters are a necessary relaxation, not discrete label proofs.
"""
import argparse,json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from fresh_h9_fifteen_facet_20261014_incidence import expand_case
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_fifteen_completion_20261014'
PAIRS=list(combinations(range(15),2));PI={p:i for i,p in enumerate(PAIRS)}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def run(case,evaluations):
 cp=f'results/{TAG}_cases/case_{case:05d}.json';doc=load(cp);r=doc['result'];assert r['verdict']=='FULL_CATEGORY_METRIC_REQUIRED'
 ip='results/fresh_h9_fifteen_facet_20261014_incidence.json';sp='results/fresh_h9_gale_frontier_20261014_geometry.json'
 inc=load(ip);source=load(sp);row=expand_case(inc['cases'][case],source)
 labels=r['label_vector'];dotted={tuple(p)for p in r['dotted_pairs']}
 variables=[(p,'dotted')for p in PAIRS if p in dotted]+[(p,'finite_ge7')for p,m in zip(PAIRS,labels)if p not in dotted and m==7]
 G=np.eye(15)
 for(i,j),m in zip(PAIRS,labels):G[i,j]=G[j,i]=0. if m==2 else-np.cos(np.pi/m)
 for(i,j),kind in variables:G[i,j]=G[j,i]=-1.2 if kind=='dotted'else-np.cos(np.pi/7)
 vertices=[tuple(i for i in range(15)if v>>i&1)for v in row['canonical_vertices']]
 vertex=max(vertices,key=lambda v:np.linalg.det(G[np.ix_(v,v)]));rest=[i for i in range(15)if i not in vertex]
 assert len(rest)==6
 lo=np.array([1+1e-10 if k=='dotted'else np.cos(np.pi/7)for p,k in variables]);hi=np.array([100 if k=='dotted'else 1-1e-10 for p,k in variables])
 def schur(x):
  for((i,j),kind),z in zip(variables,x):G[i,j]=G[j,i]=-z
  A=G[np.ix_(vertex,vertex)];B=G[np.ix_(vertex,rest)]
  return G[np.ix_(rest,rest)]-B.T@np.linalg.solve(A,B)
 def residual(x):
  S=schur(x);z=np.sqrt(np.maximum(-S.diagonal(),0))
  return np.r_[(S+np.outer(z,z))[np.triu_indices(6,1)],np.maximum(S.diagonal(),0)]
 results=[]
 for radius in([1.05,1.2,2.,5.,15.,50.]if variables else[0.]):
  x=np.array([radius if kind=='dotted'else .95 for p,kind in variables])
  if variables:
   fit=least_squares(residual,x,bounds=(lo,hi),max_nfev=evaluations,ftol=1e-12,xtol=1e-12,gtol=1e-12);x=fit.x;neval=int(fit.nfev)
  else:neval=1
  rr=residual(x);S=schur(x)
  results.append(dict(start=radius,residual=float(np.max(np.abs(rr))),parameters=x.tolist(),evaluations=neval,
   gram_eigenvalues=np.linalg.eigvalsh(G).tolist(),schur_diagonal=S.diagonal().tolist(),minimum_vertex_eigenvalue=min(float(np.linalg.eigvalsh(G[np.ix_(v,v)])[0])for v in vertices)))
  if results[-1]['residual']<1e-9 and max(S.diagonal())<0:break
 files=[cp,ip,sp,f'scripts/{TAG}_metric.py','scripts/fresh_h9_fifteen_facet_20261014_incidence.py']
 return dict(status='NUMERICAL_NECESSARY_LORENTZIAN_METRIC_SCOUT',sources={p:sha(p)for p in files},case=case,
  fixed_vertex=list(vertex),variables=[dict(pair=p,kind=k)for p,k in variables],labels=labels,dotted_upper_bound=100,max_evaluations=evaluations,
  fits=results,best_residual=min(r['residual']for r in results),scope='Numerical only. No failed fit excludes a labeling, and any successful fit requires exact rank/signature, integer labels and complete compact-volume certification.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,required=True);ap.add_argument('--evaluations',type=int,default=2000);a=ap.parse_args();d=run(a.case,a.evaluations)
 (ROOT/f'results/{TAG}_metric_{a.case}.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['case'],d['best_residual'],flush=True)
