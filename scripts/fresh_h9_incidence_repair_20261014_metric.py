#!/usr/bin/env python3
"""Rank-minor-guided completion of changed 15-facet spherical nerves.

The finite contact choices are bounded; numerical vertex tests generate
proposals, while each archived nonzero modular minor is an exact rejection.
"""
from pathlib import Path
from itertools import combinations
import argparse,json,hashlib
import numpy as np
import z3
from scipy.optimize import least_squares
from fresh_h9_incidence_repair_20261014 import source,family,finite_labels,ROOT
SOURCE=ROOT/'results/fresh_h9_incidence_repair_20261014_released_sat.json'
OUT=ROOT/'results/fresh_h9_incidence_repair_20261014_metric.json'
P=241
SQ={d:next(x for x in range(1,P)if x*x%P==d)for d in[2,3,5]}
VALUES={1:2,2:0,3:P-1,4:P-SQ[2],5:(-(1+SQ[5])*pow(2,-1,P))%P,6:P-SQ[3]}

def detmod(L,B):
 A=[[VALUES[L[i][j]]for j in B]for i in B];d=1;n=len(A)
 for k in range(n):
  r=next((i for i in range(k,n)if A[i][k]),None)
  if r is None:return 0
  if r!=k:A[r],A[k]=A[k],A[r];d=-d
  q=A[k][k];d=d*q%P;iv=pow(q,-1,P)
  for i in range(k+1,n):
   x=A[i][k]*iv%P
   for j in range(k+1,n):A[i][j]=(A[i][j]-x*A[k][j])%P
 return d%P

def matrix(L,D,x=None):
 A=np.eye(15)
 for i,j in combinations(range(15),2):A[i,j]=A[j,i]=-np.cos(np.pi/L[i][j])if L[i][j]>2 else 0
 for k,(i,j)in enumerate(D):A[i,j]=A[j,i]=-(x[k]if x is not None else 1.2)
 return A

def complete(case,model,max_rounds):
 row,G,ps=source();H,F,_=family(row,G,ps[case['pattern_index']]);L=finite_labels(H)
 for i,j in F:L[i][j]=L[j][i]=2
 for i,j,m in case['added_finite_contacts']:L[i][j]=L[j][i]=m
 D=sorted(map(tuple,model['extra_dotted_pairs']));Vs=np.array(model['vertices']);base=matrix(L,D)
 # Additional finite contacts must preserve all proposed finite vertices.
 # This numerical stage is only a proposal restriction, not an exhaustive
 # finite-label theorem. All later modular rejections are exact.
 options={}
 for i,j in combinations(range(15),2):
  if L[i][j]!=2 or(i,j)in D:continue
  opts=[2]
  for m in range(3,7):
   A=base.copy();A[i,j]=A[j,i]=-np.cos(np.pi/m)
   if np.min(np.linalg.eigvalsh(A[Vs[:,:,None],Vs[:,None,:]])[:,0])>1e-9:opts.append(m)
  if len(opts)>1:options[i,j]=opts
 X={p:{m:z3.Bool(f'x_{p[0]}_{p[1]}_{m}')for m in ms}for p,ms in options.items()};s=z3.Tactic('sat').solver();s.set(timeout=20000)
 for xs in X.values():s.add(z3.PbEq([(x,1)for x in xs.values()],1))
 Bs=[B for B in combinations(range(15),11)if not any(i in B and j in B for i,j in D)]
 Bs.sort(key=lambda B:sum(i in B and j in B for i,j in options))
 rejects=[];pd_rejects=[];metric=[];status='round_limit'
 for roundno in range(max_rounds):
  out=s.check()
  if out!=z3.sat:status=str(out);break
  M=s.model();choice={p:next(m for m,x in xs.items()if z3.is_true(M.eval(x)))for p,xs in X.items()};Q=[r[:]for r in L]
  for(i,j),m in choice.items():Q[i][j]=Q[j][i]=m
  A=matrix(Q,D);eig=np.linalg.eigvalsh(A[Vs[:,:,None],Vs[:,None,:]]);bad=np.where(eig[:,0]<=1e-9)[0]
  if len(bad):
   B=list(Vs[bad[0]]);cl=[x for p,xs in X.items()if p[0]in B and p[1]in B for m,x in xs.items()if m<choice[p]]
   assert cl;s.add(z3.Or(cl));pd_rejects.append({'round':roundno,'vertex':list(map(int,B)),'choices':[[i,j,m]for(i,j),m in choice.items()if i in B and j in B]});continue
  witness=next(((B,d)for B in Bs if(d:=detmod(Q,B))),None)
  if witness:
   B,d=witness;cl=[z3.Not(X[p][m])for p,m in choice.items()if p[0]in B and p[1]in B]
   # This exact clause is necessary: keeping every label in B fixes its
   # nonzero rank11 minor, independently of the three dotted distances.
   s.add(z3.Or(cl));rejects.append({'round':roundno,'principal_nodes':list(B),'determinant_mod241':d,'choices':[[i,j,m]for(i,j),m in choice.items()if i in B and j in B]})
   continue
  best=None
  B=next((b for b in combinations(range(15),10)if not any(i in b and j in b for i,j in D)and np.sum(np.linalg.eigvalsh(A[np.ix_(b,b)])<-1e-7)==1 and min(abs(np.linalg.eigvalsh(A[np.ix_(b,b)])))>1e-7),None)
  if B is not None:
   T=[i for i in range(15)if i not in B];ii,jj=np.triu_indices(5)
   def fun(x):
    Z=matrix(Q,D,x);C=Z[np.ix_(B,B)];R=Z[np.ix_(T,B)];return (Z[np.ix_(T,T)]-R@np.linalg.solve(C,R.T))[ii,jj]
   for h in[1.01,1.2,2,5]:
    r=least_squares(fun,np.full(len(D),h),bounds=(1.0000001,100),max_nfev=300,ftol=1e-12,xtol=1e-12,gtol=1e-12);err=float(max(abs(fun(r.x))))
    if best is None or err<best['residual']:best={'residual':err,'distances':r.x.tolist(),'basis':list(B),'eigenvalues':np.linalg.eigvalsh(matrix(Q,D,r.x)).tolist()}
  metric.append({'round':roundno,'choices':[[i,j,m]for(i,j),m in choice.items()],'metric':best})
  if best and best['residual']<1e-8:status='NUMERICAL_METRIC_CANDIDATE';break
  s.add(z3.Or([z3.Not(X[p][m])for p,m in choice.items()]))
 return {'base_case':case,'base_model':model,'allowed_new_contacts':[[i,j,ms]for(i,j),ms in options.items()],'fixed_dotted_pairs':[list(p)for p in D],'status':status,'modular_rejections':rejects,'numerical_vertex_rejections':pd_rejects,'distance_trials':metric}

def main():
 p=argparse.ArgumentParser();p.add_argument('--cases',type=int,default=4);p.add_argument('--rounds',type=int,default=500);a=p.parse_args();src=json.loads(SOURCE.read_text());rows=[]
 for case in src['cases']:
  if not case['models']:continue
  rows.append(complete(case,case['models'][0],a.rounds));r=rows[-1];print('case',len(rows),r['status'],'rank rejects',len(r['modular_rejections']),'vertex rejects',len(r['numerical_vertex_rejections']),'metric',len(r['distance_trials']),flush=True)
  d={'status':'BOUNDED_CONTACT_AND_RANK_CONSTRAINT_GENERATION','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[SOURCE,Path(__file__).resolve(),ROOT/'scripts/fresh_h9_incidence_repair_20261014.py']},'cases':rows,'scope':'Existing216-vertex nerve retained; optional new finite contacts2..6 selected by numerical individual vertex gates; exact mod241 rank nogoods; numerical simultaneous vertex and dotted-distance gates. No family completeness or compact polytope claim from solver status.'}
  OUT.write_text(json.dumps(d,indent=2)+'\n')
  if len(rows)>=a.cases:break
if __name__=='__main__':main()
