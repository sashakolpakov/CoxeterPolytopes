#!/usr/bin/env python3
"""Changed incidences: minimal affine covers split into right and dotted pairs.

All old parallel pairs become free dotted distances. Every selected cover
pair becomes either right or free dotted. At each nonsimple old cusp at
least one of its pairs must be dotted: otherwise all its old affine
components have only been weakened and the entire incident Gram is PD.
This is a fresh continuous-distance family, not a higher finite-label cutoff.
"""
import argparse,hashlib,json,time
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_h9_incidence_surgery_20261013.json'
OUT=ROOT/'results/fresh_h9_incidence_surgery_20261013_dotted.json'
def patterns(row):
 P=[tuple(p)for p in row['parallel_pairs']]
 for si,S0 in enumerate(row['finite_edge_cover']['minimum_supports']):
  S=[tuple(p)for p in S0];requirements=[]
  for c in row['ideal_vertices']:
   A=c['active']
   if any(i in A and j in A for i,j in P):continue
   requirements.append(sum(1<<k for k,(i,j)in enumerate(S)if i in A and j in A))
  masks=[m for m in range(1<<len(S))if all(m&r for r in requirements)]
  assert masks
  # Minimal number of newly dotted pairs necessary at the old cusps.
  d=min(m.bit_count()for m in masks)
  for m in masks:
   if m.bit_count()!=d:continue
   yield si,m,S,P

def search(row):
 t=time.time();name=row['id'];G=np.array(row['source_positive_gram'],float);G/=np.sqrt(np.diag(G)[:,None]*np.diag(G)[None,:]);n=len(G)
 counters={'patterns':0,'fixed_PD_or_affine_principal_gate':0,'metric_trials':0,'metric_solutions':0,'unresolved':0};witness=[];trials=[];survivors=[]
 for si,mask,S,P in patterns(row):
  counters['patterns']+=1;H=G.copy();free=set(P)
  for k,(i,j)in enumerate(S):
   if mask>>k&1:free.add((i,j))
   else:H[i,j]=H[j,i]=0
  free=sorted(free);bad=None;bases=[]
  # Only principal Grams avoiding every free pair are used for this gate.
  for B in combinations(range(n),10):
   if any(i in B and j in B for i,j in free):continue
   ev=np.linalg.eigvalsh(H[np.ix_(B,B)])
   if min(ev)>-1e-10:bad=B;break
   if min(abs(ev))>1e-8 and sum(ev<0)==1:bases.append(B)
  if bad is not None:
   counters['fixed_PD_or_affine_principal_gate']+=1;witness.append({'support_index':si,'new_dotted_mask':mask,'principal_nodes':list(bad),'principal_eigenvalues':np.linalg.eigvalsh(H[np.ix_(bad,bad)]).tolist()});continue
  if bases:B=list(bases[0])
  else:
   # A nonsingular fixed basis is convenient, not logically required.
   # Use the archived basis and solve its variable Schur equations.
   H0=H.copy()
   for i,j in free:H0[i,j]=H0[j,i]=-1.2
   best=None
   for C in combinations(range(n),10):
    ev=np.linalg.eigvalsh(H0[np.ix_(C,C)])
    if sum(ev<0)==1 and min(abs(ev))>1e-6:
     score=min(abs(ev));
     if best is None or score>best[0]:best=(score,C)
   if best is None:counters['unresolved']+=1;continue
   B=list(best[1])
  T=[i for i in range(n)if i not in B];ii,jj=np.triu_indices(len(T))
  def matrix(x):
   A=H.copy()
   for (i,j),v in zip(free,x):A[i,j]=A[j,i]=-v
   return A
  def residual(x):
   A=matrix(x)
   try:inv=np.linalg.inv(A[np.ix_(B,B)])
   except np.linalg.LinAlgError:return np.full(len(ii),1e6)
   F=A[np.ix_(T,B)];Z=A[np.ix_(T,T)]-F@inv@F.T
   return Z[ii,jj]
  best=None
  for height in[1.02,1.2,2.0]:
   r=least_squares(residual,np.full(len(free),height),bounds=(1.000001,20),max_nfev=250,ftol=1e-11,xtol=1e-11,gtol=1e-11)
   counters['metric_trials']+=1;err=max(abs(residual(r.x)))
   if best is None or err<best['residual']:best={'residual':float(err),'values':r.x.tolist(),'nfev':r.nfev}
   if err<1e-7:
    A=matrix(r.x);ev=np.linalg.eigvalsh(A)
    if sum(ev<-1e-7)==1 and sum(ev>1e-7)==9:
     counters['metric_solutions']+=1;survivors.append({'support_index':si,'new_dotted_mask':mask,'basis':B,'free_pairs':free,**best,'gram':A.tolist()});break
  trials.append({'support_index':si,'new_dotted_mask':mask,'free_pairs':free,**best})
  if counters['patterns']%20==0:print(name,counters,'best',best['residual'],flush=True)
 print(name,counters,'seconds',time.time()-t,flush=True)
 return {'id':name,'counts':counters,'fixed_principal_witnesses':witness,'metric_trials':trials,'survivors':survivors,'elapsed_seconds':time.time()-t}
def main():
 p=argparse.ArgumentParser();p.add_argument('--case');a=p.parse_args();d=json.loads(SOURCE.read_text());rows=[]
 for name in['A16_14_1_C01','A14_10_1_C02','A14_8_1_C02']:
  if a.case and name!=a.case:continue
  rows.append(search(next(r for r in d['cases']if r['id']==name)))
  OUT.write_text(json.dumps({'status':'RIGHT_DOTTED_INCIDENCE_SPLITTING_NUMERICAL_SCOUT','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'cases':rows,'scope':'Every minimal finite-edge affine cover; minimum number of its edges made dotted subject to one dotted pair per old cusp. Other cover edges made right. Existing parallel pairs also freed dotted. Numerical fixed-principal gates and bounded distance solves, not an exact classification.'},indent=2)+'\n')
if __name__=='__main__':main()
