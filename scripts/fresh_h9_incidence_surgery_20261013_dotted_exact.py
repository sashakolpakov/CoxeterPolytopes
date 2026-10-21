#!/usr/bin/env python3
"""Exact fixed-principal obstructions for151 incidence-splitting patterns.

No numerical predicates, nonlinear solver, or producer arithmetic imported.
The infinite distance parameters do not occur in any saved obstruction.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_h9_incidence_surgery_20261013.json'
PROPOSALS=ROOT/'results/fresh_h9_incidence_surgery_20261013_dotted.json'
OUT=ROOT/'results/fresh_h9_incidence_surgery_20261013_dotted_exact.json'
def rank_psd(A):
 A=[r[:]for r in A];piv=[]
 while A:
  assert all(A[i][i]>=0 for i in range(len(A)))
  k=next((i for i in range(len(A))if A[i][i]),None)
  if k is None:assert not any(any(r)for r in A);return len(piv),piv
  A[0],A[k]=A[k],A[0]
  for r in A:r[0],r[k]=r[k],r[0]
  q=A[0][0];piv.append(q)
  A=[[A[i][j]-A[i][0]*A[0][j]/q for j in range(1,len(A))]for i in range(1,len(A))]
 return len(piv),piv

def poskernel(A):
 n=len(A);M=[r[:]for r in A];p=[];k=0
 for j in range(n):
  t=next((i for i in range(k,n)if M[i][j]),None)
  if t is None:continue
  M[k],M[t]=M[t],M[k];v=M[k][j];M[k]=[x/v for x in M[k]]
  for i in range(n):
   if i!=k:
    v=M[i][j];M[i]=[x-v*y for x,y in zip(M[i],M[k])]
  p.append(j);k+=1
 assert k==n-1
 f=next(j for j in range(n)if j not in p);v=[F(0)]*n;v[f]=F(1)
 for i,j in enumerate(p):v[j]=-M[i][f]
 assert all(x>0 for x in v) and all(sum(x*y for x,y in zip(row,v))==0 for row in A)
 return v

def expected(row):
 ans=set();P=row['parallel_pairs']
 for si,S in enumerate(row['finite_edge_cover']['minimum_supports']):
  constraints=[]
  for cusp in row['ideal_vertices']:
   A=cusp['active']
   if any(i in A and j in A for i,j in P):continue
   constraints.append(sum(1<<k for k,(i,j)in enumerate(S)if i in A and j in A))
  ms=[m for m in range(1<<len(S))if all(m&c for c in constraints)];minimum=min(m.bit_count()for m in ms)
  ans.update((si,m)for m in ms if m.bit_count()==minimum)
 return ans

def verify():
 d=json.loads(SOURCE.read_text());p=json.loads(PROPOSALS.read_text());assert p['source_sha256']==hashlib.sha256(SOURCE.read_bytes()).hexdigest();rows=[]
 for row in d['cases']:
  pr=next(r for r in p['cases']if r['id']==row['id']);assert not pr['survivors'] and not pr['metric_trials'] and not pr['counts']['unresolved']
  assert {(r['support_index'],r['new_dotted_mask'])for r in pr['fixed_principal_witnesses']}==expected(row)
  cert=[];counts={'positive_definite_rank10':0,'positive_affine_component':0}
  for witness in pr['fixed_principal_witnesses']:
   si=witness['support_index'];mask=witness['new_dotted_mask'];S=row['finite_edge_cover']['minimum_supports'][si];G=[[F(x)for x in r]for r in row['source_positive_gram']];free={tuple(p)for p in row['parallel_pairs']}
   for k,(i,j)in enumerate(S):
    if mask>>k&1:free.add((i,j))
    else:G[i][j]=G[j][i]=F(0)
   B=witness['principal_nodes'];assert len(B)==10 and len(set(B))==10
   assert not any(i in B and j in B for i,j in free)
   A=[[G[i][j]for j in B]for i in B];rank,pivots=rank_psd(A)
   out={'support_index':si,'new_dotted_mask':mask,'principal_nodes':B,'rank':rank,'positive_pivots':[str(x)for x in pivots]}
   if rank==10:counts['positive_definite_rank10']+=1;out['obstruction']='positive_index_at_least10'
   else:
    unseen=set(B);affine=None
    while unseen:
     c={min(unseen)};stack=list(c)
     while stack:
      i=stack.pop()
      for j in B:
       if j not in c and G[i][j]:c.add(j);stack.append(j)
     unseen-=c;nodes=sorted(c);AA=[[G[i][j]for j in nodes]for i in nodes];rr,_=rank_psd(AA)
     if rr<len(nodes):affine=(nodes,poskernel(AA));break
    assert affine is not None
    out['obstruction']='unchanged_positive_affine_relation';out['affine_nodes']=affine[0];out['positive_kernel']=[str(x)for x in affine[1]];counts['positive_affine_component']+=1
   cert.append(out)
  rows.append({'id':row['id'],'counts':counts,'certificates':cert})
 return {'status':'EXACT_ALL_DISTANCE_OBSTRUCTION_FOR151_RIGHT_DOTTED_PATTERNS','source_hashes':{str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest()for x in[SOURCE,PROPOSALS]},'cases':rows,
 'scope':'All151 patterns selected by the stated minimal-cover/minimal-dotted rule are excluded for every real choice of their dotted weights. This is not all mutations, all incidence surgeries, or all compact H9 polytopes.'}
def main():
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=verify()
 if a.verify:assert json.loads(OUT.read_text())==d;print('PASS saved Fraction certificates')
 else:OUT.write_text(json.dumps(d,indent=2)+'\n');print('WROTE',OUT)
 for r in d['cases']:print(r['id'],r['counts'])
if __name__=='__main__':main()
