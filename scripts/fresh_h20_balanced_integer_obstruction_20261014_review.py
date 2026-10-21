#!/usr/bin/env python3
"""Separate matrix/witness arithmetic; not a second composition enumeration."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import argparse, hashlib, json, math
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_balanced_integer_obstruction_20261014_review'
SOURCE='scripts/fresh_h20_balanced_integer_obstruction_20261014.py'
DATA='results/fresh_h20_balanced_integer_obstruction_20261014.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def invert(B):
 n=len(B);A=[list(row)+[F(i==j)for j in range(n)]for i,row in enumerate(B)]
 for j in range(n):
  k=next(i for i in range(j,n)if A[i][j]);A[j],A[k]=A[k],A[j]
  d=A[j][j];A[j]=[x/d for x in A[j]]
  for i in range(n):
   if i!=j:
    d=A[i][j];A[i]=[x-d*y for x,y in zip(A[i],A[j])]
 return[row[n:]for row in A]
def run():
 d=json.loads((ROOT/DATA).read_text());assert all(sha(p)==h for p,h in d['sources'].items())
 Js=[];factorV=[]
 for n in[6,6,7]:
  B=[[F(i==j)-F(abs(i-j)==1,2)for j in range(n)]for i in range(n)]
  J=invert(B);Js.append(J);vs=[]
  for omit in range(n+1):
   rhs=[1-(n+1)*(j+1==omit)for j in range(n)]
   v=[sum(x*y for x,y in zip(row,rhs))for row in J]
   assert all(x.denominator==1 for x in v);vs.append(list(map(int,v)))
  factorV.append(vs)
 verts=[sum(parts,[])for parts in product(*factorV)];assert len(verts)==392
 for value in[F(7,4),F(7,2),F(21,4)]:
  assert math.isqrt(value.numerator)**2!=value.numerator or math.isqrt(value.denominator)**2!=value.denominator
 replayed=endpointchecks=0
 for row in d['probes']:
  lam=F(row['multiple']);seen=set()
  for r in row['records']+row['tangent_records']:
   gaps=r['gaps'];key=tuple(map(tuple,gaps));assert key not in seen;seen.add(key)
   norms=[];coeff=[];rhs=[]
   for g,J in zip(gaps,Js,strict=True):
    n=len(J);assert len(g)==n+1 and sum(g)==n+1
    assert all(isinstance(x,int)and not isinstance(x,bool)and x>=0 for x in g)
    a=[1-x for x in g[1:]];c=[sum(x*y for x,y in zip(b,a))for b in J]
    norms.append(sum(x*y for x,y in zip(a,c)));rhs+=a;coeff+=c
   norm=sum(norms);assert norm*28==r['norm28']
   power=norm-2*lam*sum(coeff)+196*lam*lam-4
   assert power*700==r['probe_power_numerator'] and power<=0
   central=(norm in[F(8),F(12)]or norm>=16);assert central==r['central_ok']
   bad=None
   # Rational scalar/vertex products reproduce all relevant endpoint distances.
   n28=int(norm*28)
   for i,v in enumerate(verts):
    dist28=n28+5488-56*sum(x*y for x,y in zip(rhs,v));a=dist28-896
    endpointchecks+=1
    if a==0 or(a>0 and a*a>=351232):continue
    bad={'vertex':i,'distance28':dist28};break
   assert bad==r['first_incompatible_endpoint'];replayed+=1
  assert len(row['records'])==555 and not row['tangent_records']
 assert F(1764,25)-2*F(98783,2800)==F(1,1400)
 assert F(1,700)>F(1,1400)>0
 paths=[SOURCE,DATA,f'scripts/{TAG}.py']
 return dict(status='PASS_INTEGER_FAMILY_SOURCE_AND_SEPARATE_WITNESS_REVIEW',
  sources={p:sha(p)for p in paths},record_checks=replayed,endpoint_pair_checks=endpointchecks,
  independent_method='Fraction inversion of literal finite A6/A7 Gram matrices, exact saved-recipe reconstruction and endpoint rejection witnesses; no producer imports.',
  source_proof_review='Complete weak compositions; factor-distance pruning is necessary; rational/sqrt7 all-angle gate is exhaustive; positive powers quantized in 1/700 give norm1/1400 strict escape.',
  primary_saved_replay='PASS',scope='Source-level completeness audit plus fresh primary replay and separate saved-record arithmetic, not a second exhaustive composition enumeration.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();f=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert json.loads(f.read_text())==d
 else:f.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['record_checks'],d['endpoint_pair_checks'])
