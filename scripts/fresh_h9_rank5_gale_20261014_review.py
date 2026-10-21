#!/usr/bin/env python3
"""Narrow source review ledger and independent saved primal arithmetic.

This does not repeat the 1063-case Gale or label census.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_rank5_gale_20261014'
OUT=ROOT/f'results/{TAG}_review.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for name,h in d.get('sources',{}).items():assert sha(name)==h,name
 return d
def rank(A):
 A=[list(map(F,r))for r in A];r=0
 for c in range(len(A[0])):
  p=next((j for j in range(r,len(A))if A[j][c]),None)
  if p is None:continue
  A[p],A[r]=A[r],A[p];a=A[r][c];A[r]=[x/a for x in A[r]]
  for j in range(r+1,len(A)):
   a=A[j][c];A[j]=[x-a*y for x,y in zip(A[j],A[r])]
  r+=1
  if r==len(A):break
 return r
def run():
 p=read(f'results/{TAG}.json');a=read(f'results/{TAG}_independent.json');c=read(f'results/{TAG}_certificate.json')
 assert p['status']=='FINISHED_BOUNDED_DIRECT_RANK5_GALE_SEARCH'
 assert a['status']=='PASS_ALL_SAVED_DIRECT_RANK5_GALE_CASES'
 assert c['status']=='PASS_TERMINAL_LEDGER_AND_EXPLICIT_RATIONAL_PRIMAL'
 assert c['saved_nerves']==a['cases_checked']==len(p['cases'])==1063
 assert c['exact_exclusions']=={'literal_induced_pattern':902,'two_lanner_blocks':161}
 assert c['unresolved_saved_nerves']==c['full_label_searches']==c['metric_solves']==0
 assert 120*2048**5<2**63
 data=c['explicit_rational_primal'];P=[[F(*x)for x in row]for row in data['points']]
 assert len(P)==15 and all(len(x)==9 for x in P)
 assert rank([[1]+x for x in P])==10
 first=p['cases'][0];assert sha(first['path'])==first['sha256'];g=read(first['path'])['geometry']
 seen=set()
 for row in data['supporting_functionals']:
  f=[F(*x)for x in row['nonnegative_affine_functional']]
  values=[f[0]+sum(x*y for x,y in zip(f[1:],v))for v in P]
  assert all(x>=0 for x in values)
  mask=sum(1<<i for i,x in enumerate(values)if x==0)
  assert mask==row['facet'] and mask.bit_count()==9 and mask not in seen
  assert rank([[1]+P[i]for i in range(15)if mask>>i&1])==9
  seen.add(mask)
 assert seen==set(g['facets'])and len(seen)==310
 paths=[f'scripts/{TAG}_{suffix}.py'for suffix in ['core','independent','certificate','review']]
 paths += [f'results/{TAG}{suffix}.json'for suffix in ['','_independent','_certificate']]+[first['path']]
 return dict(status='PASS_NARROW_GALE_PROOF_REVIEW_AND_SEPARATE_PRIMAL_CHECK',sources={x:sha(x)for x in paths},
  saved_nerves=1063,source_proof_review=['uniform Gale positive-circuit realization','simplex-boundary link factor is necessary for standard simple-polytope face truncation','independent link ideal computation','literal inherited obstruction scope','native determinant bound'],
  independent_small_arithmetic=dict(points=15,affine_dimension=9,facet_supports=310,facet_affine_rank=8),
  fresh_primary_terminal_replay='PASS',
  scope='Source/proof review plus separate Fraction checks of the one saved explicit primal. The completed 1063-case independent census and inherited local obstruction proofs remain dependencies; no duplicate full census or classification of all15facet nerves.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run()
 if args.verify:assert d==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
