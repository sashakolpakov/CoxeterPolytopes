#!/usr/bin/env python3
"""Unconditional rank29 Vinberg-weight savings targets; no new face rule."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact29_savings_frontier_20261018'
BASE='fresh_compact_dimension_bound_20261015'
BARRIER='fresh_compact29_small_faces_20261015'
TYPES={4:['A4','B4','D4','F4','H4'],7:['A7','B7','D7','E7'],8:['A8','B8','D8','E8'],9:['A9','B9','D9'],10:['A10','B10','D10']}
SAVING={4:F(3,4),7:F(3,4),8:F(3,4),9:F(3,4),10:F(1,8)}
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def weight(n):
 if n<=4:return n*(n-1)
 if n<=7:return n*(n-1)//2
 if n==8:return 24
 if n==9:return 27
 return n+19

def partitions(n,start=1):
 if not n:yield ()
 for k in range(start,n+1):
  for q in partitions(n-k,k):yield(k,*q)

def typed(p):
 groups=[list(combinations_with_replacement(TYPES[r],n))for r,n in sorted(Counter(p).items())]
 return [tuple(t for g in row for t in g)for row in product(*groups)]

def run():
 old=json.loads((ROOT/f'results/{BASE}.json').read_text());bar=json.loads((ROOT/f'results/{BARRIER}.json').read_text())
 for d in(old,bar):
  assert d['status'].startswith('PASS')
  for p,h in d['sources'].items():assert digest(p)==h
 ps=list(partitions(29));assert len(ps)==4565
 rows=[];fs=[]
 for p in ps:
  w=sum(weight(r)for r in p)
  if w<=84:continue
  assert all(r in TYPES for r in p)
  ts=typed(p);save=sum(SAVING[r]for r in p)
  assert w-save<=84
  row=dict(ranks=p,weight=w,required_saving=w-84,rank_multiplicities=dict(Counter(p)),typed_count=len(ts),conditional_saving=save,conditional_upper=w-save)
  rows.append(row)
  for t in ts:fs.append(dict(types=t,old_weight=w,required_saving=w-84,conditional_saving=save,conditional_upper=w-save))
 assert len(rows)==9 and len(fs)==2126 and len({tuple(r['types'])for r in fs})==2126
 frozen={tuple(t)for r in old['extremal_rank_partitions']for t in r['typed_forests']}
 assert frozen=={tuple(r['types'])for r in fs}
 counts=dict(sorted(Counter(r['old_weight']for r in fs).items()));assert counts=={85:18,86:660,87:1448}
 maxima=[r for r in rows if r['weight']==87];assert len(maxima)==6
 # Any rank-only savings bounded by t must satisfy3*x7+x8>=3, so t>=3/4.
 assert next(r for r in rows if r['ranks']==(7,7,7,8))['required_saving']==3
 assert max(SAVING.values())==F(3,4)
 # At that cap, x9+2*x10>=1 forces x10>=1/8; the saved vector attains both.
 assert SAVING[9]+2*SAVING[10]==1 and min(r['conditional_upper']for r in rows)>0
 obstructing=('E7','E7','E7','E8');assert sum(weight(int(t[1:]))for t in obstructing)==87
 assert max(r['conditional_upper']for r in fs)==84
 paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',f'results/{BASE}.md',f'results/{BASE}_review.md',f'scripts/{BARRIER}.py',f'results/{BARRIER}.json',f'results/{BARRIER}.md',f'results/{BARRIER}_review.md']
 return enc(dict(status='PASS_UNCONDITIONAL_COMPACT29_SAVINGS_DIAGNOSTIC',sources={p:digest(p)for p in paths},
  dimension=29,arithmetic_assumption=False,all_rank_partitions=4565,threshold='84',current_maximum='87',
  above_threshold_typed_forests=fs,rank_inequalities=rows,weight_histogram=counts,exact_maximum_rank_patterns=maxima,
  conditional_rank_savings=SAVING,conditional_maximum='84',minimal_max_rank_saving='3/4',minimal_rank10_saving_at_that_cap='1/8',
  minimax_witness_rank_pattern=[7,7,7,8],rank10_witness_pattern=[9,10,10],
  scalar_target_is_not_proved_face_rule=True,existing_local_barrier_forest=obstructing,
  current_balance_lower_bound='1050',necessary_heavy_vertex_count=350,
  new_face_theorem=False,new_compact_dimension_bound=None,new_finite_volume_example=False,
  scope='Fresh exact rank/type and savings-target diagnostic under the published unconditional Vinberg rule. Original2126 classification and79 local witnesses are inherited and explicitly not new. No arithmetic field restriction, face reweighting proof or compact29 exclusion.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/f'results/{TAG}.json'
 assert a.verify or not p.exists();out=run()
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['weight_histogram'],flush=True)
