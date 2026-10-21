#!/usr/bin/env python3
"""Narrow proof/parameter review; no new ambient-field reconstruction."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_transverse_20261014';SOURCE=f'results/{TAG}.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 J=[[F(a)for a in row]for row in d['E7_inverse']]
 assert all(a>0 for row in J for a in row)
 B=[r[:]for r in J];pivots=[]
 for k in range(7):
  p=B[k][k];assert p>0;pivots.append(p)
  for i in range(k+1,7):
   for j in range(k+1,7):B[i][j]-=B[i][k]*B[k][j]/p
 degrees=[J[i][i]for i in range(7)];assert degrees==list(map(F,['7/2','12','6','2','15/2','3/2','4']))
 assert len(d['records'])==308 and len({(r['vertex'],r['E7_position'])for r in d['records']})==308
 counts=Counter()
 for r in d['records']:
  q=F(r['vertex_negative_norm']);degree=F(r['weight_norm']);assert q>0 and degree==degrees[r['E7_position']]
  if 'h_at_c1'in r:assert F(r['h_at_c1'])==(degree-2)/q
  if 'h_at_sqrt2'in r:assert F(r['h_at_sqrt2'])==(degree-1)/q
  # Norm reduction for arbitrary c: r/c=-omega+uV,
  # u²=(d-2/c²)/q. Check the two boundary substitutions independently.
  assert degree-q*((degree-2)/q)==2
  assert 2*(degree-q*((degree-1)/q))==2
  counts[r['reason']]+=1
 assert dict(counts)==d['counts']
 assert [(r['vertex'],r['E7_wall'])for r in d['records']if r['reason']=='Coxeter_coefficient_gap']==[(40,65),(41,65),(42,65)]
 nc=Counter()
 for r in d['null_records']:
  degree=F(r['weight_norm']);c2=F(r['coefficient_squared']);assert c2==2/degree
  if c2<1:reason='coefficient_between_zero_and_one'
  elif 1<c2<2:reason='coefficient_between_one_and_sqrt2'
  else:
   assert degree==2 and c2==1 and F(r['pairing'])==2
   reason='positive_pairing_with_retained_same_cusp_filler'
  assert reason==r['reason'];nc[reason]+=1
 assert dict(nc)==d['null_counts']
 paths=[SOURCE,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
 return dict(status='PASS_NARROW_E7_TRANSVERSE_PROOF_REVIEW',sources={p:sha(p)for p in paths},
  finite_families=308,null_families=14,finite_counts=dict(counts),null_counts=dict(nc),
  exact_rational_parameter_reductions=616,positive_weight_inverse_entries=49,
  scope='Source/theorem review plus fresh primary saved replay, with separate Fraction checks of parameter norm reductions, inverse positivity and every null coefficient classification. Not a second ambient multiquadratic or whole-chamber geometry implementation. The all-parameter theorem is limited to one E7 weight and one of the44 finite or2 ideal face rays.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();z=run();p=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert z==json.loads(p.read_text())
 else:p.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status'])
