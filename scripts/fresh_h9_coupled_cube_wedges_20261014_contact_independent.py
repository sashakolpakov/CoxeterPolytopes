#!/usr/bin/env python3
"""Separate exact LDL verification of the first state's forced-right gates."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from fresh_h9_coupled_cube_wedges_20261014_independent import F,Z,O,COS,add,neg,mul,inv,sign
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_coupled_cube_wedges_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 paths=[f'results/{TAG}.json',f'results/{TAG}_csp.json',f'results/{TAG}_first_completion_probe.json',
  f'scripts/{TAG}_independent.py',f'scripts/{TAG}_contact_independent.py']
 g,s,p=[json.loads((ROOT/x).read_text())for x in paths[:3]];g=g['geometry']
 base={(a,b):m for a,b,m in s['solution']['pair_labels']};V=[z['incident']for z in g['vertices']]
 missing={tuple(m)for m in g['minimal_nonfaces']if len(m)==2}
 omitted=set(combinations(range(18),2))-set(base)-missing
 forced={tuple(q)for q,i in p['forced_right']};other={tuple(q)for q,d in p['remaining_domains']}
 assert len(forced)==111 and len(other)==7 and not forced&other and forced|other==omitted
 rows=[]
 for pair,k in p['forced_right']:
  a=V[k];assert set(pair)<=set(a);lab=base|{tuple(pair):3}
  A=[[O if i==j else neg(COS[lab.get(tuple(sorted((i,j))),2)])for j in a]for i in a]
  piv=[]
  for j in range(9):
   q=A[j][j];piv.append(q)
   if sign(q)<=0:break
   for x in range(j+1,9):
    for y in range(j+1,9):A[x][y]=add(A[x][y],neg(mul(mul(A[x][j],A[j][y]),inv(q))))
  assert sign(piv[-1])<=0
  rows.append(dict(pair=pair,vertex=k,leading_order=len(piv),last_pivot=[str(x)for x in piv[-1]],sign=sign(piv[-1])))
 return dict(status='PASS_SEPARATE_111_FORCED_RIGHT_CONTACTS',sources={x:sha(x)for x in paths},records=rows,
  forced_pairs=111,remaining_pairs=7,
  scope='Separate Q(sqrt2) Fraction LDL checks each literal actual-vertex failure at the least nonright label3. Perron-Frobenius monotonicity then excludes every larger finite label and any simultaneous additional nonright contacts. Applies only to the saved mandatory assignment. The listed seven remaining domains are not independently exhausted here.')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');a=a.parse_args();d=run();p=ROOT/f'results/{TAG}_contact_independent.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['forced_pairs'])
