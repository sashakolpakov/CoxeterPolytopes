#!/usr/bin/env python3
"""Remove two unused cut leaves from a reviewed five-family rank proof.

This is a source-bound literal-dependency certificate, not a new determinant
enumeration. No rank-monotonicity assertion is made for the full Gram matrix.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_single_right_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
BASE='results/fresh_h9_star_zero_20261014_gates.json'
OLD='results/fresh_h9_star_dense_frontier_20261014.json'
S=[3,4,6,7,9,10,12,13]
NAMES=['fresh_h9_star_zero_20261014_own_three','fresh_h9_star_zero_20261014_opposite_single',
       'fresh_h9_star_lastblock_20261014','fresh_h9_star_shared_high_20261014',
       'fresh_h9_star_opposite_high_20261014']
FREE={(6,15),(7,15),(9,16),(10,16)}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 seen={}
 def read(p):
  if p in seen:return seen[p]
  d=json.loads((ROOT/p).read_text());seen[p]=d
  for q,h in d.get('sources',{}).items():
   assert sha(q)==h,q
   if q.endswith('.json'):read(q)
  return d
 g=read(GEOM)['representatives'][0]['geometry'];V=[sorted(v['incident'])for v in g['vertices']]
 base=read(BASE);old=read(OLD)
 assert old['status']=='PASS_MONOTONE_NONCORE_OWN_CONTACT_WIDENING'
 core=sorted(S+[11,14,17]);cp=set(combinations(core,2))
 assert len(core)==11 and all(p not in cp for p in FREE)
 assert all(any(set(p)<=set(v)for v in V)for p in cp)
 assert any(set(S)<=set(v)for v in V)
 by={f['id']:f for f in base['families']}
 sources=[GEOM,BASE,OLD];audits=[];bound_rows=[];occurrences=0
 def witness(row):
  vi=row['actual_vertex'];assert 0<=vi<len(V)
  if 'support'in row:support=sorted(row['support']);mode='literal_support'
  else:
   order=row.get('prefix_order',row.get('leading_order',9))
   assert 1<=order<=9;support=V[vi][:order];mode='leading_principal_prefix'
  assert set(support)<=set(V[vi])
  assert not any(set(pair)<=set(support)for pair in FREE),(row,support)
  return dict(actual_vertex=vi,support=support,mode=mode)
 def scan(x,path,rows):
  if isinstance(x,dict):
   if 'actual_vertex'in x:rows.append(dict(path=path,**witness(x)))
   for k,v in x.items():scan(v,path+[k],rows)
  elif isinstance(x,list):
   for i,v in enumerate(x):scan(v,path+[i],rows)
 first=read('results/fresh_h9_star_cube_metric_20261014.json')
 for name in NAMES:
  p=f'results/{name}.json';r=f'results/{name}_review.json';py=f'scripts/{name}.py'
  d=read(p);review=read(r)
  assert d['status'].startswith('PASS_')and review['status'].startswith('PASS_')
  assert review['sources'][p]==sha(p)
  rows=[]
  if name.endswith('own_three'):
   assert d['zero_rank_cases']==[] and d['unresolved']==[]
   for i,row in enumerate(d['foreign_budgets']):
    if row['leaf']in(11,14,17):rows.append(dict(path=['foreign_budgets',i],**witness(row)))
   assert len(rows)==3
   needed=[row for row in first['forced_right_witnesses']if tuple(row['pair'])in cp]
   assert len(needed)==25 and any(row['pair']==[11,17]for row in needed)
   needed+=first['third_contact_upper_witnesses'];assert len(needed)==27
   for i,row in enumerate(needed):rows.append(dict(path=['inherited_needed_guard',i],**witness(row)))
   ignored_pairs=[row['pair']for row in first['forced_right_witnesses']if tuple(row['pair'])not in cp]
   assert ignored_pairs==[[2,14],[5,15],[8,16]]
   ignored_budgets=[row['leaf']for row in d['foreign_budgets']if row['leaf']not in(11,14,17)]
   assert sorted(ignored_budgets)==[2,5,8,15,16]
  else:
   scan(d,[],rows);ignored_pairs=[];ignored_budgets=[]
  counts=Counter(tuple(row['support'])for row in rows)
  occurrences+=len(rows)
  high=d.get('arbitrary_high_block_reduction',d.get('first_three_high_reductions'))
  for row in high:
   f=by[row['reviewed_gate']];assert not f['fixed']
   assert all((a,b)not in FREE for a,b,m in f['bounds'])
   bound_rows.append(dict(source=name,family=f['id'],bounds=f['bounds']))
  audits.append(dict(name=name,needed_witness_occurrences=len(rows),
    unique_witness_supports=[dict(support=list(k),occurrences=v)for k,v in sorted(counts.items())],
    literal_witnesses=rows,unused_own_pair_guards=ignored_pairs,unused_leaf_budgets=ignored_budgets,
    rank_core=core,arithmetic_change='none'))
  sources.extend([p,r,py])
 assert [a['needed_witness_occurrences']for a in audits]==[30,96,708,479,3643]
 assert occurrences==4956 and len(bound_rows)==17
 composed=read('results/fresh_h9_star_lastblock_complete_20261014.json')
 logical=[]
 for row in composed['strong_own_threshold_cases']+composed['dense_own_cases']:
  f=by[row['family']];assert not f['fixed']and f['bounds']==row['bounds']
  assert all((a,b)not in FREE for a,b,m in f['bounds'])
  logical.append(dict(family=f['id'],bounds=f['bounds']))
 assert len(logical)==6
 previous=next(f for f in old['families']if f['id']=='dense_noncore_own_last_internal_ge4')
 removed=[row for row in previous['bounds']if tuple(row[:2])in FREE]
 assert sorted(removed)==[[6,15,3],[7,15,3],[9,16,3],[10,16,3]]
 fixed=previous['fixed'];bounds=[row for row in previous['bounds']if tuple(row[:2])not in FREE]
 assert fixed==[[3,14,3],[4,14,3]]and len(bounds)==10
 families=old['families']+[dict(id='unrestricted_cut_leaves_15_16',
   origin='five_family_literal_dependency_removal',fixed=fixed,bounds=bounds)]
 sources+=['results/fresh_h9_star_cube_metric_20261014.json',
   'results/fresh_h9_star_lastblock_complete_20261014.json',f'scripts/{TAG}.py']
 return dict(status='PASS_TWO_NONCORE_CUT_LEAVES_UNRESTRICTED',case=0,
  sources={p:sha(p)for p in sources},rank_core=core,
  unrestricted_own_pairs=[list(p)for p in sorted(FREE)],fixed=fixed,bounds=bounds,
  five_family_audits=audits,witness_occurrences=occurrences,
  bound_only_arithmetic_reductions=bound_rows,bound_only_composition_reductions=logical,
  inherited_residual_cases=155680,families=families,
  proof=['All four removed contact pairs are outside the eleven-normal rank core, so every residual formula and rank determinant is unchanged.',
   'Every used actual-vertex witness has a saved principal support containing none of the four removed pairs. Its previously checked nonpositive Gram is literally unchanged when those pairs vary.',
   'The all-own-three proof needs only the three core-leaf budgets,24 cross-block right guards,the mutual anchor right guard and two third-contact upper guards. Its three unused leaf-pair guards and five unused leaf budgets are omitted; it has no zero-rank cases.',
   'Every selected high-internal or high-own reduction has lower-bound premises independent of these four pairs. The last-internal>=4 case partition is unchanged.',
   'The remaining six positive-leaf own contacts may increase from3 by the previously reviewed noncore monotonicity argument. Four internal>=4 premises remain explicit.'],
  verification_scope='New literal dependency and stored-support audit, inheriting five separately reviewed arithmetic banks. No new determinant bank or whole-incidence traversal.',
  scope='Star0: positive leaves2,5,8 each have both own labels>=3; leaf14 own labels exactly3/3; all four internal labels>=4. Own contacts of leaves15,16 and both last anchors11,17 are unrestricted finite labels; all other meeting contacts and dotted distances remain free. No whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['witness_occurrences'],len(d['families']))
