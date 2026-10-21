#!/usr/bin/env python3
"""Source-bound dependency widening on ten own contacts outside the rank core."""
import argparse,ast,json
from pathlib import Path
from hashlib import sha256
from itertools import combinations
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_dense_frontier_20261014'
COMPOSED='results/fresh_h9_star_lastblock_complete_20261014.json';BASE='results/fresh_h9_star_zero_20261014_gates.json'
DEPENDENCY_REVIEW='results/fresh_h9_star_dense_frontier_20261014_dependency_review.json'
FAMILIES=[('fresh_h9_star_zero_20261014_own_three','all_own_contacts',21952),('fresh_h9_star_zero_20261014_opposite_single','fixed_other_own_contacts',896),('fresh_h9_star_lastblock_20261014','fixed_first_three_own',77952),('fresh_h9_star_shared_high_20261014','fixed_own_contacts',7840),('fresh_h9_star_opposite_high_20261014','fixed_own_contacts',47040)]
S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 seen=set()
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p not in seen:
   seen.add(p)
   for q,h in d.get('sources',{}).items():
    assert sha(q)==h,q
    if q.endswith('.json'):read(q)
  return d
 old=read(COMPOSED);base=read(BASE);dep=read(DEPENDENCY_REVIEW)
 assert old['status']=='PASS_STAR_ZERO_COMPLETE_LAST_INTERNAL_GE4_COMPOSITION'
 assert dep['status']=='PASS_FIVE_FAMILY_NONCORE_MONOTONE_DEPENDENCY_REVIEW'
 core=sorted(S+[11,14,17]);assert len(core)==11
 cp={tuple(p)for p in combinations(core,2)}
 common={tuple(sorted((lf,r))):3 for i in range(3)for lf in(P[i],C[i])for r in S[2*i:2*i+2]}
 exact={p:m for p,m in common.items()if p in cp};free={p:m for p,m in common.items()if p not in cp}
 assert exact=={(3,14):3,(4,14):3} and len(free)==10
 assert {a for p in free for a in p if a not in core}=={2,5,8,15,16}
 assert all(len(set(p)&set(core))==1 for p in free)
 source_paths=[COMPOSED,BASE,DEPENDENCY_REVIEW];audits=[];high_rows=[]
 for name,own_key,expected in FAMILIES:
  p=f'results/{name}.json';r=f'results/{name}_review.json';py=f'scripts/{name}.py'
  d=read(p);review=read(r)
  assert d['status'].startswith('PASS_')and review['status'].startswith('PASS_')
  assert review['sources'][p]==sha(p)
  labels={tuple((a,b)):m for a,b,m in d[own_key]}
  assert all(labels.get(pair)==3 for pair in common)
  # Every primary declares this identical spatial coordinate list. The three
  # additional normals in every residual formula are11,17,14. Their precise
  # source-level use is independently reviewed, rather than inferred from a
  # rank-monotonicity assertion for the whole Gram matrix.
  tree=ast.parse((ROOT/py).read_text());Sdefs=[]
  for node in tree.body:
   if isinstance(node,ast.Assign)and any(isinstance(t,ast.Name)and t.id=='S'for t in node.targets):Sdefs.append(ast.literal_eval(node.value))
  assert Sdefs==[S]
  if name.endswith('own_three'):count=d['residual_cases'];rows=d['arbitrary_high_block_reduction']
  elif name.endswith('opposite_single'):count=len(d['cases']);rows=d['arbitrary_high_block_reduction']
  elif name=='fresh_h9_star_lastblock_20261014':count=d['residual_cases'];rows=d['first_three_high_reductions']
  else:count=len(d['cases']);rows=d['first_three_high_reductions']
  assert count==expected
  for row in rows:
   selected=next(f for f in base['families']if f['id']==row['reviewed_gate'])
   assert selected['fixed']==[]
   high_rows.append(dict(source=name,family=selected['id'],bounds=selected['bounds']))
  audits.append(dict(family=name,primary=p,independent_review=r,residual_cases=count,spatial_support=S,residual_normals=[11,17,14],core_own_equalities=sorted([list(q)+[m]for q,m in exact.items()]),relaxed_noncore_own_labels=sorted([list(q)+[m]for q,m in free.items()]),arithmetic_change='none: every changed own pair lies outside the entire eleven-normal rank core',outside_core_usage='bound-only high-label reductions and monotone nonpositive actual-vertex witnesses'))
  source_paths += [p,r,py]
 assert len(high_rows)==17
 # The six additional logical-composition reductions are also bound-only.
 logical=[]
 for row in old['strong_own_threshold_cases']+old['dense_own_cases']:
  f=next(f for f in base['families']if f['id']==row['family']);assert f['fixed']==[]and f['bounds']==row['bounds']
  logical.append(dict(family=f['id'],bounds=f['bounds']))
 assert len(logical)==6
 fixed=sorted([list(p)+[m]for p,m in exact.items()]);bounds=sorted([list(p)+[m]for p,m in free.items()]+[S[2*i:2*i+2]+[4]for i in range(4)])
 assert dep['rank_core']==core and dep['exact_own_contacts']==fixed
 assert dep['noncore_own_lower_bounds']==sorted([list(p)+[m]for p,m in free.items()])
 assert dep['internal_lower_bounds']==[S[2*i:2*i+2]+[4]for i in range(4)]
 assert sum(r['rank_cases_inherited']for r in dep['proof_dependencies'])==sum(r['residual_cases']for r in audits)==155680
 fs=old['families'][:];fs.append(dict(id='dense_noncore_own_last_internal_ge4',origin='reviewed_core_equality_and_monotone_vertex_domains',fixed=fixed,bounds=bounds))
 source_paths += [f'scripts/{TAG}.py']
 return dict(status='PASS_MONOTONE_NONCORE_OWN_CONTACT_WIDENING',case=0,sources={p:sha(p)for p in source_paths},families=fs,rank_core=core,rank_core_pair_count=len(cp),fixed_core_own=fixed,free_noncore_own_lower_bounds=sorted([list(p)+[m]for p,m in free.items()]),internal_lower_bounds=[S[2*i:2*i+2]+[4]for i in range(4)],five_family_dependency_audit=audits,arithmetic_high_reduction_rows=high_rows,composition_high_reduction_rows=logical,proof_steps=['The actual rank-core Gram and every residual arithmetic formula are literally unchanged by changing these ten pairs. No whole-Gram rank monotonicity is asserted.','For symmetric matrices with nonpositive off-diagonal entries, a nonpositive quadratic witness can be replaced by its componentwise absolute value. Decreasing off-diagonal entries therefore preserves failure of positive definiteness.','Each occurrence of a changed own contact outside the rank core is in a relaxed actual-vertex witness or a lower-bound gate. Labels>=3 only decrease its normalized off-diagonal entry from the relaxed label3 value.','The first-three internal>=4 premises are explicit; they are not inferred from changed own-label missing triples.','All five inherited profile-domain proofs and their rank tests therefore persist; the last-internal>=4 logical partition remains unchanged.'],verification_scope='Automated source/hash/core/domain inventory plus source-level dependency proof reviewed independently. This checker does not rerun the inherited determinant banks.',scope='Star0: all four internal labels>=4; leaf14 owncontacts exactly3/3; the other ten first-three-block owncontacts arbitrary finite>=3. Both last-anchor ownprofiles and all remaining finite labels/dotted distances free. Not a whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']),len(d['arithmetic_high_reduction_rows']))
