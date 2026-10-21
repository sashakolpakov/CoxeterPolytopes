#!/usr/bin/env python3
"""Separate literal-premise and generic-routing review; never executes a DFS."""
import argparse,ast,json
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_tnext_20261015_scope_review'
BASE='fresh_h9_tnext_20261015'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
TARGET='results/fresh_h9_dotted_central_20261014_next_targets.json'
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 paths=[f'scripts/{BASE}_{s}.py'for s in ['gates','search','independent']]+[
  f'results/{BASE}_gates_case_{c:02d}.json'for c in [20,11,31]]+[
  GEOM,TARGET,'results/fresh_h9_star_three_large_blocks_20261015.json',
  'results/fresh_h9_star_three_large_blocks_20261015_review.json',
  'scripts/fresh_h9_t6_frontier_20261015_independent.py',
  'results/fresh_h9_t6_frontier_20261015_gates_review.json',f'scripts/{TAG}.py']
 seen={}
 def load(p):
  if p in seen:return seen[p]
  d=json.loads((ROOT/p).read_text());seen[p]=d
  for q,h in d.get('sources',{}).items():
   assert digest(q)==h,(p,q)
   if q.endswith('.json'):load(q)
  return d
 for p in paths:
  if p.endswith('.json'):load(p)
 def funcs(p):
  return {x.name:ast.dump(x,include_attributes=False)for x in ast.walk(ast.parse((ROOT/p).read_text()))if isinstance(x,ast.FunctionDef)}
 old,new=map(funcs,['scripts/fresh_h9_t6_frontier_20261015_independent.py',f'scripts/{BASE}_independent.py'])
 assert old.keys()==new.keys()
 identical=[k for k in old if k!='run']
 for k in identical:assert old[k]==new[k],k
 # The enclosing function routes incidence, gates and saved output by the same case.
 texts={s:(ROOT/f'scripts/{BASE}_{s}.py').read_text()for s in ['gates','search','independent']}
 assert "['representatives'][case]['geometry']" in texts['gates'] and "['cases'][case]" in texts['gates']
 assert "['representatives'][case]['geometry']" in texts['search'] and 'GATE.format(case=case)' in texts['search']
 assert "['representatives'][case]['geometry']" in texts['independent'] and 'GATE.format(case=case)' in texts['independent']
 assert "assert gate['case']==case" in texts['search'] and "assert gate['case'] == case" in texts['independent']
 assert 'while stack:' in texts['search'] and "'PAUSED_EXACT_PENDING_DFS' if stack" in texts['search']
 assert "complete=True, cutoff=None, survivor=None, unexplored_frontier=None" in texts['independent']
 assert 'visit((0,)*len(pairs)' in texts['independent']
 out=[]
 expected={20:(488,155,2,157,2288),11:(488,154,2,156,2288),31:(462,149,4,153,2184)}
 for case in [20,11,31]:
  g=load(GEOM)['representatives'][case]['geometry'];missing=set(map(frozenset,g['minimal_nonfaces']))
  vertices={tuple(v['incident'])for v in g['vertices']};vlist=[set(v['incident'])for v in g['vertices']]
  def face(s):return not any(m<=set(s)for m in missing)
  assert {len(m)for m in missing}=={2,3}
  assert vertices=={v for v in combinations(range(18),9)if face(v)}
  gate=load(f'results/{BASE}_gates_case_{case:02d}.json');target=load(TARGET)['cases'][case]
  assert gate['case']==target['case']==case
  assert gate['status']=='PASS_TNEXT_LITERAL_REVIEWED_CONDITIONS'
  def witness(w):
   assert face(w['support']) and set(w['support'])<=vlist[w['actual_vertex']]
  bounds=[];occurrences=0;counts={}
  for name,rec in target['gates'].items():
   for r in rec['maps']:
    A,B,x,y=tuple(r['A']),tuple(r['B']),r['x'],r['y']
    assert len(set(A+B+(x,y)))==6
    own={frozenset(A+(x,)),frozenset(B+(y,))}
    assert own==set(map(frozenset,r['missing_own_triples'])) and own<=missing
    premises=[A+B,B+(x,),A+(y,)]+[(x,y,a,b)for a in A for b in B]
    def t(a,b,k):return (*sorted((a,b)),k)
    if name=='dense_high':thresholds=[(*sorted(A),6),(*sorted(B),7)]+[t(x,a,3)for a in A]
    elif name=='strong_contact':
     z=[b for b in B if list(t(y,b,4))in r['bounds']];assert len(z)==1
     premises.append(A+(y,z[0]));thresholds=[(*sorted(A),4),(*sorted(B),6)]+[t(x,a,3)for a in A]+[t(y,z[0],4)]
    else:
     assert name=='two_dense';premises.extend((x,y,b)for b in B)
     thresholds=[(*sorted(A),6),(*sorted(B),4)]+[t(x,a,3)for a in A]+[t(y,b,3)for b in B]
    assert sorted(tuple(sorted(s))for s in premises)==sorted(tuple(w['support'])for w in r['supports'])
    assert sorted(thresholds)==sorted(map(tuple,r['bounds']))
    assert sorted(v[:2]for v in thresholds)==sorted(tuple(w['support'])for w in r['threshold_pair_witnesses'])
    for w in r['supports']+r['threshold_pair_witnesses']:witness(w);occurrences+=1
   for c in rec['conditions']:
    assert c['bounds']==rec['maps'][c['witness_map']]['bounds'];bounds.append(c['bounds'])
   counts[name]=len(rec['conditions'])
  inherited=len(bounds)
  assert bounds==[f['bounds']for f in gate['families'][:inherited]]
  assert all(f['fixed']==[]for f in gate['families'])
  blocks=[tuple(b)for b in g['blocks']if len(b)==2];assert len(blocks)==4
  leaves=[[i for i in range(18)if frozenset(b+(i,))in missing]for b in blocks]
  assert not set(sum(leaves,[]))&set(sum(blocks,()))
  choices={};crossgood={}
  for chosen in combinations(range(4),3):
   ok=all(face(blocks[i]+(b,))and face(blocks[j]+(a,))for i,j in combinations(chosen,2)for a in blocks[i]for b in blocks[j])
   crossgood[chosen]=ok
   choices[chosen]=[ls for ls in product(*(leaves[i]for i in chosen))if ok and len(set(ls))==3 and face(ls)and all(face(blocks[j]+(ls[k],))for k,i in enumerate(chosen)for j in chosen if j!=i)]
  transfers=gate['three_large_transfers']
  assert {tuple(r['chosen_blocks'])for r in transfers}=={c for c,v in choices.items()if v}
  assert gate['uncovered_three_large_block_triples']==[list(c)for c,v in choices.items()if not v]
  for r in transfers:
   c=tuple(r['chosen_blocks']);ls=tuple(r['chosen_leaves'])
   assert ls==choices[c][0] and len(choices[c])==r['eligible_leaf_choices']
   assert r['family']['bounds']==[list(blocks[i])+[6]for i in c]
   assert {frozenset(m)for m in r['own_Lanner_triples']}=={frozenset(blocks[i]+(ls[k],))for k,i in enumerate(c)}
   assert len(r['cross_right_guards'])==len(r['foreign_right_guards'])==12
   for w in r['cross_right_guards']+r['foreign_right_guards']:
    witness(w);assert tuple(w['large_pair'])in blocks
    assert set(w['large_pair'])|set(w['trial_nonright_pair'])==set(w['support'])
    assert len(set(w['large_pair'])&set(w['trial_nonright_pair']))==1
   witness(r['leaf_triangle'])
  for take in product((0,1),repeat=4):
   perm=list(range(18))
   for bit,(a,b)in zip(take,blocks):
    if bit:perm[a],perm[b]=b,a
   assert {tuple(sorted(perm[i]for i in v))for v in vertices}==vertices
   assert {frozenset(perm[i]for i in m)for m in missing}==missing
  assert (len(vertices),inherited,len(transfers),len(gate['families']),occurrences)==expected[case]
  assert counts==gate['old_counts'] and occurrences==gate['inherited_support_occurrences']
  assert len({tuple(map(tuple,f['bounds']))for f in gate['families']})==len(gate['families'])
  out.append(dict(case=case,vertices=len(vertices),inherited_predicates=inherited,new_predicates=len(transfers),initial_total=len(gate['families']),inherited_witness_occurrences=occurrences,all_three_block_choices=[dict(blocks=list(c),cross_guards_hold=crossgood[c],leaf_choices=[list(x)for x in v])for c,v in choices.items()],simultaneous_symmetries=16))
 return dict(status='PASS_TNEXT_INITIAL_BINDERS_AND_GENERIC_ROUTING_SCOPE_REVIEW',sources={p:digest(p)for p in paths},recursive_manifests=len(seen),cases=out,unchanged_independent_function_bodies=identical,
  scope='Independent literal-premise reconstruction and source/routing review only. No DFS rerun, new analytic predicate proof, or independent whole-incidence exclusion. The initial case31 binder has153 conditions and does not close case31 by itself; its later155-condition transfer/resume and completed traversal are separate reviewed dependencies outside this note. Symbolic category7 denotes every finite label>=7; isolated rank2 components remain unrestricted. No finite-label cutoff is inferred.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],[(r['case'],r['initial_total'])for r in d['cases']])
