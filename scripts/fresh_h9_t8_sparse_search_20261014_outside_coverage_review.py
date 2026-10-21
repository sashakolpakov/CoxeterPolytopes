#!/usr/bin/env python3
"""Independent finite coverage guards and sampled certificates; source review of DFS."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_t8_sparse_search_20261014_outside_coverage_review'
SRC='results/fresh_h9_t8_sparse_search_20261014_outside_case_08.json'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def mul(a,b):return(a[0]*b[0]+5*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def inv(a):
 d=a[0]*a[0]-5*a[1]*a[1];assert d;return(a[0]/d,-a[1]/d)
def sign(a):
 x,y=a;s=lambda q:int(q>0)-int(q<0)
 if not y:return s(x)
 if not x or s(x)==s(y):return s(y)
 return s(x)*s(x*x-5*y*y)
Z=(F(0),F(0));O=(F(1),F(0))
W2={3:(F(1,4),F(0)),4:(F(1,2),F(0)),5:(F(3,8),F(1,8)),6:(F(3,4),F(0))}
def pd_tree(edges):
 """Exact forest Schur test, with a uniform high-edge three-path gate."""
 adj={}
 for (a,b),m in edges.items():
  if m<=2:continue
  adj.setdefault(a,{})[b]=m;adj.setdefault(b,{})[a]=m
 unseen=set(adj)
 while unseen:
  root=min(unseen);component={root};queue=[root];unseen.remove(root)
  for a in queue:
   for b in adj[a]:
    if b in unseen:unseen.remove(b);component.add(b);queue.append(b)
  if len(component)<=2:continue
  if sum(len(adj[a]) for a in component)!=2*(len(component)-1):return False
  if any(m>=7 for a in component for m in adj[a].values()):return False
  def visit(a,parent):
   pivot=O
   for b,m in adj[a].items():
    if b==parent:continue
    p=visit(b,a)
    if sign(p)<=0:return Z
    z=mul(W2[m],inv(p));pivot=add(pivot,(-z[0],-z[1]))
   return pivot
  if sign(visit(root,None))<=0:return False
 return True
def run():
 d=json.loads((ROOT/SRC).read_text());g=json.loads((ROOT/GEOM).read_text())['representatives'][8]['geometry']
 for obj in [d,json.loads((ROOT/GEOM).read_text())]:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 assert d['status']=='EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY'
 assert d['nodes']==59423<d['settings']['node_limit']==100000
 assert d['solution'] is None and d['unexplored_frontier'] is None
 ms=[tuple(x) for x in g['minimal_nonfaces']];ts=[x for x in ms if len(x)==3]
 vs={tuple(v['incident']) for v in g['vertices']};V=sorted(vs)
 # Reconstruct every possible rank-nine incidence from the saved clutter.
 assert vs=={v for v in combinations(range(18),9) if not any(set(m)<=set(v) for m in ms)}
 assert len(vs)==462 and len(ts)==10 and set(map(len,ms))=={2,3}
 pairs=list(combinations(range(18),2));assert pairs==[tuple(x) for x in d['pairs']]
 index={p:i for i,p in enumerate(pairs)};T=[tuple(index[p] for p in combinations(t,2)) for t in ts]
 opts={c for c in product(range(2,8),repeat=3) if sum(F(1,m) for m in c)<1}
 # Every triple with a >=7 category is hyperbolic unless the other two
 # labels are both 2. Raising that category never crosses the threshold.
 high=0
 for c in product(range(2,8),repeat=3):
  if 7 not in c:continue
  expected=not any(c[i]==c[j]==2 for i,j in combinations(range(3),2))
  assert (c in opts)==expected;high+=1
 blocks=[tuple(b) for b in g['blocks'] if len(b)==2];chosen=[];domains=[set(opts) for _ in ts]
 doubled={i for b in blocks for i in b}
 for a,b in blocks:
  perm=list(range(18));perm[a],perm[b]=b,a
  assert {tuple(sorted(perm[i] for i in v)) for v in vs}==vs
  assert {tuple(sorted(perm[i] for i in m)) for m in ms}==set(ms)
  k=next(k for k,t in enumerate(ts) if a in t and b in t);leaf=next(i for i in ts[k] if i not in (a,b))
  assert leaf not in doubled
  cp=list(combinations(ts[k],2));u=cp.index(tuple(sorted((leaf,a))));v=cp.index(tuple(sorted((leaf,b))))
  domains[k]={c for c in opts if c[u]<=c[v]};chosen.append(dict(pair=[a,b],triangle=k,leaf=leaf))
  assert all(c in domains[k] or tuple(c[v] if j==u else c[u] if j==v else c[j] for j in range(3)) in domains[k] for c in opts)
 for bits in product([0,1],repeat=4):
  perm=list(range(18))
  for bit,(a,b) in zip(bits,blocks):
   if bit:perm[a],perm[b]=b,a
  assert {tuple(sorted(perm[i] for i in v)) for v in vs}==vs
 clauses=[]
 for aa,bb in combinations(ms,2):
  A,B=set(aa),set(bb)
  if A&B or sum(set(m)<=A|B for m in ms)!=2:continue
  cross=tuple(sorted(index[tuple(sorted((a,b)))] for a in aa for b in bb))
  assert all(any(set(pairs[k])<=set(v) for v in vs) for k in cross)
  clauses.append(dict(blocks=[list(aa),list(bb)],cross_pairs=list(cross)))
 assert clauses==d['bridge_clauses'] and len(clauses)==60
 def reconstruct(path):
  arr={};selected=set()
  for row in path:
   k=row['triangle'];c=tuple(row['labels']);assert k not in selected and c in domains[k];selected.add(k)
   for i,m in zip(T[k],c):assert i not in arr or arr[i]==m;arr[i]=m
  return arr,selected
 # Recheck the saved negative-block witnesses using exact leaf Schur
 # elimination, without the producer's finite Dynkin classifier.
 vertex_order=[tuple(v['incident']) for v in g['vertices']];failurechecks=0
 for row in d['prune_samples']:
  arr,selected=reconstruct(row['path']);clause=clauses[row['clause']]
  assert [k for k,v in row['forced_cross_witnesses']]==clause['cross_pairs']
  for k,vi in row['forced_cross_witnesses']:
   if k in arr:assert arr[k]==2 and vi==-1
   else:
    support=set(vertex_order[vi]);assert set(pairs[k])<=support
    ed={pairs[i]:m for i,m in arr.items() if set(pairs[i])<=support};ed[pairs[k]]=3
    assert not pd_tree(ed);failurechecks+=1
 families={f['id']:f for f in d['proved_family_conditions']}
 allblocks=json.loads((ROOT/'results/fresh_h9_t8_sparse_search_20261014_all_blocks.json').read_text())
 twin=json.loads((ROOT/'results/fresh_h9_t8_sparse_search_20261014_twin.json').read_text())
 variants=json.loads((ROOT/'results/fresh_h9_t8_sparse_search_20261014_bridge_variants.json').read_text())
 expected_families=[('all_blocks_sparse',allblocks['fixed_labels'],allblocks['finite_label_lower_bounds']),('all_blocks_sparse_twin',twin['fixed_labels'],twin['finite_label_lower_bounds'])]
 expected_families += [(f'two_high_placement_{k}',[],bounds) for k,bounds in enumerate(variants['distinct_condition_sets'])]
 assert set(families)=={x[0] for x in expected_families}
 for name,fixed,bounds0 in expected_families:
  assert all(m<=6 for a,b,m in fixed) and all(m<=7 for a,b,m in bounds0)
  assert families[name]['fixed']==fixed and families[name]['bounds']==bounds0
  assert families[name]['conditions']==[[index[a,b],'eq',m] for a,b,m in fixed]+[[index[a,b],'ge',m] for a,b,m in bounds0]
 for row in d['sparse_partial_samples']:
  arr,selected=reconstruct(row['path']);assert len(ts)-len(selected)==row['remaining_triangle_count']
  expected=[]
  for i,op,m in families[row['family']]['conditions']:
   if i in arr:possible={arr[i]}
   else:
    pools=[]
    for j,t in enumerate(T):
     if i not in t:continue
     pool={c[t.index(i)] for c in domains[j] if all(k not in arr or arr[k]==x for k,x in zip(t,c))};pools.append(pool)
    assert pools;possible=set.intersection(*pools)
   assert possible and all(x==m if op=='eq' else x>=m for x in possible)
   expected.append(dict(pair=list(pairs[i]),relation=op,bound=m,possible_categories=sorted(possible),already_assigned=i in arr))
  assert expected==row['conditions']
 paths=[SRC,GEOM,'scripts/fresh_h9_t8_sparse_search_20261014_outside.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py',f'scripts/{TAG}.py']
 return dict(status='PASS_INDEPENDENT_T8_DFS_COVERAGE_SOURCE_AND_SMALL_CERTIFICATES',sources={p:sha(p) for p in paths},vertices=462,mandatory_triangles=10,triangle_categories=len(opts),high_category_patterns=high,twin_group_size=16,twin_ordering=chosen,negative_block_clauses=60,independent_sampled_tree_schur_failures=failurechecks,sampled_bridge_prunes=len(d['prune_samples']),sampled_family_prunes=len(d['sparse_partial_samples']),saved_completed_nodes=59423,scope='Independent incidence/category/symmetry/domain and sampled-prune arithmetic, plus full DFS source review. This does not independently repeat all59423 nodes. The complete saved traversal and prior exact metric theorems are explicitly inherited; a separate full traversal/classifier audit is owned by root.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
