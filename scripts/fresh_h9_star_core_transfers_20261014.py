#!/usr/bin/env python3
"""Literal induced-core transfers of the two single-contact theorems."""
import argparse,json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_core_transfers_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
THEOREMS=['results/fresh_h9_star_mixed_singles_20261014.json','results/fresh_h9_star_opposite_core_20261014.json']
B=[(3,4),(6,7),(9,10),(12,13)];P=[2,5,8,11];C=[14,15,16,17]
S=set(sum((list(x)for x in B),[]));CORE=S|{11,14,17}
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
 g=read(GEOM)['representatives'][0]['geometry'];theorems=[read(p)for p in THEOREMS]
 assert [d['status']for d in theorems]==['PASS_CORE_SHARED_SINGLE_ALL_INTERNAL_LABELS','PASS_CORE_OPPOSITE_SINGLE_ALL_INTERNAL_LABELS']
 V=[set(v['incident'])for v in g['vertices']];M=list(map(set,g['minimal_nonfaces']))
 def maximal(core):
  faces={tuple(sorted(v&core))for v in V}
  return sorted(f for f in faces if not any(set(f)<set(g)for g in faces))
 F=maximal(CORE);MF=sorted(tuple(sorted(m))for m in M if m<=CORE)
 assert len(F)==9 and MF==[(3,4,14),(11,12,13),(12,13,17)]
 assert all(any(set(p)<=v for v in V)for p in combinations(CORE,2))
 for d in theorems:
  assert set(d['rank_core'])==CORE and sorted(map(tuple,d['complete_maximal_core_faces']))==F
  assert len(d['families'])==2
  for f in d['families']:
   assert len(f['fixed'])==6 and len(f['bounds'])==4
   assert all(set(row[:2])<=CORE for kind in ('fixed','bounds')for row in f[kind])
   assert sorted(f['bounds'])==[[3,4,4],[6,7,4],[9,10,4],[12,13,7]]
   assert [3,14,3]in f['fixed']and[4,14,3]in f['fixed']
 placements=[];families=[];keys=set()
 for anchor in range(4):
  for owner in range(4):
   if anchor==owner:continue
   third=C[owner]if owner<anchor else P[owner]
   other=[i for i in range(4)if i not in (anchor,owner)];blockmap=[owner,*other,anchor]
   mp={x:y for i,j in enumerate(blockmap)for x,y in zip(B[i],B[j])}
   mp.update({11:P[anchor],17:C[anchor],14:third})
   assert set(mp)==CORE and len(set(mp.values()))==11
   target=set(mp.values());TF=maximal(target)
   mappedF=sorted(tuple(sorted(mp[x]for x in f))for f in F)
   mappedMF=sorted(tuple(sorted(mp[x]for x in m))for m in MF)
   assert mappedF==TF
   assert mappedMF==sorted(tuple(sorted(m))for m in M if m<=target)
   assert all(any(set(p)<=v for v in V)for p in combinations(target,2))
   witnesses=[dict(support=list(f),actual_vertex=next(i for i,v in enumerate(V)if set(f)<=v))for f in mappedF]
   fid=[]
   for path,d in zip(THEOREMS,theorems):
    for f in d['families']:
     image={k:sorted([sorted((mp[a],mp[b]))+[m]for a,b,m in f[k]])for k in ('fixed','bounds')}
     key=tuple(tuple(map(tuple,image[k]))for k in ('fixed','bounds'));assert key not in keys;keys.add(key)
     ident=f'local_anchor{anchor}_owner{owner}_'+f['id'];fid.append(ident)
     families.append(dict(id=ident,origin='reviewed_literal_core_transfer',source_theorem=path,source_family=f['id'],**image))
   placements.append(dict(anchor_block=anchor,owner_block=owner,anchors=[P[anchor],C[anchor]],third_leaf=third,
    core_map=sorted([a,b]for a,b in mp.items()),block_map=blockmap,target_core=sorted(target),
    induced_minimal_nonfaces=[list(m)for m in mappedMF],maximal_face_witnesses=witnesses,family_ids=fid))
 assert len(placements)==12 and len(families)==48
 paths=[GEOM,*THEOREMS,f'scripts/{TAG}.py']
 return dict(status='PASS_TWELVE_LITERAL_CORE_TRANSFERS',case=0,sources={p:sha(p)for p in paths},
  source_core=sorted(CORE),source_maximal_faces=[list(f)for f in F],source_induced_minimal_nonfaces=[list(m)for m in MF],
  placements=placements,families=families,placement_count=12,predicate_count=48,
  proof='Each map is injective on the eleven normal labels, preserves all induced minimal nonfaces and maximal core faces, maps the four blocks with their exact lower bounds, and transports all six fixed own contacts. Every domain witness and residual argument of the core-only theorems therefore applies literally. No full-nerve automorphism is asserted.',
  scope='Conditional on the two bound source theorems,48 label families on star0. All remaining contacts/distances free. Only mapped rank-core hypotheses are used; not a whole-incidence exclusion. A search binder must additionally require independent reviews of both source theorems and these transfers.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['placement_count'],d['predicate_count'])
