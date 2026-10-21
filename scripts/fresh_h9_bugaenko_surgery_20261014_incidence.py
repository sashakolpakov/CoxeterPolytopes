#!/usr/bin/env python3
"""Exact Bugaenko-nerve wedge, then two face truncations; no metric claim."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
import pynauty

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_surgery_20261014'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
OUT=ROOT/f'results/{TAG}_incidence.json'
BASE_BLOCKS=[(0,1,2,3),(4,),(5,6,7,8),(9,),(10,)]
PAIRS={(0,1),(0,3),(1,2),(2,4),(3,4)}
def masks(rows):return tuple(sorted(sum(1<<i for i in row)for row in rows))
def tup(mask):return tuple(i for i in range(mask.bit_length())if mask>>i&1)
def minimal(rows):
 out=[]
 for x in sorted(set(rows),key=lambda v:(v.bit_count(),v)):
  if not any(y&x==y for y in out):out.append(x)
 return tuple(sorted(out))
def stellar(M,V,F,new):
 assert all(m&F!=m for m in M)
 A=minimal([F,*[m for m in M if m&F!=F],*[(1<<new)|(m&~F)for m in M]])
 B=set()
 for v in V:
  if v&F==F:
   for i in tup(F):B.add((v&~(1<<i))|(1<<new))
  else:B.add(v)
 assert all(v.bit_count()==9 for v in B)
 return A,tuple(sorted(B))
def symmetries(weights):
 return [p for p in permutations(range(5))if all(weights[i]==weights[p[i]]for i in range(5))and{tuple(sorted((p[a],p[b])))for a,b in PAIRS}==PAIRS]
def remap_counts(counts,p):
 out=[None]*5
 for i,j in enumerate(p):out[j]=counts[i]
 return tuple(out)
def first_faces(blocks,M):
 H=symmetries(list(map(len,blocks)))
 for cc in product(*(range(len(b)+1)for b in blocks)):
  if not 2<=sum(cc)<=5 or cc!=min(remap_counts(cc,p)for p in H):continue
  F=sum(1<<x for b,k in zip(blocks,cc)for x in b[:k])
  if all(m&F!=m for m in M):yield F,cc,[p for p in H if remap_counts(cc,p)==cc]
def second_faces(blocks,F,H,M):
 parts=[(tuple(i for i in b if F>>i&1),tuple(i for i in b if not(F>>i&1)))for b in blocks]
 options=[list(product(range(len(a)+1),range(len(b)+1)))for a,b in parts]
 for cc in product(*options):
  k=sum(sum(z)for z in cc)
  for new in[0,1]:
   if not 2<=k+new<=5 or cc!=min(remap_counts(cc,p)for p in H):continue
   S=(1<<12)if new else 0
   for (a,b),(ca,cb)in zip(parts,cc):
    for i in a[:ca]+b[:cb]:S|=1<<i
   if all(m&S!=m for m in M):yield S,cc,new
def canonical(M,n):
 mm=list(M);adj={i:[]for i in range(n+len(mm))}
 for j,m in enumerate(mm,n):
  for i in tup(m):adj[i].append(j);adj[j].append(i)
 g=pynauty.Graph(n+len(mm),adjacency_dict=adj,vertex_coloring=[set(range(n)),set(range(n,n+len(mm)))])
 order=pynauty.canon_label(g)
 # Nauty's output lists old vertices in canonical order. Verify the colored
 # partition and use only the resulting explicitly relabelled clutter.
 assert set(order[:n])==set(range(n))
 mp={old:new for new,old in enumerate(order[:n])}
 C=tuple(sorted(sum(1<<mp[i]for i in tup(m))for m in mm))
 return C,tuple(order[:n])
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
 src=json.loads((ROOT/SOURCE).read_text());oldV=masks(src['geometry']['vertices'])
 BM=masks([BASE_BLOCKS[a]+BASE_BLOCKS[b]for a,b in PAIRS])
 assert oldV==tuple(sorted(sum(1<<i for i in a)for a in combinations(range(11),8)if all(sum(1<<i for i in a)&m!=m for m in BM)))
 assert len(oldV)==41
 allgroup=symmetries(list(map(len,BASE_BLOCKS)));orbits=[];left=set(range(5))
 while left:
  a=min(left);orb={p[a]for p in allgroup};orbits.append(sorted(orb));left-=orb
 assert len(orbits)==3
 unique={};recipes=[];stats=Counter();wedges=[]
 for orbit in orbits:
  wb=orbit[0];carrier=BASE_BLOCKS[wb][0];blocks=[b+((11,)if i==wb else())for i,b in enumerate(BASE_BLOCKS)]
  M=masks([blocks[a]+blocks[b]for a,b in PAIRS]);V=[]
  for v in oldV:
   if v>>carrier&1:V.append(v|(1<<11))
   else:V.extend([v|(1<<carrier),v|(1<<11)])
  V=tuple(sorted(V));assert len(V)==len(set(V))
  assert V==tuple(sorted(sum(1<<i for i in a)for a in combinations(range(12),9)if all(sum(1<<i for i in a)&m!=m for m in M)))
  first_count=second_count=eligible=0
  for F,fc,H in first_faces(blocks,M):
   first_count+=1;A,B=stellar(M,V,F,12)
   for S,sc,new in second_faces(blocks,F,H,A):
    second_count+=1;stats['ordered_recipe_orbits']+=1
    bad=[m for m in A if m.bit_count()>5]
    if any(m&S!=S for m in bad):stats['old_large_nonface_survives']+=1;continue
    C,W=stellar(A,B,S,13)
    if max(map(int.bit_count,C))>5:stats['new_large_nonface']+=1;continue
    eligible+=1;stats['passes_missing_order']+=1
    canon,order=canonical(C,14);mapping={old:new for new,old in enumerate(order)}
    CW=tuple(sorted(sum(1<<mapping[i]for i in tup(v))for v in W))
    if canon in unique:
     uid=unique[canon];assert tuple(recipes[uid]['canonical_vertices'])==CW
     recipes[uid]['recipe_orbit_multiplicity']+=1;continue
    uid=len(recipes);unique[canon]=uid
    recipes.append(dict(case=uid,wedge_block=wb,wedge_carrier=carrier,first_face=list(tup(F)),second_face=list(tup(S)),first_block_counts=list(fc),second_block_counts=[list(a)for a in sc],second_uses_new_facet=bool(new),minimal_nonfaces=list(C),actual_vertices=list(W),canonical_order=list(order),canonical_minimal_nonfaces=list(canon),canonical_vertices=list(CW),recipe_orbit_multiplicity=1))
  wedges.append(dict(wedge_block=wb,carrier_orbit_blocks=orbit,blocks=[list(b)for b in blocks],vertex_count=len(V),first_face_orbits=first_count,second_recipe_orbits=second_count,passes_missing_order=eligible))
  print('wedge',wb,'first',first_count,'second',second_count,'eligible',eligible,'unique_so_far',len(recipes),flush=True)
 for r in recipes:
  V=r['actual_vertices'];assert set(Counter(v&~(1<<i)for v in V for i in tup(v)).values())=={2}
  assert all(v.bit_count()==9 for v in V)
 files=[SOURCE,f'scripts/{TAG}_incidence.py']
 return dict(status='EXACT_BUGAENKO_WEDGE_TWO_TRUNCATION_INCIDENCE_CENSUS',sources={p:sha(p)for p in files},dimension=9,facets=14,base_minimal_nonfaces=list(BM),base_vertices=list(oldV),wedge_orbits=wedges,counts=dict(stats),distinct_admissible_incidence_types=len(recipes),cases=recipes,scope='All three combinatorial wedge-carrier orbits, then ordered truncations of every face of codimension2..5, including newly introduced facets. Within-block and block-graph stabilizers give safe operation orbits; colored incidence canonicalization merges final nerves. Passing records are explicitly polytopal by wedge and sufficiently small face truncations, not Coxeter metric realizations.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'],out['distinct_admissible_incidence_types'])
