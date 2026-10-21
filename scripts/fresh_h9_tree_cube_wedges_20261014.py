#!/usr/bin/env python3
"""All ordered ridge cuts along the two nonpath five-vertex trees.

Combinatorial stellar cuts enumerate 240 recipes. Every identification is
checked as an explicit incidence isomorphism. Each representative receives
an exact rational bounded simple realization, with four simultaneous wedges.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
import pynauty
from fresh_h9_coupled_cube_wedges_20261014 import active,dot,enc,rank,minimal
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_tree_cube_wedges_20261014'
TREES={'star':((0,1),(0,2),(0,3),(0,4)),
       'T':((0,1),(0,2),(0,3),(3,4))}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def digest(x):return sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def combinatorics(order,unwedged):
 V={tuple(2*i+s for i,s in enumerate(ss))for ss in product([0,1],repeat=5)}
 for k,(i,j)in enumerate(order):
  a,b=2*i,2*j;new=set()
  for v in V:
   if a in v and b in v:
    for r in(a,b):new.add(tuple(sorted((set(v)-{r})|{10+k})))
   else:new.add(v)
  V=new
  assert set(Counter(s for v in V for s in combinations(v,4)).values())=={2}
 MB=minimal(V,14);assert all(len(m)==2 for m in MB)
 chosen=[2*i+1 for i in range(5)if i!=unwedged]
 blocks=[];n=0
 for i in range(14):
  w=2 if i in chosen else 1;blocks.append(tuple(range(n,n+w)));n+=w
 W=set()
 for v in V:
  fixed=tuple(j for i in v for j in blocks[i]);out=[i for i in chosen if i not in v]
  for choice in product(*(blocks[i]for i in out)):W.add(tuple(sorted(fixed+choice)))
 M=sorted(tuple(sorted(j for i in m for j in blocks[i]))for m in MB)
 assert n==18 and set(map(len,W))=={9}
 assert set(Counter(s for v in W for s in combinations(v,8)).values())=={2}
 # A flag base and the wedge formula give this complete minimal clutter.
 assert all(not set(m)<=set(v)for m in M for v in W)
 return sorted(V),MB,chosen,blocks,sorted(W),M
def canonical(M):
 n=18;adj={i:[]for i in range(n+len(M))};colors=[set(range(n))]
 for size in sorted(set(map(len,M))):colors.append({n+k for k,m in enumerate(M)if len(m)==size})
 for k,m in enumerate(M):
  for i in m:adj[i].append(n+k);adj[n+k].append(i)
 g=pynauty.Graph(len(adj),adjacency_dict=adj,vertex_coloring=colors)
 order=[i for i in pynauty.canon_label(g)if i<n];pos={i:k for k,i in enumerate(order)}
 assert sorted(order)==list(range(n))
 key=tuple(sorted(tuple(sorted(pos[i]for i in m))for m in M))
 return key,order
def geometry(order,unwedged):
 rows=[];b=[];names=[]
 for i in range(5):
  for sg,label in[(1,'positive'),(-1,'negative')]:
   r=[F(0)]*5;r[i]=sg;rows.append(r);b.append(F(1));names.append(f'{label}{i}')
 V={tuple(map(F,x)):active(x,rows,b)for x in product([-1,1],repeat=5)};steps=[]
 for k,(i,j)in enumerate(order):
  support=(2*i,2*j);r=[F(0)]*5;r[i]=r[j]=1
  face={x for x,a in V.items()if set(support)<=set(a)};assert face
  eps=min(2-dot(r,x)for x in V if x not in face)/2;level=2-eps
  assert eps>0
  old=V;edges=[(x,y)for x,y in combinations(V,2)if len(set(V[x])&set(V[y]))==4]
  assert len(edges)*2==len(V)*5
  kept={x for x in V if x not in face};cross=[]
  for x,y in edges:
   sx=dot(r,x)-level;sy=dot(r,y)-level;assert sx and sy
   if sx*sy<0:
    t=-sx/(sy-sx);z=tuple(a+t*(d-a)for a,d in zip(x,y));kept.add(z)
    cross.append(dict(old_edge=[old[x],old[y]],fraction=t,point=z))
  rows.append(r);b.append(level);names.append(f'cut{i}_{j}')
  V={x:active(x,rows,b)for x in sorted(kept)}
  steps.append(dict(support=support,epsilon=eps,level=level,old_vertices=len(old),removed_vertices=len(face),new_vertices=len(V),edge_intersections=cross))
 bv,MB,chosen,blocks,W,M=combinatorics(order,unwedged)
 assert sorted(V.values())==bv and minimal(V.values(),14)==MB
 assert any(set(chosen)<=set(a)for a in V.values())
 rr=[];hh=[];nn=[]
 for i,(r,h,name)in enumerate(zip(rows,b,names)):
  for sg in([-1,1]if i in chosen else[0]):
   q=list(r)+[F(0)]*4
   if i in chosen:q[5+chosen.index(i)]=sg
   rr.append(q);hh.append(h);nn.append(name+('minus'if sg<0 else'plus'if sg>0 else''))
 vv={}
 for x,a in V.items():
  outside=[i for i in chosen if i not in a]
  for signs in product([-1,1],repeat=len(outside)):
   z=[F(0)]*4
   for i,sg in zip(outside,signs):z[chosen.index(i)]=sg*(1+x[i//2])
   point=x+tuple(z);inc=active(point,rr,hh);assert inc not in vv;vv[inc]=point
 assert sorted(vv)==W and minimal(vv,18)==M
 assert all(h>0 for h in hh)
 for i in range(18):
  pts=[x for a,x in vv.items()if i in a]
  assert rank([tuple(x[j]-pts[0][j]for j in range(9))for x in pts[1:]])==8
 hit={0}
 while True:
  grow=hit|{i for m in M if hit.intersection(m)for i in m}
  if grow==hit:break
  hit=grow
 assert hit==set(range(18))
 return dict(base_dimension=5,base_normals=rows,base_offsets=b,base_names=names,
  base_vertices=[dict(point=x,incident=a)for x,a in V.items()],truncations=steps,
  base_minimal_nonfaces=MB,wedged_facets=chosen,blocks=blocks,dimension=9,facet_count=18,
  normals=rr,offsets=hh,names=nn,vertices=[dict(point=x,incident=a)for a,x in sorted(vv.items())],
  minimal_nonfaces=M,missing_size_counts=dict(Counter(map(len,M))),missing_hypergraph_components=[sorted(hit)])
def run():
 oldpath='results/fresh_h9_coupled_cube_wedges_20261014.json';old=json.loads((ROOT/oldpath).read_text())
 for p,h in old['sources'].items():assert sha(p)==h
 oldkey,_=canonical(old['geometry']['minimal_nonfaces'])
 reps=[];bykey={};recipes=[]
 for tree,edges in TREES.items():
  for order in permutations(edges):
   for unwedged in range(5):
    _,_,_,_,V,M=combinatorics(order,unwedged);key,perm=canonical(M)
    if key not in bykey:
     k=len(reps);bykey[key]=k
     reps.append(dict(id=k,tree=tree,cut_order=order,unwedged_coordinate=unwedged,
      canonical_facet_order=perm,canonical_clutter=key,vertices=V,minimal_nonfaces=M,
      isomorphic_to_frozen_path=key==oldkey))
    k=bykey[key];rep=reps[k];to_rep=[0]*18
    for i,j in zip(perm,rep['canonical_facet_order']):to_rep[i]=j
    assert sorted(to_rep)==list(range(18))
    assert sorted(tuple(sorted(to_rep[i]for i in v))for v in V)==rep['vertices']
    assert sorted(tuple(sorted(to_rep[i]for i in m))for m in M)==rep['minimal_nonfaces']
    recipes.append(dict(tree=tree,cut_order=order,unwedged_coordinate=unwedged,
     representative=k,facet_isomorphism=to_rep,vertex_count=len(V)))
 print('deduplicated',len(recipes),'recipes to',len(reps),'incidences',flush=True)
 for r in reps:
  g=geometry(r['cut_order'],r['unwedged_coordinate'])
  assert [v['incident']for v in g['vertices']]==r.pop('vertices')
  assert g['minimal_nonfaces']==r.pop('minimal_nonfaces');r['geometry']=g
  r['recipe_multiplicity']=sum(z['representative']==r['id']for z in recipes)
  print('realized',r['id'],r['tree'],r['unwedged_coordinate'],len(g['vertices']),g['missing_size_counts'],flush=True)
 paths=[f'scripts/{TAG}.py','scripts/fresh_h9_coupled_cube_wedges_20261014.py',oldpath,
  'scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 assert len(recipes)==240
 return enc(dict(status='EXACT_REALIZABLE_NONPATH_TREE_WEDGE_FAMILY',sources={p:sha(p)for p in paths},
  trees=TREES,recipe_count=len(recipes),representative_count=len(reps),recipes=recipes,representatives=reps,
  scope='All 24 cut orders on each of the star and T trees on five positive cube facets, and all five choices of a single unwedged negative facet. Exact small cuts remove only their original ridge vertices; each of four negative facets is then wedged once. Canonicalization is used only to propose identifications, all of which are checked as literal bijections of both complete vertex sets and minimal nonfaces. Each representative has a complete bounded rational simple 9-polytope realization with 18 actual facets and connected missing-face hypergraph. These are incidence candidates, not hyperbolic metric realizations.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['representative_count'],flush=True)
