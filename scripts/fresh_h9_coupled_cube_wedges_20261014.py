#!/usr/bin/env python3
"""Exact five-cube ridge cuts, four negative-facet wedges, necessary labels.

Missing pairs keep arbitrary strictly dotted distances. Category7 represents
all finite integer labels >=7 only in triangle and spherical tests.
"""
import argparse,json
from fractions import Fraction as F
from itertools import combinations,product
from collections import Counter
from pathlib import Path
from hashlib import sha256
from fresh_h9_highbase_construction_20261014_cell import join
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_coupled_cube_wedges_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 return x
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def rank(A):
 A=[list(map(F,r))for r in A];k=0
 for j in range(len(A[0])):
  p=next((i for i in range(k,len(A))if A[i][j]),None)
  if p is None:continue
  A[k],A[p]=A[p],A[k];q=A[k][j];A[k]=[x/q for x in A[k]]
  for i in range(k+1,len(A)):
   q=A[i][j];A[i]=[x-q*y for x,y in zip(A[i],A[k])]
  k+=1
  if k==len(A):break
 return k
def active(x,rows,b):
 a=tuple(i for i,(r,h)in enumerate(zip(rows,b))if dot(r,x)==h)
 assert all(dot(r,x)<=h for r,h in zip(rows,b))
 assert len(a)==len(x)and rank([rows[i]for i in a])==len(x)
 return a
def minimal(V,n):
 faces={0}
 for v in V:
  mask=sum(1<<i for i in v);s=mask
  while s:
   faces.add(s);s=(s-1)&mask
 out=[]
 for m in range(1,1<<n):
  if m in faces:continue
  if all(m^(1<<i)in faces for i in range(n)if m>>i&1):out.append(tuple(i for i in range(n)if m>>i&1))
 return sorted(out)
def geometry():
 rows=[];b=[];names=[]
 for i in range(5):
  for s,label in [(1,'positive'),(-1,'negative')]:
   r=[F(0)]*5;r[i]=s;rows.append(r);b.append(F(1));names.append(f'{label}{i}')
 V={tuple(map(F,x)):active(x,rows,b)for x in product([-1,1],repeat=5)};steps=[]
 for k in range(4):
  support=(2*k,2*(k+1));r=[F(0)]*5;r[k]=r[k+1]=1
  face={x for x,a in V.items()if set(support)<=set(a)};assert face
  eps=min(2-dot(r,x)for x in V if x not in face)/2;level=2-eps
  assert eps>0 and all(dot(r,x)>level for x in face)
  oldV=V;E=[(x,y)for x,y in combinations(V,2)if len(set(V[x])&set(V[y]))==4]
  assert len(E)*2==len(V)*5
  kept={x for x in V if x not in face};cross=[]
  for x,y in E:
   sx=dot(r,x)-level;sy=dot(r,y)-level
   assert sx and sy
   if sx*sy<0:
    t=-sx/(sy-sx);z=tuple(a+t*(d-a)for a,d in zip(x,y));kept.add(z)
    cross.append(dict(old_edge=[list(oldV[x]),list(oldV[y])],fraction=t,point=z))
  rows.append(r);b.append(level);names.append(f'cut{k}_{k+1}')
  V={x:active(x,rows,b)for x in sorted(kept)}
  assert set(Counter(a for v in V.values()for a in combinations(v,4)).values())=={2}
  steps.append(dict(support=support,epsilon=eps,level=level,old_vertices=len(oldV),removed_vertices=len(face),new_vertices=len(V),edge_intersections=cross))
 assert len(rows)==14
 for i,r in enumerate(rows):
  verts=[x for x,a in V.items()if i in a]
  assert rank([tuple(x[j]-verts[0][j]for j in range(5))for x in verts[1:]])==4
 MB=minimal(V.values(),14);assert all(len(m)==2 for m in MB)
 chosen=[1,3,5,7];assert any(set(chosen)<=set(a)for a in V.values())
 blocks=[];rr=[];hh=[];nn=[]
 for i,(r,h,name)in enumerate(zip(rows,b,names)):
  block=[]
  for sg in([-1,1]if i in chosen else[0]):
   block.append(len(rr));q=list(r)+[F(0)]*4
   if i in chosen:q[5+chosen.index(i)]=sg
   rr.append(q);hh.append(h);nn.append(name+('minus'if sg<0 else'plus'if sg>0 else''))
  blocks.append(tuple(block))
 vv={}
 for x,a in V.items():
  outside=[i for i in chosen if i not in a]
  for signs in product([-1,1],repeat=len(outside)):
   z=[F(0)]*4
   for i,sg in zip(outside,signs):z[chosen.index(i)]=sg*(1+x[i//2])
   point=x+tuple(z);inc=active(point,rr,hh);assert inc not in vv;vv[inc]=point
 assert len(rr)==18 and all(h>0 for h in hh)
 M=sorted(tuple(sorted(t for i in m for t in blocks[i]))for m in MB)
 assert minimal(vv,18)==M
 assert set(Counter(a for v in vv for a in combinations(v,8)).values())=={2}
 assert all(rank([tuple(x[j]-pts[0][j]for j in range(9))for x in pts[1:]])==8 for i in range(18)for pts in[[x for a,x in vv.items()if i in a]])
 left=set(range(18));components=[]
 while left:
  hit={min(left)}
  while True:
   grow=hit|{i for m in M if hit.intersection(m)for i in m}
   if grow==hit:break
   hit=grow
  components.append(sorted(hit));left-=hit
 assert len(components)==1
 return dict(base_dimension=5,base_normals=rows,base_offsets=b,base_names=names,base_vertices=[dict(point=x,incident=a)for x,a in V.items()],
  truncations=steps,base_minimal_nonfaces=MB,wedged_facets=chosen,blocks=blocks,dimension=9,facet_count=18,
  normals=rr,offsets=hh,names=nn,vertices=[dict(point=x,incident=a)for a,x in sorted(vv.items())],minimal_nonfaces=M,
  missing_size_counts=dict(Counter(map(len,M))),missing_hypergraph_components=components)
def run(maxstates):
 g=geometry();print('actual rational geometry',len(g['base_vertices']),len(g['vertices']),g['missing_size_counts'],flush=True)
 label=join(g['minimal_nonfaces'],[v['incident']for v in g['vertices']],maxstates)
 paths=[f'scripts/{TAG}.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return enc(dict(status='EXACT_COUPLED_FIVE_CUBE_WEDGE_CONSTRUCTION_CHECKPOINT',sources={p:sha(p)for p in paths},
  settings=dict(max_join_states=maxstates),geometry=g,label_search=label,
  scope='One explicit simple rational H9 incidence, obtained by four small consecutive positive ridge cuts of the five-cube followed by four distinct mutually meeting negative-facet wedges. The connected minimal-nonface hypergraph excludes a nontrivial Cartesian product. Missing pairs have unspecified continuous dotted distances; only mandatory Lanner triangles and spherical actual-vertex restrictions are joined. Category7 includes every finite integer label>=7 and is never treated as a fixed Gram entry. A survivor or truncation remains open for all omitted finite contacts and global rank/signature/compactness. No fifteen-facet-specific restrictions are used.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--max-states',type=int,default=20000);p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run(a.max_states)));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['label_search']['verdict'],flush=True)
