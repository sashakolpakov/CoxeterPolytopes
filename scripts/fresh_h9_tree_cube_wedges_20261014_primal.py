#!/usr/bin/env python3
"""Separate FLINT rational primal checks of every saved tree-wedge incidence.

Does not import the geometry producer. Completeness uses the closed vertex
edge graph of the explicitly bounded full-dimensional polytope.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from flint import fmpq,fmpq_mat
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_tree_cube_wedges_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 source=f'results/{TAG}.json';data=json.loads((ROOT/source).read_text())
 for p,h in data['sources'].items():assert sha(p)==h,p
 records=[]
 for r in data['representatives']:
  g=r['geometry'];N=fmpq_mat([[fmpq(x)for x in row]for row in g['normals']]);h=[fmpq(x)for x in g['offsets']]
  X=fmpq_mat([[fmpq(x)for x in v['point']]for v in g['vertices']]);V=[tuple(v['incident'])for v in g['vertices']]
  assert N.nrows()==18 and N.ncols()==9 and len(set(V))==len(V)
  assert all(x>0 for x in h) # Origin is strict and full-dimensional.
  B=X*N.transpose();assert X.nrows()==len(V)
  for k,v in enumerate(V):
   assert len(v)==9 and all(B[k,j]<=h[j]for j in range(18))
   assert tuple(j for j in range(18)if B[k,j]==h[j])==v
   assert fmpq_mat([[N[i,j]for j in range(9)]for i in v]).rank()==9
  # Boundedness is certified directly from the original cube rows after
  # wedging: positive row e_i remains; a negative row remains or is split
  # as -e_i +/- e_(5+j). Average the split rows to recover -e_i.
  for i in range(5):
   pos=g['blocks'][2*i];neg=g['blocks'][2*i+1];assert len(pos)==1
   expected=[int(i==j)for j in range(9)]
   assert [N[pos[0],j]for j in range(9)]==expected and h[pos[0]]==1
   assert [sum(N[k,j]for k in neg)/len(neg)for j in range(9)]==[-x for x in expected]
   assert all(h[k]==1 for k in neg)
  for j in range(5,9):
   rows=[k for k in range(18)if N[k,j]];assert len(rows)==2
   assert sorted(N[k,j]for k in rows)==[-1,1]
   assert all(sum(abs(N[k,l])for l in range(5))==1 and h[k]==1 for k in rows)
   # Thus |x_i|<=1 and |z_j|<=1+x_i<=2.
  edge_vertices={}
  for k,v in enumerate(V):
   for edge in combinations(v,8):edge_vertices.setdefault(edge,[]).append(k)
  assert set(map(len,edge_vertices.values()))=={2}
  adjacency=[set()for v in V]
  for edge,(i,j)in edge_vertices.items():
   assert i!=j;adjacency[i].add(j);adjacency[j].add(i)
  assert set(map(len,adjacency))=={9}
  reached={0};todo=[0]
  for i in todo:
   for j in adjacency[i]-reached:reached.add(j);todo.append(j)
  assert len(reached)==len(V)
  facet_counts=[]
  for j in range(18):
   ids=[k for k,v in enumerate(V)if j in v];assert ids
   assert fmpq_mat([[X[k,l]-X[ids[0],l]for l in range(9)]for k in ids[1:]]).rank()==8
   facet_counts.append(len(ids))
  # Rebuild all faces from the complete vertex list, then check that the
  # advertised missing faces account for every remaining facet subset.
  faces={0}
  for v in V:
   m=sum(1<<i for i in v);s=m
   while s:faces.add(s);s=(s-1)&m
  M=[sum(1<<i for i in m)for m in g['minimal_nonfaces']]
  for m in M:
   assert m not in faces
   assert all(m^(1<<i)in faces for i in range(18)if m>>i&1)
  for m in range(1<<18):
   if m not in faces:assert any(m&q==q for q in M)
  records.append(dict(case=r['id'],vertices=len(V),edges=len(edge_vertices),facets=18,
   facet_vertex_counts=facet_counts,all_face_subsets=len(faces),missing_faces=len(M)))
  print('primal',r['id'],len(V),len(edge_vertices),flush=True)
 paths=[source,f'scripts/{TAG}.py',f'scripts/{TAG}_primal.py']
 return dict(status='PASS_SEPARATE_EXACT_RATIONAL_PRIMAL_REPLAY',sources={p:sha(p)for p in paths},
  cases=records,total_vertices=sum(r['vertices']for r in records),total_edges=sum(r['edges']for r in records),
  scope='Separate FLINT rational matrix checks: all saved vertex inequalities and exact active sets, independent rank9, all18 actual facets, origin interior, explicit bounded coordinate box, connected closed complete vertex-edge graph, and every missing face. A bounded polytope has connected graph; every edge at each saved simple vertex has its other endpoint saved, so no unlisted vertex exists. The 240-recipe canonicalization is not independently rerun here.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}_primal.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['total_vertices'],flush=True)
