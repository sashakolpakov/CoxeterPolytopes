#!/usr/bin/env python3
"""Separate Fraction primal/Gram replay; no producer/project imports."""
import argparse,json
from fractions import Fraction as F
from itertools import combinations,product
from collections import Counter
from pathlib import Path
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_coupled_cube_wedges_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d['sources'].items():assert sha(f)==h,f
 return d
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def solve(A,b):
 n=len(A);M=[list(a)+[c]for a,c in zip(A,b)]
 for k in range(n):
  p=next((j for j in range(k,n)if M[j][k]),None)
  if p is None:return None
  M[k],M[p]=M[p],M[k];q=M[k][k];M[k]=[v/q for v in M[k]]
  for i in range(n):
   if i!=k:
    q=M[i][k];M[i]=[a-q*b for a,b in zip(M[i],M[k])]
 return tuple(r[-1]for r in M)
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def mul(x,y):return(x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def inv(x):
 q=x[0]*x[0]-2*x[1]*x[1];assert q;return(x[0]/q,-x[1]/q)
def sign(x):
 a,b=x;s=lambda q:int(q>0)-int(q<0)
 if not b:return s(a)
 if not a:return s(b)
 if s(a)==s(b):return s(a)
 return s(a)*s(a*a-2*b*b)
Z=(F(0),F(0));O=(F(1),F(0));COS={2:Z,3:(F(1,2),F(0)),4:(F(0),F(1,2))}
def ldl(A):
 A=[r[:]for r in A];out=[]
 for k in range(len(A)):
  q=A[k][k];assert q!=Z;out.append(q)
  for i in range(k+1,len(A)):
   for j in range(k+1,len(A)):A[i][j]=add(A[i][j],neg(mul(mul(A[i][k],A[k][j]),inv(q))))
 return out
def run():
 source=f'results/{TAG}.json';csp=f'results/{TAG}_csp.json';p=read(source);s=read(csp);g=p['geometry']
 B=[list(map(F,r))for r in g['base_normals']];h=list(map(F,g['base_offsets']));old={}
 trials=0
 for S in combinations(range(14),5):
  trials+=1;x=solve([B[i]for i in S],[h[i]for i in S])
  if x is None or any(dot(r,x)>a for r,a in zip(B,h)):continue
  inc=tuple(i for i,(r,a)in enumerate(zip(B,h))if dot(r,x)==a);assert len(inc)==5
  old[inc]=x
 assert len(old)==70 and old=={tuple(v['incident']):tuple(map(F,v['point']))for v in g['base_vertices']}
 assert B[:10]==[[F(s if i==k else 0)for k in range(5)]for i in range(5)for s in[1,-1]]
 assert all(h[i]==1 for i in range(10))and all(a>0 for a in h)
 assert B[10:]==[[F(int(k in(i,i+1)))for k in range(5)]for i in range(4)]
 selected=g['wedged_facets'];assert selected==[1,3,5,7]
 rows=[list(map(F,r))for r in g['normals']];off=list(map(F,g['offsets']));blocks=g['blocks']
 for i,bl in enumerate(blocks):
  assert len(bl)==(2 if i in selected else 1)
  for k,j in enumerate(bl):
   r=B[i]+[F(0)]*4
   if i in selected:r[5+selected.index(i)]=[-1,1][k]
   assert rows[j]==r and off[j]==h[i]
 vertices={}
 for inc,x in old.items():
  un=[i for i in selected if i not in inc]
  for ss in product([-1,1],repeat=len(un)):
   y=list(x)+[F(0)]*4
   for i,t in zip(un,ss):y[5+selected.index(i)]=t*(1+x[i//2])
   a=tuple(j for j,(r,b)in enumerate(zip(rows,off))if dot(r,y)==b)
   assert all(dot(r,y)<=b for r,b in zip(rows,off))and len(a)==9
   assert solve([rows[j]for j in a],[off[j]for j in a])==tuple(y)
   assert a not in vertices;vertices[a]=tuple(y)
 assert len(vertices)==512 and vertices=={tuple(z['incident']):tuple(map(F,z['point']))for z in g['vertices']}
 faces={0}
 for a in vertices:
  v=sum(1<<i for i in a);q=v
  while q:faces.add(q);q=(q-1)&v
 missing=[]
 for m in range(1,1<<18):
  if m in faces:continue
  if all(m^(1<<i)in faces for i in range(18)if m>>i&1):missing.append(tuple(i for i in range(18)if m>>i&1))
 assert sorted(missing)==sorted(map(tuple,g['minimal_nonfaces']))
 assert Counter(map(len,missing))=={2:9,3:11}
 assert set(Counter(e for v in vertices for e in combinations(v,8)).values())=={2}
 connected={0}
 while True:
  nxt=connected|{i for m in missing if connected.intersection(m)for i in m}
  if nxt==connected:break
  connected=nxt
 assert connected==set(range(18))
 assert s['status']=='FULL_MANDATORY_LABEL_SURVIVOR'
 labels={(a,b):m for a,b,m in s['solution']['pair_labels']}
 assert len(labels)==26 and set(labels)=={p for t in missing if len(t)==3 for p in combinations(t,2)}
 def gram(v):return[[O if i==j else neg(COS[labels.get(tuple(sorted((i,j))),2)])for j in v]for i in v]
 triangle=[]
 for m in missing:
  if len(m)!=3:continue
  diag=ldl(gram(m));assert list(map(sign,diag))==[1,1,-1]
  assert all(all(sign(q)>0 for q in ldl(gram(e)))for e in combinations(m,2))
  triangle.append(dict(support=list(m),inertia=[2,1],pivots=[[str(a),str(b)]for a,b in diag]))
 for v in vertices:assert all(sign(q)>0 for q in ldl(gram(v)))
 paths=[source,csp,f'scripts/{TAG}_independent.py']
 return dict(status='PASS_SEPARATE_PRIMAL_AND_GRAM_REPLAY',sources={p:sha(p)for p in paths},
  base_supports_tested=trials,base_vertices=70,wedged_vertices=512,rational_inequalities_checked=512*18,
  independent_vertex_solves=512,minimal_nonface_counts={'2':9,'3':11},nonproduct=True,
  exact_lanner_triangles=triangle,exact_spherical_grams=512,unassigned_dotted_pairs=9,
  scope='Separate Fraction elimination exhausts every five-row base support, then the literal wedge inequalities independently reconstruct every higher-dimensional vertex. Full downward face comparison verifies the minimal nonfaces and connected hypergraph. Separate exact Q(sqrt2) LDL checks all11 mandatory triangles and512 relaxed vertex Grams. This validates the saved necessary assignment, not global Gram rank/signature, actual hyperbolic realization or completeness of label search.')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');a=a.parse_args();d=run();p=ROOT/f'results/{TAG}_independent.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['base_vertices'],d['wedged_vertices'])
