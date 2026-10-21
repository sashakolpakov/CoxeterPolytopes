#!/usr/bin/env python3
"""Separate sparse-radical replay of both affine E7 cusp fillers.

No producer or project imports. Fraction elimination and rational192bit
enclosures replace the primary recursive quadratic sign implementation.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_intrinsic_caps_20261014_cusps'
P=f'results/{TAG}.json';BASE='results/fresh_h20_e7_intrinsic_caps_20261014.json'
RAD=[1,2,79,158,11,22,869,1738];Z=(F(0),)*8;DEN=1<<192
BOUND=[]
for q in RAD:
 n=isqrt(q*DEN*DEN);BOUND.append((F(n,DEN),F(n+int(n*n!=q*DEN*DEN),DEN)))
def field(v):return tuple(map(F,v))
def vector(v):return list(map(field,v))
def c(x):return(F(x),)+(F(0),)*7
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def scale(x,s):return tuple(a*s for a in x)
def sub(x,y):return add(x,scale(y,-1))
def mul(x,y):
 z=[F(0)]*8
 for i,a in enumerate(x):
  for j,b in enumerate(y):
   if a and b:z[i^j]+=a*b*RAD[i&j]
 return tuple(z)
def total(values):
 v=Z
 for x in values:v=add(v,x)
 return v
def sign(x):
 if x==Z:return 0
 lo=sum(a*(u if a>=0 else v)for a,(u,v)in zip(x,BOUND));hi=sum(a*(v if a>=0 else u)for a,(u,v)in zip(x,BOUND))
 assert lo>0 or hi<0,('unresolved192bit',x,lo,hi)
 return 1 if lo>0 else -1
def inv(A):
 n=len(A);B=[list(map(F,row))+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
 for k in range(n):
  p=next(i for i in range(k,n)if B[i][k]);B[k],B[p]=B[p],B[k];a=B[k][k];B[k]=[x/a for x in B[k]]
  for i in range(n):
   if i!=k:
    a=B[i][k];B[i]=[x-a*y for x,y in zip(B[i],B[k])]
 return[r[n:]for r in B]
def pd(A):
 B=[list(map(F,r))for r in A]
 for k in range(len(B)):
  p=B[k][k];assert p>0
  for i in range(k+1,len(B)):
   for j in range(k+1,len(B)):B[i][j]-=B[i][k]*B[j][k]/p
def rational(x):assert x[1:]==(0,)*7;return x[0]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 return x
def fdiv(a,b):
 M=[[mul(b,c(1) if j==0 else tuple(F(k==j)for k in range(8)))[i]for j in range(8)]for i in range(8)]
 I=inv(M);return tuple(sum(I[i][j]*a[j]for j in range(8))for i in range(8))
def lin(vectors,coeff):
 return [total(mul(a,v[k])for a,v in zip(coeff,vectors))for k in range(21)]
def run():
 from collections import Counter
 d,p=map(read,[BASE,P]);Q=[list(map(F,r))for r in d['form']];cache={};dc={}
 sparse=[[(j,x)for j,x in enumerate(row)if x]for row in Q]
 def dot(v,w):
  key=(tuple(v),tuple(w))
  if key in dc:return dc[key]
  if key[1]not in cache:cache[key[1]]=[tuple(sum(a*w[j][k]for j,a in row)for k in range(8))for row in sparse]
  out=total(mul(a,b)for a,b in zip(v,cache[key[1]]));dc[key]=out;dc[(key[1],key[0])]=out;return out
 old=list(map(vector,d['normals']));H=vector(d['strict_center']);E=vector(d['strict_escape']['point'])
 assert len(old)==85 and p['form']==d['form'] and p['strict_center']==d['strict_center']
 S=d['E7_support'];G=[[rational(dot(old[i],old[j]))for j in S]for i in S];pd(G);I=inv(G)
 nodes=[i for i in range(7)if I[i][i]==2];assert nodes==[3]and S[3]==15
 omega=lin([old[i]for i in S],list(map(c,I[3])))
 assert enc(omega)==p['weight'] and dot(omega,omega)==c(2)
 assert all(dot(omega,old[j])==c(i==3)for i,j in enumerate(S))
 rays=[vector(r['ambient_point'])for r in d['face_ideal_vertices']];assert len(rays)==2
 t=F(p['selected_shift']);assert t==2
 new=[lin([omega,nu],[c(-1),c(t)])for nu in rays];N=old+new
 assert enc(N)==p['normals']
 norms=[rational(dot(n,n))for n in N];assert enc(norms)==p['norms']and norms[-2:]==[2,2]
 counts=Counter();pairrows=[]
 for i in range(85,87):
  for j in range(i):
   g=dot(N[i],N[j]);assert sign(g)<=0
   if g==Z:kind='right'
   elif mul(g,g)==c(norms[i]*norms[j]/4):kind='pi3'
   else:assert sign(sub(mul(g,g),c(norms[i]*norms[j])))>0;kind='dotted'
   counts[kind]+=1;pairrows.append(dict(walls=[i,j],pairing=g,relation=kind))
 assert len(pairrows)==171 and enc(pairrows)==p['new_pairs']and dict(counts)==p['new_pair_counts']
 assert sign(dot(H,H))<0 and all(sign(dot(H,n))<0 for n in N)
 fp=[]
 for i,n in enumerate(N):
  point=lin([H,n],[c(1),scale(dot(H,n),-1/norms[i])])
  assert sign(dot(point,point))<0 and sign(dot(point,H))<0
  assert all(dot(point,r)==Z if i==j else sign(dot(point,r))<0 for j,r in enumerate(N))
  fp.append(point)
 assert enc(fp)==p['facet_points']
 sq2=(F(0),F(1))+(F(0),)*6
 def rationalize(ids):
  B=[[dot(N[i],N[j])for j in ids]for i in ids];bits={}
  # Infer the scaling parity from physical Gram entries, without using saved
  # component types, source facet labels, or the producer's scale list.
  for start in range(len(ids)):
   if start in bits:continue
   bits[start]=0;todo=[start]
   while todo:
    i=todo.pop()
    for j,g in enumerate(B[i]):
     if i==j or g==Z:continue
     assert not any(g[2:])and bool(g[0])!=bool(g[1]),g
     value=bits[i]^int(bool(g[1]))
     if j in bits:assert bits[j]==value
     else:bits[j]=value;todo.append(j)
  scales=[sq2 if bits[i]else c(1)for i in range(len(ids))]
  A=[[rational(mul(mul(scales[i],scales[j]),B[i][j]))for j in range(len(ids))]for i in range(len(ids))]
  return A,scales
 links=[]
 for k,nu in enumerate(rays):
  assert dot(nu,nu)==Z and sign(dot(nu,H))<0 and dot(nu,omega)==Z
  signs=[sign(dot(nu,n))for n in N];assert max(signs)==0
  ids=[i for i,s in enumerate(signs)if s==0];assert ids==p['completed_cusp_links'][k]['incident']
  unseen=set(ids);components=[]
  while unseen:
   todo=[min(unseen)];unseen.remove(todo[0]);block=[]
   while todo:
    i=todo.pop();block.append(i)
    js=[j for j in unseen if dot(N[i],N[j])!=Z]
    for j in js:unseen.remove(j);todo.append(j)
   components.append(sorted(block))
  records=[]
  for ids0 in components:
   A,sc=rationalize(ids0);n=len(A);minor=[r[:-1]for r in A[:-1]];pd(minor);J=inv(minor)
   kernel=[-sum(J[i][j]*A[j][-1]for j in range(n-1))for i in range(n-1)]+[F(1)]
   assert min(kernel)>0 and all(sum(a*b for a,b in zip(row,kernel))==0 for row in A)
   coefficients=[scale(s,q)for s,q in zip(sc,kernel)];v=lin([N[i]for i in ids0],coefficients)
   pos=next(j for j,x in enumerate(nu)if x!=Z);mu=fdiv(v[pos],nu[pos]);assert sign(mu)>0
   assert v==[mul(mu,x)for x in nu]
   records.append(dict(walls=ids0,scaled_gram=A,scales=sc,positive_kernel=coefficients,null_multiple=mu,quotient_rank=n-1))
  assert sorted(i for r in records for i in r['walls'])==ids
  assert sum(r['quotient_rank']for r in records)==19
  expected=sorted(sorted(r['walls'])for r in p['completed_cusp_links'][k]['positive_affine_blocks'])
  assert sorted(r['walls']for r in records)==expected
  for a in range(len(records)):
   for b in range(a):assert all(dot(N[i],N[j])==Z for i in records[a]['walls']for j in records[b]['walls'])
  links.append(dict(incident=ids,all87_signs=signs,independently_derived_affine_components=records,quotient_rank=19))
 assert sign(dot(*rays))<0
 vertices=[]
 for row in d['face_finite_vertices']:
  point=vector(row['ambient_point']);assert sign(dot(point,point))<0 and sign(dot(point,H))<0
  signs=[sign(dot(point,n))for n in N];assert max(signs)==0
  ids=[i for i,s in enumerate(signs)if s==0];assert ids==row['ambient_incident']and len(ids)==20
  A,_=rationalize(ids);pd(A)
  assert signs[-2:]==[-1,-1]
  vertices.append(dict(incident=ids,new_pairings=[dot(point,n)for n in new]))
 assert len(vertices)==44 and enc(vertices)==p['preserved_face_vertices']
 ep=[dot(E,n)for n in N];assert all(sign(x)<0 for x in ep)
 assert sign(dot(E,E))>0 and sign(dot(E,H))<0
 assert enc(dict(point=E,norm=dot(E,E),all87_pairings=ep))==p['strict_escape']
 # E7 projections: r + omega = 2nu; the projected future null inequalities
 # are strict on every hyperbolic point of the unchanged H13 face.
 assert all(lin([r,omega],[c(1),c(1)])==[scale(x,t)for x in nu]for r,nu in zip(new,rays))
 paths=[BASE,P,f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
 return enc(dict(status='PASS_INDEPENDENT_SPARSE_FIELD_TWO_E7_CUSP_FILLERS',sources={x:sha(x)for x in paths},
  field='Q(sqrt2,sqrt79,sqrt11)',sign_method='outward192bit rational radical bounds',
  reconstructed_caps=new,new_pair_counts=dict(counts),all_norms_checked=87,actual_facets_checked=87,
  facet_inequalities_checked=87*87,preserved_finite_vertices=44,finite_vertex_inequalities_checked=44*87,
  finite_vertex_PD20_checks=44,completed_links=links,escape_norm=dot(E,E),escape_inequalities_checked=87,
  source_dependencies='The old85 all-pair certificate and Lorentz signature are inherited via the hash-bound predecessor; no repeat of its full pair census. All new arithmetic, physical incidence and positive-link ranks reconstructed without project imports.',
  scope='Two specified full H20 cusp links are bounded. The complete H13 face and44 finite facevertices survive. All87 facets are actual. The strict spacelike escape proves infinite H20 volume; no full ambient vertex/cusp census or general completion theorem.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}_review.json'
 if args.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
