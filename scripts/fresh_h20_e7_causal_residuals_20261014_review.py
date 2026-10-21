#!/usr/bin/env python3
"""Narrow independent recipe/new-pair replay and source-review metadata."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from collections import Counter
import fresh_h20_e8_interior_caps_20261014_exact_review as K

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_causal_residuals_20261014'
BASE='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
FACE='results/fresh_h20_e7_intrinsic_caps_20261014.json'
Z=K.Z+K.Z
def lift(x):return x+K.Z
def sc(x):return lift(K.scalar(x))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,q):return tuple(q*x for x in a)
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):return K.add(K.multiply(a[:8],b[:8]),K.scale(K.multiply(a[8:],b[8:]),43))+K.add(K.multiply(a[:8],b[8:]),K.multiply(a[8:],b[:8]))
def sign(a):
 sa,sb=K.sign(a[:8]),K.sign(a[8:])
 if not sb:return sa
 if not sa or sa==sb:return sb
 q=K.sub(K.multiply(a[:8],a[:8]),K.scale(K.multiply(a[8:],a[8:]),43));sg=K.sign(q);assert sg
 return sa*sg
def total(xs):
 y=Z
 for x in xs:y=add(y,x)
 return y
def vec(v):return tuple(tuple(map(F,x))+(F(0),)*(16-len(x))for x in v)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(y)for y in x]
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 return x

def run():
 docs=[json.loads((ROOT/p).read_text())for p in [BASE,FACE,f'results/{TAG}_exact.json',f'results/{TAG}_contact.json']]
 for d in docs:
  for p,h in d['sources'].items():assert sha(p)==h,p
 base,face,d0,d1=docs;Q=[[F(x)for x in row]for row in base['form']]
 sparse=[[(j,a)for j,a in enumerate(row)if a]for row in Q];cache={}
 def dot(v,w):
  if w not in cache:cache[w]=tuple(tuple(sum(a*w[j][k]for j,a in row)for k in range(16))for row in sparse)
  return total(mul(a,b)for a,b in zip(v,cache[w]))
 def rational(x):assert not any(x[1:]);return x[0]
 def lin(vs,cs):return tuple(total(mul(c,v[k])for c,v in zip(cs,vs))for k in range(21))
 N=list(map(vec,base['normals']));H=vec(base['strict_center']);E=vec(base['strict_escape']['point'])
 S=base['E7_support'];A=face['face_facet_walls'][:13]
 I=K.inverse([[rational(dot(N[i],N[j]))for j in S]for i in S])
 J=K.inverse([[rational(dot(N[i],N[j]))for j in A]for i in A])
 W=[lin([N[i]for i in S],list(map(sc,row)))for row in I]
 L=[lin([N[i]for i in A],list(map(sc,row)))for row in J]
 sq2=(F(0),F(1))+(F(0),)*14
 T=lin([vec(face['face_basis'][-1])],[scale(sq2,F(1,2))])
 w=lin([W[i]for i in [2,4,6]],[sc(1)]*3)
 assert dot(w,w)==sc(F(95,2))and dot(T,T)==sc(-2)
 x1=tuple(F(1,11)if i==13 else F(0)for i in range(16))
 t1=tuple(F(17,11)if i==4 else F(0)for i in range(16))
 assert mul(x1,x1)==sc(F(86,11))and mul(t1,t1)==sc(F(289,11))
 reports=[]
 for k,(doc,x,t,deleted)in enumerate([(d0,sc(2),sc(5),[61,75]),(d1,x1,t1,[61])]):
  y=lin([L[0],L[12],T],[scale(sq2,-1),scale(x,-1),t]);r=lin([w,y],[sc(-1),sc(1)])
  assert r==vec(doc['normal'])and y==vec(doc['face_projection'])
  assert dot(y,y)==sc(F(-91,2))and dot(r,r)==sc(2)and sign(dot(y,T))<0
  assert sign(dot(r,E))>0 and sign(dot(r,H))<0
  retained=[i for i in range(87)if i not in deleted];newN=[N[i]for i in retained]+[r]
  assert newN==list(map(vec,doc['normals']))
  counts=Counter();pairings=[]
  for record,i in zip(doc['new_pairs'],retained):
   g=dot(r,N[i]);prod=2*F(base['norms'][i]);assert sign(g)<=0
   if g==Z:kind='right'
   elif mul(g,g)==sc(prod/4):kind='pi3'
   elif mul(g,g)==sc(prod/2):kind='pi4'
   elif mul(g,g)==sc(3*prod/4):kind='pi6'
   elif mul(g,g)==sc(prod):kind='parallel'
   else:assert sign(sub(mul(g,g),sc(prod)))>0;kind='dotted'
   counts[kind]+=1;assert record['old_wall']==i and record['relation']==kind and vec([record['pairing']])[0]==g
   pairings.append(g)
  assert dict(counts)==doc['new_pair_counts']
  if k==0:
   for i in deleted:
    g=dot(r,N[i]);assert sign(g)<0 and sign(sub(sc(F(base['norms'][i])/2),mul(g,g)))>0
  else:
   assert dot(r,N[75])==sc(-1)and sign(dot(r,N[61]))>0
   p=lin([H,r,N[75]],[sc(1),scale(add(scale(dot(H,r),2),dot(H,N[75])),F(-1,3)),scale(add(dot(H,r),scale(dot(H,N[75]),2)),F(-1,3))])
   assert p==vec(doc['actual_pi3_ridge']['point'])and sign(dot(p,p))<0 and sign(dot(p,H))<0
   assert all(dot(p,n)==Z if i in doc['actual_pi3_ridge']['walls']else sign(dot(p,n))<0 for i,n in enumerate(newN))
  # Independently check the new facet point and all preservation signs.
  p=lin([H,r],[sc(1),scale(dot(H,r),F(-1,2))])
  assert p==vec(doc['facet_points'][-1])and sign(dot(p,p))<0 and dot(p,r)==Z
  assert all(sign(dot(p,N[i]))<0 for i in retained)
  for row in face['face_finite_vertices']:assert sign(dot(vec(row['ambient_point']),r))<0
  for row in base['completed_cusp_links']:
   nu=vec(row['point']);assert sign(dot(nu,r))<0 and all(sign(dot(nu,N[i]))<0 for i in deleted)
  EE=vec(doc['strict_escape']['point']);assert sign(dot(EE,EE))>0 and sign(dot(EE,H))<0
  assert all(sign(dot(EE,n))<0 for n in newN)
  assert dot(EE,EE)==vec([doc['strict_escape']['norm']])[0]
  reports.append(dict(facets=len(newN),new_pairs=len(pairings),pair_counts=dict(counts),old_escape_cut=dot(r,E),new_escape_norm=dot(EE,EE),preserved_finite_points=44,preserved_bounded_cusps=2))
 paths=[BASE,FACE,f'scripts/{TAG}_exact.py',f'results/{TAG}_exact.json',f'scripts/{TAG}_contact.py',f'results/{TAG}_contact.json',f'scripts/{TAG}_review.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
 return enc(dict(status='PASS_SOURCE_PROOF_REVIEW_AND_SEPARATE_RECIPE_NEW_PAIR_REPLAY',sources={p:sha(p)for p in paths},checks=reports,
                 scope='Source/proof review plus fresh producer replays and a separate nested-quadratic reconstruction of both recipes/all171 new pairs/newfacet points/refined pi3 ridge/preservation signs/escapes. No second full old-pair, old-facet or finite-neighbor subset implementation; those claims are source-reviewed and replayed in the primary checker.'))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
