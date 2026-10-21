#!/usr/bin/env python3
"""Separate nested-quadratic recipe/pair replay; no construction imports."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_causal_residuals_20261014_bridge'
BASE='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json';FACE='results/fresh_h20_e7_intrinsic_caps_20261014.json'
P=[2,79,11,3];Z=(F(0),)*16

def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,q):return tuple(x*q for x in a)
def sub(a,b):return add(a,scale(b,-1))
def sc(q):return(F(q),)+(F(0),)*15

def mul(a,b):
 if not any(a)or not any(b):return(F(0),)*len(a)
 if not any(a[1:]):return scale(b,a[0])
 if not any(b[1:]):return scale(a,b[0])
 h=len(a)//2
 if not h:return(a[0]*b[0],)
 aa,ab,ba,bb=a[:h],a[h:],b[:h],b[h:];p=P[h.bit_length()-1]
 return add(mul(aa,ba),scale(mul(ab,bb),p))+add(mul(aa,bb),mul(ab,ba))

def sign(a):
 if not any(a):return 0
 if len(a)==1:return 1 if a[0]>0 else-1
 h=len(a)//2;x,y=a[:h],a[h:];sx,sy=sign(x),sign(y)
 if not sy:return sx
 if not sx or sx==sy:return sy
 ss=sign(sub(mul(x,x),scale(mul(y,y),P[h.bit_length()-1])));assert ss
 return sx*ss

def inv(M):
 n=len(M);A=[list(row)+[F(i==j)for j in range(n)]for i,row in enumerate(M)]
 for k in range(n):
  j=next(i for i in range(k,n)if A[i][k]);A[k],A[j]=A[j],A[k];q=A[k][k];A[k]=[x/q for x in A[k]]
  for i in range(n):
   if i!=k:
    q=A[i][k];A[i]=[x-q*y for x,y in zip(A[i],A[k])]
 return[row[n:]for row in A]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(y)for y in x]
 if isinstance(x,dict):return{k:enc(y)for k,y in x.items()}
 return x

def run():
 paths=[BASE,FACE,f'results/{TAG}.json'];base,face,d=[json.loads((ROOT/p).read_text())for p in paths]
 for doc in [base,face,d]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 D=tuple(map(F,d['field']['discriminant']));assert sign(D)>0
 EZ=(Z,Z)
 def ee(x):return(x,Z)
 def ea(x,y):return(add(x[0],y[0]),add(x[1],y[1]))
 def es(x,q):return(scale(x[0],q),scale(x[1],q))
 def em(x,y):return(add(mul(x[0],y[0]),mul(D,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0])))
 def sg(x):
  a,b=x;sa,sb=sign(a),sign(b)
  if not sb:return sa
  if not sa or sa==sb:return sb
  return sa*sign(sub(mul(a,a),mul(D,mul(b,b))))
 def total(xs):
  out=EZ
  for x in xs:out=ea(out,x)
  return out
 def parse(x):return tuple(tuple(map(F,a))for a in x)
 def bvec(v):return tuple(ee(tuple(map(F,a))+(F(0),)*(16-len(a)))for a in v)
 def evec(v):return tuple(map(parse,v))
 Q=[[F(x)for x in row]for row in base['form']];sp=[[(j,q)for j,q in enumerate(row)if q]for row in Q];cache={}
 def dot(v,w):
  if w not in cache:cache[w]=tuple(total(es(w[j],q)for j,q in row)for row in sp)
  return total(em(a,b)for a,b in zip(v,cache[w]))
 def lin(vs,cs):return tuple(total(em(c,v[k])for c,v in zip(cs,vs))for k in range(21))
 def rational(x):assert x[1]==Z and not any(x[0][1:]);return x[0][0]
 N=list(map(bvec,base['normals']));H=bvec(base['strict_center']);E=bvec(base['strict_escape']['point']);S=base['E7_support'];A=face['face_facet_walls'][:13]
 I=inv([[rational(dot(N[i],N[j]))for j in S]for i in S]);J=inv([[rational(dot(N[i],N[j]))for j in A]for i in A])
 W=[lin([N[i]for i in S],[ee(sc(q))for q in row])for row in I];L=[lin([N[i]for i in A],[ee(sc(q))for q in row])for row in J]
 sqrt2=ee(tuple(F(i==1)for i in range(16)));sqrt3=ee(tuple(F(i==8)for i in range(16)))
 T=lin([bvec(face['face_basis'][-1])],[es(sqrt2,F(1,2))]);w=lin([W[4],W[6]],[sqrt3,sqrt3]);assert dot(w,w)==ee(sc(F(129,2)))and dot(T,T)==ee(sc(-2))
 aa,bb,tt=[parse(d['recipe'][k])for k in ['a','b','t']];expected_t=[F(0)]*16;expected_t[4]=F(1,22);expected_t[12]=F(1)
 assert tt==ee(tuple(expected_t))and sg(ea(aa,ee(sc(-2))))>0 and sg(ea(bb,ee(sc(-2))))>0
 qa,qb,qc=[ee(tuple(map(F,x)))for x in d['recipe']['quadratic_coefficients']]
 assert em(qb,qb)==ea(ee(D),es(em(qa,qc),4))and sg(qa)>0
 assert total([em(qa,em(bb,bb)),em(qb,bb),qc])==EZ
 assert ea(es(em(qa,bb),2),qb)==(Z,sc(-1))and d['recipe']['branch']==-1
 target=ea(es(em(tt,tt),2),ee(sc(F(-125,2))))
 assert total([es(em(aa,aa),F(11,12)),es(em(bb,bb),F(2,3)),es(target,-1)])==EZ
 y=lin([L[0],L[12],T],[es(aa,-1),es(bb,-1),tt]);r=lin([w,y],[ee(sc(-1)),ee(sc(1))]);assert r==evec(d['normals'][-1])and dot(r,r)==ee(sc(2))
 assert dot(y,y)==ee(sc(F(-125,2)))and sg(dot(y,T))<0
 assert dot(r,N[61])==dot(r,N[75])==ee(sc(-1))
 counts=Counter()
 for i,(n,rec)in enumerate(zip(N,d['new_pairs'])):
  g=dot(r,n);prod=2*F(base['norms'][i]);assert sg(g)<=0
  if sg(g)==0:kind='right'
  elif em(g,g)==ee(sc(prod/4)):kind='pi3'
  elif em(g,g)==ee(sc(prod/2)):kind='pi4'
  elif em(g,g)==ee(sc(3*prod/4)):kind='pi6'
  elif em(g,g)==ee(sc(prod)):kind='parallel'
  else:assert sg(ea(em(g,g),ee(sc(-prod))))>0;kind='dotted'
  assert rec['relation']==kind and parse(rec['pairing'])==g;counts[kind]+=1
 assert dict(counts)==d['new_pair_counts'];NN=N+[r];assert all(sg(dot(H,n))<0 for n in NN)
 for rec in d['actual_pi3_ridges']:
  j=rec['walls'][0];hr,hn=dot(H,r),dot(H,N[j]);p=lin([H,r,N[j]],[ee(sc(1)),es(ea(es(hr,2),hn),F(-1,3)),es(ea(hr,es(hn,2)),F(-1,3))])
  assert p==evec(rec['point'])and sg(dot(p,p))<0 and sg(dot(p,H))<0
  assert all(dot(p,n)==EZ if i in [j,87]else sg(dot(p,n))<0 for i,n in enumerate(NN))
 for row in face['face_finite_vertices']:assert sg(dot(bvec(row['ambient_point']),r))<0
 for row in base['completed_cusp_links']:assert sg(dot(bvec(row['point']),r))<0
 X=evec(d['strict_escape']['point']);assert dot(X,X)==parse(d['strict_escape']['norm'])and sg(dot(X,X))>0 and sg(dot(X,H))<0
 assert all(sg(dot(X,n))<0 for n in NN)and sg(dot(E,r))>0
 E42=bvec(json.loads((ROOT/'results/fresh_h20_e7_replacement_escape_20261014.json').read_text())['point']);assert sg(dot(E42,r))>0
 paths+=['results/fresh_h20_e7_replacement_escape_20261014.json',f'scripts/{TAG}.py',f'scripts/{TAG}_review.py']
 return enc(dict(status='PASS_SEPARATE_NESTED_QUADRATIC_RECIPE_NEW_PAIRS_AND_SOURCE_PROOF_REVIEW',sources={p:sha(p)for p in paths},new_pairs=dict(counts),checked_actual_ridges=2,preserved_finite_points=44,preserved_bounded_cusps=2,strict_escape_norm=dot(X,X),scope='Separate Fraction nested-quadratic reconstruction and signs, with no producer imports: recipe, all87 new pairs, two actual ridges, center,44 vertex/two cusp preservation and both escape cuts/newescape. Facet actuality and complete20-support finite-vertex exclusion source-reviewed and primary replayed; no second full old87 implementation.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_review.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
