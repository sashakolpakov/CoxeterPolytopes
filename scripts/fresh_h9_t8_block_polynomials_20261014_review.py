#!/usr/bin/env python3
"""Independent exact polynomial checks for all own-contact support patterns.

Pure Fraction and a separate Q(sqrt2) implementation. No SymPy, primary
arithmetic, or polynomial-certificate coefficients enter the construction.
"""
import argparse,json
from fractions import Fraction as F
from itertools import product,permutations
from math import comb
from pathlib import Path
from hashlib import sha256

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_block_polynomials_20261014_review'
Z=(F(0),F(0));O=(F(1),F(0));SQ2=(F(0),F(1));HALF=(F(1,2),F(0));COS4=(F(0),F(1,2))
NAMES=('e','a','b','d','c','delta');KEY=(0,)*6
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[0],-a[1]
def mul(a,b):return a[0]*b[0]+2*a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def scale(a,q):return a[0]*q,a[1]*q
def sign(x):
 a,b=x
 if not b:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if a*b>0:return (a>0)-(a<0)
 t=a*a-2*b*b;return ((a>0)-(a<0))*((t>0)-(t<0))
def const(x):
 if not isinstance(x,tuple):x=(F(x),F(0))
 return {KEY:x}if x!=Z else {}
def var(name):
 k=list(KEY);k[NAMES.index(name)]=1;return {tuple(k):O}
def plus(*ps):
 d={}
 for p in ps:
  for k,v in p.items():d[k]=add(d.get(k,Z),v)
 return {k:v for k,v in d.items()if v!=Z}
def sc(p,q):return {k:scale(v,q)for k,v in p.items()if scale(v,q)!=Z}
def times(*ps):
 d=const(1)
 for p in ps:
  z={}
  for i,a in d.items():
   for j,b in p.items():
    k=tuple(x+y for x,y in zip(i,j));z[k]=add(z.get(k,Z),mul(a,b))
  d={k:v for k,v in z.items()if v!=Z}
 return d
def sq(p):return times(p,p)
def minus(a,b):return plus(a,sc(b,-1))
def determinant(M):
 d={};n=len(M)
 for p in permutations(range(n)):
  v=times(*(M[i][p[i]]for i in range(n)))
  s=(-1)**sum(p[i]>p[j]for i in range(n)for j in range(i+1,n))
  d=plus(d,sc(v,s))
 return d
def power(x,n):
 p=O
 for _ in range(n):p=mul(p,x)
 return p
def unit_poly(p,lower):
 out={}
 for ex,co in p.items():
  terms={KEY:co}
  for i,n in enumerate(ex):
   if not n:continue
   lo=lower[NAMES[i]];width=add(O,neg(lo));factor={}
   for k in range(n+1):
    ix=list(KEY);ix[i]=k
    factor[tuple(ix)]=scale(mul(power(lo,n-k),power(width,k)),F(comb(n,k)))
   terms=times(terms,factor)
  out=plus(out,terms)
 return out
def bernstein(p,lower):
 u=unit_poly(p,lower);degree=tuple(max(k[i]for k in u)for i in range(6));values=[]
 for ix in product(*(range(d+1)for d in degree)):
  x=Z
  for k,a in u.items():
   if all(j<=i for j,i in zip(k,ix)):
    factor=F(1)
    for j,i,n in zip(k,ix,degree):factor*=F(comb(i,j),comb(n,j))
    x=add(x,scale(a,factor))
  assert sign(x)>0,(degree,ix,x)
  values.append((ix,x))
 return degree,u,values
def cases():
 patterns=[None,4,5,7,8]
 specs=[('asymmetric',('a','b','d'),COS4,False),('double',('e','a','b','d'),COS4,False),
        ('same_high',('a','d'),(F(9,10),F(0)),False),('opposite_high',('a','b'),(F(9,10),F(0)),True),
        ('same_low',('a','d'),COS4,False),('opposite_low',('a','b'),COS4,True)]
 for block,active,clower,opposite in specs:
  low=block.endswith('low');pool=[None]if low else patterns
  for ia,ib,ma,mb in product(pool,pool,[2,3],[2,3]):
   lower={name:(COS4 if low and name in active else HALF)for name in NAMES}
   lower['c']=clower;lower['delta']=Z
   own={name:var(name)if name in active else const(0)for name in['e','a','b','d']}
   if ia is not None:
    for name in ['e','a']:
     if name in active:own[name]=const(HALF)
   if ib is not None:
    for name in ['b','d']:
     if name in active:own[name]=const(HALF)
   yield block,ia,ib,ma,mb,own,lower,var('delta')if opposite else const(0)
def run():
 rows=[];count=0;minimum=None
 for block,ia,ib,ma,mb,own,lower,delta in cases():
  e,a,b,d=[own[n]for n in['e','a','b','d']];c=var('c');B=minus(const(1),sq(c))
  # Rebuild all foreign projections with the normalized B2 inverse.
  J=[[Z for _ in range(4)]for _ in range(4)]
  for i in range(4):
   J[i][i]=scale(O,2)
   J[i][i^1]=scale(SQ2,1)
  def foreign(i):return [neg(HALF)if j==i else Z for j in[4,5,7,8]]
  fx,fy=foreign(ia),foreign(ib);fz=[neg(HALF),neg(HALF),Z,Z]
  def dot(x,y):
   s=Z
   for i in range(4):
    for j in range(4):s=add(s,mul(x[i],mul(J[i][j],y[j])))
   return s
  u,v,p=map(const,[dot(fx,fx),dot(fy,fy),dot(fx,fy)])
  ac=const(add(Z if ma==2 else neg(HALF),neg(dot(fx,fz))))
  bc=const(add(Z if mb==2 else neg(HALF),neg(dot(fy,fz))))
  Q=const((F(25,19),F(1,2)));h=plus(delta,p)
  P=plus(sq(e),sc(times(c,e,a),2),sq(a));C=plus(sq(b),sc(times(c,b,d),2),sq(d))
  R=plus(times(e,b),times(a,d),times(c,plus(times(e,d),times(a,b))))
  qa=minus(times(minus(const(1),u),B),P)
  qb=minus(times(minus(const(1),v),B),C)
  ab=sc(plus(times(h,B),R),-1)
  scaled=[[qa,ab,times(ac,B)],[ab,qb,times(bc,B)],[times(ac,B),times(bc,B),sc(times(Q,B),-1)]]
  A=plus(times(minus(times(minus(const(1),u),minus(const(1),v)),sq(h)),B),
         sc(times(minus(const(1),u),C),-1),sc(times(minus(const(1),v),P),-1),sc(times(h,R),-2),sq(minus(times(e,d),times(a,b))))
  N=plus(sc(times(Q,A),-1),sc(times(B,plus(times(minus(const(1),u),sq(bc)),times(minus(const(1),v),sq(ac)),sc(times(h,ac,bc),2))),-1),
         times(C,sq(ac)),times(P,sq(bc)),sc(times(R,ac,bc),-2))
  assert determinant(scaled)==times(sq(B),N)
  assert determinant([r[:2]for r in scaled[:2]])==times(B,A)
  cert=[]
  for name,pol in [('N',N),('minus_A',sc(A,-1))]:
   deg,terms,coeff=bernstein(pol,lower);count+=len(coeff)
   for _,x in coeff:
    if minimum is None or sign(add(x,neg(minimum)))<0:minimum=x
   cert.append(dict(polynomial=name,degree=deg,unit_power_terms=[dict(index=k,value=v)for k,v in sorted(terms.items())],positive_Bernstein_coefficients=[dict(index=k,value=v)for k,v in coeff]))
  rows.append(dict(block=block,foreign=[ia,ib],third_labels=[ma,mb],lower_bounds=lower,certificates=cert))
 assert len(rows)==408
 # Exact Chebyshev lower bound for m>=7.
 f=lambda x:8*x**3+4*x*x-4*x-1
 assert f(F(31,50))==F(-559,15625) and (1+F(31,50))/2==F(81,100)
 source=f'scripts/{TAG}.py'
 return dict(status='PASS_INDEPENDENT_408_RESIDUAL_BLOCK_POLYNOMIALS',sources={source:sha256((ROOT/source).read_bytes()).hexdigest()},
  variables=NAMES,coefficient_field='Q(sqrt2), [rational,sqrt2 coefficient]',case_count=len(rows),positive_coefficient_count=count,
  minimum_positive_coefficient=minimum,cases=rows,
  scope='Independent algebraic reconstruction and exact Bernstein positivity only. Application to an incidence requires separate literal actual-vertex and mandatory-label guards. Does not by itself exclude a whole incidence or supply a Coxeter example.')
def encode(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(tuple,list)):return [encode(y)for y in x]
 if isinstance(x,dict):return {str(k):encode(v)for k,v in x.items()}
 return x
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=encode(run());out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['case_count'],d['positive_coefficient_count'],d['minimum_positive_coefficient'],flush=True)
