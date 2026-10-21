#!/usr/bin/env python3
"""Separate Fraction/quartic reconstruction for the changed middle block.

Reuses only the independent sparse-polynomial operations; their coefficient
operations are replaced in this process by Q(sqrt2,sqrt5) arithmetic.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_t8_block_polynomials_20261014_review as P
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_middle_20261014_review_algebra'
q2sign=P.sign;q2mul=P.mul
Z=(F(0),)*4;O=(F(1),F(0),F(0),F(0))
def pad(a):return tuple(a)+(F(0),)*(4-len(a))
def add(a,b):return tuple(x+y for x,y in zip(pad(a),pad(b)))
def neg(a):return tuple(-x for x in pad(a))
def scale(a,q):return tuple(x*q for x in pad(a))
def mul(a,b):
 a,b=pad(a),pad(b);z=[F(0)]*4
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i^j]+=x*y*(2 if i&j&1 else 1)*(5 if i&j&2 else 1)
 return tuple(z)
def sign(a):
 a=pad(a);u,v=a[:2],a[2:];su,sv=q2sign(u),q2sign(v)
 if not sv:return su
 if not su or su==sv:return sv
 uu,vv=q2mul(u,u),q2mul(v,v)
 return su*q2sign((uu[0]-5*vv[0],uu[1]-5*vv[1]))
def const(q):
 if not isinstance(q,tuple):q=(F(q),)
 q=pad(q);return {P.KEY:q}if q!=Z else {}
P.Z=Z;P.O=O;P.HALF=(F(1,2),F(0),F(0),F(0));P.SQ2=(F(0),F(1),F(0),F(0));P.COS4=scale(P.SQ2,F(1,2))
P.add=add;P.neg=neg;P.mul=mul;P.scale=scale;P.sign=sign;P.const=const

def run():
 rows=[];count=0;minimum=None
 for block,ia,ib,ma,mb,own,lower,delta in P.cases():
  if ia not in[7,8]and ib not in[7,8]:continue
  e,a,b,d=[own[n]for n in['e','a','b','d']];c=P.var('c');B=P.minus(const(1),P.sq(c))
  J=[[Z for _ in range(4)]for _ in range(4)]
  for i in[0,1]:J[i][i]=scale(O,2);J[i][i^1]=P.SQ2
  for i in[2,3]:J[i][i]=(F(2),F(0),F(2,5),F(0));J[i][i^1]=(F(1),F(0),F(3,5),F(0))
  cos5=(F(1,4),F(0),F(1,4),F(0))
  # Independently invert the middle normalized Gram, not the producer costs.
  for i in range(2):
   for j in range(2):
    x=add(J[2+i][2+j],neg(mul(cos5,J[2+(i^1)][2+j])))
    assert x==(O if i==j else Z)
  def foreign(w):return [neg(P.HALF)if j==w else Z for j in[4,5,7,8]]
  fx,fy=foreign(ia),foreign(ib);fz=[neg(P.HALF),neg(P.HALF),Z,Z]
  def dot(x,y):
   s=Z
   for i in range(4):
    for j in range(4):s=add(s,mul(x[i],mul(J[i][j],y[j])))
   return s
  u,v,p=map(const,[dot(fx,fx),dot(fy,fy),dot(fx,fy)])
  ac=const(add(Z if ma==2 else neg(P.HALF),neg(dot(fx,fz))))
  bc=const(add(Z if mb==2 else neg(P.HALF),neg(dot(fy,fz))))
  Q=const((F(25,19),F(1,2),F(0),F(0)));h=P.plus(delta,p)
  pc=P.plus(P.sq(e),P.sc(P.times(c,e,a),2),P.sq(a));cc=P.plus(P.sq(b),P.sc(P.times(c,b,d),2),P.sq(d))
  cross=P.plus(P.times(e,b),P.times(a,d),P.times(c,P.plus(P.times(e,d),P.times(a,b))))
  qa=P.minus(P.times(P.minus(const(1),u),B),pc);qb=P.minus(P.times(P.minus(const(1),v),B),cc)
  ab=P.sc(P.plus(P.times(h,B),cross),-1)
  M=[[qa,ab,P.times(ac,B)],[ab,qb,P.times(bc,B)],[P.times(ac,B),P.times(bc,B),P.sc(P.times(Q,B),-1)]]
  aa=P.plus(P.times(P.minus(P.times(P.minus(const(1),u),P.minus(const(1),v)),P.sq(h)),B),
   P.sc(P.times(P.minus(const(1),u),cc),-1),P.sc(P.times(P.minus(const(1),v),pc),-1),P.sc(P.times(h,cross),-2),P.sq(P.minus(P.times(e,d),P.times(a,b))))
  nn=P.plus(P.sc(P.times(Q,aa),-1),P.sc(P.times(B,P.plus(P.times(P.minus(const(1),u),P.sq(bc)),P.times(P.minus(const(1),v),P.sq(ac)),P.sc(P.times(h,ac,bc),2))),-1),
   P.times(cc,P.sq(ac)),P.times(pc,P.sq(bc)),P.sc(P.times(cross,ac,bc),-2))
  assert P.determinant(M)==P.times(P.sq(B),nn)
  assert P.determinant([r[:2]for r in M[:2]])==P.times(B,aa)
  cert=[]
  for name,pol in [('N',nn),('minus_A',P.sc(aa,-1))]:
   degree,terms,coeff=P.bernstein(pol,lower);count+=len(coeff)
   for _,x in coeff:
    if minimum is None or sign(add(x,neg(minimum)))<0:minimum=x
   cert.append(dict(polynomial=name,degree=degree,unit_power_terms=[dict(index=k,value=v)for k,v in sorted(terms.items())],positive_Bernstein_coefficients=[dict(index=k,value=v)for k,v in coeff]))
  rows.append(dict(block=block,foreign=[ia,ib],third_labels=[ma,mb],lower_bounds=lower,certificates=cert))
 assert len(rows)==256
 paths=[f'scripts/{TAG}.py','scripts/fresh_h9_t8_block_polynomials_20261014_review.py']
 return P.encode(dict(status='PASS_INDEPENDENT_256_QUARTIC_MIDDLE_BLOCK_CASES',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
  field_basis=[1,2,5,10],case_count=len(rows),positive_coefficient_count=count,minimum_positive_coefficient=minimum,cases=rows,
  scope='Independent algebra for the256 changed label5 middle-block cases; literal incidence/domain inheritance is a separate obligation. No whole-incidence exclusion.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['case_count'],d['positive_coefficient_count'],d['minimum_positive_coefficient'],flush=True)
