#!/usr/bin/env python3
"""Separate rational-interval scalar audit; no producer arithmetic imports."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_facet_cover_20261014_line_review';BITS=192;SCALE=1<<BITS
class I:
 def __init__(self,x=0,y=None):
  if isinstance(x,I):self.l,self.h=x.l,x.h
  else:self.l=F(x);self.h=self.l if y is None else F(y)
 def __add__(self,x):x=I(x);return I(self.l+x.l,self.h+x.h)
 __radd__=__add__
 def __neg__(self):return I(-self.h,-self.l)
 def __sub__(self,x):return self+-I(x)
 def __rsub__(self,x):return I(x)+-self
 def __mul__(self,x):
  x=I(x);v=[a*b for a in[self.l,self.h]for b in[x.l,x.h]];return I(min(v),max(v))
 __rmul__=__mul__
 def __truediv__(self,x):
  x=I(x);assert x.l*x.h>0;return self*I(1/x.h,1/x.l)
 def __rtruediv__(self,x):return I(x)/self
 def sqrt(self):
  assert self.l>=0
  def low(x):return F(isqrt(x.numerator*SCALE*SCALE//x.denominator),SCALE)
  return I(low(self.l),low(self.h)+F(1,SCALE))
 def enc(self):return [str(self.l),str(self.h)]
s5=I(5).sqrt();rad=((16+4*s5)/11).sqrt()
def dec(x):return I(x[0][0])+I(x[0][1])*s5+(I(x[1][0])+I(x[1][1])*s5)*rad
def dv(v):return [dec(x)for x in v]
def pair(x,y,Q):return sum((x[i]*Q[i][j]*y[j]for i in range(len(x))for j in range(len(y))),I())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 base='results/fresh_h9_mixed_d4_coupling_20261014.json';cap='results/fresh_h9_mixed_d4_global_caps_20261014.json';param='results/fresh_h9_h8_facet_cover_20261014_parameter_gate.json';line='results/fresh_h9_h8_facet_cover_20261014_line_exclusion.json'
 b,c,p,l=[json.loads((ROOT/x).read_text())for x in[base,cap,param,line]]
 for doc in[p,l]:
  for path,h in doc['sources'].items():assert sha(path)==h
 Q=[[dec(x)for x in row]for row in b['form']];N=[dv(v)for v in b['final_normals']];H=dv(b['center']);E=[I(x)for x in c['literal_escape']]
 hh,he,ee=[pair(x,y,Q)for x,y in[(H,H),(H,E),(E,E)]];a=hh+ee-2*he;bb=2*(he-hh);cc=hh
 assert a.h<0 and (2*a+bb).l>0
 def q(t):return a*t*t+bb*t+cc
 left,right=map(F,p['null_bracket']);assert q(left).h<0 and q(right).l>0
 gh=[pair(H,n,Q)for n in N];ge=[pair(E,n,Q)for n in N];nn=[pair(n,n,Q)for n in N]
 deriv=[]
 for i in range(12):
  assert gh[i].h<0 and ge[i].h<0 and nn[i].l>0
  d=ge[i]-gh[i];lin=d*bb-2*a*gh[i];con=2*d*cc-gh[i]*bb
  vals=[lin*left+con,lin+con];assert min(v.l for v in vals)>0
  deriv.append([x.enc()for x in vals])
 def kappa(t):
  z=[h*(1-t)+e*t for h,e in zip(H,E)];qq=pair(z,z,Q);assert qq.l>0
  gs=[pair(z,n,Q)for n in N];assert all(g.h<0 for g in gs)
  out=[g*g/qq-v for g,v in zip(gs,nn)];assert all(x.l>0 for x in out);return out
 low,up=map(F,p['necessary_open_t_interval']);T=F(l['split_parameter']);assert low<T<up
 ka=kappa(low);kt=kappa(T);small=(s5-1)/4;bound=4*small
 margins=[2*small-ka[7],kt[7]-small]+[bound-kt[i]for i in range(4)];assert all(x.l>0 for x in margins)
 # Independent symbolic normalization/derivative identities.
 t,g0,g1,q0,q1,q2=s.symbols('t g0 g1 q0 q1 q2');g=g0+g1*t;qq=q0+q1*t+q2*t*t
 assert s.expand(2*s.diff(g,t)*qq-g*s.diff(qq,t)-(g1*q1-2*g0*q2)*t-(2*g1*q0-g0*q1))==0
 x,y,z=s.symbols('x y z',positive=True)
 # derivative of (xy+z)/sqrt((x²-1)(y²-1)) has negative numerator.
 f=(x*y+z)/s.sqrt((x*x-1)*(y*y-1))
 expected=-(y+x*z)/((x*x-1)*s.sqrt((x*x-1)*(y*y-1)))
 assert s.simplify(s.diff(f,x)-expected)==0
 paths=[base,cap,param,line,'scripts/fresh_h9_h8_facet_cover_20261014_parameter_gate.py','scripts/fresh_h9_h8_facet_cover_20261014_line_exclusion.py','results/fresh_h9_h8_facet_cover_20261014_gap_gate.json',f'scripts/{TAG}.py']
 return dict(status='PASS_INDEPENDENT_PARAMETER_LINE_SCALAR_REVIEW',sources={p:sha(p)for p in paths},interval_bits=BITS,
  independent_null_norm_bounds=[q(left).enc(),q(right).enc()],independent_gamma_derivative_endpoints=deriv,
  independent_split_margins=[x.enc()for x in margins],symbolic_identity_count=2,
  scope='Independent physical-coordinate rational interval scalar checks and two symbolic identities. Source/theorem review covers the exact diameter argument and all-label line split; saved primary replays reuse producer arithmetic. Does not independently repeat820source distance comparisons or3657gap vertices.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
