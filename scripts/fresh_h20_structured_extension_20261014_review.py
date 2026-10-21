#!/usr/bin/env python3
"""Narrow independent frame/identity audit, not a duplicate tail census."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_structured_extension_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def sign(x):return(x>0)-(x<0)
def quadratic_sign(a,b):
 if not a:return sign(b)
 if not b or sign(a)==sign(b):return sign(a)
 return sign(a)*sign(a*a-2*b*b)
def run():
 paths=[f'scripts/{TAG}.py',f'results/{TAG}.json','results/fresh_h20_coupled_lift_20261014.json',
  'results/fresh_dimension_construction_routes_20261014.json',f'scripts/{TAG}_review.py']
 proof=load(paths[1]);d=load(paths[2]);old=load(paths[3]);Q=[[F(x)for x in r]for r in d['form']]
 def pair(x,y):return sum(x[i]*sum(Q[i][j]*y[j]for j in range(21))for i in range(21))
 H=list(map(F,d['strict_center']));X=list(map(F,d['strict_escape']['point']))
 assert pair(H,H)==-4296 and pair(X,X)==F(1411253767129583933,200000000000000000)and pair(X,H)<0
 assert all(F(r['d'])==2 for r in d['normals'][:60])
 for r in d['normals'][:60]:
  a=list(map(F,r['a']));b=list(map(F,r['b']))
  assert quadratic_sign(pair(a,H),pair(b,H))<0
  assert quadratic_sign(pair(a,X),pair(b,X))<0
 E=[F(x['a'])/5 for x in old['strict_escape']['boundary_point']]
 assert all(F(x['b'])==0 for x in old['strict_escape']['boundary_point'])
 assert pair(E,E)==1 and pair(E,H)<0
 frames=[]
 for r in d['base_recipes']:
  h=list(map(F,r['h0']));w=list(map(F,r['w']));z=[a-pair(h,E)*b for a,b in zip(h,E)]
  assert pair(z,z)<0<pair(w,w)and pair(z,H)<0 and pair(w,H)==0
  assert pair(E,z)==pair(E,w)==pair(z,w)==0
  rows=[]
  for j in r['support']:
   a=list(map(F,d['normals'][j]['a']));assert all(F(b)==0 for b in d['normals'][j]['b'])
   assert pair(a,E)==pair(a,z)==pair(a,w)==0;rows.append(a)
  assert len(rows)==18 and S.Matrix(rows).rank()==18
  frames.append(dict(support=r['support'],T_raw_square=str(pair(z,z)),W_raw_square=str(pair(w,w))))
 # Independently verify the identities used for the entire real tails.
 x,y=S.symbols('x y',positive=True);c=lambda t:(t+1/t)/2;s=lambda t:(t-1/t)/2
 assert S.expand(c(x)**2-s(x)**2)==1
 a,b,c0,w=S.symbols('a b c0 w',real=True)
 mean=a*c(x)+b*s(x)
 assert S.simplify(2*x*(mean+w)-((a+b)*x*x+2*w*x+a-b))==0
 assert S.simplify(S.diff(mean,x)-(a*(1-x**-2)+b*(1+x**-2))/2)==0
 A,B,C,D,E0,F0,G,H0,K,eps,eta=S.symbols('A B C D E F G H K eps eta',real=True)
 ip=c(x)*c(y)*A+s(x)*c(y)*B+c(x)*s(y)*C+s(x)*s(y)*D+eta*(c(x)*E0+s(x)*F0)+eps*(c(y)*G+s(y)*H0)+eps*eta*K
 dx=s(x)*(c(y)*A+s(y)*C+eta*E0)+c(x)*(c(y)*B+s(y)*D+eta*F0)
 dy=s(y)*(c(x)*A+s(x)*B+eps*G)+c(y)*(c(x)*C+s(x)*D+eps*H0)
 assert S.simplify(x*S.diff(ip,x)-dx)==0 and S.simplify(y*S.diff(ip,y)-dy)==0
 assert len(proof['range_certificates'])==4 and len(proof['transported_families'])==16 and len(proof['cross_family_dotted'])==120
 for r in proof['transported_families']:
  assert F(r['uniform_upper'][1])<0 and F(r['derivative_pairing'][1])<0
 for r in proof['cross_family_dotted']:
  assert all(F(a[1])<0 for a in r['derivative_gates'])
  assert all(F(a['pair'][1])<-2 for a in r['corner_pairs'])
 return dict(status='PASS_NARROW_SOURCE_PROOF_AND_SEPARATE_FRAME_IDENTITIES',sources={p:sha(p)for p in paths},
  frames=frames,old_center_and_escape_signs=120,independent_symbolic_identities=5,tail_families=16,cross_family_gates=120,
  scope='Source/theorem review plus separate Fraction rational-frame and Qsqrt2 sign checks, and independent SymPy Laurent identities. The complete rigorous interval tail gates are producer-replayed and source-reviewed, not recomputed by a second interval implementation. Scope requires both caps at each parameter, fixed W, the four supports/height flips and all60 predecessor walls.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['old_center_and_escape_signs'],d['tail_families'],flush=True)
