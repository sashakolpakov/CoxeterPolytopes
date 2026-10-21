#!/usr/bin/env python3
"""Independent frame and alternative-coordinate quadratic check of the S61 gate."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import sympy as sp
from fresh_h20_bridge_edge_closure_20261014_review import Field,parse
ROOT=Path(__file__).resolve().parents[1];BASE='fresh_h20_s61_single_cap_20261014';TAG=BASE+'_review';BITS=256;DEN=1<<BITS
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json','results/fresh_h20_e7_bridge_face_20261014_exact.json','results/fresh_h20_bridge_edge_closure_20261014.json',f'results/{BASE}.json']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def out(a):return(F(floor(a[0]*DEN),DEN),F(ceil(a[1]*DEN),DEN))
def pt(x):return(F(x),F(x))
def plus(a,b):return out((a[0]+b[0],a[1]+b[1]))
def neg(a):return(-a[1],-a[0])
def minus(a,b):return plus(a,neg(b))
def times(a,b):
 p=[x*y for x in a for y in b];return out((min(p),max(p)))
def div(a,b):
 assert b[0]>0 or b[1]<0;return times(a,(1/b[1],1/b[0]))
def root(a):
 assert a[0]>0
 lo=isqrt(floor(a[0]*DEN*DEN));hi=isqrt(floor(a[1]*DEN*DEN))+1
 z=F(lo,DEN),F(hi,DEN);assert z[0]**2<=a[0] and z[1]**2>a[1];return z
def overlap(a,b):return a[0]<=b[1] and b[0]<=a[1]
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(tuple,list)):return [enc(v) for v in x]
 if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
 return x
def run():
 old,pre,edge,d=[json.loads((ROOT/p).read_text()) for p in FILES]
 for obj in [old,pre,edge,d]:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 assert d['status']=='PASS_COMPLETE_ONE_CAP_S61_QUADRILATERAL_GATE'
 ff=Field([parse(old['field']['discriminant'],0),parse(pre['field']['new_positive_sqrt_radical'],1)])
 add,mul,sc,total,c,sign=ff.add,ff.mul,ff.scale,ff.total,ff.const,ff.sign;zero=ff.zeros[2]
 Q=[list(map(F,row)) for row in old['form']];rows=[[(j,q) for j,q in enumerate(row) if q] for row in Q];cache={}
 def vp(raw,k):return tuple(ff.lift(parse(x,k),k) for x in raw)
 def dot(v,w):
  if w not in cache:cache[w]=tuple(total(sc(w[j],q) for j,q in row) for row in rows)
  return total(mul(x,y) for x,y in zip(v,cache[w]))
 rec=edge['recipe'];assert int(rec['boost_integer'])==1 and rec['support']==d['support']==pre['recipe']['right_support']
 X=vp(rec['X'],2);W=vp(rec['W'],2);r=vp(pre['new_normal'],2);omega=vp(rec['dual_weight'],2);N0=vp(old['normals'][0],1)
 assert X==vp(d['frame']['X'],2) and W==vp(d['frame']['W'],2) and r==vp(d['frame']['r95'],2)
 Av=sc(dot(X,X),-1);wv=dot(W,W)
 assert Av==parse(rec['A'],2)==parse(d['frame']['A'],2) and wv==parse(rec['W_norm'],2)==parse(d['frame']['W_norm'],2)
 assert sign(Av)>0 and sign(add(wv,c(-2)))>0
 assert dot(r,r)==c(2) and dot(r,N0)==c(-1) and dot(omega,omega)==c(F(2,3))
 assert all(dot(v,u)==zero for v,u in [(X,W),(X,r),(W,r),(omega,X),(omega,W),(omega,r)])
 assert all(dot(v,vp(old['normals'][i],1))==zero for v in [X,W,r] for i in rec['support'])
 assert parse(edge['field']['new_positive_sqrt_radical'],2)==mul(Av,add(wv,c(F(-5,6))))
 part=tuple(add(sc(a,-1),sc(b,F(-1,2))) for a,b in zip(omega,r));assert part==vp(rec['particular'],2)
 for side,raw in zip([1,-1],edge['new_normals']):
  coords=[parse(x,3) for x in raw]
  assert tuple(x[0] for x in coords)==tuple(mul(Av,add(p,sc(w,side))) for p,w in zip(part,W))
  assert tuple(x[1] for x in coords)==X
 bounds=lambda x:out(ff.interval(x,2,BITS))
 A=bounds(Av);w=bounds(wv);v=minus(w,pt(F(5,6)));L=minus(v,pt(F(1,2)))
 assert L[0]>0
 rt2=root(pt(2));rt5=root(pt(5));rA=root(A);rw=root(w);rv=root(v)
 gx=div(bounds(dot(X,N0)),rA);ge=neg(div(pt(1),rt2));gf=div(bounds(dot(W,N0)),rw)
 # Reconstruct the frame; a cap's norm equation is solved in b, rather
 # than in a as in the producer. This avoids copying its root formula.
 b,s,dd,vv,ww=sp.symbols('b s d v w',real=True)
 norm=-(s-b/sp.sqrt(2))**2/vv+b*b+dd*dd/ww-2
 polynomial=(vv-sp.Rational(1,2))*b*b+sp.sqrt(2)*s*b+vv*dd*dd/ww-s*s-2*vv
 assert sp.expand(vv*norm-polynomial)==0
 discr=sp.expand(2*s*s-4*(vv-sp.Rational(1,2))*(vv*dd*dd/ww-s*s-2*vv))
 assert sp.simplify(discr-4*vv*(s*s+(vv-sp.Rational(1,2))*(2-dd*dd/ww)))==0
 cos={2:pt(0),3:pt(F(1,2)),4:div(rt2,pt(2)),5:div(plus(pt(1),rt5),pt(4))}
 assert d['complete_finite_contact_labels']==list(cos) and F(d['excluded_m_ge6_pivot_upper'])==F(-1,4)
 for m,p in cos.items():assert minus(pt(2),times(pt(3),times(p,p)))[0]>0
 assert F(2)-3*F(3,4)==F(-1,4)
 saved={(tuple(row['labels']),row['branch']):row for row in d['candidates']};assert len(saved)==32
 result=[];counts=Counter()
 for m,p in cos.items():
  for n,q in cos.items():
   ss=plus(p,q);difference=minus(q,p)
   disc=times(v,plus(times(ss,ss),times(L,minus(pt(2),div(times(difference,difference),w)))))
   radical=root(disc)
   for side in [-1,1]:
    bb=div(plus(neg(div(ss,rt2)),times(pt(side),radical)),L)
    aa=div(minus(ss,div(bb,rt2)),rv);cc=div(difference,rw)
    g95=times(rt2,bb);g0=plus(plus(times(aa,gx),times(bb,ge)),times(cc,gf))
    blocker=95 if side==1 else 0;g=g95 if blocker==95 else g0
    assert g[0]>F(99,100);counts[blocker]+=1
    row=saved[(m,n),-side];assert row['positive_blocker']==blocker
    assert all(overlap(z,tuple(map(F,raw))) for z,raw in zip([aa,bb,cc],row['frame_coeff_bounds']))
    assert overlap(g0,tuple(map(F,row['old0_pair_bounds']))) and overlap(g95,tuple(map(F,row['old95_pair_bounds'])))
    # Outward-rounded scalar arithmetic certifies the asserted positive
    # pairing. Exact quadratic algebra, not interval overlap, proves roots.
    result.append(dict(labels=[m,n],producer_branch=-side,positive_blocker=blocker,uniform_pairing_lower_bound='99/100'))
 assert dict(counts)=={0:16,95:16}
 paths=FILES+[f'scripts/{BASE}.py',f'scripts/{TAG}.py','scripts/fresh_h20_bridge_edge_closure_20261014_review.py','scripts/fresh_h20_bridge_faces_20261014_review.py','results/fresh_h20_bridge_edge_closure_20261014_review.json']
 return enc(dict(status='PASS_INDEPENDENT_S61_SINGLE_CAP_FRAME_AND_32_BRANCHES',sources={p:sha(p) for p in paths},implementation='Independent earlier sparse-radical arithmetic, new dyadic interval operations, and symbolic elimination in the other coordinate; no primary gate/field imports.',interval_bits=BITS,candidates=result,positive_blocker_counts=counts,complete_labels=[2,3,4,5],scope='One norm-two cap perpendicular toS61 and meeting both retained96/97 walls at finite endpoints. All labels and both real quadratic branches are covered; multiple caps, nonright support contacts and changed oldwalls remain outside.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
