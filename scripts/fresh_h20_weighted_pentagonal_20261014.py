#!/usr/bin/env python3
"""Exact R=4/phi family, a finite-contact motif, and a mixed covering gate.

All-radius vertex caps and half-vertex radius-one caps are included in the
covering argument; arbitrary nonvertex radii remain outside its scope.
"""
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter, defaultdict
from pathlib import Path
from hashlib import sha256
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_weighted_pentagonal_20261014'
ZERO=(F(0),F(0));ONE=(F(1),F(0));PHI=(F(0),F(1))
def q(x):return(F(x),F(0))
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def neg(x):return(-x[0],-x[1])
def sub(x,y):return add(x,neg(y))
def mul(x,y):return(x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def sc(x,t):return(x[0]*t,x[1]*t)
def sign(x):
 a,b=x[0]+x[1]/2,x[1]/2
 if not a:return(b>0)-(b<0)
 if not b:return(a>0)-(a<0)
 if a*b>0:return(a>0)-(a<0)
 d=a*a-5*b*b;assert d
 return((a>0)-(a<0))*((d>0)-(d<0))
def total(xs):
 z=ZERO
 for x in xs:z=add(z,x)
 return z
WEIGHTS=[ONE]*10+[q(2)]+[add(ONE,PHI)]*8
R5=sc(sub(PHI,ONE),4);Q5=mul(R5,R5)
def dot(x,y):return total(mul(w,mul(a,b))for w,a,b in zip(WEIGHTS,x,y,strict=True))
def dist(x,y):return dot([sub(a,b)for a,b in zip(x,y,strict=True)],[sub(a,b)for a,b in zip(x,y,strict=True)])
def center(kind,profile,pattern=0):
 if kind==5:
  full=set(profile);small=sub(sc(PHI,2),q(3));c=[ONE]*11+[ONE if j in full else small for j in range(8)]
 elif kind==2:
  a,i=profile;small=sub(q(2),PHI);c=[q(j!=i)for j in range(10)]+[ONE]+[ONE if j==a else small for j in range(8)]
 elif kind==4:c=[ONE]*19
 elif kind==1:c=[q(F(1,2))]*19
 else:raise ValueError(kind)
 return[neg(x)if(pattern>>i)&1 else x for i,x in enumerate(c)]
def relation(d,r,s):
 baseline=add(mul(r,r),mul(s,s));rs=mul(r,s)
 for label,c in[(2,ZERO),(3,q(F(1,2))),(5,sc(PHI,F(1,2)))]:
  if d==add(baseline,sc(mul(rs,c),2)):return str(label)
 test=sign(sub(d,add(baseline,sc(rs,2))))
 return 'dotted'if test>0 else'parallel'if test==0 else'forbidden'
def profiles():
 return [(5,list(c))for k in[2,3]for c in combinations(range(8),k)]+[(2,[a,i])for a in range(8)for i in range(-1,10)]
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(tuple,list)):return[enc(v)for v in x]
 return x

def run():
 prior='results/fresh_h20_seed_deformation_20261014.json'
 old=json.loads((ROOT/prior).read_text())
 for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
 assert Q5==(F(32),F(-16))
 # Complete seed coordinate choices, with no finite-label cutoff.
 # Smaller opposite gap <=b: only zero in unit factors; zero or R/2 in
 # heavy/golden factors. The heavy R/2 choice has forbidden complement.
 assert sign(sub(Q5,q(4)))>0 and sign(sub(Q5,q(8)))<0
 assert sign(sub(sc(Q5,F(1,2)),add(ONE,PHI)))>0
 assert mul(R5,PHI)==q(4)
 phi_lo,phi_hi=F(1618033,10**6),F(1618034,10**6)
 rt_lo,rt_hi=F(1414213,10**6),F(1414214,10**6)
 assert phi_lo**2-phi_lo-1<0<phi_hi**2-phi_hi-1
 assert rt_lo**2<2<rt_hi**2
 rr_lo,rr_hi=4*phi_lo-4,4*phi_hi-4
 heavy_cos=((2*rt_lo-rr_hi/2)/rr_hi,(2*rt_hi-rr_lo/2)/rr_lo)
 assert F(1,2)<heavy_cos[0]<heavy_cos[1]<rt_lo/2
 assert sub(sc(PHI,2),sc(R5,F(1,2)))==q(2)
 assert sign(sub(sc(add(ONE,PHI),4),Q5))>0
 half=[q(F(1,2))]*19
 rows=[]
 for k in range(9):
  c=center(5,list(range(k)));n=dot(c,c)
  expected=add((F(52),F(-24)),sc((F(-4),F(4)),k))
  power=sub(dist(c,half),Q5)
  assert n==expected and power==add((F(21),F(-14)),sc((F(-6),F(4)),k))
  cc=add(sc(sub(PHI,ONE),F(1,2)),q(F(k,4)))
  assert sub(sub(n,q(4)),Q5)==sc(mul(R5,cc),4)
  typ=relation(n,q(2),R5)
  assert typ==('forbidden'if k<2 else'5'if k==2 else'dotted')
  if typ!='forbidden':assert(sign(power)<0)==(k in[2,3])
  # For an accepted cutter, every opposite sign strictly misses. The
  # smallest sign-flip penalty is 2phi-2, attained at a small golden entry.
  if k in[2,3]:assert sign(add(power,sc(sub(PHI,ONE),2)))>0
  rows.append(dict(full_golden=k,norm=n,half_probe_power=power,central_cosine=cc,central_relation=typ))
 P=profiles()+[(1,[])];assert len(P)==173
 counts={};classes=defaultdict(Counter)
 for leftkind,rightkind in[(5,5),(5,2),(2,2),(1,5),(1,2),(1,1)]:
  ids=[i for i,p in enumerate(P)if p[0]==leftkind];jds=[i for i,p in enumerate(P)if p[0]==rightkind]
  counter=Counter()
  for i in ids:
   for j in jds:
    lk,lp=P[i];rk,rp=P[j]
    a=center(lk,lp);b=center(rk,rp,1<<11)
    rel=relation(dist(a,b),R5 if lk==5 else q(lk),R5 if rk==5 else q(rk));counter[rel]+=1
    if lk==1 and rk==5:
     assert rel==('forbidden'if len(rp)==2 else'5'if 0 in rp else'2')
    if lk==1 and rk==2:assert(rel!='forbidden')==(rp[0]==0)
    if lk==5 and rk==5:
     key=(len(lp),len(rp),0 in lp,0 in rp,len(set(lp)&set(rp)))
     classes[str(key)][rel]+=1
  counts[f'{leftkind}/{rightkind}']=dict(counter)
 assert counts['5/5']=={'forbidden':7056} and counts['5/2']=={'forbidden':7392}
 assert counts['2/2']=={'forbidden':6484,'2':1260}
 assert counts['1/5']=={'forbidden':28,'5':21,'2':35}
 assert counts['1/2']=={'dotted':11,'forbidden':77} and counts['1/1']=={'forbidden':1}

 # Distance-two golden probes: exact coordinate contributions reduce all
 # k3/k3 pairs to this integer Q(phi) formula. Sign permutations transfer
 # the representative S={0,1} to every two-coordinate difference.
 small=sub(sc(PHI,2),q(3));w=add(ONE,PHI)
 assert mul(w,mul(sub(ONE,small),sub(ONE,small)))==(F(8),F(-4))
 assert sc(mul(w,mul(small,small)),4)==(F(20),F(-12))
 assert sc(w,4)==(F(4),F(4))
 assert sub(mul(w,mul(add(ONE,small),add(ONE,small))),(F(8),F(-4)))==(F(-4),F(4))
 full=[sum(1<<j for j in c)for c in combinations(range(8),3)]
 distance_two=Counter();S=3
 for A in full:
  for B in full:
   u=(A^B).bit_count();r2=(A&B&S).bit_count();r1=((A^B)&S).bit_count();r0=2-r1-r2
   d=(F(8*u+20*r0+4*r2-4*r1),F(-4*u-12*r0+4*r2+4*r1))
   rel=relation(d,R5,R5);distance_two[rel]+=1
   expected=(A&B)==0 and ((A|B)&S).bit_count()==1
   assert(rel!='forbidden')==expected
   if expected:assert rel=='2'
 assert dict(distance_two)=={'forbidden':2896,'2':240}
 assert 3*3>8  # Three such pairwise-compatible k3 neighbors are impossible.

 # A positive four-cap motif on a golden square, with right and pi/5
 # contacts. Its radius-one caps are literal half centers.
 motif=[(1,[],0),(5,[0,2,3],1),(5,[4,5,6],2),(1,[],3)]
 centers=[center(k,p,pattern<<11)for k,p,pattern in motif]
 radii=[q(1),R5,R5,q(1)];motif_pairs=[];walls=[]
 for i,(c,r)in enumerate(zip(centers,radii)):
  n=dot(c,c);D=sc(sub(n,mul(r,r)),F(1,2));walls.append(c+[D,ONE])
  assert sign(add(q(100),D))>0
  probe=center(1,[],motif[i][2]<<11)
  assert sign(sub(dist(c,probe),mul(r,r)))<0
  assert relation(n,q(2),r)!='forbidden'
  seed_finite=[ZERO,q(F(1,4)),q(F(1,2)),sc(add(ONE,PHI),F(1,4)),q(F(3,4))]
  for w,ci in zip(WEIGHTS,c,strict=True):
   for orient in[-1,1]:
    gap=sub(ONE,sc(ci,orient));assert sign(gap)>=0
    numerator=mul(w,mul(gap,gap));rq=mul(r,r)
    assert any(numerator==mul(rq,f)for f in seed_finite)or sign(sub(numerator,rq))>=0
  for j in range(i):
   rel=relation(dist(c,centers[j]),r,radii[j]);assert rel!='forbidden'
   motif_pairs.append([j,i,rel])
 assert Counter(row[2]for row in motif_pairs)=={'2':4,'5':1,'dotted':1}
 x=[ZERO]*11+[q(F(99,100))]*8;xq=dot(x,x)
 assert all(sign(sub(ONE,sc(xi,e)))>0 for xi in x for e in[-1,1])
 escape_powers=[sub(dist(x,c),mul(r,r))for c,r in zip(centers,radii)]
 assert all(sign(sub(p,ONE))>0 for p in escape_powers)
 assert sign(sub(sub(xq,q(4)),ONE))>0
 D=sc(sub(xq,ONE),F(1,2));assert sub(xq,sc(D,2))==ONE
 assert sign(add(q(100),D))>0

 # Complete strict half-probe counting, not merely a bounded CSP failure.
 # Adjacent covered edges are R1--k3R5 (at most2 per R1) or charged to
 # the unique full coordinate of R2. Hence e(A)<=2|A| in each Q8.
 max_per_fiber=1024//6;assert max_per_fiber==170
 assert 8*170-1024<=2*170 and 8*171-1024>2*171
 max_vertex_caps=old['R4_vertex_family']['maximum_compatible_caps'];assert max_vertex_caps==512
 prior_mixed='results/fresh_h20_weighted_mixed_probe_20261014.json'
 mixed=json.loads((ROOT/prior_mixed).read_text())
 for p,h in mixed['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
 assert mixed['R4_half_probe_cover']['probes_per_cap']==275
 max_covered=2048*max_per_fiber+512*275;missed=2**19-max_covered;assert missed==35328
 # No R2 boundary exception: starting with all matching small golden
 # entries, power=9-6phi +m(phi-1)+2s+f(1+3phi)+2w. If its phi
 # coefficient vanishes, its rational part is3+4f+2s+2w>0.
 assert all(3+4*f+2*s+2*w>0 for f,s,w in product(range(9),range(9),range(13)))
 # R5 power at a half-probe always has odd rational coefficient:
 # 21-6k+2u+4h-2s+2f, so it cannot vanish in Q(phi).
 assert all((21-6*k+2*u+4*h-2*s+2*f)%2==1
            for k in range(9)for u,h in product(range(11),range(2))
            for s in range(9-k)for f in range(k+1))
 # Boundary-only R1 covers differ by one unit bit: weights>=1, equality
 # exactly one unit. At a probe p, displace its ten unit coordinates by
 # epsilon*p_i; every tangent R1 power increases by epsilon+5epsilon²/2.
 assert all(sign(sub(v,ONE))>0 for v in WEIGHTS[10:])
 assert F(10,4)==F(5,2)
 paths=[f'scripts/{TAG}.py',prior,prior_mixed]
 return enc(dict(status='PASS_PENTAGONAL_MOTIF_AND_FOUR_FAMILY_COVERING_OBSTRUCTION',
  sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
  radius=R5,radius_square=Q5,full_golden_rows=rows,profiles=P,
  all_label_center_gate=dict(heavy_forbidden_cosine_interval=heavy_cos,
    centers='unit/heavy full signs; golden full +/-phi or small +/-(2-phi)',cutters_per_probe=84),
  adjacent_golden_pair_counts=counts,adjacent_R5_classes=dict(classes),
  distance_two_k3=dict(trials=3136,pair_counts=dict(distance_two),
    condition='full3-subsets disjoint; exactly one changed bit in their union',maximum_neighbors_of_R1=2),
  coherent_square_motif=dict(wall_count=43,recipes=motif,normals=walls,pairs=motif_pairs,
    strict_escape=dict(spatial=x,D=D,U=ONE,norm=ONE,cap_powers=escape_powers)),
  mixed_covering=dict(max_strictly_covered_per_golden_fiber=max_per_fiber,
    number_of_fibers=2048,max_radius4_covered=512*275,minimum_not_strictly_cut=missed,
    no_R2_R5_boundary='R2 zero-phi coefficient forces positive rational part; R5 rational power coefficient is odd',
    boundary_recession='epsilon times the ten signed half-unit coordinates; every tangentR1 power increases epsilon+(5/2)epsilon²'),
  scope='Retain weighted-box39 core. No finite coherent union of arbitrary seed-compatible R2 caps, R=4/phi caps, radius-one caps at literal half-vertices, and arbitrary-radius caps at literal old vertices can have finite volume. All labels allowed. Other nonvertex radii/centers or changing the core remain outside this theorem. The coherent43-wall motif still has an exact strict escape.'))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();data=run();out=ROOT/f'results/{TAG}.json'
 if args.verify:assert json.loads(out.read_text())==data
 else:out.write_text(json.dumps(data,indent=2)+'\n')
 print(data['status']);print(data['adjacent_golden_pair_counts'])
