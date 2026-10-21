#!/usr/bin/env python3
"""Independent compressed counts and Q(sqrt5) geometry; no producer imports."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import comb
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_outer_unit_code_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def K(a=0,b=0):return F(a),F(b)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def mul(x,y):return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def scale(x,a):return x[0]*a,x[1]*a
def sg(x):
 a,b=x
 if a==0:return (b>0)-(b<0)
 if b==0 or a*b>0:return (a>0)-(a<0)
 d=a*a-5*b*b
 return ((a>0)-(a<0))*((d>0)-(d<0))
def decode(x):return K(F(x[0])+F(x[1])/2,F(x[1])/2)
def out(x):return [str(t)for t in x]
PHI=K(F(1,2),F(1,2));W=[K(1)]*10+[K(2)]+[add(K(1),PHI)]*8
def dot(x,y):
 r=K()
 for a,b,w in zip(x,y,W,strict=True):r=add(r,scale(w,a*b))
 return r
def dist(x,y):return dot([a-b for a,b in zip(x,y)],[a-b for a,b in zip(x,y)])
def smallest(xs):
 it=iter(xs);r=next(it)
 for x in it:
  if sg(sub(x,r))<0:r=x
 return r

def unit_histogram(X,Y):
 # Each common nonzero coordinate contributes 2+2t^4, each single
 # nonzero coordinate contributes 2t. Zero/zero contributes 1.
 out=Counter()
 for x in X:
  for y in Y:
   single=(x^y).bit_count();both=10-(x|y).bit_count()
   for r in range(both+1):out[single+4*r]+=2**(both+single)*comb(both,r)
 return out

def family_pairs(X,Y,nx,ny,kind):
 h=unit_histogram(X,Y);counts=Counter();same=kind in ('CC','HH')
 for d,n in h.items():
  if d==0:
   assert same;label='identity'
  elif kind=='CC':
   assert d>=4 and d%2==0;label={4:'right',6:'pi3',8:'parallel'}.get(d,'dotted')
  elif kind=='HH':
   assert d>=2 and d%2==0;label={2:'right',4:'parallel'}.get(d,'dotted')
  else:
   assert d>=3 and d%2==1;label={3:'right',5:'pi4'}.get(d,'dotted')
  counts[label]+=512*n
 assert sum(counts.values())==nx*ny//512
 off=nx*ny-sum(counts.values())
 if kind=='CC':counts['parallel']+=nx;off-=nx
 counts['dotted']+=off
 if same:
  assert counts.pop('identity')==nx
  assert all(n%2==0 for n in counts.values())
  counts=Counter({k:v//2 for k,v in counts.items()})
 assert sum(counts.values())==(nx*(nx-1)//2 if same else nx*ny)
 return dict(sorted(counts.items()))

def small_pairs(X,nx,kind):
 c=Counter();ns=2**19
 for m in X:
  k=10-m.bit_count()
  for r in range(k+1):
   label=({0:'right',1:'parallel'}.get(r,'dotted')if kind=='SC' else ('pi4'if r==0 else 'dotted'))
   c[label]+=512*1024*comb(k,r)
 assert sum(c.values())==nx*ns//512
 c['dotted']+=nx*ns-sum(c.values())
 return dict(sorted(c.items()))

def run():
 src=f'results/{TAG}.json';saved=json.loads((ROOT/src).read_text())
 for p,h in saved['sources'].items():assert sha(p)==h,p
 cols=(0,1,2,3,4,5,8,10,12,15)
 def paritycheck(m):
  return tuple(sum((m>>j&1)*(cols[j]>>i&1)for j in range(10))%2 for i in range(4))
 C=[m for m in range(1024)if m.bit_count()%2==0 and paritycheck(m)==(0,0,0,0)]
 H=[m for m in range(1024)if m.bit_count()%2 and min((m^c).bit_count()for c in C)>=3]
 assert C==saved['code']and H==saved['holes']and len(C)==32 and len(H)==192
 assert all(a^b in C for a in C for b in C)
 assert min(m.bit_count()for m in C if m)==4
 assert all(min((h^c).bit_count()for c in C)==3 for h in H)
 nc=512*sum(2**(10-m.bit_count())for m in C);nh=512*sum(2**(10-m.bit_count())for m in H);ns=2**19
 counts=dict(code=nc,hole=nh,small=ns,binary=ns,core=39)
 assert counts==saved['counts']and sum(counts.values())==saved['actual_facets']==7545383
 pairs=dict(code_code=family_pairs(C,C,nc,nc,'CC'),hole_hole=family_pairs(H,H,nh,nh,'HH'),
  code_hole=family_pairs(C,H,nc,nh,'CH'),small_code=small_pairs(C,nc,'SC'),small_hole=small_pairs(H,nh,'SH'))
 ss={'right':ns*10//2,'parallel':ns*comb(10,2)//2}
 ss['dotted']=ns*(ns-1)//2-sum(ss.values());pairs['small_small']=dict(sorted(ss.items()))
 assert pairs==saved['new_pair_counts']
 assert sg(sub(scale(add(K(1),PHI),4),K(8)))>0
 # Seed, central, binary gates are formulas for all mask/sign choices.
 def vector(m):return [F(not(m>>i&1))for i in range(10)]+[F(1)]*9
 seed_rows=[]
 for family,masks,q in [('code',C,F(2)),('hole',H,F(1))]:
  for k in sorted({m.bit_count()for m in masks}):
   n=add(K(20-k),scale(PHI,8));num=sub(n,K(q+4))
   assert sg(num)>0 and sg(sub(mul(num,num),K(16*q)))>0
   labels={'right':19-k,'dotted':19-k,'pi4'if q==2 else 'parallel':2*k}
   assert sum(labels.values())==38
   seed_rows.append([family,k,labels,n])
 n=add(K(F(25,2)),scale(PHI,8));num=sub(n,K(F(9,2)))
 assert sg(num)>0 and sg(sub(mul(num,num),K(8)))>0
 seed_rows.append(['small',None,{'right':9,'pi4':10,'dotted':19},n])
 for got,expect in zip(seed_rows,saved['seed_counts'],strict=True):
  assert got[:3]==expect[:3]and got[3]==decode(expect[3])
 assert F(8)<9 # sqrt2<3/2; strict thresholds below 6 and 3.
 bounds={'code':add(K(F(17,4)),scale(PHI,2)),'hole':add(K(F(17,4)),scale(PHI,2)),
  'small':add(K(F(13,4)),scale(PHI,2))}
 for name,b in bounds.items():
  assert b==decode(saved['binary_minimum_distance_squares'][name])
  assert sg(sub(b,K({'code':6,'hole':4,'small':3}[name])))>0
 for full in ({0,1,2},{0,3,4}):
  b=[F(1)if i in full else F(1,2)for i in range(19)]
  for m in C+H:assert dist(vector(m),b)==add(bounds['code'],K(sum(m>>i&1 for i in full)))
  assert dist([F(1,2)]*10+[F(1)]*9,b)==bounds['small']
 assert saved['binary_new_dotted_pairs']==ns*(nc+nh+ns)
 # Original finite endpoint margins and a common strict center.
 assert 0 in C and 1023 in C
 assert min(k.bit_count()for k in C if k)>0
 minD=scale(sub(add(K(10),scale(PHI,8)),K(2)),F(1,2))
 assert minD==decode(saved['minimum_new_D'])and sg(minD)>0
 assert saved['preserved_old_finite_endpoints']==2**19
 assert saved['endpoint_inequality_positive_margins']=={'other_code':'4','hole':'1','small':'3/2'}
 oldx=[F(0)]*10+[F(99,100)]*9
 oldpower=sub(dist(oldx,vector(1023)),K(2));assert oldpower==decode(saved['old_escape_repair']['power'])and sg(oldpower)<0
 x=[F(t)for t in saved['strict_escape']['spatial']]
 assert x==[F(1,2),F(0),F(1,2),F(99,100),F(1,2),F(99,100),F(1,2),F(1,2),F(1,2),F(99,100)]+[F(99,100)]*9
 powers={'code':smallest(sub(dist(x,vector(m)),K(2))for m in C),
 'hole':smallest(sub(dist(x,vector(m)),K(1))for m in H),
 'small':sub(dist(x,[F(1,2)]*10+[F(1)]*9),K(F(1,2))),
 'binary_lower_bound':smallest(sub(dist(x,[F(1)if i in f else F(1,2)for i in range(19)]),K(1))for f in ({0,1,2},{0,3,4})),
 'central':sub(dot(x,x),K(4))}
 assert all(abs(t)<1 for t in x)
 assert min(1-abs(t)for t in x)==F(1,100)
 for name,v in powers.items():assert v==decode(saved['strict_escape']['powers'][name])and sg(sub(v,K(F(1,4))))>0
 D=scale(sub(dot(x,x),K(F(1,4))),F(1,2))
 assert D==decode(saved['strict_escape']['D'])and sub(dot(x,x),scale(D,2))==K(F(1,4))and sg(add(D,K(100)))>0
 # Literal actual pi5 ridge, checked without the binary producer.
 r2=K(1);w2=add(K(1),PHI);rw=scale(w2,F(-1,2));hr=neg(add(K(F(825,8)),PHI));hw=neg(w2)
 ar=add(K(167),scale(PHI,F(173,2)));aw=add(K(F(169,2)),scale(PHI,F(173,4)))
 assert add(hr,add(mul(ar,r2),mul(aw,rw)))==K()
 assert add(hw,add(mul(ar,rw),mul(aw,w2)))==K()
 det=sub(mul(r2,w2),mul(rw,rw));assert det==scale(add(PHI,K(2)),F(1,4))and sg(det)>0
 p2=add(K(-200),add(mul(ar,hr),mul(aw,hw)))
 assert p2==neg(add(K(F(141089,8)),scale(PHI,F(149517,16))))and sg(p2)<0 and sg(ar)>0 and sg(aw)>0
 # An even simpler retained obstruction is the original central facet.
 N=add(K(20),scale(PHI,8));corner_q=sub(N,K(4))
 assert corner_q==add(K(16),scale(PHI,8))and sg(corner_q)>0
 central_ridge_q=sub(add(K(1),PHI),K(4))
 assert central_ridge_q==sub(PHI,K(3))and sg(central_ridge_q)<0
 assert mul(PHI,PHI)==add(K(1),PHI)
 # Central norm4, golden seed normphi^2 and pair-phi^2 give cosphi/2.
 assert mul(neg(add(K(1),PHI)),neg(add(K(1),PHI)))==mul(add(K(1),PHI),mul(PHI,PHI))
 paths=[f'scripts/{TAG}_review.py',f'scripts/{TAG}.py',src,
 'results/fresh_h20_weighted_binary_20261014_local_review.json',
 'results/fresh_h20_weighted_binary_20261014_dotted_completion_review.md']
 return dict(status='PASS_INDEPENDENT_COMPRESSED_CODE_COUNTS_AND_ESCAPE',sources={p:sha(p)for p in paths},
  counts=counts,total_facets=sum(counts.values()),new_pair_counts=pairs,
  arithmetic='Independent stdlib Fraction basis 1,sqrt5; no producer or project arithmetic imports.',
  minimum_powers_Qsqrt5={k:out(v)for k,v in powers.items()},old_escape_cut_power_Qsqrt5=out(oldpower),
  strict_escape_norm='1/4',literal_pi5_ridge_projection_norm_Qsqrt5=out(p2),
  central_facet_obstruction=dict(all_coordinate_corner_norm_Qsqrt5=out(corner_q),actual_pi5_ridge_point_norm_Qsqrt5=out(central_ridge_q)),
  scope='Separate combinatorial generating-polynomial counts for all new-new pairs; exact new-seed/central/binary bounds, endpoint preservation margins, center and old/new escapes. Inherited binary old-old gates and local no-finite-vertex theorem are hash-bound, not rerun. Actual-facet and dotted-only-completion arguments are source/proof reviewed, not a new full vertex or cusp census.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['total_facets'])
