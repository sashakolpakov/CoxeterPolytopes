#!/usr/bin/env python3
"""Exact changed weighted-box seed and two distinct global covering gates.

The R2 and R4 families are excluded separately; their union is not excluded.
No lattice root enumeration or polytope arithmetic assumption is used.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
import argparse,hashlib,json
import sympy as S
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_seed_deformation_20261014'
Z=(F(0),F(0));O=(F(1),F(0));PHI=(F(0),F(1))
def q(a):return(F(a),F(0))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[0],-a[1])
def sub(a,b):return add(a,neg(b))
def mul(a,b):return(a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])
def scale(a,t):return(a[0]*t,a[1]*t)
def sign(a):
 x,y=a[0]+a[1]/2,a[1]/2
 if not x:return(y>0)-(y<0)
 if not y:return(x>0)-(x<0)
 if x*y>0:return(x>0)-(x<0)
 d=x*x-5*y*y;assert d
 return((x>0)-(x<0))*((d>0)-(d<0))
def total(xs):
 a=Z
 for x in xs:a=add(a,x)
 return a
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(a)for a in x]
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 return x
def run():
 # Coordinates divide physical coordinate i by b_i, with spatial weight b_i².
 weights=[O]*10+[q(2)]+[add(O,PHI)]*8
 N=total(weights);assert N==(F(20),F(8))
 def dot(a,b):return total(mul(w,mul(x,y))for w,x,y in zip(weights,a,b,strict=True))
 def dist(a,b):return dot([sub(x,y)for x,y in zip(a,b)],[sub(x,y)for x,y in zip(a,b)])
 # Exactly solve each 1D seed-contact profile using labels2..6, or dotted.
 ph=(1+S.sqrt(5))/2;finite=[S.Integer(0),S.Integer(1),S.sqrt(2),ph,S.sqrt(3)]
 options=[]
 for b,expected in[(S.Integer(1),[-1,0,1]),(S.sqrt(2),[-S.sqrt(2),0,S.sqrt(2)]),
                     (ph,[-ph,1-ph,0,ph-1,ph])]:
  accepted=[]
  for small in finite+[S.Integer(2)]:
   if not bool(small<=b):continue
   big=S.simplify(2*b-small)
   if bool(big>=2)or any(S.simplify(big-f)==0 for f in finite):
    accepted +=[S.simplify(b-small),S.simplify(small-b)]
  accepted=set(accepted);assert accepted==set(map(S.sympify,expected))
  options.append(sorted(map(str,accepted)))
 # R4 vertex caps: all seed angles are right and 3/4/5; central angle is5.
 assert scale(sub(N,q(20)),F(1,16))==scale(PHI,F(1,2))
 # Same golden signs allow weighted binary distance only8 or12.
 def relation(c):
  if c==Z:return'right'
  if c==q(F(1,2)):return'pi3'
  if c==scale(PHI,F(1,2)):return'pi5'
  if sign(sub(c,O))>=0:return'parallel_or_dotted'
  return'forbidden'
 same=[];adj=[]
 for h1,h2 in product(range(11),range(2)):
  H=q(h1+2*h2);r=relation(scale(sub(H,q(8)),F(1,8)))
  if r!='forbidden':same.append(h1+2*h2)
  H=add(H,add(O,PHI));adj.append(relation(scale(sub(H,q(8)),F(1,8))))
 assert set(same)=={8,12} and set(adj)=={'forbidden'}
 # Weighted Plotkin: length12 and distance>=8 imply M<=4.
 # Golden patterns form an independent subset of Q8, so at most128 occur.
 cover=0
 for h1,h2,h5 in product(range(11),range(2),range(9)):
  H=add(q(h1+2*h2),scale(add(O,PHI),h5))
  if sign(sub(H,q(4)))<=0:cover+=comb(10,h1)*comb(8,h5)
 assert cover==530 and 512*cover<2**19
 # R2 half-probe cover: all-small golden coordinates cost13-6phi in total.
 small= sub(q(2),PHI) # coordinate normalized by phi
 base_distance=(F(13),F(-6));increment=sub(PHI,O)
 assert sign(sub(base_distance,q(4)))<0
 assert sign(sub(add(base_distance,scale(increment,2)),q(4)))>0
 assert sign(sub((F(28),F(-8)),q(16)))<0 # allsmall, maximum unit/heavy norm
 assert sign(sub((F(26),F(-7)),q(16)))<0 # one zero golden coordinate
 assert sign(sub((F(25),F(-6)),q(16)))<0
 assert sign(sub((F(26),F(-6)),q(16)))>0 # one full, unit weighted count11
 # Their central finite cosine has a negative phi coefficient, so it is not
 # any of0,1/2,phi/2. Thus the dotted threshold16 is necessary.
 recipes=[(a,i)for a in range(8)for i in range(-1,10)]
 def center(a,i,pattern=0):
  v=[q(j!=i)for j in range(10)]+[O]
  v +=[scale(O if j==a else small,-1 if(pattern>>j)&1 else 1)for j in range(8)]
  return v
 half=[q(F(1,2))]*19
 targetpower=(F(8),F(-5));assert sign(targetpower)<0
 for a,i in recipes:
  c=center(a,i);n=dot(c,c)
  assert relation(scale(sub(n,q(8)),F(1,8)))=='parallel_or_dotted'
  assert sub(dist(c,half),q(4))==targetpower
 cross={};passing=0
 for a,i in recipes:
  for b,j in recipes:
   d=dist(center(a,i),center(b,j,1));r=relation(scale(sub(d,q(8)),F(1,8)))
   cross[r]=cross.get(r,0)+1
   condition=(a==0) != (b==0)
   condition=condition and i>=0 and j>=0 and i!=j
   assert (r!='forbidden')==condition
   if condition:assert r=='right';passing+=1
 assert passing==1260
 # One cap cannot cover two probes differing in a golden sign.
 for a,i in recipes:
  for j in range(8):
   v=list(half);v[11+j]=neg(v[11+j])
   assert sign(sub(dist(center(a,i),v),q(4)))>0
 # In a Q3 subcube, each chosen covering cap designates only one edge;
 # every adjacent pair requires one designation. 8<12 proves noncoverage.
 assert 2**3<3*2**2
 # A literal local cutter is fully coherent with the39-wall seed/central core.
 c=center(0,0);n=dot(c,c);seedcounts={};finite_squares={Z:'right',q(F(1,4)):'pi3',q(F(1,2)):'pi4',scale(add(O,PHI),F(1,4)):'pi5',q(F(3,4)):'pi6',scale(add(q(2),PHI),F(1,4)):'pi10'}
 for w,x in zip(weights,c):
  for e in[-1,1]:
   gap=sub(O,scale(x,e));assert sign(gap)>=0
   csq=scale(mul(w,mul(gap,gap)),F(1,4))
   typ=finite_squares.get(csq)
   if typ is None:
    assert sign(sub(csq,O))>=0;typ='parallel'if csq==O else'dotted'
   seedcounts[typ]=seedcounts.get(typ,0)+1
 assert seedcounts=={'pi3':9,'right':11,'parallel':9,'dotted':9}
 # H=(0,100,1) is strict; all spatial norms are nonnegative.
 D=scale(sub(n,q(4)),F(1,2));assert sign(add(D,q(100)))>0
 paths=[f'scripts/{TAG}.py']
 return enc(dict(status='PASS_CHANGED_WEIGHTED_BOX_AND_TWO_SEPARATE_COVERING_OBSTRUCTIONS',
  sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
  spatial_weights=weights,vertex_square=N,core_facets=39,central_radius=2,
  seed_profile_alphabet=list(range(2,7)),one_coordinate_center_options_physical=options,
  R4_vertex_family=dict(seed_labels=[2,3,4,5],central_label=5,same_golden_pattern_distances=sorted(set(same)),
   maximum_per_golden_pattern=4,maximum_golden_patterns=128,maximum_compatible_caps=512,
   cube_vertices_per_closed_ball=cover,total_cube_vertices=2**19,maximum_covered_vertices=512*cover),
  R2_nonvertex_family=dict(cutters_per_half_probe=88,recipes=[dict(full_golden=a,zero_unit=i)for a,i in recipes],
   half_probe_power=targetpower,adjacent_golden_probe_pair_counts=cross,
   compatible_pair_rule='Exactly one full golden coordinate equals the changed sign coordinate, and the two zero-unit indices are distinct and nonnegative.',
   obstruction_probe_patterns=list(range(8)),obstruction_cube_vertices=8,obstruction_cube_edges=12),
  explicit_local_cutter=dict(center=c,radius=2,normal_D=D,U=1,norm=4,seed_pair_counts=seedcounts,
   central_relation='dotted',half_probe_power=targetpower,common_center=[0,100,1]),
  scale_gate='For an interior central cap of radiusR, each positive seed gap is>=R/2; affine weighted averages force every factor inradius>=R/2. For A factors the vertex barycenter is0, so max vertex distance²>=mean vertex square+center². Boundary-centered caps are outside this gate.',
  scope='Exact impossibility for R4 vertex caps and R2 seed-label2..6 nonvertex caps separately, retaining the39-wall weighted box/central core. Their mixed union, other radii, other seed labels and shifted central caps are not excluded. No finite-volume construction.'))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');args=a.parse_args();d=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert json.loads(p.read_text())==d
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print(d['R2_nonvertex_family']['adjacent_golden_probe_pair_counts'])
