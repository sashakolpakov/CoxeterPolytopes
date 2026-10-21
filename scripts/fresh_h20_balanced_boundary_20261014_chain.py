#!/usr/bin/env python3
"""Exact Q(sqrt65) seven inner caps, two full outer orbits, and a bank escape.

The 415-wall baseline omits the previous three-cap motif. The selected 791 new
caps are mutually compatible; the larger 1176-cap escape bank need not be.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter
from pathlib import Path
import json,hashlib,argparse
from fresh_h20_balanced_20261014 import ROOT,exact_data
TAG='fresh_h20_balanced_boundary_20261014_chain'
BASE='results/fresh_h20_balanced_20261014_bridge_exact.json'
ZERO=(F(0),F(0));ONE=(F(1),F(0))
def q(x):return (F(x),F(0))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[0],-a[1])
def sub(a,b):return add(a,neg(b))
def mul(a,b):return(a[0]*b[0]+65*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def sc(a,b):return(a[0]*b,a[1]*b)
def sg(a):
 x,y=a
 if not y:return (x>0)-(x<0)
 if not x:return (y>0)-(y<0)
 if x*y>=0:return (x>0)-(x<0)
 d=x*x-65*y*y;assert d
 return ((d>0)-(d<0))*((x>0)-(x<0))
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 return x
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 base=json.loads((ROOT/BASE).read_text());assert all(sha(p)==h for p,h in base['sources'].items())
 E=exact_data();V=E['factor_vertices'];seed=E['seed_vectors']
 verts=[[sum(V[k][a[k]][i]for k in range(3))for i in range(19)]for a in product(range(7),range(7),range(8))]
 def dot(a,b):return sum(x*y for x,y in zip(a,b))-sum((a[i]*b[i+1]+a[i+1]*b[i])/2 for lo,hi in[(0,6),(6,12),(12,19)]for i in range(lo,hi-1))
 assert all(dot(v,v)==196 for v in verts)
 k=(F(33,32),F(-1,32));assert sg(k)>0 and sg(sub(ONE,k))>0
 assert mul(k,k)==sub(sc(k,F(33,16)),ONE)
 s=[ONE]
 for j in range(1,5):s.append(mul(s[-1],k))
 rays=[a*56+a*8+a for a in range(7)];assert rays==[65*a for a in range(7)]
 assert min(392-2*dot(verts[i],verts[j])for n,i in enumerate(rays)for j in rays[:n])==280
 caps=[];counts=Counter()
 assert min(map(F,base['vertex_distance_distribution']))==84
 # Every literal product vertex has the same19 right and3 omitted seed gaps.
 for z in verts:
  gaps=[1-dot(z,w)for w in seed]
  assert gaps.count(0)==19 and sorted(x for x in gaps if x)==[7,7,8]
 gs={j:sub(ONE,s[j])for j in[1,2,3]}
 # Pair classes use squared old-vertex distance d. All coefficients of d below
 # are positive, so the exact lower-distance gates cover the entire orbits.
 gates=[]
 for i,j,dmin,label in[(1,1,280,'dotted'),(1,2,0,'pi3'),(1,2,84,'dotted'),(1,3,0,'dotted'),(2,2,84,'dotted'),(2,3,0,'pi3'),(2,3,84,'dotted'),(3,3,84,'dotted')]:
  ti,tj=s[i],s[j];D=add(sc(mul(sub(ti,tj),sub(ti,tj)),196),sc(mul(gs[i],gs[j]),dmin));qi=sc(mul(ti,ti),4);qj=sc(mul(tj,tj),4)
  p=sc(sub(add(qi,qj),D),F(1,2));rprod=sc(mul(ti,tj),4)
  assert sg(gs[i])>0 and sg(gs[j])>0
  if label=='pi3':assert p==sc(rprod,F(-1,2))
  else:assert sg(add(p,rprod))<0
  gates.append(dict(levels=[i,j],old_distance_lower_bound=dmin,relation=label,pairing=p))
 for j in[1,2,3]:
  t=s[j];g=gs[j];r=sc(t,2);norm=mul(r,r);N=sc(mul(g,g),196)
  # Retained seed contacts have cosine1/2; the omitted gaps are strictlydotted.
  for h in[7,8]:assert sg(sub(sub(q(h),sc(t,h-1)),r))>0
  p=sc(sub(add(q(4),norm),N),F(1,2))
  if j==1:assert p==neg(r)
  else:assert sg(add(p,sc(r,2)))<0
  # The own old vertex is the nearest old vertex to this center; every other
  # old-vertex squared distance adds (1-t)*d with d>=0.
  ownD=sc(mul(t,t),196);p=sc(sub(add(norm,q(28)),ownD),F(1,2))
  assert sg(p)<0 and sg(sub(mul(p,p),sc(norm,28)))>0
  assert sg(sub(norm,add(ownD,q(28))))<0
  assert sg(add(sub(N,norm),q(400)))>0
  for ray in rays if j==1 else range(392):
   z=verts[ray];center=[sc(g,c)for c in z]
   caps.append(dict(ray=ray,level=j,gamma=g,radius=r,norm=norm,center=center,D=sc(sub(N,norm),F(1,2)),U=1))
 counts.update(pi3=791*19+7+7+392,dotted=791*3+784+791*392+21+(7*392-7)+7*392+2*(392*391//2)+(392*392-392))
 assert len(caps)==791 and sum(counts.values())==640710
 totalcounts=Counter(pi3=828,right=7657,parallel=392,dotted=77028)
 assert sum(totalcounts.values())==415*414//2
 totalcounts.update(counts)
 assert sum(totalcounts.values())==1206*1205//2
 # Level2 strictly cuts the prior universal .4-spatial norm1/4 escape on every
 # old ray; all392 such caps belong to the coherent selected chamber.
 lam=F(2,5);g=sub(ONE,s[2]);P=sub(sc(mul(sub(q(lam),g),sub(q(lam),g)),196),sc(mul(s[2],s[2]),4))
 assert sg(P)<0
 # The common .6 probe survives every radius layer on every ray, by the exact
 # identity P(lambda*z0,g*zv)=P(lambda*z0,g*z0)+lambda*g*|z0-zv|².
 lam=F(3,5);powers=[]
 for j in[1,2,3]:
  g=sub(ONE,s[j]);p=sub(sc(mul(sub(q(lam),g),sub(q(lam),g)),196),sc(mul(s[j],s[j]),4))
  assert sg(sub(p,q(F(1,8))))>0 and sg(g)>0;powers.append(p)
 oldvertexpower=196*(1-lam)**2-28;assert oldvertexpower==F(84,25)>F(1,8)
 x2=196*lam*lam;D=(x2-F(1,8))/2;assert D==F(14087,400)and -200-D<0
 assert min(1-lam*dot(verts[0],w)for w in seed)==F(2,5)
 # The fourth cap on any ray conflicts with its own retained old vertex cap:
 # its normal pairing is strictly positive.
 t=s[4];N4=sc(mul(t,t),192);a=sub(N4,q(28));assert sg(a)<0 # in fact positive normal pairing
 cycle=mul(sub(sc((F(33,32),F(1,32)),7),q(6)),sub(sc((F(33,32),F(1,32)),7),q(6)))
 assert cycle==(F(2353,512),F(273,512))
 # Source-independent exact normal identities in the two null coordinates.
 for cap in caps:
  g=cap['gamma'];assert sub(sc(mul(g,g),196),sc(cap['D'],2))==cap['norm']
 paths=[BASE,'scripts/fresh_h20_balanced_20261014.py',f'scripts/{TAG}.py']
 return enc(dict(status='PASS_EXACT_UNIT_FULL_OUTER_CHAIN_ORBITS_WITH_FULL_BANK_ESCAPE',sources={p:sha(p)for p in paths},field='Q(sqrt65), coefficient order [1,sqrt65]',k=k,k_relation='k²-(33/16)k+1=0, 0<k<1',baseline_wall_count=415,previous_three_cap_motif_retained=False,selected_rays=rays,selected_ray_minimum_distance_squared=280,new_caps=caps,new_pair_counts=dict(counts),total_pair_counts=dict(totalcounts),pair_class_gates=gates,coherent_wall_count=1206,preserved_finite_endpoints=392,actual_facets='The common center (0,200,1) is strict; acute distinct normals imply strict facet points by orthogonal projection.',old_gap_cut=dict(spatial_fraction='2/5',chain_level=2,ball_power=P,number_of_selected_probes=392),relaxed_all_ray_bank=dict(rays=392,levels=[1,2,3],caps=1176,mutual_compatibility_not_assumed=True),universal_escape=dict(spatial_fraction='3/5',reference_vertex=0,spatial_square=x2,D=D,U=1,norm=F(1,8),seed_minimum_gap=F(2,5),central_ball_power=x2-4,vertex_layer_minimum_power=oldvertexpower,chain_layer_minimum_powers=powers,future_center_pairing=-200-D),level4_own_vertex_pairing_sign='positive',irrational_cartan_two_cycle=cycle,rational_indefinite_restriction=[[0,-7],[-7,4]],scope='Exact 791-cap coherent contact-chain construction retaining the 415-wall baseline, without the prior three-cap motif. The larger 1176-cap relaxed bank still has a strict escape. No complete new cusp or finite-vertex census and no finite-volume claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';d=run()
 if a.verify:assert json.loads(p.read_text())==d
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'])
