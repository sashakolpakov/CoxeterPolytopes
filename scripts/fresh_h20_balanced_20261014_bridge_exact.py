#!/usr/bin/env python3
"""Exact central-contact motif and all-label interior-center radius-five gate.

This is a construction/necessary-cover certificate, not a finite-volume claim.
All rational arithmetic below is Python Fraction; no numerical proposal is trusted.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter
import argparse,json,hashlib
from pathlib import Path
from fresh_h20_balanced_20261014 import ROOT,exact_data,sha
TAG='fresh_h20_balanced_20261014_bridge_exact'
SOURCES=['scripts/fresh_h20_balanced_20261014.py',f'scripts/{TAG}.py']
def add(a,b):return [x+y for x,y in zip(a,b)]
def scale(a,t):return [t*x for x in a]
def sub(a,b):return [x-y for x,y in zip(a,b)]
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def run():
 E=exact_data();B=E['spatial_Gram'];V=E['factor_vertices'];seed=E['seed_vectors']
 def dot(a,b):return sum(a[i]*b[i]for i in range(19))-sum((a[i]*b[i+1]+a[i+1]*b[i])/2 for lo,hi in[(0,6),(6,12),(12,19)]for i in range(lo,hi-1))
 assert all(B[i][j]==(F(i==j)-F(abs(i-j)==1 and (i//6==j//6 if max(i,j)<12 else min(i,j)>=12),2))for i in range(19)for j in range(19))
 # Sylvester pivots for finite A_n unit Grams are (i+1)/(2i).
 pivots=[F(i+1,2*i)for n in E['ranks']for i in range(1,n+1)];assert all(x>0 for x in pivots)
 verts=[add(add(V[0][a],V[1][b]),V[2][c])for a,b,c in product(range(7),range(7),range(8))]
 assert len(set(tuple(v)for v in verts))==392
 assert all(dot(v,v)==196 for v in verts)
 assert all(all(x.denominator==1 for x in v)for v in verts)
 def dint(a,b):
  return sum(a[i]*b[i]for i in range(19))-sum((a[i]*b[i+1]+a[i+1]*b[i])//2 for lo,hi in[(0,6),(6,12),(12,19)]for i in range(lo,hi-1))
 # Work with doubled dot products to avoid assuming each summand even.
 def twiceint(a,b):return 2*sum(x*y for x,y in zip(a,b))-sum(a[i]*b[i+1]+a[i+1]*b[i]for lo,hi in[(0,6),(6,12),(12,19)]for i in range(lo,hi-1))
 vi=[[int(x)for x in v]for v in verts];ds=Counter()
 for i,v in enumerate(vi):
  for w in vi[:i]:
   D=F(392-twiceint(v,w));ds[str(D)]+=1
   assert D==84 or D>=112
 assert min(map(F,ds))==84 and sum(ds.values())==392*391//2
 # Seed plus centralR2 plus complete vertexR=sqrt28 is coherent.
 seedpairs=Counter()
 for i,w in enumerate(seed):
  assert dot(w,w)==1
  for r in seed[:i]:
   p=dot(w,r);assert p in[0,F(-1,2)];seedpairs['right'if p==0 else'pi3']+=1
 for z in verts:
  gaps=[1-dot(z,w)for w in seed]
  assert gaps.count(0)==19 and sorted(x for x in gaps if x)==[7,7,8]
  assert all(x==0 or x*x>28 for x in gaps)
 assert (196-4-28)**2>4*4*28
 z0=verts[0];z2=verts[2]
 bary=add(add(V[0][0],V[1][0]),scale(add(add(V[2][0],V[2][2]),V[2][5]),F(1,3)))
 new=[dict(name='face_bridge',center=scale(bary,F(3,11)),radius=F(16,11)),dict(name='radial_000',center=scale(z0,F(1,4)),radius=F(3,2)),dict(name='radial_002',center=scale(z2,F(1,4)),radius=F(3,2))]
 assert dot(bary,bary)==F(364,3)
 newcounts=Counter();records=[]
 for k,row in enumerate(new):
  c,R=row['center'],row['radius'];N=dot(c,c);gaps=[1-dot(c,w)for w in seed]
  labels=[]
  for g in gaps:
   co=g/R
   if co==0:typ='right'
   elif co==F(1,2):typ='pi3'
   elif co==1:typ='parallel'
   else:assert co>1;typ='dotted'
   labels.append(typ);newcounts[typ]+=1
  central=(N-R*R-4)/(4*R)
  assert central==[F(1,2),F(1),F(1)][k]
  newcounts['pi3'if central==F(1,2)else'parallel']+=1
  for z in verts:
   d=dot(sub(c,z),sub(c,z));p=d-R*R-28
   assert p>0 and p*p>4*R*R*28;newcounts['dotted']+=1
   # Old finite vertex (z,(z²+28)/2,1) remains strict.
   assert R*R-d-28<0
  for old in new[:k]:
   d=dot(sub(c,old['center']),sub(c,old['center']));co=(d-R*R-old['radius']**2)/(2*R*old['radius'])
   typ='right'if co==0 else'dotted'
   assert co==0 or co>1;newcounts[typ]+=1
  records.append(dict(**row,norm=R*R,spatial_square=N,seed_labels=labels,central_cosine=central))
 cycle=4*(1-dot(new[0]['center'],seed[0]))**2/new[0]['radius']**2
 assert cycle==F(841,64) and cycle.denominator!=1
 # A strict common center H=(0,200,1), H²=-400.
 assert all(400+dot(row['center'],row['center'])-row['radius']**2>0 for row in new)
 # E=(2z0/5,(|2z0/5|²-1)/2,1) has square+1.
 x=scale(z0,F(2,5));x2=dot(x,x);assert x2==F(784,25)
 assert all(1-dot(x,w)>0 for w in seed)
 powers=[x2-4]+[dot(sub(x,z),sub(x,z))-28 for z in verts]+[dot(sub(x,r['center']),sub(x,r['center']))-r['radius']**2 for r in new]
 assert min(powers)>1
 # The radius-five lemma uses only these exact constants.
 # min(14g+2(1-g),2+4(1-g)) <=5, crossing atg=1/4.
 assert 12*F(1,4)+2==6-4*F(1,4)==5
 assert F(2,5)*14-5==F(3,5)
 out=dict(status='PASS_EXACT_BALANCED_CENTRAL_BRIDGE_MOTIF_AND_RADIUS_FIVE_GATE',sources={p:sha(p)for p in SOURCES},signature=[20,1],spatial_positive_pivots=pivots,old_vertex_count=392,old_vertex_square=196,vertex_layer_radius_squared=28,vertex_distance_distribution=dict(ds),seed_pair_counts=dict(seedpairs),new_caps=records,new_pair_counts=dict(newcounts),coherent_wall_count=418,common_center=dict(spatial=[0]*19,D=200,U=1,norm=-400),strict_escape=dict(spatial=x,D=(x2-1)/2,U=1,norm=1,minimum_ball_power=min(powers),center_pairing=-200-(x2-1)/2),preserved_finite_endpoints=392,nonintegral_cartan_two_cycle=dict(seed_index=0,cap='face_bridge',value=cycle),actual_facets='All418: a strict timelike center and pairwise nonpositive distinct spacelike normals give a strict facet point by orthogonal projection.',radius_five_gate=dict(central_radius=2,vertex_norm_bound=14,scope='Any Coxeter-compatible hemisphere with center strictly inside the unit A6^2 x A7 alcove and finite or parallel contact with centralR2.',proof='Let d>0 be the least of22 seed distances. Factor sums force d<=1; center/(1-d) belongs toC. Coxeter necessity gives radius<=2d. Putg=1-d: norm(center)+radius<=min(12g+2,6-4g)<=5. Thed=1 center0 cannot have acute central contact.',missed_boundary_point=x,uniform_ball_power_lower_bound=F(9,25)),scope='Exact coherent three-cap contact motif on the central+fullvertex seed, still infinite volume. No global endpoint/cusp census or arbitrary boundary-center cap obstruction. Numerical scouts are separate.')
 return enc(out)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';out=run()
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['new_pair_counts'],'minimum escape power',out['strict_escape']['minimum_ball_power'])
