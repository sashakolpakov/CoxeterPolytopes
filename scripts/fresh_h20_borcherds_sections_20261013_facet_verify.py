#!/usr/bin/env python3
"""Optimizer-free exact exclusion of all center/facet section pencils."""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_borcherds_sections_20261013_facet_verify'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
SYMS='data/enumerations/fresh_h20_global_interior_sections_20261013.json'
REC='results/fresh_h20_borcherds_sections_20261013_facet_recipes.json'
OUT=ROOT/'results'/f'{TAG}.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(A,v):return [dot(r,v)for r in A]
def pair(a,b,Q):return dot(a,mv(Q,b))
def add(a,b):return [x+y for x,y in zip(a,b,strict=True)]
def scale(c,a):return [c*x for x in a]
def enc(v):return list(map(str,v))
def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(v)for v in r]for r in raw['basis_gram']];W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 C=[mv(Q,r)for r in W];G=[[dot(a,c)for c in C]for a in W]
 assert all(G[i][j]<=0 for i in range(210)for j in range(i))
 v=[sum(r[j]for r in W)for j in range(22)];gg=math.gcd(*v);h=[x//gg for x in v];hp=[dot(h,c)for c in C];assert pair(h,h,Q)==-14 and set(hp)=={-1,-8}
 gens=json.loads((ROOT/SYMS).read_text())['ambient_symmetry_generators'];perms=[]
 for rec in gens:
  A=rec['matrix'];p=rec['wall_permutation'];assert sorted(p)==list(range(210))
  assert mv(A,h)==h and all(mv(A,r)==W[p[i]]for i,r in enumerate(W))
  # Literal ambient form, not only a combinatorial diagram action.
  assert [[sum(A[k][i]*Q[k][l]*A[l][j]for k in range(22)for l in range(22))for j in range(22)]for i in range(22)]==Q
  perms.append(p)
 cases=[];covered=set()
 for rec in json.loads((ROOT/REC).read_text()):
  i=rec['wall'];j,k=rec['ridge_pair'];T=F(rec['T']);x=list(map(F,rec['maximum_point']));weights={int(i):F(w)for i,w in rec['dual_weights']}
  n=G[i][i];d=-hp[i]
  assert [n,d,T]in([2,1,F(1,3)],[4,8,F(1)])
  assert all(w>0 for w in weights.values())
  assert all(T*h[a]-W[i][a]==sum(w*W[l][a]for l,w in weights.items())for a in range(22))
  assert pair(x,h,Q)==-1 and pair(x,x,Q)<=0 and dot(x,C[i])==-T
  xp=[dot(x,c)for c in C];assert max(xp)<=0
  assert G[j][j]==G[k][k]==2 and hp[j]==hp[k]==-1
  assert G[i][j]==G[i][k]==G[j][k]==xp[j]==xp[k]==0
  # All three old walls are orthogonal, so these are exact relative-interior
  # points of the facet and triple face. The common-center projection also
  # proves all other inequalities strict, directly checked here.
  zi=add(h,scale(F(d,n),W[i]));z=add(add(zi,scale(F(1,2),W[j])),scale(F(1,2),W[k]))
  for p,zeros in((zi,[i]),(z,sorted([i,j,k]))):
   pp=[dot(p,c)for c in C]
   assert pair(p,p,Q)<0 and pair(p,h,Q)<0 and [l for l,y in enumerate(pp)if y==0]==zeros
   assert all(y<0 for l,y in enumerate(pp)if l not in zeros)
  zi=scale(-1/pair(zi,h,Q),zi);z=scale(-1/pair(z,h,Q),z)
  # p(t)=(T-t)z+t*x is timelike for 0<t<T, vanishes only at j,k,
  # and lies in (r_i-t*h)^perp. The analogous zi combination is strict.
  for p in (zi,z):
   assert pair(p,h,Q)==-1 and dot(p,C[i])==0
  e0=W[i];e1=scale(-1,h);p0=scale(T,z);p1=add(x,scale(-1,z))
  assert pair(p0,e0,Q)==0 and pair(p0,e1,Q)+pair(p1,e0,Q)==0 and pair(p1,e1,Q)==0
  assert all(dot(p0,C[l])==dot(p1,C[l])==0 for l in(j,k))
  # q=n+2dt-14t², projected pair: diagonal (2q-t²)/q,
  # offdiagonal -t²/q. Concavity and these endpoint bounds imply
  # 2q-3t²>0 throughout [0,T], hence 0<cos<1/2 on (0,T).
  endpoints=[F(2*n),F(2*n)+4*d*T-31*T*T];assert min(endpoints)>0
  orbit={i};queue=[i]
  for a in queue:
   for perm in perms:
    b=perm[a]
    if b not in orbit:orbit.add(b);queue.append(b)
  assert all(G[a][a]==n for a in orbit) and covered.isdisjoint(orbit);covered|=orbit
  cases.append({'representative':i,'orbit_size':len(orbit),'normal':'r_i-t*h','complete_interior_parameter_interval':['0',str(T)],
                'upper_bound_dual_weights':[[l,str(w)]for l,w in sorted(weights.items())],
                'maximizing_point':enc(x),'maximizing_point_norm':str(pair(x,x,Q)),
                'ridge_walls':[j,k],'triple_face_point':enc(z),'facet_point':enc(zi),
                'ridge_point_polynomial':[enc(p0),enc(p1)],'section_normal_norm_polynomial':[n,2*d,-14],
                'ridge_cosine':'t^2/(2*n+4*d*t-29*t^2)',
                'strict_cosine_upper_bound_concave_polynomial':[2*n,4*d,-31],
                'positive_endpoint_bounds':enc(endpoints)})
 assert covered==set(range(210))
 return {'status':'PASS_ALL_REAL_210_CENTER_FACET_SECTION_PENCILS_EXCLUDED',
   'sources':{p:sha(p)for p in(RAW,SYMS,REC,'scripts/'+Path(__file__).name)},'cases':cases,
   'scope':'All genuine interior H20 sections of the frozen H21 chamber whose normal belongs to span(h,r_i) for any of210 original walls. Every such section has an actual ridge with acute cosine strictly between0 and1/2. All real projective coefficients are covered; no field/label cutoff. Other sections and general H20 polytopes are not excluded.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run()
 if a.verify:assert json.loads(OUT.read_text())==d
 else:OUT.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
