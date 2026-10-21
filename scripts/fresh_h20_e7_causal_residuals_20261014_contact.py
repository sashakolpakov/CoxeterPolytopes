#!/usr/bin/env python3
"""Exact pi3 contact refinement, deleting only old wall61."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from itertools import combinations
from collections import Counter
from hashlib import sha256
import fresh_h20_e7_causal_residuals_20261014_exact as P
K=P.K;ROOT=P.ROOT;TAG=P.TAG;RAD=K.RAD+[43*x for x in K.RAD];Z=(F(0),)*16;DEN=1<<192
BOUNDS=[]
for q in RAD:
 k=isqrt(q*DEN*DEN);BOUNDS.append((F(k,DEN),F(k+int(k*k!=q*DEN*DEN),DEN)))
def ext(x):return tuple(x)+(F(0),)*8
def scalar(x):return(F(x),)+(F(0),)*15
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,q):return tuple(q*x for x in a)
def mul(a,b):
 c=[F(0)]*16
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if x and y:c[i^j]+=x*y*RAD[i&j]
 return tuple(c)
def total(values):
 a=Z
 for b in values:a=add(a,b)
 return a
def sign(a):
 if a==Z:return 0
 lo=sum(x*(u if x>=0 else v)for x,(u,v)in zip(a,BOUNDS));hi=sum(x*(v if x>=0 else u)for x,(u,v)in zip(a,BOUNDS))
 assert lo>0 or hi<0,('unresolved',a,lo,hi)
 return 1 if lo>0 else -1
def lin(vs,cs):return[total(mul(a,v[k])for a,v in zip(cs,vs))for k in range(21)]
def rational(a):assert not any(a[1:]);return a[0]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d,f,Q,old0,H0,E0,S,A,J,I,omega0,lam0,T0,w0,y0,r0,_,_,_,_=P.setup()
 old=[[ext(x)for x in v]for v in old0];H=list(map(ext,H0));oldE=list(map(ext,E0));T=list(map(ext,T0));w=list(map(ext,w0));lam=[[ext(x)for x in v]for v in lam0]
 cache={};sp=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
 def dot(v,u):
  key=tuple(u)
  if key not in cache:cache[key]=[tuple(sum(q*u[j][k]for j,q in row)for k in range(16))for row in sp]
  return total(mul(a,b)for a,b in zip(v,cache[key]))
 s2=ext((F(0),F(1))+(F(0),)*6);x=tuple(F(1,11)if i==13 else F(0)for i in range(16));t=tuple(F(17,11)if i==4 else F(0)for i in range(16))
 assert mul(x,x)==scalar(F(86,11))and mul(t,t)==scalar(F(289,11))and sign(x)>0 and sign(t)>0
 assert mul(t,t)==scale(add(scalar(71),mul(x,x)),F(1,3))and sign(add(x,scalar(-2)))>0
 y=lin([lam[0],lam[12],T],[scale(s2,-1),scale(x,-1),t]);r=lin([w,y],[scalar(-1),scalar(1)])
 assert dot(r,r)==scalar(2)and dot(y,y)==scalar(F(-91,2))and sign(dot(y,T))<0
 assert dot(r,old[75])==scalar(-1)and sign(dot(r,old[61]))>0
 retained=[i for i in range(87)if i!=61];N=[old[i]for i in retained]+[r];norms=[rational(dot(n,n))for n in N]
 counts=Counter();pairs=[]
 for j,i in enumerate(retained):
  g=dot(r,old[i]);assert sign(g)<=0;prod=2*norms[j]
  if g==Z:kind='right'
  elif mul(g,g)==scalar(prod/4):kind='pi3'
  elif mul(g,g)==scalar(prod/2):kind='pi4'
  elif mul(g,g)==scalar(3*prod/4):kind='pi6'
  elif mul(g,g)==scalar(prod):kind='parallel'
  else:assert sign(add(mul(g,g),scalar(-prod)))>0;kind='dotted'
  counts[kind]+=1;pairs.append(dict(walls=[86,j],old_wall=i,pairing=g,relation=kind))
 assert sign(dot(H,H))<0 and all(sign(dot(H,n))<0 for n in N)
 fp=[]
 for i,n in enumerate(N):
  p=lin([H,n],[scalar(1),scale(dot(H,n),-1/norms[i])]);assert sign(dot(p,p))<0 and sign(dot(p,H))<0
  assert all(dot(p,v)==Z if i==j else sign(dot(p,v))<0 for j,v in enumerate(N));fp.append(p)
 # An actual pi3 ridge to retained old wall75, with an explicit strict
 # relative-interior point obtained by orthogonal projection.
 hr=dot(H,r);hn=dot(H,old[75]);ridge=lin([H,r,old[75]],[scalar(1),scale(add(scale(hr,2),hn),F(-1,3)),scale(add(hr,scale(hn,2)),F(-1,3))])
 rid=[retained.index(75),86];assert sign(dot(ridge,ridge))<0 and sign(dot(ridge,H))<0
 assert all(dot(ridge,n)==Z if i in rid else sign(dot(ridge,n))<0 for i,n in enumerate(N))
 # E7 projection of deleted61 is timelike; all fifteen source facets
 # remain. The new projection y is future timelike and redundant on F.
 op=[dot(old[61],old[i])for i in S];coeff=[total(scale(op[j],J[i][j])for j in range(7))for i in range(7)]
 projection=lin([old[61]]+[old[i]for i in S],[scalar(1)]+[scale(z,-1)for z in coeff])
 assert sign(dot(projection,projection))<0 and sign(dot(projection,T))<0
 assert 61 not in S+f['face_facet_walls']
 vertices=[]
 for row in f['face_finite_vertices']:
  point=[ext(z)for z in K.vector(row['ambient_point'])];assert sign(dot(point,r))<0
  vals=[dot(point,n)for n in N];assert all(sign(z)<=0 for z in vals)
  assert [retained[i]for i in range(86)if vals[i]==Z]==row['ambient_incident'];vertices.append(dict(old_incident=row['ambient_incident'],new_pairing=dot(point,r)))
 links=[]
 for row in d['completed_cusp_links']:
  nu=[ext(z)for z in K.vector(row['point'])];assert sign(dot(nu,old[61]))<0 and sign(dot(nu,r))<0
  assert all(sign(dot(nu,n))<=0 for n in N)
  assert [retained[i]for i in range(86)if dot(nu,N[i])==Z]==row['incident']
  links.append(dict(point=nu,old_incident=row['incident'],new_pairing=dot(nu,r),inherited_positive_affine_quotient_rank=19))
 # Complete possible finite-vertex supports on the new facet. Reject
 # an old parallel/dotted pair first; the remaining old19 supports have
 # rational Gram and their exact Schur complement is negative.
 finite=[row['old_wall']for row in pairs if row['relation']in ['right','pi3','pi4','pi6']];assert len(finite)==20
 supports=[]
 for ids in combinations(finite,19):
  blocker=None
  for i,j in combinations(ids,2):
   g=dot(old[i],old[j]);prod=F(d['norms'][i])*F(d['norms'][j])
   if sign(add(mul(g,g),scalar(-prod)))>=0:blocker=[i,j];break
  if blocker:supports.append(dict(old_support=list(ids),nonfinite_old_pair=blocker));continue
  B=[[rational(dot(old[i],old[j]))for j in ids]for i in ids]
  test=[row[:]for row in B];obstruction=None
  for k in range(19):
   pivot=test[k][k]
   if pivot<=0:obstruction=dict(index=k,nonpositive_pivot=pivot);break
   for i in range(k+1,19):
    for j in range(k+1,19):test[i][j]-=test[i][k]*test[j][k]/pivot
  if obstruction:supports.append(dict(old_support=list(ids),old_Gram_not_PD=obstruction));continue
  BI=K.inv(B);b=[dot(r,old[i])for i in ids]
  residual=add(scalar(2),scale(total(scale(mul(b[i],b[j]),BI[i][j])for i in range(19)for j in range(19)),-1))
  assert sign(residual)<0;supports.append(dict(old_support=list(ids),negative_schur_residual=residual))
 assert len(supports)==20 and sign(dot(r,oldE))>0
 escape_path='results/fresh_h20_e7_replacement_escape_20261014.json';ed=K.read(escape_path);E=[ext(z)for z in K.vector(ed['point'])];ep=[dot(E,n)for n in N]
 assert sign(dot(E,E))>0 and sign(dot(E,H))<0 and all(sign(z)<0 for z in ep)
 paths=[P.BASE,P.FACE,escape_path,f'scripts/{TAG}_exact.py',f'scripts/{TAG}_contact.py','scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps_review.py']
 return K.enc(dict(status='PASS_REFINED_87_ONE_DELETION_ACTUAL_PI3_RIDGE_GLOBAL_ESCAPE',sources={p:sha(p)for p in paths},field='Q(sqrt2,sqrt79,sqrt11,sqrt43)',square_root_basis=RAD,form=Q,normals=N,norms=norms,removed_old87_wall=61,retained_old87_walls=retained,recipe=dict(E7_positions=[2,4,6],source_contact0=s2,source_contact12=x,time_coefficient=t),normal=r,face_projection=y,face_projection_norm=dot(y,y),deleted_wall_pairing=dot(r,old[61]),new_pair_counts=dict(counts),new_pairs=pairs,strict_center=H,facet_points=fp,actual_pi3_ridge=dict(walls=rid,old_neighbor=75,point=ridge,norm=dot(ridge,ridge)),deleted_face_projection_norm=dot(projection,projection),preserved_face_vertices=vertices,preserved_bounded_cusps=links,new_facet_finite_vertex_exclusion=supports,old_escape_cut=dot(r,oldE),strict_escape=dict(point=E,norm=dot(E,E),all87_pairings=ep),scope='Exactly one old off-face wall removed. New norm2 cap has an actual pi3 ridge with retained75; all pairs/facets exact. Complete H13 face,44 finitevertices and2boundedcusps unchanged. Exhaustive20 possible new finite-vertex supports fail. Saved old escape is cut but a new strict escape survives. No finite-volume completion claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_contact.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'])
