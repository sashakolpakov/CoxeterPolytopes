#!/usr/bin/env python3
"""Mixed-label E7 cap joining both retained off-face anchors."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from itertools import combinations
from collections import Counter
from hashlib import sha256
import fresh_h20_e7_causal_residuals_20261014_exact as P
K=P.K;ROOT=P.ROOT;TAG=P.TAG;RAD=K.RAD+[3*x for x in K.RAD];Z=(F(0),)*16;DEN=1<<192
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
def div(a,b):
 M=[[mul(b,tuple(F(k==j)for k in range(16)))[i]for j in range(16)]for i in range(16)];I=K.inv(M)
 return tuple(sum(I[i][j]*a[j]for j in range(16))for i in range(16))
def run():
 d,f,Q,old0,H0,E0,S,A,J,I,omega0,lam0,T0,_,_,_,_,_,_,dot0=P.setup()
 old=[[ext(z)for z in v]for v in old0];H=list(map(ext,H0));T=list(map(ext,T0));lam=[[ext(z)for z in v]for v in lam0];cache={};sp=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
 def bdot(v,u):
  key=tuple(u)
  if key not in cache:cache[key]=[tuple(sum(q*u[j][k]for j,q in row)for k in range(16))for row in sp]
  return total(mul(a,b)for a,b in zip(v,cache[key]))
 s3=tuple(F(i==8)for i in range(16));w=lin([[ext(z)for z in omega0[i]]for i in [4,6]],[s3,s3]);wq=bdot(w,w)
 assert all(bdot(lam[j],old[75])==Z for j in [0,12])
 tt=div(add(bdot(w,old[75]),scalar(-1)),bdot(T,old[75]));assert sign(tt)>0
 ma=scale(bdot(lam[0],old[61]),-1);mb=scale(bdot(lam[12],old[61]),-1)
 rhs=add(add(scalar(-1),bdot(w,old[61])),scale(mul(tt,bdot(T,old[61])),-1))
 u=div(rhs,ma);v=div(scale(mb,-1),ma);qA=F(11,12);qB=F(2,3)
 target=add(add(scale(mul(tt,tt),2),scalar(2)),scale(wq,-1))
 qa=add(scale(mul(v,v),qA),scalar(qB));qb=scale(mul(u,v),2*qA);qc=add(scale(mul(u,u),qA),scale(target,-1));disc=add(mul(qb,qb),scale(mul(qa,qc),-4))
 assert sign(qa)>0 and sign(disc)>0
 EZ=(Z,Z)
 def ee(x):return(x,Z)
 def ea(x,y):return(add(x[0],y[0]),add(x[1],y[1]))
 def es(x,q):return(scale(x[0],q),scale(x[1],q))
 def em(x,y):return(add(mul(x[0],y[0]),mul(disc,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0])))
 def sg(x):
  a,b=x;sa,sb=sign(a),sign(b)
  if not sb:return sa
  if not sa or sa==sb:return sb
  return sa*sign(add(mul(a,a),scale(mul(disc,mul(b,b)),-1)))
 def et(xs):
  out=EZ
  for x in xs:out=ea(out,x)
  return out
 def elin(vs,cs):return[et(em(a,v[k])for a,v in zip(cs,vs))for k in range(21)]
 ecache={}
 def dot(v,u):
  key=tuple(u)
  if key not in ecache:ecache[key]=[et(es(u[j],q)for j,q in row)for row in sp]
  return et(em(a,b)for a,b in zip(v,ecache[key]))
 olde=[[ee(z)for z in v]for v in old];he=[ee(z)for z in H]
 # Select the quadratic branch by the exact inequalities a>2,b>2 and
 # positivity at the old escape. All acceptance gates below are exact.
 recipes=[]
 for branch in [-1,1]:
  bb=(div(scale(qb,-1),scale(qa,2)),div(scalar(branch),scale(qa,2)))
  aa=ea(ee(u),em(ee(v),bb));t=ee(tt)
  assert sg(ea(ea(es(em(aa,aa),qA),es(em(bb,bb),qB)),es(ee(target),-1)))==0
  if sg(ea(aa,ee(scalar(-2))))<=0 or sg(ea(bb,ee(scalar(-2))))<=0:continue
  y=elin([[ee(z)for z in lam[0]],[ee(z)for z in lam[12]],[ee(z)for z in T]],[es(aa,-1),es(bb,-1),t])
  r=elin([[ee(z)for z in w],y],[ee(scalar(-1)),ee(scalar(1))])
  ep=dot(r,[ee(ext(z))for z in E0])
  if sg(ep)>0:recipes.append((branch,aa,bb,t,y,r))
 assert len(recipes)==1
 branch,aa,bb,t,y,r=recipes[0];assert sg(ea(dot(r,r),ee(scalar(-2))))==0
 assert sg(dot(y,y))<0 and sg(dot(y,[ee(z)for z in T]))<0
 assert sg(ea(dot(r,olde[61]),ee(scalar(1))))==sg(ea(dot(r,olde[75]),ee(scalar(1))))==0
 N=olde+[r];norms=list(map(F,d['norms']))+[F(2)];counts=Counter();pairs=[]
 for j,n in enumerate(olde):
  g=dot(r,n);assert sg(g)<=0;prod=2*norms[j]
  if sg(g)==0:kind='right'
  elif sg(ea(em(g,g),ee(scalar(-prod/4))))==0:kind='pi3'
  elif sg(ea(em(g,g),ee(scalar(-prod/2))))==0:kind='pi4'
  elif sg(ea(em(g,g),ee(scalar(-3*prod/4))))==0:kind='pi6'
  elif sg(ea(em(g,g),ee(scalar(-prod))))==0:kind='parallel'
  else:assert sg(ea(em(g,g),ee(scalar(-prod))))>0;kind='dotted'
  counts[kind]+=1;pairs.append(dict(walls=[87,j],pairing=g,relation=kind))
 assert sg(dot(he,he))<0 and all(sg(dot(he,n))<0 for n in N)
 # Full old87 acute-pair certificate plus these87 new gates and a strict
 # center prove every facet actual: its center projection is timelike,
 # and every other inequality becomes strictly more negative.
 fp=[]
 for i,n in enumerate(N):
  h=dot(he,n);point=elin([he,n],[ee(scalar(1)),es(h,-1/norms[i])]);assert sg(dot(point,n))==0 and sg(dot(point,point))<0 and sg(dot(point,he))<0
  fp.append(point)
 ridges=[]
 for j in [61,75]:
  assert F(d['norms'][j])==2
  hr=dot(he,r);hn=dot(he,olde[j]);point=elin([he,r,olde[j]],[ee(scalar(1)),es(ea(es(hr,2),hn),F(-1,3)),es(ea(hr,es(hn,2)),F(-1,3))])
  assert sg(dot(point,point))<0 and sg(dot(point,he))<0
  assert all(sg(dot(point,n))==0 if i in [j,87]else sg(dot(point,n))<0 for i,n in enumerate(N))
  ridges.append(dict(walls=[j,87],point=point,norm=dot(point,point)))
 vertices=[]
 for row in f['face_finite_vertices']:
  point=[ee(ext(z))for z in K.vector(row['ambient_point'])];g=dot(point,r);assert sg(g)<0;vertices.append(dict(old_incident=row['ambient_incident'],new_pairing=g))
 links=[]
 for row in d['completed_cusp_links']:
  point=[ee(ext(z))for z in K.vector(row['point'])];g=dot(point,r);assert sg(g)<0;links.append(dict(old_incident=row['incident'],new_pairing=g,inherited_quotient_rank=19))
 # Complete twenty possible vertex supports on the new facet.
 finite=[j for j,row in enumerate(pairs)if row['relation']in ['right','pi3','pi4','pi6']];assert len(finite)==20
 supports=[]
 for ids in combinations(finite,19):
  blocker=None
  for i,j in combinations(ids,2):
   g=dot0(old0[i],old0[j]);prod=F(d['norms'][i])*F(d['norms'][j])
   if K.sign(K.sub(K.mul(g,g),K.c(prod)))>=0:blocker=[i,j];break
  if blocker:supports.append(dict(old_support=list(ids),nonfinite_old_pair=blocker));continue
  B=[[K.rational(dot0(old0[i],old0[j]))for j in ids]for i in ids];test=[row[:]for row in B];obstruction=None
  for k in range(19):
   pivot=test[k][k]
   if pivot<=0:obstruction=dict(index=k,nonpositive_pivot=pivot);break
   for i in range(k+1,19):
    for j in range(k+1,19):test[i][j]-=test[i][k]*test[j][k]/pivot
  if obstruction:supports.append(dict(old_support=list(ids),old_Gram_not_PD=obstruction));continue
  BI=K.inv(B);g=[dot(r,olde[i])for i in ids];residual=ea(ee(scalar(2)),es(et(es(em(g[i],g[j]),BI[i][j])for i in range(19)for j in range(19)),-1));assert sg(residual)<0
  supports.append(dict(old_support=list(ids),negative_schur_residual=residual))
 escape_path='results/fresh_h20_e7_two_anchor_escape_20261014_scout.json';ed=json.loads((ROOT/escape_path).read_text())
 # Rational point is a proposal only; recheck every original and new sign.
 raw=ed.get('point',ed.get('escape_point'))
 if raw is None:raise AssertionError(list(ed))
 X0=[K.c(F(round(float(z)*10**9),10**9))for z in raw];X=[ee(ext(z))for z in X0];ep=[ee(ext(dot0(X0,n)))for n in old0]+[dot(X,r)]
 assert all(sg(z)<0 for z in ep)and K.sign(dot0(X0,X0))>0 and K.sign(dot0(X0,H0))<0
 originalE=[ee(ext(z))for z in E0];e42=K.read('results/fresh_h20_e7_replacement_escape_20261014.json');E42=[ee(ext(z))for z in K.vector(e42['point'])]
 assert sg(dot(originalE,r))>0 and sg(dot(E42,r))>0
 paths=[P.BASE,P.FACE,escape_path,'results/fresh_h20_e7_replacement_escape_20261014.json',f'scripts/{TAG}_exact.py',f'scripts/{TAG}_bridge.py','scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps_review.py']
 return K.enc(dict(status='PASS_88_WALL_MIXED_E7_TWO_ANCHOR_BRIDGE_STILL_INFINITE',sources={p:sha(p)for p in paths},field=dict(base_square_root_basis=RAD,extension='positive sqrt(discriminant)',discriminant=disc),form=Q,normals=N,norms=norms,recipe=dict(E7_positions=[4,6],E7_labels=[6,6],anchor_walls=[61,75],anchor_labels=[3,3],source_spatial_support=[0,12],a=aa,b=bb,t=t,quadratic_coefficients=[qa,qb,qc],branch=branch),new_pair_counts=dict(counts),new_pairs=pairs,strict_center=he,facet_points=fp,actual_pi3_ridges=ridges,preserved_face_vertices=vertices,preserved_bounded_cusps=links,new_facet_finite_vertex_exclusion=supports,old_escape_cut=dot(originalE,r),common42_escape_cut=dot(E42,r),strict_escape=dict(point=X,norm=ee(ext(dot0(X0,X0))),all88_pairings=ep),scope='Retains all87. New mixed-pi6 E7 cap makes actual pi3 ridges with both dotted-to-each-other old61/75; all87 new pairs exact.44 facevertices/2boundedcusps and completeH13face preserved. All88 facets actual by center-projection lemma. No finitevertices on newfacet after20support checks. Both earlier escapes cut, but newstrictescape remains.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_bridge.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'])
