#!/usr/bin/env python3
"""Exact mixed-E7 interior cap after deleting two off-face caps."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from hashlib import sha256
import fresh_h20_e7_intrinsic_caps_20261014_cusps_review as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_causal_residuals_20261014'
BASE='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json';FACE='results/fresh_h20_e7_intrinsic_caps_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def setup():
 d,f=map(K.read,[BASE,FACE]);Q=[list(map(F,r))for r in d['form']];old=list(map(K.vector,d['normals']));H=K.vector(d['strict_center']);E=K.vector(d['strict_escape']['point']);cache={};sp=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
 def dot(v,w):
  key=tuple(w)
  if key not in cache:cache[key]=[tuple(sum(q*w[j][k]for j,q in row)for k in range(8))for row in sp]
  return K.total(K.mul(a,b)for a,b in zip(v,cache[key]))
 S=d['E7_support'];A=f['face_facet_walls'][:13]
 G=[[K.rational(dot(old[i],old[j]))for j in S]for i in S];J=K.inv(G)
 C=[[K.rational(dot(old[i],old[j]))for j in A]for i in A];I=K.inv(C)
 omega=[K.lin([old[i]for i in S],list(map(K.c,row)))for row in J]
 lam=[K.lin([old[i]for i in A],list(map(K.c,row)))for row in I]
 sqrt2=(F(0),F(1))+(F(0),)*6
 T=[K.mul(x,K.scale(sqrt2,F(1,2)))for x in K.vector(f['face_basis'][-1])]
 assert dot(T,T)==K.c(-2)and all(dot(T,old[i])==K.Z for i in S+A)
 w=K.lin([omega[i]for i in [2,4,6]],[K.c(1)]*3)
 y=K.lin([lam[0],lam[12],T],[K.scale(sqrt2,-1),K.c(-2),K.c(5)])
 r=K.lin([w,y],[K.c(-1),K.c(1)])
 assert dot(w,w)==K.c(F(95,2))and dot(y,y)==K.c(F(-91,2))and dot(r,r)==K.c(2)
 assert K.sign(dot(y,T))<0
 deleted=[61,75];retained=[i for i in range(87)if i not in deleted];N=[old[i]for i in retained]+[r]
 return d,f,Q,old,H,E,S,A,J,I,omega,lam,T,w,y,r,deleted,retained,N,dot
def run():
 d,f,Q,old,H,E,S,A,J,I,omega,lam,T,w,y,r,deleted,retained,N,dot=setup();norms=[K.rational(dot(n,n))for n in N]
 bad=[]
 for i in deleted:
  g=dot(r,old[i]);assert K.sign(g)<0 and K.sign(K.sub(K.c(2*F(d['norms'][i])/4),K.mul(g,g)))>0
  bad.append(dict(old_wall=i,pairing=g,squared_cosine=K.scale(K.mul(g,g),1/(2*F(d['norms'][i]))),reason='strictly_between_right_and_pi3'))
 counts=Counter();pairs=[]
 for j in range(85):
  g=dot(r,N[j]);prod=2*norms[j];assert K.sign(g)<=0
  if g==K.Z:kind='right'
  elif K.mul(g,g)==K.c(prod/4):kind='pi3'
  elif K.mul(g,g)==K.c(prod/2):kind='pi4'
  elif K.mul(g,g)==K.c(3*prod/4):kind='pi6'
  elif K.mul(g,g)==K.c(prod):kind='parallel'
  else:assert K.sign(K.sub(K.mul(g,g),K.c(prod)))>0;kind='dotted'
  counts[kind]+=1;pairs.append(dict(walls=[85,j],old_wall=retained[j],pairing=g,relation=kind))
 assert K.sign(dot(H,H))<0 and all(K.sign(dot(H,n))<0 for n in N)
 facetpoints=[]
 for i,n in enumerate(N):
  p=K.lin([H,n],[K.c(1),K.scale(dot(H,n),-1/norms[i])])
  assert K.sign(dot(p,p))<0 and K.sign(dot(p,H))<0
  assert all(dot(p,v)==K.Z if j==i else K.sign(dot(p,v))<0 for j,v in enumerate(N));facetpoints.append(p)
 assert K.sign(dot(r,E))>0
 # Each deleted normal had strict timelike face projection: check this
 # directly. The new projection y is future timelike too.
 projections=[]
 for i in deleted:
  pair=[dot(old[i],old[j])for j in S]
  coeff=[K.total(K.scale(pair[b],J[a][b])for b in range(7))for a in range(7)]
  p=K.lin([old[i]]+[old[j]for j in S],[K.c(1)]+[K.scale(x,-1)for x in coeff])
  assert K.sign(dot(p,p))<0 and K.sign(dot(p,T))<0;projections.append(dict(old_wall=i,projection=p,norm=dot(p,p)))
 assert not(set(deleted)&set(f['face_facet_walls']+S))
 vertices=[]
 for row in f['face_finite_vertices']:
  p=K.vector(row['ambient_point']);assert K.sign(dot(p,r))<0
  signs=[K.sign(dot(p,n))for n in N];assert max(signs)==0
  assert [retained[i]for i in range(85)if signs[i]==0]==row['ambient_incident']
  vertices.append(dict(old_incident=row['ambient_incident'],new_pairing=dot(p,r)))
 links=[]
 for row in d['completed_cusp_links']:
  nu=K.vector(row['point']);assert K.sign(dot(nu,r))<0
  assert all(K.sign(dot(nu,old[i]))<0 for i in deleted)
  inc=[retained[i]for i in range(85)if dot(nu,N[i])==K.Z];assert inc==row['incident']
  assert all(K.sign(dot(nu,n))<=0 for n in N)
  links.append(dict(point=nu,old_incident=inc,new_pairing=dot(nu,r),inherited_quotient_rank=19))
 # Every finite neighbor is one of the twenty old corner walls. Exhaust
 # the possible 20-wall vertex supports containing this cap by omitting
 # one old corner wall, with exact Schur complements.
 finite=[retained[j]for j,row in enumerate(pairs)if row['relation']in ['right','pi3','pi4','pi6']]
 assert set(finite).issubset(set(S+A))
 from itertools import combinations
 supports=[]
 for ids in combinations(finite,19):
  B=[[K.rational(dot(old[i],old[j]))for j in ids]for i in ids];K.pd(B);BI=K.inv(B)
  b=[dot(r,old[i])for i in ids];cost=K.total(K.scale(K.mul(b[i],b[j]),BI[i][j])for i in range(19)for j in range(19))
  residual=K.sub(K.c(2),cost);assert K.sign(residual)<0
  supports.append(dict(old_support=list(ids),schur_residual=residual))
 assert supports
 paths=[BASE,FACE,f'scripts/{TAG}_exact.py','scripts/fresh_h20_e7_intrinsic_caps_20261014_cusps_review.py']
 result=dict(status='PASS_86_WALL_CHANGED_CORE_MIXED_E7_CAP_OLD_ESCAPE_CUT',sources={x:sha(x)for x in paths},form=Q,normals=N,norms=norms,deleted_old87_walls=deleted,retained_old87_walls=retained,E7_support_old87=S,source_corner_old87=A,recipe=dict(E7_positions=[2,4,6],source_contacts={'0':'sqrt2','12':'2'},time_coefficient=5),normal=r,face_projection=y,face_projection_norm=dot(y,y),deleted_wall_failures=bad,new_pair_counts=dict(counts),new_pairs=pairs,strict_center=H,facet_points=facetpoints,old_escape_cut=dot(r,E),deleted_face_projections=projections,preserved_face_vertices=vertices,preserved_bounded_cusps=links,no_finite_vertex_on_new_cap=supports,scope='Deletes two specified off-face caps and adds one norm2 mixed-E7 cap. Exact Coxeter pairs, actual facets, complete unchanged H13 face,44 finite facevertices and2 bounded cusps. Cuts prior global escape; does not establish finite volume. New cap has no finite vertex, by exhaustive finite-neighbor support check.')
 escape_path='results/fresh_h20_e7_replacement_escape_20261014.json'
 ed=K.read(escape_path);X=K.vector(ed['point']);values=[dot(X,n)for n in N]
 assert K.sign(dot(X,X))>0 and K.sign(dot(X,H))<0 and all(K.sign(x)<0 for x in values)
 result['strict_escape']=dict(point=X,norm=dot(X,X),all86_pairings=values);result['sources'][escape_path]=sha(escape_path)
 result['status']='PASS_86_WALL_CHANGED_CORE_OLD_ESCAPE_CUT_NEW_ESCAPE_REMAINS'
 return K.enc(result)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_exact.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'])
