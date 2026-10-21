#!/usr/bin/env python3
"""Exact H12 E8-face reduction for the frozen76-wall replacement.
No new caps or search. The universal Schur bound covers every Coxeter label.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import fresh_h20_contact_closure_20261014_replacement as b
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_contact_closure_20261014_e8_face';SOURCE='results/fresh_h20_contact_closure_20261014_replacement.json';Z=b.ZZ

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h
 Q=[[F(x)for x in row]for row in d['form']];N=[[tuple(map(F,x))for x in row]for row in d['normals']];H=[b.c(F(x))for x in d['strict_center']]
 cache={}
 def dot(x,y):
  key=tuple(y)
  if key not in cache:cache[key]=[tuple(sum(Q[i][j]*y[j][k]for j in range(21))for k in range(4))for i in range(21)]
  g=Z
  for a,c in zip(x,cache[key]):g=b.add(g,b.mul(a,c))
  return g
 def linear(vectors,coeff):
  out=[Z]*21
  for a,v in zip(coeff,vectors):out=[b.add(x,b.mul(a,y))for x,y in zip(out,v)]
  return out
 E=[13,14,15,16,19,27,34,68];A=[64,1,2,3,4,5,6,7,8,9,10,11]
 B=[]
 for i in E:
  row=[]
  for j in E:
   g=dot(N[i],N[j]);assert g[1:]==(0,0,0);row.append(g[0])
  B.append(row)
 b.rat.pd(B);I=b.rat.inv(B);diag=[I[i][i]for i in range(8)];assert diag==list(map(F,[2,6,12,20,8,30,14,4]));assert all(x>0 for row in I for x in row)
 def project(v):
  pairs=[dot(N[i],v)for i in E];coefs=[tuple(sum(I[a][j]*pairs[j][k]for j in range(8))for k in range(4))for a in range(8)];shift=linear([N[i]for i in E],coefs)
  return[b.add(x,b.neg(y))for x,y in zip(v,shift)]
 P=project(H);assert b.sign(dot(P,P))<0 and b.sign(dot(P,H))<0
 assert all(dot(P,r)==Z if i in E else b.sign(dot(P,r))<0 for i,r in enumerate(N))
 orth=[i for i,r in enumerate(N)if all(dot(r,N[j])==Z for j in E)];assert orth==sorted(A)
 projections=[project(r)for r in N];remaining=[i for i in range(76)if i not in E+A]
 rows=[]
 for i in remaining:
  q=dot(projections[i],projections[i]);orientation=dot(projections[i],P);assert b.sign(q)<0 and b.sign(orientation)<0;rows.append(dict(wall=i,projection_norm=q,orientation=orientation))
 assert len(rows)==56
 # Literal12-node A12 chain, unchanged by projection.
 C=[]
 for i in A:
  row=[]
  for j in A:
   g=dot(N[i],N[j]);assert g[1:]==(0,0,0);row.append(g[0])
  C.append(row)
 assert C==[[F(2)if i==j else F(-1)if abs(i-j)==1 else F(0)for j in range(12)]for i in range(12)]
 b.rat.pd(C);CI=b.rat.inv(C);weights=[-sum(row)for row in CI];Y=linear([N[i]for i in A],[b.c(t)for t in weights]);V=[[tuple(map(F,x))for x in row['point']]for row in d['finite_vertices']if set(row['incident'])==set(E+A)];assert len(V)==1;V=V[0]
 assert all(dot(Y,N[i])==b.c(-1)for i in A);assert all(dot(Y,N[i])==Z for i in E)
 assert all(dot(V,N[i])==Z for i in E+A);assert dot(V,Y)==Z and b.sign(dot(V,V))<0 and b.sign(dot(V,H))<0 and b.sign(dot(Y,Y))>0
 # Normalize the timelike line: the resulting face form is rational.
 aa=F(4628);bb=F(2303);den=aa*aa-79*bb*bb
 k=(aa,F(0),bb,F(0));assert dot(V,V)==b.scale(b.mul(k,k),F(-1,117))
 factor=(3*aa/den,F(0),-3*bb/den,F(0));assert b.sign(factor)>0
 T=[b.mul(factor,x)for x in V];assert dot(T,T)==b.c(F(-1,13))
 assert all(dot(T,N[i])==Z for i in E+A)and b.sign(dot(T,H))<0
 reduced=[]
 for w,proj in zip(N,projections):
  pairs=[dot(N[i],w)for i in A]
  coefs=[tuple(sum(CI[a][j]*pairs[j][k]for j in range(12))for k in range(4))for a in range(12)]
  coefs.append(b.scale(dot(T,w),-13))
  assert linear([N[i]for i in A]+[T],coefs)==proj
  reduced.append(coefs)
 form=[row+[F(0)]for row in C]+[[F(0)]*12+[F(-1,13)]]
 # V+sqrt(-V²/Y²)Y is future null and strict on every non-E8 wall.
 # Its tiny outward continuation is an intrinsic strict spacelike escape.
 paths=[SOURCE,f'scripts/{TAG}.py','scripts/fresh_h20_contact_closure_20261014_replacement.py','scripts/fresh_h20_coupled_lift_20261014.py']
 return b.enc(dict(status='PASS_E8_FACE_IS_A12_H12_CONE_NEEDS_NEW_E8_ORTHOGONAL_WALLS',sources={p:sha(p)for p in paths},E8_walls=E,gram=B,inverse=I,min_inverse_diagonal=min(diag),positive_inverse_entries=True,orthogonal_face_walls=A,A12_gram=C,nonorthogonal_projected_walls=rows,face_strict_interior=P,finite_apex=V,apex_norm=dot(V,V),normalized_timelike_basis=T,timelike_scaling_from_apex=factor,induced_form=form,induced_form_determinant=F(-1),ambient_embedding_field='Q(sqrt79)',all76_reduced_projection_coordinates=reduced,strict_outgoing_spatial_direction=Y,direction_norm=dot(Y,Y),universal_schur_bound='For any acute Coxeter cap r of normq>0 not orthogonal to E8, a nonzero coefficient c_i=−r.e_i has c_i²>=q/2. Positive E8inverse and min diagonal2 give projection²>=q and residualnorm<=0. It induces no facet of the H12 face.',scope='Exact current H12 face has exactly12A12 facets and is infinitevolume. A finitevolume completion retaining this E8 support must add walls orthogonal to all8 E8 roots. Such additions remain possible; no unrestricted noncompletion theorem and no compactness claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['nonorthogonal_projected_walls']))
