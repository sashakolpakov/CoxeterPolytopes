#!/usr/bin/env python3
"""Exact H13 face of the new A11+A2+E7 corner; all-label extension gate."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import fresh_h20_e8_interior_caps_20261014_exact as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_face_20261014';SOURCE='results/fresh_h20_e8_interior_caps_20261014_exact.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 Q=[[F(x)for x in row]for row in d['form']];N=[[tuple(map(F,x))for x in row]for row in d['normals']];H=[tuple(map(F,x))for x in d['strict_center']];S=[12,13,14,15,19,33,67];A=list(range(11))+[18,83]
 cache={};dc={}
 def dot(x,y):
  pair=(tuple(x),tuple(y))
  if pair in dc:return dc[pair]
  if pair[1]not in cache:cache[pair[1]]=[tuple(sum(Q[i][j]*y[j][k]for j in range(21))for k in range(8))for i in range(21)]
  out=K.Z
  for a,b in zip(x,cache[pair[1]]):out=K.add(out,K.mul(a,b))
  dc[pair]=out;dc[(pair[1],pair[0])]=out;return out
 def linear(vectors,coefficients):
  out=[K.Z]*21
  for a,v in zip(coefficients,vectors):out=[K.add(x,K.mul(a,y))for x,y in zip(out,v)]
  return out
 B=[]
 for i in S:
  row=[]
  for j in S:
   g=dot(N[i],N[j]);assert g[1:]==(0,)*7;row.append(g[0])
  B.append(row)
 K.rat.pd(B);I=K.rat.inv(B);diag=[I[i][i]for i in range(7)]
 assert diag==[F(7,2),12,6,2,F(15,2),F(3,2),4]and all(x>0 for row in I for x in row)
 assert min(diag)==F(3,2)and S[diag.index(F(3,2))]==33
 assert sorted(diag)[:2]==[F(3,2),F(2)]
 def project(v):
  p=[dot(v,N[i])for i in S];coef=[tuple(sum(I[a][j]*p[j][k]for j in range(7))for k in range(8))for a in range(7)]
  shift=linear([N[i]for i in S],coef);return[K.add(x,K.scale(y,-1))for x,y in zip(v,shift)]
 P=project(H);assert K.sign(dot(P,P))<0 and K.sign(dot(P,H))<0
 projections=list(map(project,N));positive=[];negative=[];zero=[]
 for i,p in enumerate(projections):
  q=dot(p,p);sg=K.sign(q)
  if sg>0:positive.append(i);assert p==N[i]
  elif sg<0:negative.append(dict(wall=i,projection_norm=q,orientation=dot(p,P)));assert K.sign(dot(p,P))<0
  else:assert p==[K.Z]*21;zero.append(i)
 assert positive==A and zero==S and len(negative)==67
 G=[]
 for i in A:
  row=[]
  for j in A:
   g=dot(N[i],N[j]);assert g[1:]==(0,)*7;row.append(g[0])
  G.append(row)
 assert G==[[F(2)if i==j else F(-1)if ((i<11 and j<11)or(i>=11 and j>=11))and abs(i-j)==1 else F(0)for j in range(13)]for i in range(13)]
 K.rat.pd(G);GI=K.rat.inv(G)
 vertex=next(v for v in d['new_finite_vertices']if set(v['incident'])==set(S+A));V=[tuple(map(F,x))for x in vertex['point']]
 Y=linear([N[i]for i in A],[K.c(-sum(row))for row in GI])
 assert all(dot(V,N[i])==K.Z for i in S+A)and K.sign(dot(V,V))<0 and K.sign(dot(V,H))<0
 assert dot(Y,Y)==K.c(145)and dot(Y,V)==K.Z and all(dot(Y,N[i])==K.c(-1)for i in A)
 assert all(dot(Y,N[i])==K.Z for i in S)
 # Exact denesting of the timelike apex scale. The resulting intrinsic
 # basis has rational form A11+A2+<-2>, despite its nonrational embedding.
 s=tuple(map(F,[F(-3204,24),0,F(-153,24),0,F(5537,24),0,F(356,24),0]))
 assert K.sign(s)>0 and dot(V,V)==K.scale(K.mul(s,s),-2)
 conjugates=[tuple(-a if (i&mask).bit_count()%2 else a for i,a in enumerate(s))for mask in [2,4,6]]
 numerator=K.c(1)
 for z in conjugates:numerator=K.mul(numerator,z)
 ns=K.mul(s,numerator);assert ns[1:]==(0,)*7 and ns[0]!=0
 invs=K.scale(numerator,1/ns[0]);assert K.mul(s,invs)==K.c(1)
 T=[K.mul(invs,x)for x in V]
 assert dot(T,T)==K.c(-2)and K.sign(dot(T,H))<0
 assert all(dot(T,N[i])==K.Z for i in S+A)
 assert all(x[j]==0 for v in [T]+[N[i]for i in A]for x in v for j in [1,3,5,7])
 intrinsic_form=[row+[F(0)]for row in G]+[[F(0)]*13+[F(-2)]]
 determinant=F(-2)
 for pivot in K.rat.pd(G):determinant*=pivot
 assert determinant==-72
 reduced=[]
 for w,p in zip(N,projections):
  pair=[dot(N[i],w)for i in A]
  coefficients=[tuple(sum(GI[i][j]*pair[j][k]for j in range(13))for k in range(8))for i in range(13)]
  coefficients.append(K.scale(dot(T,w),F(-1,2)))
  assert linear([N[i]for i in A]+[T],coefficients)==p
  reduced.append(coefficients)
 # The unique nonright meeting profile has c_i²=q/2 at inverse diagonal3/2.
 residual_fraction=F(1)-F(1,2)*min(diag);assert residual_fraction==F(1,4)
 omega=linear([N[i]for i in S],[K.c(x)for x in I[5]])
 assert dot(omega,omega)==K.c(F(3,2))and all(dot(omega,N[i])==K.c(int(j==5))for j,i in enumerate(S))
 paths=[SOURCE,f'scripts/{TAG}.py','scripts/fresh_h20_e8_interior_caps_20261014_exact.py','scripts/fresh_h20_coupled_lift_20261014.py']
 return K.enc(dict(status='PASS_H13_E7_FACE_A11_A2_CONE_WITH_UNIQUE_NONRIGHT_EXTENSION_PROFILE',sources={p:sha(p)for p in paths},E7_support=S,gram=B,inverse=I,inverse_diagonal=diag,face_walls=A,face_gram=G,projected_walls=projections,causal_nonfacets=negative,face_strict_center=P,finite_apex=V,apex_norm=dot(V,V),strict_outgoing_spatial_direction=Y,direction_norm=dot(Y,Y),intrinsic_normalization=dict(apex_scale=s,timelike_basis=T,form=intrinsic_form,determinant=determinant,spatial_basis_walls=A,embedding_field='Q(sqrt79,sqrt11)',ambient_coefficient_basis=K.RAD,projected_wall_coordinates=reduced),unique_nonright_meeting_profile=dict(wall=33,label=3,other_six_right=True,projected_norm_fraction=residual_fraction,fundamental_weight=omega,weight_norm=F(3,2)),necessary_lift_recipes='For intrinsic face normal y of square qf>0, orthogonal lift is r=y. The only nonright lift is r=y-sqrt(2qf)*omega33, square4qf, pi3 atwall33 and right toother6; all remaining ambient pairs must still be gated.',scope='Exact current E7 face is an infinite-volume H13 A11+A2 cone. Any added actual facefacet in an acute Coxeter completion is orthogonal toE7 or has exactlyonepi3 atwall33. No arbitrarycompletion exclusion; no arithmeticity claim and no finitevolume claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
