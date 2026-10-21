#!/usr/bin/env python3
"""Complete15wall H13 source: the E8 face at Borcherds corner243.

Standalone Fraction arithmetic. All105 projective active sets are checked.
The source has an A11+A2 corner; transfer is allowed to adjoin sqrt2.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h13_a11a2_completion_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def rr(A):
 A=[list(map(F,r))for r in A];p=[];k=0
 for j in range(len(A[0])):
  z=next((i for i in range(k,len(A))if A[i][j]),None)
  if z is None:continue
  A[k],A[z]=A[z],A[k];a=A[k][j];A[k]=[x/a for x in A[k]]
  for i in range(len(A)):
   if i!=k:
    a=A[i][j];A[i]=[x-a*y for x,y in zip(A[i],A[k])]
  p.append(j);k+=1
  if k==len(A):break
 return A,p
def inv(A):
 n=len(A);R,p=rr([list(r)+[int(i==j)for j in range(n)]for i,r in enumerate(A)]);assert p[:n]==list(range(n));return[r[n:]for r in R]
def pd(A):
 A=[list(map(F,r))for r in A];ps=[]
 for i in range(len(A)):
  x=A[i][i];assert x>0;ps.append(x)
  for j in range(i+1,len(A)):
   for k in range(i+1,len(A)):A[j][k]-=A[j][i]*A[i][k]/x
 return ps
def mv(A,x):return[sum(a*b for a,b in zip(r,x))for r in A]
def lin(V,a):return[sum(c*v[i]for c,v in zip(a,V))for i in range(len(V[0]))]
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return{str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return[enc(x)for x in v]
 return v
def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-F(x)for x in row]for row in raw['basis_gram']];W=[list(map(F,r))for r in raw['wall_coordinates_as_columns_transposed']]
 walk=json.loads((ROOT/WALK).read_text());row=walk['corner_orbits'][243];assert row['id']==243 and row['type']=='finite';V=list(map(F,row['primitive_ray']))
 def dot(a,b):return sum(x*y for x,y in zip(a,mv(Q,b)))
 S=[8,11,22,25,26,28,33,39];A=[6,34,27,23,17,24,46,13,14,48,4];C=[36,162];corner=A+C
 B=[[dot(W[i],W[j])for j in S]for i in S];BI=inv(B);pd(B)
 assert sorted(sum(x!=0 for j,x in enumerate(r)if j!=i)for i,r in enumerate(B))==[1,1,1,2,2,2,2,3]
 # The tree is E8: arm lengths1,2,4.
 adj={i:[j for j in range(8)if i!=j and B[i][j]]for i in range(8)};center=next(i for i in adj if len(adj[i])==3);arms=[]
 for j in adj[center]:
  prev=center;length=0
  while True:
   length+=1;nxt=[k for k in adj[j]if k!=prev]
   if not nxt:break
   assert len(nxt)==1;prev,j=j,nxt[0]
  arms.append(length)
 assert sorted(arms)==[1,2,4]and all(B[i][i]==2 for i in range(8))and all(B[i][j]in[0,-1]for i in range(8)for j in range(i))
 G0=[[dot(W[i],W[j])for j in corner]for i in corner]
 target=[[F(2)if i==j else F(-1)if abs(i-j)==1 else F(0)for j in range(11)]for i in range(11)]
 assert [r[:11]for r in G0[:11]]==target
 assert [r[11:]for r in G0[11:]]==[[4,-2],[-2,4]]and all(G0[i][j]==0 for i in range(11)for j in range(11,13))
 GI=inv(G0);pd(G0);assert dot(V,V)==-4 and all(dot(V,W[i])==0 for i in S+corner)
 assert set(row['active_walls'])==set(S+corner)
 assert all(dot(W[i],W[j])==0 for i in S for j in corner)
 assert len(rr([W[i]for i in S+corner]+[V])[1])==22
 def project(w):
  coeff=mv(BI,[dot(w,W[i])for i in S]);shift=lin([W[i]for i in S],coeff);return[x-y for x,y in zip(w,shift)]
 projections=list(map(project,W));positive=[];null=[];negative=[];zero=[]
 for i,p in enumerate(projections):
  q=dot(p,p)
  if q>0:positive.append(i);assert p==W[i]
  elif p==[0]*22:zero.append(i)
  else:
   assert dot(p,V)<0
   (null if q==0 else negative).append(dict(wall=i,norm=q,orientation=dot(p,V)))
 assert positive==[3,4,6,9,13,14,17,23,24,27,34,36,46,48,162]and zero==S and len(null)==2 and len(negative)==185
 # Order13 corner facets, then sourcewall3 andsourcewall9.
 ids=corner+[3,9];N=[W[i]for i in ids];G=[[dot(r,s)for s in N]for r in N];counts=Counter()
 for i,j in combinations(range(15),2):
  a=G[i][j];assert a<=0;q=a*a/(G[i][i]*G[j][j])
  assert q in[0,F(1,4),F(1,2),F(1)]or q>1
  kind={F(0):'right',F(1,4):'pi3',F(1,2):'pi4',F(1):'parallel'}.get(q,'dotted');counts[kind]+=1
 weights=[lin(N[:13],r)for r in GI];Y=lin(weights,[-1]*13);H=[x+y/F(100)for x,y in zip(V,Y)]
 assert dot(H,H)<0 and dot(H,V)<0 and all(dot(H,r)<0 for r in N)
 facetpoints=[]
 for i,r in enumerate(N):
  p=[x-dot(H,r)/G[i][i]*y for x,y in zip(H,r)]
  assert dot(p,p)<0 and dot(p,V)<0 and all(dot(p,s)==0 if i==j else dot(p,s)<0 for j,s in enumerate(N));facetpoints.append(p)
 # With projective time coefficient V=1, gaps z_i=-(p,N_i)>=0
 # give p=V-sum z_i*weight_i. Two positive linear bounds close the slice.
 bounds=[[-dot(w,r)for w in weights]for r in N[13:]]
 assert [dot(V,r)for r in N[13:]]==[-1,-1]
 assert all(x>=0 for r in bounds for x in r)and all(sum(r[i]for r in bounds)>0 for i in range(13))
 inequalities=[[-F(int(i==j))for j in range(13)]for i in range(13)]+bounds;rhs=[F(0)]*13+[F(1)]*2
 vertices={};tests=0
 for J in combinations(range(15),13):
  tests+=1;R,piv=rr([inequalities[i]+[rhs[i]]for i in J])
  if piv!=list(range(13)):continue
  z=[r[-1]for r in R]
  if any(sum(x*y for x,y in zip(a,z))>b for a,b in zip(inequalities,rhs)):continue
  p=[x-y for x,y in zip(V,lin(weights,z))];q=dot(p,p);assert q<=0 and dot(p,V)==-4
  active=[i for i,r in enumerate(N)if dot(p,r)==0]
  assert all(dot(p,r)<=0 for r in W)
  if q<0:assert len(active)==13;pd([[G[i][j]for j in active]for i in active])
  vertices[tuple(z)]=dict(gaps=z,point=p,norm=q,incident=active)
 assert tests==105
 cusps=[]
 for data in vertices.values():
  if data['norm']:continue
  point=data['point'];J=data['incident'];un=set(J);blocks=[];rank=0
  while un:
   todo=[min(un)];C0=set()
   while todo:
    i=todo.pop()
    if i in C0:continue
    C0.add(i);todo.extend(j for j in J if j!=i and G[i][j])
   un-=C0;K=sorted(C0);B0=[[G[i][j]for j in K]for i in K];RR,pp=rr(B0);assert len(pp)==len(K)-1
   free=next(i for i in range(len(K))if i not in pp);c=[F(0)]*len(K);c[free]=1
   for r,i in zip(RR,pp):c[i]=-r[free]
   assert min(c)>0;pd([[B0[i][j]for j in range(len(K))if j!=free]for i in range(len(K))if i!=free])
   v=lin([N[i]for i in K],c);k=next(i for i,x in enumerate(point)if x);mu=v[k]/point[k]
   assert mu>0 and v==[mu*x for x in point];rank+=len(K)-1;blocks.append(dict(incident=K,positive_kernel=c,null_multiple=mu))
  assert rank==12;cusps.append(dict(point=point,incident=J,blocks=blocks,quotient_rank=rank))
 # Transfer to normalizedA11+A2 and T²=-2. Coefficients are a+b*sqrt2.
 normalcoords=[]
 for r in N:
  ca=mv(GI,[dot(r,s)for s in N[:13]]);beta=-dot(r,V)/4
  normalcoords.append([[x,0]for x in ca[:11]]+[[0,x]for x in ca[11:]]+[[0,beta]])
 # Independently evaluate every saved reduced pairing in Q(sqrt2).
 transfer_form=[[F(0)]*14 for _ in range(14)]
 for i in range(11):
  for j in range(11):transfer_form[i][j]=target[i][j]
 transfer_form[11][11]=transfer_form[12][12]=2
 transfer_form[11][12]=transfer_form[12][11]=-1;transfer_form[13][13]=-2
 for i,a in enumerate(normalcoords):
  for j,b in enumerate(normalcoords):
   value=[F(0),F(0)]
   for k in range(14):
    for l in range(14):
     t=transfer_form[k][l]
     if t:
      value[0]+=t*(a[k][0]*b[l][0]+2*a[k][1]*b[l][1])
      value[1]+=t*(a[k][0]*b[l][1]+a[k][1]*b[l][0])
   assert value==[G[i][j],0]
 # Formula weight coefficients use inverseA11/A2, not root coordinates.
 assert [dot(N[13],w)for w in N[:11]]==[0]*10+[-1]
 assert [dot(N[13],w)for w in N[11:13]]==[0,-2]
 assert [dot(N[14],w)for w in N[:11]]==[0,0,-1]+[0]*8
 assert all(dot(N[14],w)==0 for w in N[11:13])and G[13][14]==0
 paths=[RAW,WALK,f'scripts/{TAG}.py']
 return enc(dict(status='PASS_COMPLETE_15_WALL_H13_SOURCE_WITH_A11_A2_CORNER',sources={p:sha(p)for p in paths},
  ambient_form=Q,E8_support=S,E8_gram=B,E8_inverse=BI,source_corner_orbit=243,source_corner_ray=V,
  source_wall_ids=ids,normals=N,gram=G,pair_counts=dict(counts),projections=projections,zero_projection_walls=zero,causal_null_nonfacets=null,causal_timelike_nonfacets=negative,
  strict_center=H,strict_facet_points=facetpoints,corner_gram=G0,corner_inverse=GI,projective_positive_bounds=bounds,projective_active_sets=tests,
  projective_vertices=list(vertices.values()),finite_vertex_count=sum(x['norm']<0 for x in vertices.values()),ideal_vertex_count=len(cusps),bounded_cusp_links=cusps,
  normalized_transfer=dict(field='Q(sqrt2)',coefficient_encoding='a+b*sqrt2',A11_source_order=A,A2_source_order=C,A2_scale='1/sqrt2',time='T=V/sqrt2; T²=-2',
   normal_coordinates=normalcoords,form=transfer_form,form_name='diag(C_A11,C_A2,-2)',normal_coordinate_scope='Original15 source normals in the normalized frame; the two originalA2 walls still have norm4 and coordinates sqrt2 times the corresponding basis root.',extra_formulas=['R=-lambda11-sqrt2*mu2+(sqrt2/4)*T','S=-lambda3+(sqrt2/4)*T'],extra_norms=[2,2],mutual_pairing=0),
  scope='Complete finite-volume intrinsicH13 source with15actualfacets, independently proved by bounded projective105active-set census and positive full-rank cusp links. Extracted fromE8face of archivedH21 chamber; no ambientH20 compatibility or finitevolume assertion. This geometric transfer adjoins sqrt2 and is not restricted to the integral latticeA11+A2+<-2>.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['finite_vertex_count'],d['ideal_vertex_count'],d['pair_counts'])
