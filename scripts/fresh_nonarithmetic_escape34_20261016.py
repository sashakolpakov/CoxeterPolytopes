#!/usr/bin/env python3
"""Exact escape cutter closing a tetrahedral three-face of the 34-wall chamber."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_nonarithmetic_wall3_cap_20261016 as P
R,C,K,Z,O,X=P.R,P.C,P.K,P.Z,P.O,P.X
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_escape34_20261016'
PREVIOUS='results/fresh_nonarithmetic_wall3_cap_20261016.json'
REVIEW='results/fresh_nonarithmetic_wall3_cap_20261016_review.json'
D=24*O+72*R.BASIS[1]/5
SUPPORT=[0,1,2,3,5,6,7,12]
FACE=[0,1,2,3,5,12]
THETA=K.from_sympy(K.ext.as_expr())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d['sources'].items():assert sha(f)==h,f
 return d
@lru_cache(None)
def decode(a):
 z=Z
 for c in a:z=z*THETA+K.convert(F(c))
 return z
def db(a):return decode(tuple(a))
def dx(a):return X(db(a[0]),db(a[1]))
class Y:
 """a+b sqrt(D), over the previous formal quadratic extension."""
 __slots__=('a','b')
 def __init__(self,a=0,b=0):
  if isinstance(a,Y):self.a,self.b=a.a,a.b
  else:self.a,self.b=X(a),X(b)
 def __add__(self,o):
  o=Y(o);return Y(self.a+o.a,self.b+o.b)
 __radd__=__add__
 def __neg__(self):return Y(-self.a,-self.b)
 def __sub__(self,o):return self+-Y(o)
 def __rsub__(self,o):return Y(o)+-self
 def __mul__(self,o):
  o=Y(o);return Y(self.a*o.a+D*self.b*o.b,self.a*o.b+self.b*o.a)
 __rmul__=__mul__
 def __truediv__(self,o):
  o=Y(o);assert o.b==X() and not o.a.b and o.a.a
  return Y(self.a/o.a,self.b/o.a)
 def __eq__(self,o):
  o=Y(o);return self.a==o.a and self.b==o.b
 def sign(self):
  a,b=self.a.sign(),self.b.sign()
  if not a:return b
  if not b or a==b:return a
  return a*(self.a*self.a-D*self.b*self.b).sign()
def enc(x):
 if isinstance(x,Y):return [P.enc(x.a),P.enc(x.b)]
 if isinstance(x,X):return P.enc(x)
 if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 if hasattr(x,'to_list'):return R.enc(x)
 return x
def dot(a,b):return sum((Y(x)*Y(y)for x,y in zip(a,b)),Y())
def mv(M,v):return [dot(row,v)for row in M]
def pair(a,b,Q):return dot(a,mv(Q,b))
def base(x):
 x=Y(x);assert x.b==X() and not x.a.b;return x.a.a

def run():
 old,rev=map(read,(PREVIOUS,REVIEW))
 assert rev['sources'][PREVIOUS]==sha(PREVIOUS)
 assert old['finite_vertex_count']==95 and old['finite_edge_count']==257
 Q=[[db(x)for x in row]for row in old['form']]
 H=list(map(db,old['center']));oldE=list(map(db,old['strict_escape']['point']))
 N=[[dx(x)for x in row]for row in old['normals']]
 assert len(N)==34
 # The new recipe uses only original base-field normals.
 ids=SUPPORT+[8]
 oldcov=[P.mv(Q,n)for n in N]
 rows=[[base(x)for x in oldcov[i]]for i in ids]
 _,pivots=C.eliminate(rows);assert len(pivots)==9
 k=C.null_one(rows)
 p=C.mv(C.inverse(rows+[C.mv(Q,k)]),[Z]*8+[-O,Z])
 pp,kk,pk=C.pair(p,p,Q),C.pair(k,k,Q),C.pair(p,k,Q)
 assert pp==14*O+6*R.BASIS[1] and kk==-(5*O+R.BASIS[1])/16 and pk==Z
 assert pp+D*kk==2*O and R.sg(D)>0
 branches=[]
 for sign in[-1,1]:
  r=[Y(a,sign*b)for a,b in zip(p,k)]
  vals=[dot(r,c)for c in oldcov]
  branches.append(dict(sign=sign,norm=pair(r,r,Q),old_pairing_signs=[v.sign()for v in vals],center_sign=pair(r,H,Q).sign(),escape_sign=pair(r,oldE,Q).sign()))
 assert branches[0]['center_sign']==-1 and branches[0]['escape_sign']==1
 assert max(branches[0]['old_pairing_signs'])==0 and 1 in branches[1]['old_pairing_signs']
 r=[Y(a,-b)for a,b in zip(p,k)]
 normals=[[Y(x)for x in n]for n in N]+[r]
 cov=[mv(Q,n)for n in normals];norms=[dot(n,c)for n,c in zip(normals,cov)]
 gram=[[dot(n,c)for c in cov]for n in normals]
 hp=[dot(H,c)for c in cov]
 assert all(x.sign()<0 for x in hp) and pair(H,H,Q).sign()<0
 relations=[];neighbors=[];counts=Counter()
 for i in range(34):
  g=gram[34][i];assert g.sign()<=0
  cycle=4*g*g/(norms[34]*norms[i])
  if g==Y():kind,label='right',2
  elif i==8:
   assert g==Y(-1) and norms[i]==Y(2) and cycle==Y(1)
   kind,label='finite',3
  else:
   assert (cycle-4).sign()>0
   kind,label='dotted',None
  if kind!='dotted':neighbors.append(i)
  counts[kind]+=1;relations.append(dict(wall=i,pairing=g,cartan_cycle=cycle,kind=kind,label=label))
 assert counts==dict(right=8,finite=1,dotted=25) and neighbors==sorted(ids)
 assert all(gram[i][j].sign()<=0 for i in range(35)for j in range(i))
 facets=[]
 for i,n in enumerate(normals):
  V=[Y(h)-(hp[i]/norms[i])*a for h,a in zip(H,n)]
  vals=[dot(V,c)for c in cov];norm=pair(V,V,Q)
  assert norm.sign()<0 and pair(V,H,Q).sign()<0
  assert [v.sign()for v in vals]==[0 if j==i else -1 for j in range(35)]
  facets.append(dict(wall=i,point=V,norm=norm,pairings=vals))
 print('35 actual facets;34 new exact Coxeter pairs;renewed escape cut PASS',flush=True)
 vertices=[];preservation=[]
 for v in old['finite_vertices']:
  V=[Y(dx(x))for x in v['point']];g=dot(V,cov[34]);assert g.sign()<0
  preservation.append(g);vertices.append(dict(incident=v['incident'],point=V,norm=pair(V,V,Q)))
 candidates=[];new=[]
 for sub in combinations(neighbors,8):
  support=list(sub)+[34];G=[[base(gram[i][j])for j in support]for i in support]
  pd,piv=R.pd(G);candidates.append(dict(support=support,positive_definite=pd,pivots=piv))
  if not pd:continue
  inv=C.inverse(G);coeff=[sum((Y(c)*hp[j]for c,j in zip(row,support)),Y())for row in inv]
  assert all(c.sign()<0 for c in coeff)
  V=[Y(H[t])-sum((c*normals[i][t]for c,i in zip(coeff,support)),Y())for t in range(10)]
  vals=[dot(V,c)for c in cov];norm=pair(V,V,Q)
  assert norm.sign()<0 and pair(V,H,Q)==norm
  assert [a.sign()for a in vals]==[0 if j in support else -1 for j in range(35)]
  new.append(dict(incident=support,point=V,norm=norm,pairings=vals,active_gram=G,inverse_active_gram=inv,projection_coefficients=coeff))
  vertices.append(dict(incident=support,point=V,norm=norm))
 assert len(new)==3 and len(vertices)==98 and len(candidates)==9
 edges=[]
 for i,j in combinations(range(98),2):
  common=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
  if len(common)==8:edges.append(dict(vertices=[i,j],support=common))
 new_edges=[e for e in edges if max(e['vertices'])>=95]
 assert [e for e in edges if max(e['vertices'])<95]==old['finite_edges']
 assert len(edges)==263 and len(new_edges)==6
 assert [e['vertices']for e in new_edges]==[[1,95],[1,96],[1,97],[95,96],[95,97],[96,97]]
 print('Complete3 new finite vertices;6 new edges;compact tetrahedron candidates PASS',flush=True)
 # Exact complete three-face: six white roots and four boundary walls.
 FG=[[base(gram[i][j])for j in FACE]for i in FACE];assert R.pd(FG)[0]
 FI=C.inverse(FG);boundary=[6,7,8,34]
 projected=[]
 for i in boundary:
  cs=C.mv(FI,[base(gram[i][j])for j in FACE])
  projected.append([normals[i][t]-sum((Y(c)*normals[j][t]for c,j in zip(cs,FACE)),Y())for t in range(10)])
 BG=[[base(pair(a,b,Q))for b in projected]for a in projected]
 BI=C.inverse(BG)
 assert all(R.pd([[BG[i][j]for j in sub]for i in sub])[0]for sub in combinations(range(4),3))
 assert all(R.sg(x)<0 for row in BI for x in row)
 dual_vertices=[]
 for j in range(4):
  Vj=[-sum((Y(BI[i][j])*projected[i][t]for i in range(4)),Y())for t in range(10)]
  vals=[dot(Vj,c)for c in cov]
  incident=sorted(FACE+[boundary[i]for i in range(4)if i!=j])
  assert pair(Vj,Vj,Q)==Y(BI[j][j]) and pair(Vj,H,Q).sign()<0
  assert [a.sign()for a in vals]==[0 if i in incident else -1 for i in range(35)]
  assert [vals[i]for i in boundary]==[Y(-1 if i==j else 0)for i in range(4)]
  dual_vertices.append(dict(incident=incident,point=Vj,norm=pair(Vj,Vj,Q),pairings=vals))
 face_ids=[1,95,96,97]
 assert set.intersection(*(set(vertices[i]['incident'])for i in face_ids))==set(FACE)
 for i in face_ids:
  assert len(set(vertices[i]['incident'])&set(boundary))==3
 fc=[sum((vertices[i]['point'][t]for i in face_ids),Y())for t in range(10)]
 fcp=[dot(fc,c)for c in cov]
 assert pair(fc,fc,Q).sign()<0 and pair(fc,H,Q).sign()<0
 assert [a.sign()for a in fcp]==[0 if j in FACE else -1 for j in range(35)]
 # A different old exterior ray survives; no complete exterior census is claimed.
 source=read(P.BASE);escape_index=8
 escape_ids=source['complete_census']['hyperideal_rays'][escape_index]
 assert escape_ids==[0,1,2,4,5,6,7,9,17]
 V=C.null_one([[base(x)for x in cov[i]]for i in escape_ids]);scale=-C.pair(V,H,Q);V=[a/scale for a in V]
 raypairs=[dot(V,c)for c in cov];assert R.sg(C.pair(V,V,Q))>0
 assert all(a.sign()<=0 for a in raypairs)
 E=[v+K.convert(F(1,8192))*h for v,h in zip(V,H)]
 ep=[dot(E,c)for c in cov];en=C.pair(E,E,Q)
 assert R.sg(en)>0 and R.sg(C.pair(E,H,Q))<0 and all(a.sign()<0 for a in ep)
 # Its projection lies strictly inside the newly compact three-face.
 ec=C.mv(FI,[base(dot(E,cov[i]))for i in FACE])
 EP=[Y(E[t])-sum((Y(c)*normals[i][t]for c,i in zip(ec,FACE)),Y())for t in range(10)]
 epp=[dot(EP,c)for c in cov]
 assert pair(EP,EP,Q).sign()<0 and pair(EP,H,Q).sign()<0
 assert [a.sign()for a in epp]==[0 if j in FACE else -1 for j in range(35)]
 barycentric=[-epp[i]for i in boundary]
 assert all(a.sign()>0 for a in barycentric)
 assert EP==[sum((a*v['point'][t]for a,v in zip(barycentric,dual_vertices)),Y())for t in range(10)]
 paths=[f'scripts/{TAG}.py',PREVIOUS,REVIEW,'scripts/fresh_nonarithmetic_wall3_cap_20261016.py','scripts/fresh_h9_joint_corner_20261015_bridge_review.py',P.BASE]
 return enc(dict(status='PASS_EXACT_ESCAPE_CUTTER_AND_COMPACT_TETRAHEDRAL_FACE',sources={f:sha(f)for f in paths},
  field=dict(base='Previous formal extension over the frozen octic engine',extra_radical_square=D,extra_radical_expression='24+72*sqrt(5)/5',encoding='[a,b] means a+b*positive_sqrt(extra_radical_square), with a,b encoded as in the previous wall3 package; no minimal-degree claim.'),
  form=Q,center=H,normals=normals,norms=norms,gram=gram,
  cap=dict(index=34,right_support=SUPPORT,anchor=8,label=3,linear_pivots=pivots,particular=p,kernel=k,particular_norm=pp,kernel_norm=kk,cross_pairing=pk,normal=r,norm=norms[34],branch=-1,branches=branches,center_pairing=hp[34],old_escape_pairing=dot(oldE,cov[34]),relations=relations),
  pair_counts=dict(counts),facet_witnesses=facets,actual_facets=35,
  old_vertex_preservation_pairings=preservation,complete_new_cap_finite_census=dict(finite_neighbors=neighbors,candidate_supports=candidates,vertices=new),
  finite_vertices=vertices,finite_vertex_incidence=[v['incident']for v in vertices],finite_vertex_count=98,finite_edges=edges,finite_edge_count=263,new_finite_edges=new_edges,open_finite_germ_count=356,finite_graph_components=P.components(98,edges),
  compact_three_face=dict(support=FACE,support_gram=FG,inverse_support_gram=FI,boundary=boundary,projected_normals=projected,projected_boundary_gram=BG,inverse_projected_boundary_gram=BI,vertices=face_ids,dual_vertex_rays=dual_vertices,strict_center=fc,center_pairings=fcp,proof='All four proper rank-three boundary Grams are PD and the full inverse is strictly negative. The four finite dual rays form a complete compact tetrahedron in the support perpendicular. Their actual old/new signs prove the whole retained face equals this tetrahedron. No intrinsic Coxeter-angle assertion is made.'),
  old_escape=dict(point=oldE,norm=C.pair(oldE,oldE,Q),new_cap_pairing=dot(oldE,cov[34])),
  strict_escape=dict(point=E,norm=en,pairings=ep,ray=V,ray_active=escape_ids,ray_index=escape_index,ray_pairings=raypairs,center_epsilon='1/8192',compact_face_projection=EP,projection_norm=pair(EP,EP,Q),projection_pairings=epp,positive_dual_vertex_coefficients=barycentric),
  scope='All35 facets are actual. Complete finite-corner update adds3 vertices and6 edges, closes3 old germs, and certifies a compact tetrahedral three-face. A new strict spacelike escape proves infinite volume. No complete exterior-ray census or finite-volume example is claimed. Any further coherent acute wall perpendicular to the six-root compact-face support cannot cut this escape, by finite-vertex preservation and its interior face projection.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text());print('PASS saved exact tetrahedral cap replay')
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['finite_vertex_count'],out['finite_edge_count'],out['open_finite_germ_count'])
