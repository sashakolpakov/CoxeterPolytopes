#!/usr/bin/env python3
"""Exact new escape cutter closing a compact Lambert quadrilateral in H9."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_nonarithmetic_escape34_20261016 as B
R,C,K,Z,O,X=B.R,B.C,B.K,B.Z,B.O,B.Y
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_escape35_20261016'
PREVIOUS='results/fresh_nonarithmetic_escape34_20261016.json'
REVIEW='results/fresh_nonarithmetic_escape34_20261016_review.json'
BASE32=B.P.BASE
D=10*O+68*R.BASIS[1]/15
SUPPORT=[0,1,2,4,5,6,7,9];FACE=[0,1,2,4,5,6,9]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d['sources'].items():assert sha(f)==h,f
 return d
def dx(a):return X(B.dx(a[0]),B.dx(a[1]))
class Y:
 __slots__=('a','b')
 def __init__(self,a=0,b=0):
  if isinstance(a,Y):self.a,self.b=a.a,a.b
  else:self.a,self.b=X(a),X(b)
 def __add__(self,o):o=Y(o);return Y(self.a+o.a,self.b+o.b)
 __radd__=__add__
 def __neg__(self):return Y(-self.a,-self.b)
 def __sub__(self,o):return self+-Y(o)
 def __rsub__(self,o):return Y(o)+-self
 def __mul__(self,o):
  o=Y(o);return Y(self.a*o.a+D*self.b*o.b,self.a*o.b+self.b*o.a)
 __rmul__=__mul__
 def __truediv__(self,o):
  o=Y(o);assert o.b==X();q=B.base(o.a);assert q
  return Y(self.a/q,self.b/q)
 def __eq__(self,o):o=Y(o);return self.a==o.a and self.b==o.b
 def sign(self):
  a,b=self.a.sign(),self.b.sign()
  if not a:return b
  if not b or a==b:return a
  return a*(self.a*self.a-D*self.b*self.b).sign()
def enc(x):
 if isinstance(x,Y):return [B.enc(x.a),B.enc(x.b)]
 if isinstance(x,X):return B.enc(x)
 if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 if hasattr(x,'to_list'):return R.enc(x)
 return x
def dot(a,b):return sum((Y(x)*Y(y)for x,y in zip(a,b)),Y())
def mv(M,v):return [dot(row,v)for row in M]
def pair(a,b,Q):return dot(a,mv(Q,b))
def base(x):x=Y(x);assert x.b==X();return B.base(x.a)
def cross(a,b):return[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def load_geometry():
 d=read(PREVIOUS)
 return d,[[B.db(x)for x in row]for row in d['form']],[[dx(x)for x in row]for row in d['normals']],list(map(B.db,d['center']))
def run():
 old,Q,N,H=load_geometry();review=read(REVIEW);assert review['sources'][PREVIOUS]==sha(PREVIOUS)
 assert len(N)==35 and old['finite_vertex_count']==98 and old['finite_edge_count']==263
 oldE=list(map(B.db,old['strict_escape']['point']));oldcov=[B.mv(Q,n)for n in N]
 ids=SUPPORT+[8];rows=[[B.base(x)for x in oldcov[i]]for i in ids]
 _,pivots=C.eliminate(rows);assert len(pivots)==9
 k=C.null_one(rows);assert k==[Z]*9+[O]
 p=C.mv(C.inverse(rows+[C.mv(Q,k)]),[Z]*8+[-O,Z])
 pp,kk,pk=C.pair(p,p,Q),C.pair(k,k,Q),C.pair(p,k,Q)
 assert pp==-5*O-5*R.BASIS[1]/2 and kk==R.BASIS[1]/4 and pk==Z
 assert pp+D*kk==2*O/3 and R.sg(D)>0
 branches=[]
 for sign in[-1,1]:
  r=[Y(a,sign*b)for a,b in zip(p,k)];vals=[dot(r,c)for c in oldcov]
  branches.append(dict(sign=sign,norm=pair(r,r,Q),old_pairing_signs=[a.sign()for a in vals],center_sign=pair(r,H,Q).sign(),escape_sign=pair(r,oldE,Q).sign()))
 assert branches[0]['center_sign']==-1 and branches[0]['escape_sign']==1 and max(branches[0]['old_pairing_signs'])==0
 assert 1 in branches[1]['old_pairing_signs']
 r=[Y(a,-b)for a,b in zip(p,k)];normals=[[Y(a)for a in n]for n in N]+[r]
 cov=[mv(Q,n)for n in normals];norms=[dot(n,c)for n,c in zip(normals,cov)]
 gram=[[dot(n,c)for c in cov]for n in normals];hp=[dot(H,c)for c in cov]
 assert all(x.sign()<0 for x in hp) and pair(H,H,Q).sign()<0
 assert all(gram[i][j].sign()<=0 for i in range(36)for j in range(i))
 relations=[];neighbors=[];counts=Counter()
 for i in range(35):
  g=gram[35][i];cycle=4*g*g/(norms[35]*norms[i]);assert g.sign()<=0
  if g==Y():kind,label='right',2
  elif i==8:
   assert g==Y(-1) and norms[i]==Y(2) and cycle==Y(3);kind,label='finite',6
  else:assert (cycle-4).sign()>0;kind,label='dotted',None
  if kind!='dotted':neighbors.append(i)
  counts[kind]+=1;relations.append(dict(wall=i,pairing=g,cartan_cycle=cycle,kind=kind,label=label))
 assert counts==dict(right=8,finite=1,dotted=26) and neighbors==sorted(ids)
 facets=[]
 for i,n in enumerate(normals):
  V=[Y(h)-(hp[i]/norms[i])*a for h,a in zip(H,n)];vals=[dot(V,c)for c in cov]
  assert pair(V,V,Q).sign()<0 and pair(V,H,Q).sign()<0
  assert [a.sign()for a in vals]==[0 if j==i else -1 for j in range(36)]
  facets.append(dict(wall=i,point=V,norm=pair(V,V,Q),pairings=vals))
 print('36 actual facets;35 new exact Coxeter pairs;current escape cut PASS',flush=True)
 vertices=[];preservation=[]
 for v in old['finite_vertices']:
  V=[Y(dx(x))for x in v['point']];g=dot(V,cov[35]);assert g.sign()<0
  preservation.append(g);vertices.append(dict(incident=v['incident'],point=V,norm=pair(V,V,Q)))
 candidates=[];new=[]
 for sub in combinations(neighbors,8):
  support=list(sub)+[35];G=[[base(gram[i][j])for j in support]for i in support]
  positive,piv=R.pd(G);candidates.append(dict(support=support,positive_definite=positive,pivots=piv))
  if not positive:continue
  inv=C.inverse(G);coeff=[sum((Y(c)*hp[j]for c,j in zip(row,support)),Y())for row in inv]
  assert all(a.sign()<0 for a in coeff)
  V=[Y(H[t])-sum((a*normals[i][t]for a,i in zip(coeff,support)),Y())for t in range(10)]
  vals=[dot(V,c)for c in cov];norm=pair(V,V,Q)
  assert norm.sign()<0 and pair(V,H,Q)==norm
  assert [a.sign()for a in vals]==[0 if j in support else -1 for j in range(36)]
  new.append(dict(incident=support,point=V,norm=norm,pairings=vals,active_gram=G,inverse_active_gram=inv,projection_coefficients=coeff));vertices.append(dict(incident=support,point=V,norm=norm))
 assert len(candidates)==9 and len(new)==2 and len(vertices)==100
 edges=[]
 for i,j in combinations(range(100),2):
  S=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
  if len(S)==8:edges.append(dict(vertices=[i,j],support=S))
 new_edges=[e for e in edges if max(e['vertices'])>=98]
 assert [e for e in edges if max(e['vertices'])<98]==old['finite_edges']
 assert [e['vertices']for e in new_edges]==[[6,98],[8,99],[98,99]] and len(edges)==266
 print('Complete2 new finite vertices;100 total;old6--8 cycle closed PASS',flush=True)
 # Complete six-pair line census in the support-perpendicular Lorentz3 space.
 FG=[[base(gram[i][j])for j in FACE]for i in FACE];assert R.pd(FG)[0];FI=C.inverse(FG)
 boundary=[7,8,11,35];projected=[]
 for i in boundary:
  cc=C.mv(FI,[base(gram[i][j])for j in FACE])
  projected.append([normals[i][t]-sum((Y(a)*normals[j][t]for a,j in zip(cc,FACE)),Y())for t in range(10)])
 BG=[[pair(a,b,Q)for b in projected]for a in projected]
 frame=projected[:3];framegram=[[base(pair(a,b,Q))for b in frame]for a in frame];assert len(C.eliminate(framegram)[1])==3
 fcov=[[pair(n,b,Q)for b in frame]for n in projected]
 line_tests=[];finite=[]
 for i,j in combinations(range(4),2):
  cc=cross(fcov[i],fcov[j]);V=[sum((a*v[t]for a,v in zip(cc,frame)),Y())for t in range(10)]
  assert all(pair(V,n,Q)==Y()for n in[projected[i],projected[j]])
  assert any(a.sign()for a in V)
  signs=[]
  for side in[-1,1]:
   W=[side*a for a in V];vals=[pair(W,n,Q)for n in projected]
   signs.append(dict(side=side,boundary_pairing_signs=[a.sign()for a in vals]))
  good=[a for a in signs if max(a['boundary_pairing_signs'])<=0]
  if good:
   assert len(good)==1;W=[good[0]['side']*a for a in V];vals=[dot(W,c)for c in cov]
   incident=sorted(FACE+[boundary[i],boundary[j]])
   assert pair(W,W,Q).sign()<0 and pair(W,H,Q).sign()<0
   assert [a.sign()for a in vals]==[0 if t in incident else -1 for t in range(36)]
   finite.append(dict(boundary_pair=[boundary[i],boundary[j]],incident=incident,point=W,norm=pair(W,W,Q),pairings=vals))
  else:assert pair(V,V,Q).sign()>0
  line_tests.append(dict(pair=[boundary[i],boundary[j]],norm_sign=pair(V,V,Q).sign(),orientations=signs,feasible=bool(good)))
 assert len(finite)==4
 face_ids=[6,8,98,99];assert {tuple(v['incident'])for v in finite}=={tuple(vertices[i]['incident'])for i in face_ids}
 facecenter=[sum((v['point'][t]for v in finite),Y())for t in range(10)];fp=[dot(facecenter,c)for c in cov]
 assert pair(facecenter,facecenter,Q).sign()<0 and [a.sign()for a in fp]==[0 if i in FACE else -1 for i in range(36)]
 # A literal surviving exterior ray gives the next strict escape.
 source=read(BASE32);escape_index=9;active=source['complete_census']['hyperideal_rays'][escape_index]
 assert active==[0,1,2,4,5,6,7,10,27]
 V=C.null_one([[base(x)for x in cov[i]]for i in active]);scale=-C.pair(V,H,Q);V=[a/scale for a in V]
 rp=[dot(V,c)for c in cov];assert R.sg(C.pair(V,V,Q))>0 and all(a.sign()<=0 for a in rp)
 E=[v+K.convert(F(1,8192))*h for v,h in zip(V,H)];ep=[dot(E,c)for c in cov];en=C.pair(E,E,Q)
 assert R.sg(en)>0 and R.sg(C.pair(E,H,Q))<0 and all(a.sign()<0 for a in ep)
 cc=C.mv(FI,[base(dot(E,cov[i]))for i in FACE]);EP=[Y(E[t])-sum((Y(a)*normals[i][t]for a,i in zip(cc,FACE)),Y())for t in range(10)]
 epp=[dot(EP,c)for c in cov];assert pair(EP,EP,Q).sign()<0 and pair(EP,H,Q).sign()<0
 assert [a.sign()for a in epp]==[0 if i in FACE else -1 for i in range(36)]
 paths=[f'scripts/{TAG}.py',PREVIOUS,REVIEW,'scripts/fresh_nonarithmetic_escape34_20261016.py',BASE32]
 return enc(dict(status='PASS_EXACT_ESCAPE_CUTTER_AND_COMPACT_QUADRILATERAL_FACE',sources={p:sha(p)for p in paths},
  field=dict(extra_radical_square=D,extra_radical_expression='10+68*sqrt(5)/15',encoding='[a,b] means a+b*positive_sqrt(D), with a,b in the previous formal tower; no minimal degree assertion'),form=Q,center=H,normals=normals,norms=norms,gram=gram,
  cap=dict(index=35,right_support=SUPPORT,anchor=8,label=6,linear_pivots=pivots,particular=p,kernel=k,particular_norm=pp,kernel_norm=kk,cross_pairing=pk,normal=r,norm=norms[35],branch=-1,branches=branches,relations=relations,old_escape_pairing=dot(oldE,cov[35]),center_pairing=hp[35]),
  pair_counts=dict(counts),actual_facets=36,facet_witnesses=facets,old_vertex_preservation_pairings=preservation,
  complete_new_cap_finite_census=dict(finite_neighbors=neighbors,candidate_supports=candidates,vertices=new),finite_vertices=vertices,finite_vertex_incidence=[v['incident']for v in vertices],finite_vertex_count=100,finite_edges=edges,finite_edge_count=266,new_finite_edges=new_edges,open_finite_germ_count=368,finite_graph_components=B.P.components(100,edges),
  compact_two_face=dict(support=FACE,support_gram=FG,inverse_support_gram=FI,boundary=boundary,projected_normals=projected,projected_boundary_gram=BG,frame=frame,frame_gram=framegram,complete_six_line_tests=line_tests,finite_rays=finite,vertices=face_ids,strict_center=facecenter,center_pairings=fp),
  old_escape=dict(point=oldE,norm=C.pair(oldE,oldE,Q),new_cap_pairing=dot(oldE,cov[35])),strict_escape=dict(point=E,norm=en,pairings=ep,ray=V,ray_active=active,ray_index=escape_index,ray_pairings=rp,center_epsilon='1/8192',compact_face_projection=EP,projection_norm=pair(EP,EP,Q),projection_pairings=epp),
  scope='Exact cutter of35-wall escape, with dotted contacts to tetrahedron whitewalls3/12. Complete100finitevertices/266edges;closes two old germs and a compact Lambert quadrilateral. New strict escape proves infinite volume. Any further coherent acute wall right to this quadrilateral support cannot cut the escape. Full global cone is separate.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();path=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(path.read_text());print('PASS saved exact quadrilateral cap replay')
 else:path.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['finite_vertex_count'],d['finite_edge_count'],d['open_finite_germ_count'])
