#!/usr/bin/env python3
"""Exact 275-facet intrinsic H20 golden cusp construction; still infinite volume.

Standalone Fraction arithmetic in Q(phi), phi^2=phi+1. The spatial basis is
one PD19 unit-root edge Gram; two null coordinates satisfy D.U=-1.
"""
import argparse,json,hashlib
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_intrinsic_caps_20261014_golden_seed'
SRC='results/fresh_h20_intrinsic_caps_20261014_endpoint.json'
ZERO=(F(0),F(0));ONE=(F(1),F(0));PHI=(F(0),F(1))
def q(a,b=0):return F(a),F(b)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[0],-a[1]
def sub(a,b):return add(a,neg(b))
def mul(a,b):return a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1]
def scale(a,b):return a[0]*b,a[1]*b
def sign(v):
 a,b=v[0]+v[1]/2,v[1]/2
 if not b:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if (a>0)==(b>0):return (a>0)-(a<0)
 d=a*a-5*b*b;assert d
 return ((a>0)-(a<0))if d>0 else((b>0)-(b<0))
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,tuple):return [enc(v)for v in x]
 if isinstance(x,list):return [enc(v)for v in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def inv(A):
 n=len(A);M=[list(r)+[F(i==j)for j in range(n)]for i,r in enumerate(A)]
 for j in range(n):
  p=next(i for i in range(j,n)if M[i][j]);M[j],M[p]=M[p],M[j]
  d=M[j][j];M[j]=[x/d for x in M[j]]
  for i in range(n):
   if i!=j and M[i][j]:
    d=M[i][j];M[i]=[x-d*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def pd(A):
 A=[r[:]for r in A];p=[]
 while A:
  d=A[0][0];assert d>0;p.append(d)
  A=[[A[i][j]-A[i][0]*A[0][j]/d for j in range(1,len(A))]for i in range(1,len(A))]
 return p

def run():
 src=json.loads((ROOT/SRC).read_text())
 for p,h in src['sources'].items():assert sha(p)==h,p
 ids=src['seed_order'];at={w:i for i,w in enumerate(ids)};marks=dict(zip(ids,map(F,src['positive_affine_marks'])))
 H=[[F(x,2)for x in row]for row in src['normalized_seed_Gram']]
 edge=src['types'][0]['edge_support'];ix=[at[w]for w in edge];B=[[H[i][j]for j in ix]for i in ix];I=inv(B);piv=pd(B)
 assert len(edge)==19 and all(B[i][i]==1 for i in range(19))
 A=[w for w in ids if w in[4,12,14,21,25,28,30,31,39,48,49,63]]
 T=[23,50,200];E=[22,24,26,33,34,40,45];factors=[A,T,E]
 assert set(A+T+E)==set(ids) and list(map(len,factors))==[12,3,7]
 sizes=[ONE,PHI,ONE]
 spatial={w:[dot(row,[H[j][at[w]]for j in ix])for row in I]for w in ids}
 rdot=lambda v,w:dot(v,[dot(row,w)for row in B])
 assert all(rdot(spatial[i],spatial[j])==H[at[i]][at[j]]for i in ids for j in ids)
 assert all([sum(marks[w]*spatial[w][j]for w in f)for j in range(19)]==[F(0)]*19 for f in factors)
 # Unit-inradius vertices of each affine simplex in the common spatial basis.
 fv=[];tables=[]
 for f in factors:
  total=sum(marks[w]for w in f);vv={}
  for omitted in f:
   rhs=[F(0)if w not in f else F(1)-(total/marks[w]if w==omitted else 0)for w in edge]
   v=[dot(row,rhs)for row in I];vv[omitted]=v
   assert all(rdot(spatial[w],v)==(F(1)-(total/marks[w]if w==omitted else 0))for w in f)
   assert all(rdot(spatial[w],v)==0 for w in ids if w not in f)
  fv.append(vv);tables.append({(i,j):rdot(vv[i],vv[j])for i in f for j in f})
 assert all(rdot(v,w)==0 for aa,bb in combinations(fv,2)for v in aa.values()for w in bb.values())
 combos=list(product(*factors));assert len(combos)==252
 def vertexdot(a,b):return q(tables[0][a[0],b[0]]+tables[2][a[2],b[2]]+tables[1][a[1],b[1]],tables[1][a[1],b[1]])
 vertices=[]
 for om in combos:
  z=[q(fv[0][om[0]][j]+fv[2][om[2]][j],fv[1][om[1]][j])for j in range(19)]
  vertices.append(z)
 znorm=[vertexdot(o,o)for o in combos]
 def ldot(a,b):
  out=ZERO
  for i,row in enumerate(B):
   for j,x in enumerate(row):
    if x:out=add(out,scale(mul(a[i],b[j]),x))
  return sub(sub(out,mul(a[19],b[20])),mul(a[20],b[19]))
 walls=[]
 for w in ids:
  f=next(i for i,fac in enumerate(factors)if w in fac)
  walls.append([q(x)for x in spatial[w]]+[sizes[f],ZERO])
 walls.append([ZERO]*19+[q(-2),ONE])
 walls.extend([z+[scale(sub(n,ONE),F(1,2)),ONE]for z,n in zip(vertices,znorm)])
 assert len(walls)==275
 assert [ldot(w,w)for w in walls]==[ONE]*22+[q(4)]+[ONE]*252
 count=Counter();gates=[]
 for i,j in combinations(range(22),2):
  p=H[i][j];assert p in(F(0),F(-1,2));count['right'if p==0 else'label3']+=1
 for i,w in enumerate(ids):
  p=ldot(walls[i],walls[22]);expect=neg(PHI if w in T else ONE)
  assert p==expect;count['label5'if w in T else'label3']+=1
 min_vertex_norm=znorm[0]
 for n in znorm[1:]:
  if sign(sub(n,min_vertex_norm))<0:min_vertex_norm=n
 assert all(sign(sub(n,min_vertex_norm))>=0 for n in znorm)
 min_sep=None;min_sep_pair=None
 for v,om in enumerate(combos):
  n=znorm[v]
  assert sign(sub(n,q(9)))>0
  p=scale(sub(q(5),n),F(1,2));assert sign(add(p,q(2)))<0;count['dotted']+=1
  for i,w in enumerate(ids):
   f=next(j for j,fac in enumerate(factors)if w in fac)
   expected=ZERO if w!=om[f]else scale(sizes[f],-sum(marks[x]for x in factors[f])/marks[w])
   assert ldot(walls[i],walls[23+v])==expected
   if expected==ZERO:count['right']+=1
   else:assert sign(add(expected,ONE))<0;count['dotted']+=1
  for u in range(v):
   dist=sub(add(n,znorm[u]),scale(vertexdot(om,combos[u]),2))
   assert sign(sub(dist,q(4)))>0
   if min_sep is None or sign(sub(dist,min_sep))<0:min_sep=dist;min_sep_pair=[u,v]
   count['dotted']+=1
 assert dict(count)=={'right':4998,'label3':40,'label5':3,'dotted':32634}
 assert sum(count.values())==275*274//2
 center=[ZERO]*19+[q(F(9,2)),ONE]
 assert ldot(center,center)==q(-9)
 cp=[ldot(center,w)for w in walls];assert all(sign(x)<0 for x in cp)
 facet_point=[ZERO]*19+[q(2),ONE]
 assert ldot(facet_point,facet_point)==q(-4)
 fp=[ldot(facet_point,w)for w in walls]
 assert [i for i,p in enumerate(fp)if p==ZERO]==[22]and all(sign(p)<=0 for p in fp)
 endpoints=[]
 for v,(z,n,om)in enumerate(zip(vertices,znorm,combos)):
  x=z+[scale(add(n,ONE),F(1,2)),ONE]
  assert ldot(x,x)==q(-1)and sign(ldot(x,center))<0
  zeros=[]
  # Reconstruct all275 endpoint signs from independently checked dot tables.
  for i,w in enumerate(ids):
   p=ldot(x,walls[i]);assert sign(p)<=0
   if p==ZERO:zeros.append(i)
  central=scale(sub(q(3),n),F(1,2));assert sign(central)<0
  assert ldot(x,walls[22])==central
  for j,other in enumerate(combos):
   dist=sub(add(n,znorm[j]),scale(vertexdot(om,other),2));p=scale(dist,F(-1,2))
   if j==v:assert p==ZERO;zeros.append(23+j)
   else:assert sign(p)<0
  expected=[i for i,w in enumerate(ids)if w not in om]+[23+v]
  assert zeros==expected and len(zeros)==20
  # Old support has a PD19 Gram, and the new norm-one cap is right to it.
  tr=next(r for r in src['literal_edge_transports']if r['omitted']==list(om))
  base=src['types'][tr['type']]
  assert [[H[at[i]][at[j]]for j in tr['edge_support_in_profile_order']]for i in tr['edge_support_in_profile_order']]==[[F(x,2)for x in row]for row in base['PD19_Gram']]
  endpoints.append({'omitted_source_walls':list(om),'incident_wall_indices':zeros,'point':x,'norm':q(-1)})
 assert len({tuple(r['incident_wall_indices'])for r in endpoints})==252
 assert all(any(i in r['incident_wall_indices']for r in endpoints)for i in range(22))
 oldescape=[ZERO]*19+[q(F(-1,2)),ONE]
 assert ldot(oldescape,oldescape)==ONE and ldot(oldescape,walls[22])==q(F(5,2))
 assert all(sign(ldot(oldescape,w))<0 for i,w in enumerate(walls)if i!=22)
 z=[scale(v,F(1,2))for v in vertices[0]];n=scale(znorm[0],F(1,4));escape=z+[scale(sub(n,ONE),F(1,2)),ONE]
 assert ldot(escape,escape)==ONE and sign(ldot(escape,center))<0
 es=[ldot(escape,w)for w in walls];assert all(sign(p)<0 for p in es)
 # An indefinite rational block survives alongside an irrational cycle.
 # The claim of non-quasi-arithmeticity is conditional on a cofinite completion.
 assert len(A)==12 and all(ldot(walls[at[w]],walls[22])==q(-1)for w in A)
 assert mul(PHI,PHI)==q(1,1)
 paths=[SRC,'scripts/fresh_h20_intrinsic_caps_20261014_endpoint.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return enc({'status':'PASS_EXACT_275_FACET_INTRINSIC_H20_GOLDEN_SEED_WITH_STRICT_ESCAPE','sources':{p:sha(p)for p in paths},'field':'Q(phi), phi^2=phi+1, phi>1','coordinates':'19 spatial basis coordinates followed by D,U; D.U=-1, D^2=U^2=0','spatial_basis_source_walls':edge,'spatial_Gram':B,'spatial_LDL_pivots':piv,'Lorentz_signature':[20,1],'source_wall_order':ids,'positive_affine_marks':src['positive_affine_marks'],'cusp_link':'affine A11 x affine A2 x affine E6, inradii1,phi,1','walls':walls,'wall_norms':[ONE]*22+[q(4)]+[ONE]*252,'pair_counts':dict(count),'minimum_vertex_center_square':min_vertex_norm,'minimum_vertex_center_separation_square':min_sep,'minimum_separation_pair':min_sep_pair,'common_strict_center':{'point':center,'norm':q(-9),'all_pairings':cp},'central_facet_point':{'point':facet_point,'norm':q(-4),'all_pairings':fp},'finite_endpoints':endpoints,'all252_old_cusp_edges_end_finitely':True,'old_incenter_escape_cut':{'point':oldescape,'norm':ONE,'central_pairing':q(F(5,2))},'strict_remaining_escape':{'point':escape,'norm':ONE,'center_pairing':ldot(escape,center),'all275_pairings':es},'conditional_field_obstruction':{'central_triangle_two_cycle':q(1,1),'rational_block':'12 affineA11 unit roots plus central cap','restricted_Gram_on_sum_and_cap':[[0,-12],[-12,4]],'restricted_determinant':-144},'scope':'Explicit275 actual facets with exact Coxeter pair gates and252 actual finite endpoints, retaining the source cusp diagram but choosing intrinsic triangle inradiusphi. Still infinite volume: a strict spacelike escape survives all275 walls. The252 endpoints are a certified subset, not a complete global corner census.'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['pair_counts'])
