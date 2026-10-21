#!/usr/bin/env python3
"""Exact rational pi/3 finite-edge bridges in the nonarithmetic H18 lift.

The old60 is retained; frozen92 is a frame/escape source, not retained walls.
All acceptance uses Q(sqrt2), rational Gram elimination and exact signs.
"""
import argparse,json,sys
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import product
from collections import Counter
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from coxeter_polytopes.exact.silver import Silver as K,silver_inertia
import fresh_dimension_construction_routes_20261014 as b
TAG='fresh_h20_contact_closure_20261014';SOURCE='results/fresh_dimension_construction_routes_20261014.json';PRIOR='results/fresh_h20_coupled_lift_20261014.json';SCOUT=f'results/{TAG}_scout.json';Z=K(0)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):return b.enc(x)
def null(A):
 R,p=b.rref(A);cols=[]
 for j in range(len(A[0])):
  if j in p:continue
  x=[Z]*len(A[0]);x[j]=K(1)
  for row,k in zip(R,p):x[k]=-row[j]
  cols.append(x)
 return cols

def pd(A):
 A=[r[:]for r in A];out=[]
 for i in range(len(A)):
  p=A[i][i];assert p>0;out.append(p)
  for j in range(i+1,len(A)):
   for k in range(i+1,len(A)):A[j][k]-=A[j][i]*A[i][k]/p
 return out

def run():
 src=json.loads((ROOT/SOURCE).read_text());prior=json.loads((ROOT/PRIOR).read_text());sc=json.loads((ROOT/SCOUT).read_text())
 for d in[src,prior]:
  for p,h in d['sources'].items():assert sha(p)==h,p
 Q=[[K.from_json(x)for x in row]for row in src['form']];old=[list(map(K.from_json,row))for row in src['normals']];assert len(old)==60
 @lru_cache(None)
 def cov(t):return tuple(b.matvec(Q,t))
 def dot(x,y):return sum((a*c for a,c in zip(x,cov(tuple(y)))),Z)
 H=[K(F(x))for x in prior['strict_center']];E=[K(F(x))for x in prior['strict_escape']['point']];v=list(map(K.from_json,src['strict_escape']['boundary_point']))
 assert dot(H,H)<0 and dot(E,E)>0 and dot(E,H)<0
 supports=[r['support']for r in prior['base_recipes']]
 S=supports[1];X=null([b.matvec(Q,old[i])for i in S]);assert len(X)==3
 assert [[x[i]for x in X]for i in[11,19,20]]==[[K(i==j)for j in range(3)]for i in range(3)]
 Q3=[[dot(x,y)for y in X]for x in X];assert Q3==[[K(F(x))for x in row]for row in[['143/2','-99/2','165'],['-99/2','71/2','-115'],['165','-115','381']]]
 def comb(cols,c):return[sum((x[i]*t for x,t in zip(cols,c)),Z)for i in range(21)]
 r=comb(X,[K(F(-19,2)),K(F(13,2)),K(6)]);s=comb(X,list(map(K,[-191,99,110])))
 assert dot(r,r)==2 and dot(s,s)==418 and dot(r,s)==0 and dot(r,old[1])==-1
 frames=[]
 for si in[1,3]:
  rec=prior['base_recipes'][si];hs=[K(F(x))for x in rec['h0']];w=[K(F(x))for x in rec['w']];P=[v,hs,w]
  assert all(dot(x,old[i])==0 for x in P for i in supports[si])
  frames.append(P)
 assert [[dot(x,y)for y in frames[0]]for x in frames[0]]==[[dot(x,y)for y in frames[1]]for x in frames[1]]
 assert [dot(x,old[1])for x in frames[0]]==[dot(x,old[20])for x in frames[1]]
 I=b.inverse([[dot(x,y)for y in frames[0]]for x in frames[0]])
 base=[r,s]+[comb(frames[1],b.matvec(I,[dot(x,z)for x in frames[0]]))for z in[r,s]]
 assert [dot(x,x)for x in base]==list(map(K,[2,418,2,418]))
 normals=old[:];vertices=[];edges=[];groups=[];perms=[]
 oldkeys={tuple(v):i for i,v in enumerate(old)};assert len(oldkeys)==60
 for signs in product([1,-1],repeat=2):
  def flip(x):
   x=x[:];x[19]*=signs[0];x[20]*=signs[1];return x
  perm=[oldkeys[tuple(flip(x))]for x in old];assert sorted(perm)==list(range(60));perms.append(perm)
  begin=len(normals);normals.extend(flip(x)for x in base)
  for k,si in enumerate([1,3]):
   rr=begin+2*k;ss=rr+1;a=perm[1 if si==1 else 20];support=[perm[i]for i in supports[si]]
   inc1=support+[rr,ss];inc2=support+[rr,a]
   groups.append(dict(template=si,signs=signs,new_walls=[rr,ss],anchor=a,support=support))
   edges.append(dict(support=support+[rr],endpoints=[len(vertices),len(vertices)+1]))
   vertices.extend(dict(incident=inc)for inc in[inc1,inc2])
 assert len(normals)==76
 norms=[dot(x,x)for x in normals];G=[[dot(x,y)for y in normals]for x in normals];counts=Counter()
 for i in range(60,76):
  for j in range(i):counts[b.pairgate(G[i][j],norms[i],norms[j])]+=1
 assert counts=={'right':296,'pi3':8,'dotted':776},counts
 assert all(dot(x,H)<0 for x in normals)
 # Strict orthogonal projection gives every retained/new hyperplane an actual facet point.
 facetpoints=[]
 for i,r0 in enumerate(normals):
  p=b.add(H,b.scale(r0,-dot(H,r0)/norms[i]));assert dot(p,p)<0 and dot(p,H)<0
  assert all(dot(p,w)==0 if j==i else dot(p,w)<0 for j,w in enumerate(normals))
  facetpoints.append(p)
 for row in vertices:
  J=row['incident'];B=[[G[i][j]for j in J]for i in J];ps=pd(B);c=b.matvec(b.inverse(B),[dot(normals[i],H)for i in J]);p=H[:]
  for a,i in zip(c,J):p=b.add(p,b.scale(normals[i],-a))
  assert dot(p,p)<0 and dot(p,H)<0
  assert all(dot(p,w)==0 if j in J else dot(p,w)<0 for j,w in enumerate(normals))
  row.update(point=p,norm=dot(p,p),pd20_pivots=ps)
 assert len({tuple(sorted(v['incident']))for v in vertices})==16
 for edge in edges:
  a,c=[vertices[i]for i in edge['endpoints']];assert set(a['incident'])&set(c['incident'])==set(edge['support']) and len(edge['support'])==19
  assert a['point']!=c['point'];pd([[G[i][j]for j in edge['support']]for i in edge['support']])
 # All four old escapes are cut, but the frozen92 escape remains a strict escape.
 oldE=list(map(K.from_json,src['strict_escape']['point']));targetcuts=[]
 for gi,signs in enumerate(product([1,-1],repeat=2)):
  T=oldE[:];T[19]*=signs[0];T[20]*=signs[1];values=[dot(normals[i],T)for i in range(60+4*gi,64+4*gi)];assert all(x>0 for x in values);targetcuts.append(values)
 assert all(dot(x,E)<0 for x in normals)
 cuspchecks=[]
 for cusp in src['repaired_cusps']:
  delta=list(map(K.from_json,cusp['delta']));vals=[dot(x,delta)for x in normals[60:]];assert all(x<0 for x in vals);cuspchecks.append(dict(pair=cusp['pair'],new_pairings=vals))
 # Exact positive-dual certificates on ONLY the eight unreflected-height19 supports.
 # No numerical LP feasibility claim is used for the other eight supports.
 dual=[]
 for row in sc['dual_proposals']:
  if not row['dual_feasible']:continue
  signs=row['signs'];perm=perms[list(product([1,-1],repeat=2)).index(tuple(signs))];S=[perm[i]for i in supports[row['template']]];X0=null([b.matvec(Q,old[i])for i in S]);ids=row['dual_walls'];M=[[dot(old[i],x)for i in ids]for x in X0];weights=b.matvec(b.inverse(M),[dot(E,x)for x in X0]);assert all(t>0 for t in weights)
  residual=E[:]
  for t,i in zip(weights,ids):residual=b.add(residual,b.scale(old[i],-t))
  assert all(dot(residual,x)==0 for x in X0)
  assert len(b.rref([old[i]for i in S+ids])[1])==21
  dual.append(dict(template=row['template'],signs=signs,support=S,witness_walls=ids,positive_weights=weights))
 assert len(dual)==8 and all(r['signs'][0]==1 for r in dual)
 paths=[SOURCE,PRIOR,SCOUT,f'scripts/{TAG}.py',f'scripts/{TAG}_scout.py','scripts/fresh_dimension_construction_routes_20261014.py','src/coxeter_polytopes/exact/silver.py','src/coxeter_polytopes/exact/matrix.py']
 return enc(dict(status='PASS_76_WALL_RATIONAL_PI3_BRIDGES_8_FINITE_EDGES_WITH_STRICT_ESCAPE',sources={p:sha(p)for p in paths},form=Q,normals=normals,norms=norms,strict_center=H,center_norm=dot(H,H),new_pair_counts=dict(counts),base_cap_free_coordinates=[['-19/2','13/2','6'],['-191','99','110']],base_caps=base,groups=groups,actual_facet_points=facetpoints,finite_vertices=vertices,finite_edges=edges,old_escape_cuts=targetcuts,preserved_old_cusps=cuspchecks,strict_escape=dict(point=E,norm=dot(E,E),center_pairing=dot(E,H),all76_strict=True),positive_dual_support_obstructions=dual,dual_scope='For each of these eight supports, every nonzero old60-acute normal perpendicular to all18 support walls is strict negative on the saved escape, with no field/norm/angle restriction. Height19-reflected supports are not included.',conditional_nonquasi=src['conditional_nonquasi'],scope='Retains old60 and replaces optional32 caps of frozen92. Exact76 actual facets,16 specified finite vertices and8 actual finite edges. No full vertex/edge census. Strict escape proves infinite volume. Eight positive-dual support certificates are all-parameter only in their stated supports.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],len(d['positive_dual_support_obstructions']))
