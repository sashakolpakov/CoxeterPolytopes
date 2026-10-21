#!/usr/bin/env python3
"""Exact degree-eight cap joining two components of the finite-edge graph."""
import argparse,json
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
import fresh_h9_octic_escape_caps_20261015 as O
A,D,C,ROOT=O.A,O.D,O.C,O.ROOT
X=O.X
TAG='fresh_h9_joint_corner_20261015_bridge'
BASE='results/fresh_h9_joint_pair_20261015_two_extension.json'
SCOUT='results/fresh_h9_joint_corner_20261015_scout.json'
def inverse(M):
 n=len(M);B=[[X(z)for z in row]+[X(i==j)for j in range(n)]for i,row in enumerate(M)]
 for k in range(n):
  i=next(i for i in range(k,n)if B[i][k]);B[k],B[i]=B[i],B[k];z=B[k][k];B[k]=[a/z for a in B[k]]
  for i in range(n):
   if i!=k:
    z=B[i][k]
    if z:B[i]=[a-z*b for a,b in zip(B[i],B[k])]
 return [r[n:]for r in B]
def pd(M):
 B=[list(r)for r in M]
 for k in range(len(B)):
  if B[k][k].sign()<=0:return False
  for i in range(k+1,len(B)):
   for j in range(k+1,len(B)):B[i][j]-=B[i][k]*B[k][j]/B[k][k]
 return True
def components(count,edges):
 ad=[set()for _ in range(count)]
 for i,j in edges:ad[i].add(j);ad[j].add(i)
 seen=set();out=[]
 for i in range(count):
  if i in seen:continue
  todo=[i];v=[]
  while todo:
   j=todo.pop()
   if j in seen:continue
   seen.add(j);v.append(j);todo.extend(ad[j]-seen)
  out.append(sorted(v))
 return sorted(out)
def rational_square(q):
 if q<0:return False
 return isqrt(q.numerator)**2==q.numerator and isqrt(q.denominator)**2==q.denominator
def run():
 d=O.read(BASE)
 for p,h in d['sources'].items():assert O.sha(p)==h,p
 A.RAD=A.A.Q(*d['field']['new_positive_radical_square']);C.Golden=A.L;C.ZERO=A.Z
 Q0,N0,H0,E0=D.dm(d['form']),D.dm(d['normals']),D.dv(d['center']),D.dv(d['strict_escape']['point'])
 S=[0,1,2,5,6,10,29,30];ids=S+[8];cov0=[A.mv(Q0,n)for n in N0]
 v=C._full_rank_kernel([cov0[i]for i in ids]);assert v is not None and v[9];v=[x/v[9]for x in v]
 phi=A.L(A.A.Q(F(1,2),F(1,2)));rhs=[A.L(-1)if i==30 else-phi if i==8 else A.Z for i in ids]
 p=A.mv(A.inv([cov0[i]for i in ids]+[A.mv(Q0,v)]),rhs+[A.Z]);assert not A.pair(p,v,Q0)
 delta=(A.L(2)-A.pair(p,p,Q0))/A.pair(v,v,Q0)
 assert delta==A.L(A.A.Q(F(1988340828,1990805),F(4448872824,9954025)))
 nonsquares=[]
 for z in(delta.a,delta.a/A.RAD):
  norm=z.a*z.a-5*z.b*z.b;assert not rational_square(norm);nonsquares.append(norm)
 O.DELTA=delta;Q=[[X(a)for a in row]for row in Q0];N=[[X(a)for a in row]for row in N0];H=[X(a)for a in H0];E=[X(a)for a in E0]
 rows=[]
 for branch in(-1,1):
  n=[X(a,branch*b)for a,b in zip(p,v)];assert O.pair(n,n,Q)==2
  if O.pair(n,H,Q).sign()>=0:continue
  vals=[O.pair(n,m,Q)for m in N];relations=[]
  for i,g in enumerate(vals):
   kind,m,c=O.relation(g,X(2),O.pair(N[i],N[i],Q),5)
   if kind in('obtuse','unclassified'):break
   relations.append(dict(wall=i,kind=kind,label=m,invariant=c))
  else:rows.append((branch,n,vals,relations))
 assert len(rows)==1;branch,n,vals,relations=rows[0];N.append(n);cov=[O.mv(Q,x)for x in N]
 norm=[O.pair(x,x,Q)for x in N];G=[[O.dot(x,c)for c in cov]for x in N]
 assert O.pair(n,E,Q).sign()<0 # previous strict escape survives exactly
 facets=[]
 for i,x in enumerate(N):
  t=O.pair(H,x,Q)/norm[i];V=[a-t*b for a,b in zip(H,x)];vs=[O.dot(V,c).sign()for c in cov]
  assert O.pair(V,V,Q).sign()<0 and vs==[-1]*i+[0]+[-1]*(32-i)
  facets.append(dict(wall=i,point=V,norm=O.pair(V,V,Q)))
 old_preservation=[]
 for vertex in d['finite_vertices']:
  value=O.pair([X(a)for a in D.dv(vertex['point'])],n,Q);assert value.sign()<0;old_preservation.append(value)
 neighbors=[r['wall']for r in relations if r['kind']in('right','finite')]
 assert neighbors==sorted(ids);new=[];candidates=[]
 for SS in combinations(neighbors,8):
  II=list(SS)+[32];positive=pd([[G[i][j]for j in II]for i in II]);candidates.append(dict(support=II,positive_definite=positive))
  if not positive:continue
  V=O.mv(inverse([cov[i]for i in II]+[O.mv(Q,H)]),[X()]*9+[X(-1)])
  vp=[O.dot(V,c)for c in cov];assert O.pair(V,V,Q).sign()<0
  assert all(x.sign()==(0 if i in II else -1)for i,x in enumerate(vp))
  new.append(dict(incident=II,point=V,norm=O.pair(V,V,Q),pairings=vp))
 assert [r['incident']for r in new]==[[0,1,2,5,6,8,10,29,32],[0,1,2,5,6,10,29,30,32]]
 incident=[v['incident']for v in d['finite_vertices']]+[v['incident']for v in new]
 edges=[dict(vertices=[i,j],support=sorted(set(a)&set(b)))for i,a in enumerate(incident)for j,b in enumerate(incident[i+1:],i+1)if len(set(a)&set(b))==8]
 old_edges=[r['vertices']for r in d['finite_edges']];comps=components(len(incident),[r['vertices']for r in edges]);oldcomps=components(90,old_edges)
 assert len(incident)==92 and len(edges)==253 and len(oldcomps)==8 and len(comps)==7
 joint_index=next(i for i,r in enumerate(d['finite_vertices'])if 30 in r['incident']and 31 in r['incident'])
 merged=next(x for x in comps if joint_index in x);assert merged==[19,20,21,90,91]
 ep=[O.dot(E,c)for c in cov];assert O.pair(E,E,Q).sign()>0 and all(x.sign()<0 for x in ep)
 files=[BASE,SCOUT,f'scripts/{TAG}.py','scripts/fresh_h9_octic_escape_caps_20261015.py',
        'scripts/fresh_h9_quartet_polar_second_20261014.py','scripts/fresh_h9_quartet_polar_20261014.py',
        'scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py']
 return O.enc(dict(status='PASS_EXACT_OCTIC_CAP_JOINING_FINITE_EDGE_COMPONENTS',sources={f:O.sha(f)for f in files},
   field=dict(base=d['field'],extra_positive_radical_square=delta,encoding='[a,b] denotes a+b sqrt(delta), each coefficient in oldK',relative_degree=2,absolute_degree=8,nonsquare_norm_tests=nonsquares),
   cap=dict(normal=n,norm=X(2),branch=branch,reconstruction=dict(ids=ids,particular=p,kernel=v,delta=delta),relations=relations,pairings=vals,center_pairing=O.pair(n,H,Q)),
   actual_facets=33,facet_witnesses=facets,old_finite_vertex_preservation_pairings=old_preservation,
   complete_new_cap_finite_census=dict(finite_neighbors=neighbors,candidate_supports=candidates,vertices=new),
   finite_vertex_incidence=incident,finite_edges=edges,finite_vertex_count=92,finite_edge_count=253,open_finite_germ_count=9*92-2*253,
   old_finite_graph_components=oldcomps,finite_graph_components=comps,joint_corner_joined_component=merged,
   retained_strict_escape=dict(point=E,norm=O.pair(E,E,Q),pairings=ep,cap_pairing=ep[-1]),
   scope='Exact33actualfacets;32newpairgates and allold496 bound. All90oldfinitevertices persist; exhaustive9 finite-neighbor support tests give2newvertices and merge two old finite-edge components. Full finite-vertex/edge census follows, but no full exterior-ray census is claimed. The same strict escape proves infinite volume.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['finite_vertex_count'],d['finite_edge_count'],d['open_finite_germ_count'])
