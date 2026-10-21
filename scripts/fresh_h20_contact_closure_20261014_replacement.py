#!/usr/bin/env python3
"""Exact changed-support replacement:76 walls, new Q(sqrt79) edge caps.

Standalone Fraction arithmetic in basis1,sqrt2,sqrt79,sqrt158; no floating
acceptance. Numerical recipes/escape are proposals only.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from hashlib import sha256
from math import isqrt
import fresh_h20_coupled_lift_20261014 as rat
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_contact_closure_20261014_replacement';BASE='results/fresh_h20_contact_closure_20261014.json';REC='results/fresh_h20_contact_closure_20261014_backtracking_scout.json';SCOUT=f'results/{TAG}_scout.json';ZZ=(F(0),)*4

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def c(x):return(F(x),F(0),F(0),F(0))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def neg(x):return tuple(-a for a in x)
def scale(x,s):return tuple(a*s for a in x)
def mul(x,y):
 z=[F(0)]*4
 for i,a in enumerate(x):
  if not a:continue
  for j,b in enumerate(y):
   if not b:continue
   common=i&j;z[i^j]+=a*b*(2 if common&1 else 1)*(79 if common&2 else 1)
 return tuple(z)
SCALE=1<<120
BOUNDS=[]
for d in[1,2,79,158]:
 n=isqrt(d*SCALE*SCALE);BOUNDS.append((F(n,SCALE),F(n+(n*n!=d*SCALE*SCALE),SCALE)))
def sign(x):
 if x==ZZ:return 0
 lo=hi=F(0)
 for a,(l,u)in zip(x,BOUNDS):
  lo+=a*(l if a>=0 else u);hi+=a*(u if a>=0 else l)
 assert lo>0 or hi<0,('unresolved sign',x,lo,hi)
 return 1 if lo>0 else -1

def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 return x

def run():
 base=json.loads((ROOT/BASE).read_text());rec=json.loads((ROOT/REC).read_text());sc=json.loads((ROOT/SCOUT).read_text())
 for f,h in base['sources'].items():assert sha(f)==h,f
 Q=[[F(x['a'])for x in row]for row in base['form']];assert all(F(x['b'])==0 for row in base['form']for x in row)
 old=[[(F(x['a']),F(x['b']),F(0),F(0))for x in row]for row in base['normals']];H=[F(x['a'])for x in base['strict_center']];assert all(F(x['b'])==0 for x in base['strict_center'])
 def cov(v):return[tuple(sum((Q[i][j]*v[j][k]for j in range(21)),F(0))for k in range(4))for i in range(21)]
 cache={}
 def dot(x,y):
  key=tuple(y)
  if key not in cache:cache[key]=cov(y)
  v=ZZ
  for a,b in zip(x,cache[key]):v=add(v,mul(a,b))
  return v
 new=[];recipes=[]
 for ci in[1,3,5,7]:
  row=rec['candidates'][ci];J=row['support']+[row['anchor']];assert len(J)==20
  assert all(all(x[k]==0 for k in[1,2,3])for j in J for x in old[j])
  A=[[sum(Q[k][l]*old[j][l][0]for l in range(21))for k in range(21)]+[F(-int(j in row['profile']))if j in row['support']else F(-1)]for j in J]
  R,p=rat.rref(A);assert len(p)==20 and 20 not in p
  u=[F(0)]*21;w=[F(0)]*21;w[20]=1
  for rr,j in zip(R,p):u[j]=rr[21];w[j]=-rr[20]
  aa=rat.dot(w,w,Q);bb=2*rat.dot(u,w,Q);cc=rat.dot(u,u,Q)-2
  sg=1 if ci in[1,3]else-1
  assert [aa,bb,cc]==[F(54,242),F(-214*sg,242),F(35,242)]
  ta=F(107*sg,54);tb=F(11*sg,54)
  root=[(u[j]+w[j]*ta,F(0),w[j]*tb,F(0))for j in range(21)]
  assert dot(root,root)==c(2)
  for j in J:assert dot(root,old[j])==c(-int(j in row['profile']or j==row['anchor']))
  new.append(root);recipes.append(dict(candidate_index=ci,support_prior_indices=row['support'],profile=row['profile'],anchor=row['anchor'],parameter_a=ta,parameter_b=tb,radicand=79,rational_offset=u,rational_direction=w,norm_polynomial=[aa,bb,cc]))
 # Complete both height signs, replacing all eight norm418 mates together.
 for v in new[:]:
  w=v[:];w[19]=neg(w[19]);new.append(w)
 removed=list(range(61,76,2));keep=[i for i in range(76)if i not in removed];N=[old[i]for i in keep]+new;assert len(N)==76
 norms=[]
 for r in N:
  q=dot(r,r);assert q[1:]==(0,0,0)and q[0]>0;norms.append(q[0])
 assert norms[68:]==[2]*8
 G=[[dot(r,s)for s in N]for r in N];counts=Counter()
 for i in range(68,76):
  for j in range(i):
   g=G[i][j];sg=sign(g);assert sg<=0
   if not sg:kind='right'
   elif mul(g,g)==c(norms[i]*norms[j]/4):kind='pi3'
   else:assert sign(add(mul(g,g),c(-norms[i]*norms[j])))>0;kind='dotted'
   counts[kind]+=1
 assert sum(counts.values())==8*68+28
 HH=list(map(c,H));assert sign(dot(HH,HH))<0 and all(sign(dot(r,HH))<0 for r in N)
 facetpoints=[]
 for i,r in enumerate(N):
  coef=scale(dot(HH,r),-1/norms[i]);p=[add(x,mul(coef,y))for x,y in zip(HH,r)];assert sign(dot(p,p))<0 and sign(dot(p,HH))<0
  assert all(dot(p,s)==ZZ if j==i else sign(dot(p,s))<0 for j,s in enumerate(N));facetpoints.append(p)
 # Existing D13+D7 endpoints survive; enumerate every finite vertex involving a new cap.
 vertices=[]
 for row in base['finite_vertices']:
  if any(i in removed for i in row['incident']):continue
  point=[(F(x['a']),F(x['b']),F(0),F(0))for x in row['point']];J=[keep.index(i)for i in row['incident']]
  assert all(dot(point,r)==ZZ if i in J else sign(dot(point,r))<0 for i,r in enumerate(N));vertices.append(dict(kind='retained',incident=J,point=point,norm=dot(point,point)))
 def rational_matrix(J):
  A=[]
  for i in J:
   row=[]
   for j in J:assert G[i][j][1:]==(0,0,0);row.append(G[i][j][0])
   A.append(row)
  return A
 candidate_checks=[]
 for newi in range(68,76):
  finite=[j for j in range(76)if j!=newi and sign(add(mul(G[newi][j],G[newi][j]),c(-norms[newi]*norms[j])))<0];assert len(finite)==20
  accepted=[]
  for J0 in combinations(finite,19):
   J=list(J0)+[newi]
   # Any dotted pair cannot form a spherical or affine endpoint.
   if any(sign(add(mul(G[i][j],G[i][j]),c(-norms[i]*norms[j])))>0 for i,j in combinations(J,2)):continue
   A=rational_matrix(J)
   try:piv=rat.pd(A)
   except AssertionError:continue
   C=rat.inv(A);bvec=[dot(N[i],HH)for i in J];coefs=[tuple(sum((C[a][j]*bvec[j][k]for j in range(20)),F(0))for k in range(4))for a in range(20)];point=HH[:]
   for a,i in zip(coefs,J):point=[add(x,neg(mul(a,y)))for x,y in zip(point,N[i])]
   assert sign(dot(point,point))<0 and sign(dot(point,HH))<0
   assert all(dot(point,r)==ZZ if j in J else sign(dot(point,r))<0 for j,r in enumerate(N))
   accepted.append(len(vertices));vertices.append(dict(kind='new_A12_E8',incident=J,point=point,norm=dot(point,point),pd20_pivots=piv))
  assert len(accepted)==1;candidate_checks.append(dict(new_wall=newi,finite_neighbors=finite,actual_vertex=accepted[0]))
 assert len(vertices)==16
 types=[]
 for vertex in vertices:
  J=vertex['incident'];adj={i:{j for j in J if j!=i and G[i][j]!=ZZ}for i in J};unseen=set(J);components=[]
  while unseen:
   stack=[min(unseen)];nodes=set()
   while stack:
    a=stack.pop()
    if a in nodes:continue
    nodes.add(a);stack.extend(adj[a]-nodes)
   unseen-=nodes
   assert all(norms[i]==2 for i in nodes)
   assert all(G[i][j]==c(-1) for i in nodes for j in adj[i])
   assert sum(len(adj[i])for i in nodes)==2*(len(nodes)-1)
   branches=[i for i in nodes if len(adj[i])==3]
   if not branches:
    assert max(len(adj[i])for i in nodes)<=2;typ='A'+str(len(nodes))
   else:
    assert len(branches)==1 and max(len(adj[i])for i in nodes)==3
    arms=[];center=branches[0]
    for neighbor in adj[center]:
     prev=center;cur=neighbor;n=1
     while len(adj[cur])==2:
      nxt=next(i for i in adj[cur]if i!=prev);prev,cur=cur,nxt;n+=1
     arms.append(n)
    arms.sort()
    if arms[:2]==[1,1]:typ='D'+str(len(nodes))
    else:assert len(nodes)==8 and arms==[1,2,4];typ='E8'
   components.append(typ)
  vertex['spherical_type']=sorted(components);types.append(tuple(sorted(components)))
 assert Counter(types)=={('D13','D7'):8,('A12','E8'):8}
 edges=[]
 for i,j in combinations(range(16),2):
  S=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
  if len(S)==19:rat.pd(rational_matrix(S));edges.append(dict(vertices=[i,j],support=S))
 assert len(edges)==8
 # A new common strict escape is only proposed by the numerical LP.
 pp=sc['escapes'][0]['point'];E=[F(round(x*10**9),10**9)+h/1000 for x,h in zip(pp,H)];EE=list(map(c,E));assert sign(dot(EE,EE))>0 and sign(dot(EE,HH))<0 and all(sign(dot(EE,r))<0 for r in N)
 oldE=[(F(x['a']),F(x['b']),F(0),F(0))for x in base['strict_escape']['point']];cuts=[]
 for s1 in[1,-1]:
  for s2 in[1,-1]:
   t=oldE[:];t[19]=scale(t[19],s1);t[20]=scale(t[20],s2);owners=[i for i in range(68,76)if sign(dot(N[i],t))>0];assert owners;cuts.append(dict(signs=[s1,s2],owners=owners))
 source=json.loads((ROOT/'results/fresh_dimension_construction_routes_20261014.json').read_text());cusps=[]
 for row in source['repaired_cusps']:
  delta=[(F(x['a']),F(x['b']),F(0),F(0))for x in row['delta']];values=[dot(r,delta)for r in new];assert all(sign(g)<0 for g in values);cusps.append(dict(pair=row['pair'],new_pairings=values))
 paths=[BASE,REC,SCOUT,'scripts/fresh_h20_contact_closure_20261014_backtracking_scout.py',f'scripts/{TAG}.py',f'scripts/{TAG}_scout.py','scripts/fresh_h20_coupled_lift_20261014.py']
 return enc(dict(status='PASS_76_WALL_CHANGED_SUPPORT_REPLACEMENT_OLD_GAP_CUT_NEW_STRICT_ESCAPE',sources={f:sha(f)for f in paths},field_basis=['1','sqrt2','sqrt79','sqrt158'],form=Q,removed_predecessor_walls=removed,retained_predecessor_walls=keep,new_recipes=recipes,normals=N,norms=norms,new_pair_counts=dict(counts),strict_center=H,facet_points=facetpoints,finite_vertices=vertices,finite_edges=edges,new_cap_vertex_census=candidate_checks,previous_escape_sign_images_cut=cuts,preserved_three_cusp_links=cusps,strict_escape=dict(point=E,norm=dot(EE,EE),center_pairing=dot(EE,HH),all76_strict=True),conditional_nonquasi=base['conditional_nonquasi'],scope='Replaces all8 norm418 mates by8 norm2 Qsqrt79 caps; retains allold60 and8 rational pi3 caps. All76 actual facets. Eight old finite vertices survive and8 new vertices yield8 actual finite edges, crossing changed old supports. Each new cap has exactlyone finite vertex by its complete finite-neighbor clique census. Allfour sign images of prior92 escape are cut, but a new exactstrictescape remains. No whole chamber corner/ideal census.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],d['strict_escape']['norm'])
