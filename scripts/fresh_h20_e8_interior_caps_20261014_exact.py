#!/usr/bin/env python3
"""Exact coupled interior-E8 cap pair, retaining all frozen79 walls.
Field basis sqrt2^i sqrt79^j sqrt11^k; rational interval signs, no floats.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from collections import Counter
from math import isqrt
import fresh_h20_coupled_lift_20261014 as rat
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_interior_caps_20261014_exact'
BASE='results/fresh_h20_e8_cusp_completion_20261014.json';FACE='results/fresh_h20_contact_closure_20261014_h12_replacement.json';SOURCE='results/fresh_h12_a12_source_20261014.json';SCOUT='results/fresh_h20_e8_orbit_caps_20261014_volume_scout.json'
Z=(F(0),)*8;RAD=[1,2,79,158,11,22,869,1738]
def c(x):return(F(x),)+(F(0),)*7
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def scale(x,s):return tuple(a*s for a in x)
def mul(x,y):
 out=[F(0)]*8
 for i,a in enumerate(x):
  if not a:continue
  for j,b in enumerate(y):
   if b:out[i^j]+=a*b*RAD[i&j]
 return tuple(out)
DEN=1<<160;BOUNDS=[]
for q in RAD:
 n=isqrt(q*DEN*DEN);BOUNDS.append((F(n,DEN),F(n+int(n*n!=q*DEN*DEN),DEN)))
def sign(x):
 if x==Z:return 0
 lo=sum(a*(lo if a>=0 else hi)for a,(lo,hi)in zip(x,BOUNDS));hi=sum(a*(hi if a>=0 else lo)for a,(lo,hi)in zip(x,BOUNDS))
 assert lo>0 or hi<0,('unresolved',x,lo,hi)
 return 1 if lo>0 else -1
def embed(x):return tuple(map(F,x))+(F(0),)*4
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d,f,s,sc=[json.loads((ROOT/p).read_text())for p in [BASE,FACE,SOURCE,SCOUT]]
 for data in [d,f,s,sc]:
  for p,h in data['sources'].items():assert sha(p)==h,p
 Q=[[F(x)for x in row]for row in d['form']];old=[[embed(x)for x in row]for row in d['normals']];H=[embed(x)for x in d['strict_center']];oldescape=[embed(x)for x in d['strict_escape']['point']]
 cache={};dotcache={}
 def dot(x,y):
  pairkey=(tuple(x),tuple(y))
  if pairkey in dotcache:return dotcache[pairkey]
  key=tuple(y)
  if key not in cache:cache[key]=[tuple(sum(Q[i][j]*y[j][k]for j in range(21))for k in range(8))for i in range(21)]
  out=Z
  for a,b in zip(x,cache[key]):out=add(out,mul(a,b))
  dotcache[pairkey]=out;dotcache[(pairkey[1],pairkey[0])]=out
  return out
 def linear(vectors,coefficients):
  out=[Z]*21
  for a,v in zip(coefficients,vectors):out=[add(x,mul(a,y))for x,y in zip(out,v)]
  return out
 E=f['E8_support'];I=[[F(x)for x in row]for row in f['E8_inverse']]
 B=[[dot(old[i],old[j])for j in E]for i in E];assert all(x[1:]==(0,)*7 for row in B for x in row);B=[[x[0]for x in row]for row in B]
 rat.pd(B);assert rat.inv(B)==I
 vrows=f['H12_face_finite_vertices']+f['H12_face_ideal_vertices'];by_inc={tuple(row['source_incident']):row for row in vrows}
 V0,V1=[[embed(x)for x in by_inc[tuple(row['incident'])]['point']]for row in s['projective_vertices'][:2]]
 assert dot(V0,V0)==c(-1)and dot(V1,V1)==c(F(-13,121))and dot(V0,V1)==c(-1)
 coefficients=[
  ((F(17,8),0,0,0,0,0,F(-1,24),0),(F(-2057,104),0,F(187,52),0,F(-99,52),0,F(121,312),0)),
  ((F(33,16),0,0,0,0,0,F(-1,48),0),(F(-3993,208),0,F(363,104),0,F(-99,104),0,F(121,624),0))]
 profiles=[[3,6],[4,5]];new=[];recipes=[]
 for k,(profile,(a,b))in enumerate(zip(profiles,coefficients)):
  a=tuple(map(F,a));b=tuple(map(F,b));assert sign(a)>0 and sign(b)>0
  weights=[sum(I[i][j]for i in profile)for j in range(8)];w=linear([old[i]for i in E],[c(x)for x in weights])
  y=linear([V0,V1],[a,b]);r=[add(scale(x,-1),v)for x,v in zip(w,y)]
  assert dot(r,r)==c(2)and dot(w,y)==Z and sign(dot(y,y))<0
  assert all(dot(r,old[i])==c(-int(j in profile))for j,i in enumerate(E))
  assert all(dot(r,old[i])==Z for i in range(11))and dot(r,old[19])==c(-int(k==0))
  assert sign(dot(r,oldescape))>0
  new.append(r);recipes.append(dict(E8_positions=profile,E8_walls=[E[i]for i in profile],weight_coefficients=weights,weight_norm=dot(w,w),source_vertex_indices=[0,1],a=a,b=b,anchor_wall=19,anchor_pairing=dot(r,old[19]),old_escape_pairing=dot(r,oldescape)))
 assert dot(*new)==Z
 primitive=new;new=[];orbit=[]
 for ri,r in enumerate(primitive):
  for h19,h20 in [(1,1),(1,-1),(-1,1),(-1,-1)]:
   v=r[:];v[19]=scale(v[19],h19);v[20]=scale(v[20],h20);new.append(v);orbit.append(dict(primitive=ri,height_signs=[h19,h20]))
 assert len(set(tuple(r)for r in new))==8
 N=old+new;norms=[]
 for r in N:
  q=dot(r,r);assert q[1:]==(0,)*7 and q[0]>0;norms.append(q[0])
 counts=Counter();pair_rows=[]
 for i in range(79,87):
  for j in range(i):
   g=dot(N[i],N[j]);assert sign(g)<=0;sq=mul(g,g);prod=norms[i]*norms[j]
   if g==Z:kind='right'
   elif sq==c(prod/4):kind='pi3'
   else:assert sign(add(sq,c(-prod)))>0;kind='dotted'
   counts[kind]+=1;pair_rows.append(dict(walls=[i,j],pairing=g,relation=kind))
 assert counts==Counter(right=144,pi3=20,dotted=496)and len(pair_rows)==660
 assert sign(dot(H,H))<0 and all(sign(dot(H,r))<0 for r in N)
 facetpoints=[]
 for i,r in enumerate(N):
  p=[add(x,mul(scale(dot(H,r),-1/norms[i]),y))for x,y in zip(H,r)]
  assert sign(dot(p,p))<0 and sign(dot(p,H))<0
  assert all(dot(p,w)==Z if j==i else sign(dot(p,w))<0 for j,w in enumerate(N));facetpoints.append(p)
 # Complete finite-corner census only for corners incident to either newcap.
 vertices={};candidate_checks=[]
 for i in range(79,87):
  neighbors=[j for j in range(87)if j!=i and sign(add(mul(dot(N[i],N[j]),dot(N[i],N[j])),c(-norms[i]*norms[j])))<0]
  assert len(neighbors)==21;accepted=[];tested=0
  for part in combinations(neighbors,19):
   J=sorted(list(part)+[i]);tested+=1;GG=[[dot(N[a],N[b])for b in J]for a in J]
   if any(sign(add(mul(GG[a][b],GG[a][b]),c(-norms[J[a]]*norms[J[b]])))>0 for a in range(20)for b in range(a)):continue
   assert all(x[1:]==(0,)*7 for row in GG for x in row);G=[[x[0]for x in row]for row in GG]
   try:rat.pd(G)
   except AssertionError:continue
   GI=rat.inv(G);pH=[dot(N[j],H)for j in J];coef=[tuple(sum(GI[a][b]*pH[b][k]for b in range(20))for k in range(8))for a in range(20)]
   proj=linear([N[j]for j in J],coef);p=[add(x,scale(y,-1))for x,y in zip(H,proj)]
   assert sign(dot(p,p))<0 and sign(dot(p,H))<0
   assert all(dot(p,w)==Z if j in J else sign(dot(p,w))<0 for j,w in enumerate(N))
   vertices[tuple(J)]=dict(incident=J,point=p,norm=dot(p,p),gram=G);accepted.append(J)
  assert tested==210 and len(accepted)==(1 if i<83 else 2),(i,tested,len(accepted));candidate_checks.append(dict(wall=i,finite_neighbors=neighbors,subsets_tested=tested,accepted=accepted))
 assert len(vertices)==8
 for J,vertex in vertices.items():
  unseen=set(J);components=[]
  while unseen:
   todo=[min(unseen)];C=set()
   while todo:
    i=todo.pop()
    if i in C:continue
    C.add(i);todo.extend(j for j in J if j!=i and dot(N[i],N[j])!=Z)
   unseen-=C;adj={i:[j for j in C if j!=i and dot(N[i],N[j])!=Z]for i in C}
   assert all(norms[i]==2 for i in C)and all(dot(N[i],N[j])==c(-1)for i in C for j in adj[i])
   assert sum(len(v)for v in adj.values())==2*(len(C)-1)
   branch=[i for i in C if len(adj[i])==3]
   if not branch:assert max(map(len,adj.values()))<=2;kind='A'+str(len(C))
   else:
    assert len(branch)==1 and len(C)==7;center=branch[0];arms=[]
    for n in adj[center]:
     previous=center;current=n;length=1
     while len(adj[current])==2:
      nxt=next(i for i in adj[current]if i!=previous);previous,current=current,nxt;length+=1
     arms.append(length)
    assert sorted(arms)==[1,2,3];kind='E7'
   components.append(kind)
  vertex['spherical_type']=sorted(components)
 assert Counter(tuple(row['spherical_type'])for row in vertices.values())=={('A11','A2','A7'):4,('A11','A2','E7'):4}
 vlist=list(vertices.values());edges=[]
 for i,j in combinations(range(8),2):
  S=sorted(set(vlist[i]['incident'])&set(vlist[j]['incident']))
  if len(S)==19:
   G=[[dot(N[a],N[b])for b in S]for a in S];assert all(x[1:]==(0,)*7 for row in G for x in row);rat.pd([[x[0]for x in row]for row in G]);edges.append(dict(vertices=[i,j],support=S))
 assert len(edges)==4
 preservation=[]
 for row in f['H12_face_finite_vertices']:
  p=[embed(x)for x in row['point']];values=[dot(p,r)for r in new];assert all(sign(g)<0 for g in values);preservation.append(dict(incident=row['incident'],new_pairings=values))
 # The two already completed cusp links are preserved: timelike y is strict.
 cusp_signs=[]
 for row in d['completed_cusp_links']:
  nu=[embed(x)for x in row['point']];g=[dot(nu,r)for r in new];assert all(sign(x)<0 for x in g);cusp_signs.append(g)
 # A numerical proposal is rounded, then moved inward by a rational amount.
 v=sc['escapes'][0]['point'];escape=[add(c(F(round(x*10**9),10**9)),scale(h,F(1,10000)))for x,h in zip(v,H)]
 assert sign(dot(escape,escape))>0 and sign(dot(escape,H))<0 and all(sign(dot(escape,r))<0 for r in N)
 paths=[BASE,FACE,SOURCE,SCOUT,f'scripts/{TAG}.py','scripts/fresh_h20_coupled_lift_20261014.py']
 return enc(dict(status='PASS_87_WALL_INTERIOR_E8_ORBIT_CUTS_OLD_ESCAPE_8_VERTICES_4_EDGES_STILL_ESCAPES',sources={p:sha(p)for p in paths},field_radicands=RAD,form=Q,normals=N,norms=norms,new_recipes=recipes,orbit=orbit,new_pair_counts=dict(counts),new_pairs=pair_rows,strict_center=H,facet_points=facetpoints,new_finite_vertices=vlist,new_finite_edges=edges,new_cap_finite_vertex_census=candidate_checks,preserved_H12_finite_vertices=preservation,preserved_two_cusp_links=cusp_signs,previous_escape=[dot(oldescape,r)for r in new],strict_escape=dict(point=escape,norm=dot(escape,escape),center_pairing=dot(escape,H),all87_strict=True),conditional_nonquasi=d['conditional_nonquasi'],scope='Retains all79 and adds all4heightflips of an orthogonal norm2 pair with nonzero two-node E8 contacts. Original pair cuts saved79escape. All87actualfacets,660newpairs,complete new-cap finitecorner census gives8vertices/4finiteedges,35facevertices and2boundedcusps preserved. A different strict escape remains, so infinitevolume. No full ambient corner/cusp census.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],d['strict_escape']['norm'])
