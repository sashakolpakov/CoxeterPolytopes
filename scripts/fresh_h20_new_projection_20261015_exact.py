#!/usr/bin/env python3
"""Exact new source projection, finite endpoint and renewed H20 escape."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h20_joint_route_20261015_exact as B
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_new_projection_20261015_exact'
SOURCE='results/fresh_h20_joint_route_20261015_exact.json';HELPER='scripts/fresh_h20_joint_route_20261015_exact.py'
P=B.P;K=B.K
class Z:
 def __init__(self,a=0,b=0):
  if isinstance(a,Z):self.a,self.b=a.a,a.b
  else:self.a,self.b=P(a),P(b)
 def __add__(self,x):x=Z(x);return Z(self.a+x.a,self.b+x.b)
 __radd__=__add__
 def __neg__(self):return Z(-self.a,-self.b)
 def __sub__(self,x):return self+-Z(x)
 def __rsub__(self,x):return Z(x)+-self
 def __mul__(self,x):x=Z(x);return Z(self.a*x.a+19*self.b*x.b,self.a*x.b+self.b*x.a)
 __rmul__=__mul__
 def __truediv__(self,x):
  x=Z(x);d=x.a*x.a-19*x.b*x.b;assert d;return self*Z(x.a/d,-x.b/d)
 def __rtruediv__(self,x):return Z(x)/self
 def __bool__(self):return bool(self.a or self.b)
 def __eq__(self,x):x=Z(x);return self.a==x.a and self.b==x.b
 def sign(self):
  a,b=self.a.sign(),self.b.sign()
  if not b:return a
  if not a:return b
  if a==b:return a
  return a*(self.a*self.a-19*self.b*self.b).sign()
 def enc(self):return self.a.enc()+self.b.enc()
def enc(x):
 if isinstance(x,Z):return x.enc()
 if isinstance(x,P):return Z(x).enc()
 if isinstance(x,K):return Z(P(x)).enc()
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def dec4(x):return Z(P(K(F(x[0]),F(x[1])),K(F(x[2]),F(x[3]))))
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Z())
def mv(M,v):return [dot(row,v)for row in M]
def inverse(M):
 n=len(M);a=[[Z(x)for x in r]+[Z(i==j)for j in range(n)]for i,r in enumerate(M)]
 for k in range(n):
  j=next(j for j in range(k,n)if a[j][k]);a[k],a[j]=a[j],a[k];v=a[k][k];a[k]=[x/v for x in a[k]]
  for j in range(n):
   if j!=k and a[j][k]:
    v=a[j][k];a[j]=[x-v*y for x,y in zip(a[j],a[k])]
 return [r[n:]for r in a]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 Q=[[dec4(x)for x in row]for row in d['form']];N=[[dec4(x)for x in row]for row in d['normals']]
 assert len(N)==43
 C=list(map(dec4,d['strict_center']));old_cov=[mv(Q,n)for n in N]
 shared=[1,2,3,4,5,7,8,9,10,11,13,14,15,16,18,20,21];anchors=[12,37]
 rows=[old_cov[i][:19]for i in shared+anchors]
 Y=mv(inverse(rows),[Z()]*17+[Z(-1),Z(-1)])
 Q0=[row[:19]for row in Q[:19]]
 assert dot(Y,mv(Q0,Y))==F(-15,2) and dot(Y,mv(Q0,C[:19])).sign()<0
 R=Y+[Z(0,F(-1,2)),Z()];N.append(R);cov=[mv(Q,n)for n in N]
 G=[[dot(n,c)for c in cov]for n in N];norm=[G[i][i]for i in range(44)];assert norm[-1]==2
 counts=Counter();all_counts=Counter();labels={}
 for i,j in combinations(range(44),2):
  g=G[i][j];assert g.sign()<=0
  t=4*g*g/(norm[i]*norm[j])
  if not g:label='right'
  elif t==1:label='pi3'
  elif t==2:label='pi4'
  elif t==3:label='pi6'
  elif t==4:label='parallel'
  else:assert(t-4).sign()>0;label='dotted'
  labels[i,j]=label;all_counts[label]+=1
  if j==43:counts[label]+=1
 assert counts==dict(right=18,pi3=2,parallel=1,dotted=22)
 cp=[dot(C,c)for c in cov];qC=dot(C,mv(Q,C));assert qC.sign()<0 and all(x.sign()<0 for x in cp)
 facets=[]
 for i,n in enumerate(N):
  p=[a-cp[i]/norm[i]*b for a,b in zip(C,n)]
  assert dot(p,mv(Q,p)).sign()<0
  assert [dot(p,c).sign()for c in cov]==[-1]*i+[0]+[-1]*(43-i)
  facets.append(p)
 weight={k:F({'right':0,'pi3':1,'pi4':2,'pi6':3}[v],4)for k,v in labels.items()if v in('right','pi3','pi4','pi6')}
 neighbors=[i for i in range(43)if(i,43)in weight]
 def pd(S):
  adj={i:{}for i in S};diag={i:F(1)for i in S}
  for i,j in combinations(S,2):
   p=tuple(sorted((i,j)))
   if p not in weight:return False
   if weight[p]:adj[i][j]=weight[p];adj[j][i]=weight[p]
  while adj:
   leaves=[i for i in adj if len(adj[i])<=1]
   if not leaves:return False
   i=min(leaves)
   if diag[i]<=0:return False
   if adj[i]:
    j,w=next(iter(adj[i].items()));diag[j]-=w/diag[i];del adj[j][i]
   del adj[i]
  return True
 candidates=[list(s)+[43]for s in combinations(neighbors,19)];positive=[s for s in candidates if pd(s)]
 assert len(neighbors)==20 and len(candidates)==20
 assert positive==[[1,2,3,4,5,7,8,9,10,11,13,14,15,16,18,20,21,37,41,43]]
 S=positive[0];S0=sorted(S[:-1]+[40]);assert pd(S0)
 def vertex(s):
  v=mv(inverse([cov[i]for i in s]+[mv(Q,C)]),[Z()]*20+[Z(-1)])
  q=dot(v,mv(Q,v));p=[dot(v,c)for c in cov]
  assert q.sign()<0 and dot(v,mv(Q,C))==-1
  assert all(g.sign()==(0 if i in s else -1)for i,g in enumerate(p))
  return dict(support=s,point=v,norm=q,pairings=p)
 endpoints=[vertex(S0),vertex(S)];edge=sorted(set(S0)&set(S));assert len(edge)==19 and pd(edge)
 assert endpoints[0]['pairings'][43].sign()<0 and endpoints[1]['pairings'][40].sign()<0
 retained=[]
 for row in d['finite_edge']['endpoints']:
  v=list(map(dec4,row['point']));p=[dot(v,c)for c in cov]
  assert dot(v,mv(Q,v)).sign()<0 and all(x.sign()==(0 if i in row['support']else -1)for i,x in enumerate(p))
  retained.append(dict(support=row['support'],point=v,pairings=p))
 assert len({tuple(sorted(r['support']))for r in endpoints+retained})==4
 oldE=list(map(dec4,d['strict_escape']));cut=dot(oldE,cov[-1]);assert cut.sign()>0
 E=oldE[:19]+[Z(F(-7,10)),Z(F(-1,10))];qE=dot(E,mv(Q,E));ep=[dot(E,c)for c in cov]
 assert qE.sign()>0 and all(p.sign()<0 for p in ep)and dot(E,mv(Q,C)).sign()<0
 face=list(map(dec4,d['full_H18_source_face']['point']));assert [dot(face,c).sign()for c in cov]==[-1]*40+[0,0,-1,-1]
 cycle=4*G[0][36]*G[0][36]/(norm[0]*norm[36]);minor=norm[38]*norm[39]-G[38][39]*G[38][39]
 assert cycle==Z(P(K(12,8)))and minor==-32
 paths=[f'scripts/{TAG}.py',SOURCE,HELPER]
 return enc(dict(status='PASS_EXACT44_WALL_NEW_PROJECTION_FINITE_EDGE_AND_STRICT_ESCAPE',sources={p:sha(p)for p in paths},
  field_basis=['1','sqrt2','sqrt3','sqrt6','sqrt19','sqrt38','sqrt57','sqrt114'],form=Q,normals=N,gram=G,
  cap_recipe=dict(source_rights=shared,source_anchors=anchors,anchor_pairings=[-1,-1],source_projection=Y,source_projection_norm=F(-15,2),transverse_coordinates=[Z(0,F(-1,2)),Z()]),
  dimension=20,facets=44,all_pair_counts=dict(all_counts),new_pair_counts=dict(counts),strict_center=C,center_norm=qC,center_pairings=cp,facet_points=facets,
  complete_new_cap_vertex_census=dict(finite_neighbors=neighbors,tested_supports=len(candidates),positive_supports=positive,new_vertices=1),
  finite_edge=dict(shared_support=edge,endpoints=endpoints,old_released_wall=40,new_cap=43),retained_previous_edge_vertices=retained,
  full_H18_source_face=dict(point=face,new_projection_future_timelike=True),old_escape_cut=cut,strict_escape=E,escape_norm=qE,escape_pairings=ep,escape_center_pairing=dot(E,mv(Q,C)),
  conditional_nonquasi=dict(witness_pair=[0,36],cartan_two_cycle=cycle,rational_indefinite_pair=[38,39],determinant=minor),
  scope='One exact noninducing cap with a different source projection retains all43 old walls, cuts their saved escape, and creates exactly one finite vertex on the new cap plus an actual finite edge. Four distinct finite vertices are explicitly checked, not a complete ambient census. The complete H18 face is unchanged. A renewed exact strict spacelike escape proves infinite volume; no arithmetic lattice or finite-volume claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],d['escape_norm'])
