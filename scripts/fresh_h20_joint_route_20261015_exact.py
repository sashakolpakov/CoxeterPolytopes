#!/usr/bin/env python3
"""Exact noninducing A1^2/H18 cap, complete new-cap vertex census and escape."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h20_b2_h18_coupling_20261014_review as A
K=A.K;ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_joint_route_20261015_exact'
SOURCE='results/fresh_h20_a1a1_h18_coupling_20261015.json';HELPER='scripts/fresh_h20_b2_h18_coupling_20261014_review.py'
class P:
 def __init__(self,a=0,b=0):
  if isinstance(a,P):self.a,self.b=a.a,a.b
  else:self.a,self.b=K(a),K(b)
 def __add__(self,x):x=P(x);return P(self.a+x.a,self.b+x.b)
 __radd__=__add__
 def __neg__(self):return P(-self.a,-self.b)
 def __sub__(self,x):return self+-P(x)
 def __rsub__(self,x):return P(x)+-self
 def __mul__(self,x):x=P(x);return P(self.a*x.a+3*self.b*x.b,self.a*x.b+self.b*x.a)
 __rmul__=__mul__
 def __truediv__(self,x):
  x=P(x);d=x.a*x.a-3*x.b*x.b;assert d;return self*P(x.a/d,-x.b/d)
 def __rtruediv__(self,x):return P(x)/self
 def __bool__(self):return bool(self.a or self.b)
 def __eq__(self,x):x=P(x);return self.a==x.a and self.b==x.b
 def sign(self):
  a,b=self.a.sign(),self.b.sign()
  if not b:return a
  if not a:return b
  if a==b:return a
  return a*(self.a*self.a-3*self.b*self.b).sign()
 def enc(self):return list(map(str,[self.a.a,self.a.b,self.b.a,self.b.b]))
def enc(x):
 if isinstance(x,P):return x.enc()
 if isinstance(x,K):return P(x).enc()
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def dot(x,y):return sum((a*b for a,b in zip(x,y)),P())
def mv(M,v):return [dot(r,v)for r in M]
def inverse(M):
 n=len(M);a=[[P(x)for x in r]+[P(i==j)for j in range(n)]for i,r in enumerate(M)]
 for k in range(n):
  j=next(j for j in range(k,n)if a[j][k]);a[k],a[j]=a[j],a[k];q=a[k][k];a[k]=[v/q for v in a[k]]
  for j in range(n):
   if j!=k:
    q=a[j][k]
    if q:a[j]=[v-q*w for v,w in zip(a[j],a[k])]
 return [r[n:]for r in a]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 Q0=[[A.dec(x)for x in row]for row in d['source_form']];old0=[[A.dec(x)for x in row]for row in d['source_normals']]
 H0=list(map(A.dec,d['source_center']));assert A.inertia(Q0)==(18,1,0)
 shared=[1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,21,36];anchors=[12,18]
 oldcov=[A.mv(Q0,n)for n in old0]
 Y=A.mv(A.inverse([oldcov[i]for i in shared+anchors]),[K()]*17+[K(-1),K(-1)])
 assert A.dot(Y,A.mv(Q0,Y))==K(F(-3,2)) and A.dot(Y,A.mv(Q0,H0)).sign()<0
 Q=[[P(x)for x in row]+[P(),P()]for row in Q0]+[[P()]*19+[P(2),P()],[P()]*19+[P(),P(2)]]
 N=[[P(x)for x in row]+[P(),P()]for row in old0]+[[P()]*19+[P(1),P()],[P()]*19+[P(),P(1)]]
 R=[P(x)for x in Y]+[P(0,F(-1,2)),P(-1)];N.append(R)
 assert len(N)==43
 cov=[mv(Q,n)for n in N];G=[[dot(n,c)for c in cov]for n in N];norm=[G[i][i]for i in range(43)]
 assert norm[42]==2
 kinds={};counts=Counter();all_counts=Counter()
 for i,j in combinations(range(43),2):
  g=G[i][j];assert g.sign()<=0
  t=4*g*g/(norm[i]*norm[j])
  if not g:label='right'
  elif t==1:label='pi3'
  elif t==2:label='pi4'
  elif t==3:label='pi6'
  elif t==4:label='parallel'
  else:assert(t-4).sign()>0;label='dotted'
  kinds[i,j]=label;all_counts[label]+=1
  if j==42:counts[label]+=1
 assert counts==dict(right=18,pi3=3,pi6=1,parallel=1,dotted=19)
 assert [i for i in range(42)if not G[i][42]]==shared[:15]+[20,21,36]
 C=[P(x)for x in map(A.dec,d['strict_center'])]
 qC=dot(C,mv(Q,C));cp=[dot(C,c)for c in cov]
 assert qC.sign()<0 and all(p.sign()<0 for p in cp)
 facets=[]
 for i,n in enumerate(N):
  v=[a-cp[i]/norm[i]*b for a,b in zip(C,n)]
  assert dot(v,mv(Q,v)).sign()<0
  signs=[dot(v,c).sign()for c in cov];assert signs==[-1]*i+[0]+[-1]*(42-i)
  facets.append(v)
 # Exact finite-neighbor and forest-Schur census; no floating eigenvalue gates.
 finite=lambda i,j:kinds[min(i,j),max(i,j)]in('right','pi3','pi4','pi6')
 neighbors=[i for i in range(42)if finite(i,42)]
 assert len(neighbors)==22
 weight={pair:F({'right':0,'pi3':1,'pi4':2,'pi6':3}[label],4)for pair,label in kinds.items()if label in('right','pi3','pi4','pi6')}
 def pd(ids):
  adjacency={i:{}for i in ids}
  for i,j in combinations(ids,2):
   pair=min(i,j),max(i,j)
   if pair not in weight:return False
   if weight[pair]:adjacency[i][j]=weight[pair];adjacency[j][i]=weight[pair]
  diag={i:F(1)for i in ids}
  while adjacency:
   choices=[i for i,row in adjacency.items()if len(row)<=1]
   if not choices:return False # nonright cycle has a nonpositive quadratic witness
   i=min(choices)
   if diag[i]<=0:return False
   if adjacency[i]:
    j,w=next(iter(adjacency[i].items()));diag[j]-=w/diag[i];del adjacency[j][i]
   del adjacency[i]
  return True
 positive=[];tested=0
 for S in combinations(neighbors,19):
  tested+=1
  if pd(list(S)+[42]):positive.append(list(S)+[42])
 assert tested==1540 and positive==[[1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,20,21,36,40,42]]
 S=positive[0];S0=S[:-1]+[41]
 assert pd(S0)
 def vertex(support):
  # The last equation fixes a future scale relative to the common center.
  equations=[cov[i]for i in support]+[mv(Q,C)]
  v=mv(inverse(equations),[P()]*20+[P(-1)])
  assert dot(v,mv(Q,v)).sign()<0 and dot(v,mv(Q,C))==-1
  pairings=[dot(v,c)for c in cov]
  assert all(p.sign()==(0 if i in support else -1)for i,p in enumerate(pairings))
  return dict(support=support,point=v,norm=dot(v,mv(Q,v)),pairings=pairings)
 vertices=[vertex(S0),vertex(S)]
 edge=sorted(set(S)&set(S0));assert len(edge)==19 and pd(edge)
 assert vertices[0]['pairings'][42].sign()<0 and vertices[1]['pairings'][41].sign()<0
 Eold=[P(x)for x in map(A.dec,d['strict_escape'])]
 assert dot(Eold,cov[42]).sign()>0
 E=Eold[:19]+[P(F(-3,4)),P(F(-1,1000))]
 qE=dot(E,mv(Q,E));ep=[dot(E,c)for c in cov]
 assert qE.sign()>0 and all(p.sign()<0 for p in ep) and dot(E,mv(Q,C)).sign()<0
 source_face=[P(x)for x in H0]+[P(),P()]
 assert [dot(source_face,c).sign()for c in cov]==[-1]*40+[0,0,-1]
 # Conditional non-quasi-arithmetic witness is unchanged in the retained source.
 cycle=4*G[0][36]*G[0][36]/(norm[0]*norm[36]);minor=norm[38]*norm[39]-G[38][39]*G[38][39]
 assert cycle==P(K(12,8)) and minor==-32
 paths=[f'scripts/{TAG}.py',SOURCE,HELPER]
 return enc(dict(status='PASS_EXACT_NONINDUCING_H20_CAP_AND_FINITE_EDGE_WITH_RENEWED_ESCAPE',sources={p:sha(p)for p in paths},
  field_basis=['1','sqrt2','sqrt3','sqrt6'],form=Q,normals=N,gram=G,
  cap_recipe=dict(source_rights=shared,source_anchors=anchors,anchor_pairings=[-1,-1],source_projection=Y,source_projection_norm=F(-3,2),transverse_coordinates=[P(0,F(-1,2)),P(-1)]),
  dimension=20,facets=43,all_pair_counts=dict(all_counts),new_pair_counts=dict(counts),strict_center=C,center_norm=qC,center_pairings=cp,facet_points=facets,
  full_H18_source_face=dict(point=source_face,reason='The new source projection is strictly future timelike, so its inequality is strict on every future timelike or null source ray.'),
  complete_new_cap_vertex_census=dict(finite_neighbors=neighbors,tested_supports=tested,positive_supports=positive,new_vertices=1),
  finite_edge=dict(shared_support=edge,endpoints=vertices,old_released_wall=41,new_cap=42),
  old_escape_cut=dot(Eold,cov[42]),strict_escape=E,escape_norm=qE,escape_pairings=ep,escape_center_pairing=dot(E,mv(Q,C)),
  conditional_nonquasi=dict(witness_pair=[0,36],cartan_two_cycle=cycle,rational_indefinite_pair=[38,39],determinant=minor),
  scope='Exact 43-wall geometric branch retaining the complete 40-wall H18 source, with a conditional nonquasi-arithmetic completion witness. One noninducing cap cuts the old direct-coupling escape and closes a finite edge. Complete census only for vertices on the new cap; no complete ambient vertex census. A new exact strict spacelike escape proves this chamber still has infinite volume. The current infinite-volume chamber is not asserted to be a lattice. No arithmetic root search or finite-volume claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],d['escape_norm'])
