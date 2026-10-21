#!/usr/bin/env python3
"""Exact coherent quartet from four forward endpoint recipes; no scout import."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_disjoint_walk_20261014_exact'
RAD=None
class L:
 __slots__=('a','b')
 def __init__(self,a=0,b=0):
  if isinstance(a,L):self.a,self.b=a.a,a.b
  else:self.a,self.b=A.Q(a),A.Q(b)
 def __add__(self,y):
  y=L(y);return L(self.a+y.a,self.b+y.b)
 __radd__=__add__
 def __neg__(self):return L(-self.a,-self.b)
 def __sub__(self,y):return self+-L(y)
 def __rsub__(self,y):return L(y)+-self
 def __mul__(self,y):
  y=L(y);return L(self.a*y.a+RAD*self.b*y.b,self.a*y.b+self.b*y.a)
 __rmul__=__mul__
 def __truediv__(self,y):
  y=L(y);d=y.a*y.a-RAD*y.b*y.b;assert d
  return self*L(y.a/d,-y.b/d)
 def __rtruediv__(self,y):return L(y)/self
 def __bool__(self):return bool(self.a or self.b)
 def __eq__(self,y):
  y=L(y);return self.a==y.a and self.b==y.b
 def __hash__(self):return hash((self.a,self.b))
 def sign(self):
  a,b=self.a.sign(),self.b.sign()
  if not a:return b
  if not b or a==b:return a
  return a*(self.a*self.a-RAD*self.b*self.b).sign()
 def enc(self):return [self.a.enc(),self.b.enc()]
Z=L();ONE=L(1);PHI=A.Q(F(1,2),F(1,2))
def enc(x):
 if isinstance(x,(L,A.Q)):return x.enc()
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(y)for y in x]
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 return x
def old(x):
 a=A.dec(x);assert not a.b;return L(a.a)
def dot(x,y):return sum((a*b for a,b in zip(x,y)),Z)
def mv(M,v):return[dot(r,v)for r in M]
def pair(x,y,Q):return dot(x,mv(Q,y))
def inv(M):
 n=len(M);B=[r[:]+[L(i==j)for j in range(n)]for i,r in enumerate(M)]
 for j in range(n):
  p=next(i for i in range(j,n)if B[i][j]);B[p],B[j]=B[j],B[p]
  d=B[j][j];B[j]=[x/d for x in B[j]]
  for i in range(n):
   if i!=j:
    c=B[i][j];B[i]=[a-c*b for a,b in zip(B[i],B[j])]
 return[r[n:]for r in B]
def pd(M):
 B=[r[:]for r in M]
 for i in range(len(B)):
  d=B[i][i]
  if d.sign()<=0:return False
  for j in range(i+1,len(B)):
   for k in range(i+1,len(B)):B[j][k]-=B[j][i]*B[i][k]/d
 return True
def relation(g,q,r):
 if not g:return 'right',2
 if g.sign()>0:return 'obtuse',None
 c=g*g/(q*r)
 if(c-1).sign()>=0:return('parallel'if c==1 else'dotted'),None
 for m,v in {3:L(F(1,4)),4:L(F(1,2)),5:L(A.Q(F(3,8),F(1,8))),6:L(F(3,4))}.items():
  if c==v:return'finite',m
 return'forbidden',None
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def cliques(adj,k):
 def rec(prefix,cands):
  if len(prefix)==k:yield prefix;return
  for n,i in enumerate(cands):
   tail=[j for j in cands[n+1:]if j in adj[i]]
   if len(prefix)+1+len(tail)>=k:yield from rec(prefix+[i],tail)
 yield from rec([],list(range(len(adj))))
def run():
 global RAD
 base='fresh_h9_h8_disjoint_cap_20261014';front='fresh_h9_h8_disjoint_walk_20261014'
 d=json.loads((ROOT/f'results/{base}.json').read_text());f=json.loads((ROOT/f'results/{front}.json').read_text())
 for doc in[d,f]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 # Solve old spatial coordinates from nine prescribed source pairings.
 Q0=A.dm(d['form']);N0=A.dm(d['normals']);Q8=[r[:9]for r in Q0[:9]]
 recipes=[(13,[(9,3)],2,5),(27,[(0,3),(10,3)],1,2),(9,[(10,5)],3,5),(19,[(9,5)],2,3)]
 parts=[];records=[]
 for eid,contacts,anchor,m in recipes:
  S=f['queue'][eid]['support'];ids=S+[anchor];assert len(set(ids))==9 and max(ids)<11
  targets={i:A.Z for i in S};targets.update({i:-A.E(1 if m0==3 else PHI)for i,m0 in contacts})
  targets[anchor]={2:A.Z,3:A.E(-1),5:-A.E(PHI)}[m]
  rows=[A.mv(Q8,N0[i][:9])for i in ids];p=A.mv(A.inverse(rows),[targets[i]for i in ids])
  assert all(A.pair(p,N0[i][:9],Q8)==targets[i]for i in ids)
  rad=(A.E(2)-A.pair(p,p,Q8))/Q0[9][9];assert not rad.b and rad.sign()>0
  if RAD is None:RAD=rad.a
  else:assert RAD==rad.a
  parts.append([L(x.a)for x in p]);records.append(dict(edge_id=eid,contacts=contacts,anchor=anchor,anchor_label=m,spatial_part=parts[-1]))
 assert RAD==A.Q(6,F(16,5)) and RAD.sign()>0
 assert RAD.a*RAD.a-5*RAD.b*RAD.b==F(-76,5) # not a square in Q(sqrt5)
 Q=[[old(x)for x in row]for row in d['form']];N=[[old(x)for x in row]for row in d['normals']]
 N += [p+[L(0,-1)]for p in parts]
 H=[old(x)for x in d['center']];E=[old(x)for x in d['strict_escape']['point']]
 norms=[pair(n,n,Q)for n in N];assert norms[15:]==[L(2)]*4
 G=[[pair(n,m,Q)for m in N]for n in N];pairs=[];counts=Counter();adj=[set()for _ in N]
 for i,j in combinations(range(19),2):
  typ,m=relation(G[i][j],norms[i],norms[j]);assert typ in['right','finite','dotted']
  if typ!='dotted':adj[i].add(j);adj[j].add(i)
  if j>=15:pairs.append(dict(pair=[i,j],type=typ,label=m,product=G[i][j]));counts['pi'+str(m)if typ=='finite'else typ]+=1
 assert len(pairs)==66
 assert dict(counts)==dict(right=31,pi3=6,pi5=8,dotted=21)
 assert [[G[i][j]for j in range(15,19)]for i in range(15,19)]==[[L(2),L(-2*PHI),L(-PHI),Z],[L(-2*PHI),L(2),Z,L(-PHI)],[L(-PHI),Z,L(2),Z],[Z,L(-PHI),Z,L(2)]]
 assert pair(H,H,Q).sign()<0 and all(pair(H,n,Q).sign()<0 for n in N)
 assert pair(E,E,Q).sign()>0 and pair(E,H,Q).sign()<0 and all(pair(E,n,Q).sign()<0 for n in N)
 facets=[]
 for i,n in enumerate(N):
  c=pair(H,n,Q)/norms[i];p=[h-c*x for h,x in zip(H,n)];vals=[pair(p,x,Q)for x in N]
  assert pair(p,p,Q).sign()<0 and not vals[i]and all(x.sign()<0 for j,x in enumerate(vals)if j!=i)
  facets.append(dict(facet=i,point=p,norm=pair(p,p,Q),pairings=vals))
 vertices=[];trials=bad=0
 for ids in cliques(adj,9):
  trials+=1;M=[[G[i][j]for j in ids]for i in ids]
  assert all(not x.b for row in M for x in row)
  if not pd([[x.a for x in row]for row in M]):bad+=1;continue
  J=inv(M);coef=mv(J,[pair(H,N[i],Q)for i in ids])
  p=[H[b]-sum((c*N[i][b]for c,i in zip(coef,ids)),Z)for b in range(10)]
  vals=[pair(p,n,Q)for n in N]
  assert pair(p,p,Q).sign()<0 and pair(p,H,Q).sign()<0
  assert all(x.sign()==(0 if j in ids else-1)for j,x in enumerate(vals))
  vertices.append(dict(incident=ids,point=p,norm=pair(p,p,Q),pairings=vals))
 oldsets={tuple(sorted(r['incident']))for r in d['complete_finite_vertices']}
 assert oldsets<={tuple(r['incident'])for r in vertices}
 new=[i for i,r in enumerate(vertices)if tuple(r['incident'])not in oldsets]
 edges=[]
 for i,j in combinations(range(len(vertices)),2):
  S=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
  if len(S)==8:edges.append(dict(vertices=[i,j],support=S))
 closed=[]
 for row in f['queue']:
  S=set(row['support']);hits=[i for i in new if S<=set(vertices[i]['incident'])]
  assert len(hits)<=1
  if not hits:continue
  vi=hits[0];p=vertices[vi]['point'];v=[old(x)for x in row['start_point']];D=[old(x)for x in row['direction']]
  alpha=pair(p,v,Q)/pair(v,v,Q);beta=pair(p,D,Q)/pair(D,D,Q)
  assert alpha.sign()>0 and beta.sign()>0
  assert all(x==alpha*y+beta*z for x,y,z in zip(p,v,D))
  closed.append(dict(germ=row['id'],start_vertex=row['start_vertex'],new_vertex=vi,forward_parameter=beta/alpha,support=row['support']))
 assert len(vertices)==53 and len(new)==9 and len(edges)==175
 assert [r['germ']for r in closed]==[9,13,15,19,27,34]
 paths=[f'results/{base}.json',f'results/{front}.json',f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
 return enc(dict(status='PASS_EXACT_COHERENT_FOUR_CAP_ENDPOINT_EXTENSION',sources={p:sha(p)for p in paths},
  field=dict(base='Q(sqrt5)',new_positive_radical_square=RAD,encoding='[[a,b],[c,d]] means a+b sqrt5+(c+d sqrt5)sqrt(RAD)'),
  recipes=records,form=Q,normals=N,norms=norms,gram=G,new_pairs=pairs,new_pair_counts=dict(counts),
  center=H,facet_points=facets,finite_neighbors={i:sorted(a)for i,a in enumerate(adj)},
  complete_finite_vertex_census=dict(candidate_cliques=trials,nonpositive=bad,vertices=vertices,count=len(vertices),new_indices=new),
  complete_finite_edges=edges,closed_old_germs=closed,strict_escape=dict(point=E,norm=pair(E,E,Q),pairings=[pair(E,n,Q)for n in N]),
  scope='All19facets actual. Complete finite vertices/finite edges via all9-cliques of complete finite-neighbor graph, with direct physical signs. Listed old62germs closed at distinct future finite endpoints. Existing escape survives, so global volume remains infinite; no full exterior/ideal frontier census.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['complete_finite_vertex_census']['count'],'vertices',len(out['complete_finite_edges']),'finite edges',len(out['closed_old_germs']),'closed germs',flush=True)
