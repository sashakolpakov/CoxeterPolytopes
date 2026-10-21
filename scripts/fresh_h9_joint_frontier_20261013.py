#!/usr/bin/env python3
"""Exact finite-corner and edge census of the frozen 27-wall Makarov cone.

Uses the separate Fraction Q(sqrt(5)) arithmetic, then sparse products of
positive square roots of that field. No numerical incidence is authoritative.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
from math import isqrt
import argparse, hashlib, json
import fresh_h9_makarov_joint_20261013_exact_independent as I

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_joint_frontier_20261013'
SOURCE=ROOT/'results/fresh_h9_makarov_joint_20261013_exact.json'
SCOUT=ROOT/'results/fresh_h9_makarov_joint_20261013_new_frontier_scout.json'
OUT=ROOT/'results'/f'{TAG}.json'
K=I.k; Z=I.Z; O=K(1); PHI=K(F(1,2),F(1,2)); DS=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sqrtf(x):
 if x<0:return None
 a,b=isqrt(x.numerator),isqrt(x.denominator)
 return F(a,b) if a*a==x.numerator and b*b==x.denominator else None
def sqrtk(x):
 n=sqrtf(x[0]*x[0]-5*x[1]*x[1])
 if n is not None:
  for s in [n,-n]:
   a=sqrtf((x[0]+s)/2)
   if a:
    y=(a,x[1]/(2*a))
    if I.mul(y,y)==x:return y if I.sign(y)>0 else I.neg(y)
 if not x[1]:
  b=sqrtf(x[0]/5)
  if b is not None:return K(0,b)
 return None
def C(x):return {} if x==Z else {0:x}
def add(a,b):
 o=dict(a)
 for i,x in b.items():
  y=I.add(o.get(i,Z),x)
  if y!=Z:o[i]=y
  elif i in o:del o[i]
 return o
def neg(a):return {i:I.neg(x)for i,x in a.items()}
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 o={}
 for i,x in a.items():
  for j,y in b.items():
   z=I.mul(x,y);bits=i&j
   while bits:
    bit=bits&-bits;z=I.mul(z,DS[bit.bit_length()-1]);bits-=bit
   key=i^j;z=I.add(o.get(key,Z),z)
   if z!=Z:o[key]=z
   elif key in o:del o[key]
 return o
def scale(q,a):return {i:I.mul(q,x)for i,x in a.items()if I.mul(q,x)!=Z}
def total(xs):
 o={}
 for x in xs:o=add(o,x)
 return o
@lru_cache(maxsize=100000)
def sgn_cached(items):
 a=dict(items)
 if not a:return 0
 bit=1<<(max(a).bit_length()-1) if max(a) else 0
 if not bit:return I.sign(a[0])
 low={i:x for i,x in a.items()if not i&bit};high={i^bit:x for i,x in a.items()if i&bit}
 sa,sb=sgn(low),sgn(high)
 if not sa:return sb
 if not sb:return sa
 if sa==sb:return sa
 return sa*sgn(sub(mul(low,low),scale(DS[bit.bit_length()-1],mul(high,high))))
def sgn(a):return sgn_cached(tuple(sorted(a.items())))
def root(d):
 assert I.sign(d)>0
 sq=sqrtk(d)
 if sq is not None:return C(sq)
 for i,e in enumerate(DS):
  sq=sqrtk(I.mul(d,I.reciprocal(e)))
  if sq is not None:return {1<<i:sq}
 DS.append(d);return {1<<(len(DS)-1):O}
def enc(a):return [[i,list(map(str,x))]for i,x in sorted(a.items())]
def encv(v):return list(map(enc,v))
def encm(M):return [[list(map(str,x))for x in row]for row in M]
def vadd(a,b):return list(map(add,a,b))
def vs(q,v):return [mul(q,x)for x in v]
def dot(a,b):return total(mul(x,y)for x,y in zip(a,b,strict=True))
def cov(v,Q):return [total(scale(q,x)for q,x in zip(row,v,strict=True))for row in Q]
def pair(a,b,Q):return dot(a,cov(b,Q))
def ki(M):
 n=len(M);A=[list(r)+[O if i==j else Z for j in range(n)]for i,r in enumerate(M)]
 for j in range(n):
  p=next(i for i in range(j,n)if A[i][j]!=Z);A[p],A[j]=A[j],A[p]
  a=I.reciprocal(A[j][j]);A[j]=[I.mul(a,x)for x in A[j]]
  for i in range(n):
   if i!=j:
    a=A[i][j];A[i]=[I.sub(x,I.mul(a,y))for x,y in zip(A[i],A[j],strict=True)]
 return [r[n:]for r in A]
def pd(M):
 A=[list(r)for r in M]
 while A:
  a=A[0][0]
  if I.sign(a)<=0:return False
  a=I.reciprocal(a)
  A=[[I.sub(A[i][j],I.mul(I.mul(A[i][0],A[0][j]),a))for j in range(1,len(A))]for i in range(1,len(A))]
 return True
def clique_list(adj,size):
 out=[]
 def visit(v,rest):
  if len(v)==size:out.append(tuple(v));return
  if len(v)+len(rest)<size:return
  for i in sorted(rest):visit(v+[i],{j for j in rest&adj[i]if j>i})
 visit([],set(range(len(adj))));return out

def spherical_profiles(M,ids):
 B=ki(M);left=set(range(len(M)));components=[]
 while left:
  seen={min(left)};queue=list(seen)
  for i in queue:
   for j in sorted(left-seen):
    if M[i][j]!=Z:seen.add(j);queue.append(j)
  components.append(sorted(seen));left-=seen
 which={i:k for k,row in enumerate(components)for i in row}
 rows=[{'contacts':[],'projection_square':list(map(str,Z)),'residual_norm':list(map(str,K(2)))}];tails=[]
 squares={3:O,4:K(2),5:K(F(3,2),F(1,2))}
 def append(inds,labels):
  p=sumk(I.mul(squares[m],B[i][i])for i,m in zip(inds,labels,strict=True));q=I.sub(K(2),p)
  if I.sign(q)>0:rows.append({'contacts':[[ids[i],m]for i,m in zip(inds,labels,strict=True)],'projection_square':list(map(str,p)),'residual_norm':list(map(str,q))})
 for i in range(len(M)):
  if len(components[which[i]])==1:
   assert B[i][i]==K(F(1,2));tails.append({'wall':ids[i],'labels':'every integer m >= 3','residual_norm':'2 sin(pi/m)^2'});continue
  for m in squares:append((i,),(m,))
 for n in (2,3):
  for inds in combinations(range(len(M)),n):
   if len({which[i]for i in inds})!=n:continue
   for labels in product(squares,repeat=n)if n==2 else [(3,3,3)]:append(inds,labels)
 return {'components':[[ids[i]for i in row]for row in components],'finite_profiles':rows,'dihedral_tail_profiles':tails}

def setup():
 DS.clear();sgn_cached.cache_clear()
 s=json.loads(SOURCE.read_text());Q=[[K(*x)for x in r]for r in s['form']]
 assert I.inertia(Q)==[9,1]
 W=[[C(K(*x))for x in r]for r in s['seed_normals']]
 for row in s['caps']:
  d=root(K(*row['radicand']));W.append([add(C(K(*a)),scale(K(*b),d))for a,b in row['normal']])
 # Normalize the isolated polar wall as well; every finite Gram entry is
 # then in K and equals 0,-1,-phi. The physical hyperplanes do not change.
 polar=K(F(7,11),F(-1,11));W[6]=vs(root(I.mul(K(2),I.reciprocal(polar))),W[6])
 c=[C(K(*x))for x in s['strict_center']];CV=[cov(w,Q)for w in W]
 assert len(W)==27 and sgn(pair(c,c,Q))<0
 for w,cw in zip(W,CV,strict=True):assert sgn(sub(dot(w,cw),C(K(2))))==0 and sgn(dot(c,cw))<0
 labels={tuple(row['pair']):row['type']for row in s['all_pair_gates']}
 adj=[set()for _ in W];G=[[None]*27 for _ in W];counts=Counter()
 for i in range(27):
  G[i][i]=K(2)
  for j in range(i):
   x=dot(W[i],CV[j]);tag=labels[j,i];counts[tag]+=1
   if tag=='dotted':assert sgn(add(x,C(K(2))))<0
   else:
    target={'m2':Z,'m3':K(-1),'m5':I.neg(PHI)}[tag]
    assert sgn(sub(x,C(target)))==0
    G[i][j]=G[j][i]=target;adj[i].add(j);adj[j].add(i)
 assert dict(counts)==s['pair_counts']
 return s,Q,W,CV,c,G,adj

def run():
 s,Q,W,CV,c,G,adj=setup();cliques=clique_list(adj,9);assert len(cliques)==351
 cp=[dot(c,x)for x in CV];vertices=[];rejected=Counter();proposals={tuple(v['active']):v['id']for v in json.loads(SCOUT.read_text())['vertices']if v['type']=='finite'}
 for ids in cliques:
  gram=[[G[i][j]for j in ids]for i in ids]
  if not pd(gram):rejected['not_positive_definite']+=1;continue
  B=ki(gram);coeff=[total(scale(q,cp[j])for q,j in zip(row,ids,strict=True))for row in B]
  v=[sub(c[a],total(mul(t,W[j][a])for t,j in zip(coeff,ids,strict=True)))for a in range(10)]
  signs=[sgn(dot(v,x))for x in CV]
  if max(signs)>0:rejected['infeasible']+=1;continue
  active=[i for i,x in enumerate(signs)if not x];assert active==list(ids)
  q=pair(v,v,Q);assert sgn(q)<0 and sgn(pair(v,c,Q))<0
  vertices.append({'id':len(vertices),'active':list(ids),'point':encv(v),'norm':enc(q),'all_wall_signs':signs,'numerical_proposal_id':proposals.get(ids),'_point':v,'_inverse':B})
 print('finite corners',len(vertices),'rejected',dict(rejected),flush=True)
 assert len(vertices)==65 and {tuple(v['active'])for v in vertices}==set(proposals)
 index={tuple(v['active']):v['id']for v in vertices};edges=[];seen_support={};profiles={}
 for vertex in vertices:
  ids=vertex['active'];v=vertex['_point'];B=vertex['_inverse']
  for local,released in enumerate(ids):
   support=tuple(i for i in ids if i!=released)
   if support in seen_support:
    edges[seen_support[support]]['finite_vertices'].append(vertex['id']);continue
   e={'id':len(edges),'support':list(support),'finite_vertices':[vertex['id']],'released_at_start':released}
   seen_support[support]=len(edges)
   gram=[[G[i][j]for j in support]for i in support];assert pd(gram);BI=ki(gram)
   extensions=[]
   for owner in range(27):
    if owner in support or any(owner not in adj[i]for i in support):continue
    b=[G[i][owner]for i in support];projection=sumk(I.mul(x,I.mul(BI[i][j],y))for i,x in enumerate(b)for j,y in enumerate(b))
    schur=I.sub(K(2),projection);ss=I.sign(schur)
    if ss>0:
     key=tuple(sorted(support+(owner,)))
     if key in index:extensions.append({'wall':owner,'kind':'finite','vertex':index[key]})
    elif ss==0:
     # Independent affine extension yields the null vector r-proj_U(r).
     aa=[sumk(I.mul(x,y)for x,y in zip(row,b,strict=True))for row in BI]
     z=[sub(W[owner][a],total(scale(t,W[j][a])for t,j in zip(aa,support,strict=True)))for a in range(10)]
     assert sgn(pair(z,z,Q))==0
     if not sgn(pair(z,c,Q)):raise AssertionError('zero affine residual')
     if sgn(pair(z,c,Q))>0:z=[neg(x)for x in z]
     signs=[sgn(dot(z,x))for x in CV]
     if max(signs)<=0:extensions.append({'wall':owner,'kind':'ideal','point':encv(z),'all_wall_signs':signs})
   assert any(x.get('vertex')==vertex['id']for x in extensions)
   e['all_feasible_finite_or_ideal_extensions']=extensions
   other=[x for x in extensions if x.get('vertex')!=vertex['id']]
   if other:
    assert len(other)==1;e['kind']=other[0]['kind']+'_endpoint'
   else:
    e['kind']='open_ideal_boundary'
    profiles[str(e['id'])]=spherical_profiles(gram,support)
    # A spacelike point beyond the light cone can be obtained without
    # algebraic division by using an exact rational step along the dual.
    dual=[total(scale(a,W[j][k])for a,j in zip(B[local],ids,strict=True))for k in range(10)]
    assert sgn(pair(v,dual,Q))==0 and sgn(pair(dual,dual,Q))>0
    # Numerical values propose a rational; every accepted fact is exact.
    tv=(-approx(pair(v,v,Q))/approx(pair(dual,dual,Q)))**.5
    found=None
    for den in [10**n for n in range(2,15)]:
     t=F(int(tv*den)+1,den)
     y=[sub(x,scale(K(t),d))for x,d in zip(v,dual,strict=True)]
     ps=[sgn(dot(y,x))for x in CV]
     if sgn(pair(y,y,Q))>0 and max(ps)<=0 and [i for i,x in enumerate(ps)if not x]==list(support):found=(t,y,ps);break
    assert found is not None,(e['id'],tv)
    t,y,ps=found
    e['spacelike_continuation']={'parameter':str(t),'point':encv(y),'norm':enc(pair(y,y,Q)),'all_wall_signs':ps,'center_pairing':enc(pair(y,c,Q))}
    assert sgn(pair(y,c,Q))<0
   edges.append(e)
 for e in edges:
  expected={x['vertex']for x in e['all_feasible_finite_or_ideal_extensions']if x['kind']=='finite'}
  assert set(e['finite_vertices'])==expected
 counts=Counter(e['kind']for e in edges)
 print('edges',len(edges),dict(counts),flush=True)
 vg=[set()for _ in vertices]
 for e in edges:
  if len(e['finite_vertices'])==2:
   a,b=e['finite_vertices'];vg[a].add(b);vg[b].add(a)
 reached={0};queue=[0]
 for i in queue:
  for j in sorted(vg[i]-reached):reached.add(j);queue.append(j)
 assert len(reached)==len(vertices)
 profilecounts={'finite_profiles':sum(len(x['finite_profiles'])for x in profiles.values()),'dihedral_tails':sum(len(x['dihedral_tail_profiles'])for x in profiles.values())}
 known={tuple(row['incident'])for row in s['new_finite_endpoints']}
 seed=json.loads(I.SEED.read_text());known|={tuple(row+list(range(7,11)))for row in seed['vertex_active_sets']}
 additional=[v['id']for v in vertices if tuple(v['active'])not in known];assert len(additional)==15
 print('profiles',profilecounts,flush=True)
 for v in vertices:v.pop('_point');v.pop('_inverse')
 return {'status':'PASS_COMPLETE_FINITE_CORNERS_AND_INCIDENT_EDGES','sources':{str(p.relative_to(ROOT)):sha(p)for p in [SOURCE,SCOUT,I.SEED,Path(__file__),Path(I.__file__)]},'square_root_basis_over_Qsqrt5':[list(map(str,d))for d in DS],'normalized_wall_norm':2,'normalized_normals':[encv(w)for w in W],'finite_pair_nine_cliques':len(cliques),'rejected_nine_cliques':dict(rejected),'finite_vertex_count':len(vertices),'finite_vertices':vertices,'additional_finite_vertex_ids':additional,'finite_edge_graph_connected':True,'edge_counts':dict(counts),'edges':edges,'open_edge_spherical_profiles':profiles,'profile_counts':profilecounts,'profile_completeness':'A finite rank-nine Coxeter diagram is a forest with valence at most three. The new node has at most one neighbor in each existing connected component. Valence-three links all have label three; links with label at least six occur only in an isolated I2(m) component. All remaining labels 3,4,5 and all strictly positive Schur complements are exhausted. These are necessary spherical local profiles, not globally realized new caps; affine new endpoints are outside this profile census.','scope':'All finite vertices and all edges incident to a finite vertex of the frozen 27-wall cone. A finite Coxeter vertex is simple; all 351 finite-pair nine-cliques are exhausted. Every outgoing rank-eight support is tested against every existing wall for a positive-definite or affine extension; open edges have exact feasible spacelike continuations. No census of edges with no finite endpoint or arbitrary new mirrors is claimed.'}
def sumk(xs):
 z=Z
 for x in xs:z=I.add(z,x)
 return z
def approx(x):
 from math import sqrt,prod
 return sum((float(a)+sqrt(5)*float(b))*prod(sqrt(float(DS[j][0])+sqrt(5)*float(DS[j][1]))for j in range(len(DS))if mask>>j&1)for mask,(a,b)in x.items())
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run()
 if a.verify:assert json.loads(OUT.read_text())==out
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],flush=True)
