#!/usr/bin/env python3
"""Producer-free replay: Fraction polynomial identities + rational intervals.

In contrast to the producer's recursive radical sign comparison, every
strict sign here has a rational interval certificate. Formal identities
are checked coefficient by coefficient; no radical independence is assumed.
"""
from fractions import Fraction as F
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from math import isqrt
import argparse, hashlib, json
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_joint_frontier_20261013_independent'
SOURCE=ROOT/'results/fresh_h9_joint_frontier_20261013.json'
WALLS=ROOT/'results/fresh_h9_makarov_joint_20261013_exact.json'
OUT=ROOT/'results'/f'{TAG}.json';Z=(F(0),F(0));ONE=(F(1),F(0));DS=[];MON={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ka(x,y):return x[0]+y[0],x[1]+y[1]
def kn(x):return -x[0],-x[1]
def km(x,y):return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def kr(x):
 d=x[0]*x[0]-5*x[1]*x[1];assert d
 return x[0]/d,-x[1]/d
def ks(x):
 a,b=x
 if not b:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if a*b>0:return (a>0)-(a<0)
 d=a*a-5*b*b;return ((a>0)-(a<0))*((d>0)-(d<0))
def K(x):return tuple(map(F,x))
def con(x):return {0:x}if x!=Z else {}
def dec(x):return {m:K(q)for m,q in x}
def add(a,b):
 out=dict(a)
 for m,x in b.items():out[m]=ka(out.get(m,Z),x)
 return {m:x for m,x in out.items()if x!=Z}
def neg(x):return {m:kn(a)for m,a in x.items()}
def sub(x,y):return add(x,neg(y))
def mul(a,b):
 out={}
 for m,x in a.items():
  for n,y in b.items():
   q=km(x,y)
   for j,d in enumerate(DS):
    if (m&n)>>j&1:q=km(q,d)
   out[m^n]=ka(out.get(m^n,Z),q)
 return {m:x for m,x in out.items()if x!=Z}
def scale(q,x):return {m:km(q,a)for m,a in x.items()if km(q,a)!=Z}
def sm(xs):
 a={}
 for x in xs:a=add(a,x)
 return a
def dot(a,b):return sm(mul(x,y)for x,y in zip(a,b,strict=True))
def cov(v,Q):return [sm(scale(q,x)for q,x in zip(row,v,strict=True))for row in Q]
def pair(a,b,Q):return dot(a,cov(b,Q))
DEN=2**160
def sqrt_bounds(a,b):
 assert a>0
 lo=isqrt(a.numerator*DEN*DEN//a.denominator)
 hi=isqrt(b.numerator*DEN*DEN//b.denominator)+1
 return F(lo,DEN),F(hi,DEN)
SQ5=sqrt_bounds(F(5),F(5))
def interval_k(q):
 a,b=q
 return (a+b*SQ5[0],a+b*SQ5[1])if b>=0 else(a+b*SQ5[1],a+b*SQ5[0])
def interval_mul(x,y):
 a=[u*v for u in x for v in y];return min(a),max(a)
def interval(x):
 a=b=F(0)
 for m,q in x.items():
  u,v=interval_mul(interval_k(q),MON[m]);a+=u;b+=v
 return a,b
def sign(x):
 if not x:return 0
 a,b=interval(x)
 assert a>0 or b<0,('unresolved strict interval',a,b)
 return 1 if a>0 else -1
def inverse(M):
 n=len(M);A=[list(row)+[ONE if i==j else Z for j in range(n)]for i,row in enumerate(M)]
 for j in range(n):
  i=next(i for i in range(j,n)if A[i][j]!=Z);A[i],A[j]=A[j],A[i]
  u=kr(A[j][j]);A[j]=[km(u,x)for x in A[j]]
  for i in range(n):
   if i!=j:
    u=A[i][j];A[i]=[ka(x,kn(km(u,y)))for x,y in zip(A[i],A[j],strict=True)]
 return [row[n:]for row in A]
def positive(M):
 A=[list(row)for row in M]
 while A:
  if ks(A[0][0])<=0:return False
  u=kr(A[0][0]);A=[[ka(A[i][j],kn(km(u,km(A[i][0],A[0][j]))))for j in range(1,len(A))]for i in range(1,len(A))]
 return True
def ksum(xs):
 z=Z
 for x in xs:z=ka(z,x)
 return z
def run():
 d=json.loads(SOURCE.read_text());s=json.loads(WALLS.read_text())
 for p,h in d['sources'].items():assert sha(ROOT/p)==h
 DS[:]=[K(x)for x in d['square_root_basis_over_Qsqrt5']];assert all(ks(x)>0 for x in DS)
 roots=[sqrt_bounds(*interval_k(q))for q in DS];MON.clear();MON[0]=(F(1),F(1))
 for m in range(1,2**len(DS)):
  bit=m&-m;MON[m]=interval_mul(MON[m^bit],roots[bit.bit_length()-1])
 Q=[[K(x)for x in row]for row in s['form']];W=[[dec(x)for x in row]for row in d['normalized_normals']]
 c=[con(K(x))for x in s['strict_center']];CV=[cov(w,Q)for w in W]
 # Reconstruct each physical source mirror, verifying every saved radical
 # root by positive sign and exact square rather than trusting its label.
 for i,row in enumerate(s['seed_normals']):
  old=[K(x)for x in row]
  if i!=6:assert W[i]==[con(x)for x in old]
  else:
   a=next(j for j,x in enumerate(old)if x!=Z);u=scale(kr(old[a]),W[i][a]);assert sign(u)>0
   assert W[i]==[scale(x,u)for x in old]
   q=(F(7,11),F(-1,11));assert mul(u,u)==con(km((F(2),F(0)),kr(q)))
 for i,row in enumerate(s['caps'],11):
  A=[K(x[0])for x in row['normal']];B=[K(x[1])for x in row['normal']]
  nz=[j for j,x in enumerate(B)if x!=Z]
  if not nz:assert W[i]==[con(x)for x in A];continue
  j=nz[0];u=scale(kr(B[j]),sub(W[i][j],con(A[j])))
  assert sign(u)>0 and mul(u,u)==con(K(row['radicand']))
  assert W[i]==[add(con(a),scale(b,u))for a,b in zip(A,B,strict=True)]
 two=(F(2),F(0));phi=(F(1,2),F(1,2));G=[[None]*27 for _ in W];adj=[set()for _ in W]
 for i in range(27):assert dot(W[i],CV[i])==con(two) and sign(dot(c,CV[i]))<0;G[i][i]=two
 for r in s['all_pair_gates']:
  i,j=r['pair'];p=dot(W[i],CV[j])
  if r['type']=='dotted':assert sign(add(p,con(two)))<0
  else:
   q={'m2':Z,'m3':kn(ONE),'m5':kn(phi)}[r['type']];assert p==con(q);G[i][j]=G[j][i]=q;adj[i].add(j);adj[j].add(i)
 assert sign(pair(c,c,Q))<0
 escape=[con((F(x),F(0)))for x in s['escape']['point']]
 assert sign(pair(escape,escape,Q))>0 and sign(pair(escape,c,Q))<0 and all(sign(dot(escape,x))<0 for x in CV)
 # Exhaust all finite-pair nine-cliques independently of saved proposals.
 stack=[((),set(range(27)))];cliques=[];pdsets=[]
 while stack:
  ids,candidates=stack.pop()
  if len(ids)==9:
   cliques.append(ids)
   if positive([[G[i][j]for j in ids]for i in ids]):pdsets.append(ids)
   continue
  if len(ids)+len(candidates)<9:continue
  for i in candidates:stack.append((ids+(i,),{j for j in candidates&adj[i]if j>i}))
 assert len(cliques)==351 and len(pdsets)==65
 assert set(pdsets)=={tuple(v['active'])for v in d['finite_vertices']}
 V=[]
 for row in d['finite_vertices']:
  v=[dec(x)for x in row['point']];q=pair(v,v,Q)
  assert q==dec(row['norm']) and sign(q)<0 and sign(pair(v,c,Q))<0
  signs=[sign(dot(v,x))for x in CV];assert signs==row['all_wall_signs'] and max(signs)<=0
  assert [i for i,x in enumerate(signs)if not x]==row['active'];V.append(v)
 # Every vertex release occurs exactly once as a distinct support; all
 # existing endpoint extensions are compared with the archived incidence.
 expected=Counter(tuple(i for i in v['active']if i!=j)for v in d['finite_vertices']for j in v['active'])
 assert len(expected)==441 and set(expected)=={tuple(e['support'])for e in d['edges']}
 affine_candidates=0;counts=Counter();profile_count=tail_count=0
 vmap={tuple(v['active']):v['id']for v in d['finite_vertices']}
 for e in d['edges']:
  ids=e['support'];assert expected[tuple(ids)]==len(e['finite_vertices']);B=inverse([[G[i][j]for j in ids]for i in ids])
  ends=[]
  for j in range(27):
   if j in ids or any(j not in adj[i]for i in ids):continue
   b=[G[i][j]for i in ids];schur=ka(two,kn(ksum(km(km(b[a],B[a][bidx]),b[bidx])for a in range(8)for bidx in range(8))))
   if ks(schur)>0:
    key=tuple(sorted(ids+[j]));assert key in vmap;ends.append(vmap[key])
   elif schur==Z:
    affine_candidates+=1;coeff=[ksum(km(x,y)for x,y in zip(row,b,strict=True))for row in B]
    z=[sub(W[j][a],sm(scale(t,W[i][a])for t,i in zip(coeff,ids,strict=True)))for a in range(10)]
    assert not pair(z,z,Q)
    if sign(pair(z,c,Q))>0:z=[neg(x)for x in z]
    assert any(sign(dot(z,x))>0 for x in CV)
  assert sorted(ends)==sorted(e['finite_vertices']);counts[e['kind']]+=1
  if e['kind']=='finite_endpoint':assert len(ends)==2;continue
  assert len(ends)==1 and e['kind']=='open_ideal_boundary'
  r=e['spacelike_continuation'];y=[dec(x)for x in r['point']];q=pair(y,y,Q)
  assert q==dec(r['norm']) and sign(q)>0 and sign(pair(y,c,Q))<0
  signs=[sign(dot(y,x))for x in CV];assert signs==r['all_wall_signs'] and max(signs)<=0
  assert [i for i,x in enumerate(signs)if not x]==ids
  # Independent forest profile enumeration, with each isolated dihedral
  # family kept symbolically rather than truncated at a label bound.
  profile=d['open_edge_spherical_profiles'][str(e['id'])];comp=profile['components'];assert sorted(i for row in comp for i in row)==ids
  which={i:k for k,row in enumerate(comp)for i in row}
  assert all((G[i][j]==Z)for i in ids for j in ids if which[i]!=which[j])
  for row in comp:
   reach={row[0]}
   while True:
    new=reach|{j for j in row if any(G[i][j]!=Z for i in reach)}
    if new==reach:break
    reach=new
   assert reach==set(row)
  squares={3:ONE,4:two,5:(F(3,2),F(1,2))};found={():Z};tails=[]
  for n in [1,2,3]:
   for indices in combinations(range(8),n):
    if len({which[ids[i]]for i in indices})<n:continue
    if n==1 and len(comp[which[ids[indices[0]]]])==1:
     assert B[indices[0]][indices[0]]==(F(1,2),F(0));tails.append(ids[indices[0]]);continue
    for labs in product([3,4,5],repeat=n)if n<3 else [(3,3,3)]:
     p=ksum(km(squares[m],B[i][i])for i,m in zip(indices,labs,strict=True))
     if ks(ka(two,kn(p)))>0:found[tuple((ids[i],m)for i,m in zip(indices,labs,strict=True))]=p
  saved={tuple(map(tuple,r['contacts'])):K(r['projection_square'])for r in profile['finite_profiles']};assert saved==found and len(saved)==len(profile['finite_profiles'])
  assert sorted(tails)==sorted(x['wall']for x in profile['dihedral_tail_profiles'])
  assert all(x['labels']=='every integer m >= 3' and x['residual_norm']=='2 sin(pi/m)^2'for x in profile['dihedral_tail_profiles'])
  for r in profile['finite_profiles']:assert K(r['residual_norm'])==ka(two,kn(K(r['projection_square'])))
  profile_count+=len(saved);tail_count+=len(tails)
 assert counts==d['edge_counts'] and [profile_count,tail_count]==[10057,418]
 reached={0}
 while True:
  new=reached|{i for e in d['edges']if reached.intersection(e['finite_vertices'])for i in e['finite_vertices']}
  if new==reached:break
  reached=new
 assert len(reached)==65 and d['finite_edge_graph_connected']
 return {'status':'PASS_PRODUCER_FREE_EXACT_FRONTIER_REPLAY','sources':{str(p.relative_to(ROOT)):sha(p)for p in [SOURCE,WALLS,Path(__file__)]},'arithmetic':'Independent standard-library Fraction polynomial identities and 160-bit outward rational square-root intervals; no producer imports or numerical signs.','wall_count':27,'all_wall_pairs':351,'nine_cliques':351,'finite_vertices':65,'edges':441,'edge_counts':dict(counts),'infeasible_affine_extension_candidates':affine_candidates,'finite_profile_records':profile_count,'isolated_dihedral_tails':tail_count,'source_strict_spacelike_escape_rechecked':True,'scope':'Complete finite-corner/incident-edge census and necessary spherical finite-endpoint profiles. No census of vertices or edges outside the hyperbolic ball, and no new root-search completeness claim.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run()
 if a.verify:assert json.loads(OUT.read_text())==out
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['finite_vertices'],out['edge_counts'],flush=True)
