#!/usr/bin/env python3
"""Exact rank>=4 one-omitted affine-factor direction census.

This is a local direction classification. It is not a global H20 existence
or nonexistence certificate; rank-three factors are explicitly separate.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from collections import Counter
from pathlib import Path
import argparse, hashlib, json, math

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_omitted_factor_classification_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
ZERO=(F(0),)*8
RAD=[1,2,3,6,5,10,15,30]
def ff(a):return (F(a),)+(F(0),)*7
ONE=ff(1)
def add(a,b):return tuple(x+y for x,y in zip(a,b,strict=True))
def neg(a):return tuple(-x for x in a)
def sub(a,b):return add(a,neg(b))
def scale(a,c):return tuple(x*c for x in a)
MUL=[[(i^j,math.prod(p for k,p in enumerate((2,3,5))if (i&j)&(1<<k)))for j in range(8)]for i in range(8)]
def mul(a,b):
 z=[F(0)]*8
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:
     k,c=MUL[i][j];z[k]+=c*x*y
 return tuple(z)
def sqrtq(a):
 a=F(a);assert a>=0
 if not a:return ZERO
 n=a.numerator*a.denominator;c=1;d=1;p=2
 while p*p<=n:
  e=0
  while n%p==0:n//=p;e+=1
  c*=p**(e//2)
  if e%2:d*=p
  p+=1
 if n>1:d*=n
 assert d in RAD,(a,d)
 z=[F(0)]*8;z[RAD.index(d)]=F(c,a.denominator);return tuple(z)
@lru_cache(None)
def bounds(bits):
 den=1<<bits;out=[]
 for d in RAD:
  n=math.isqrt(d*den*den);lo=F(n,den);hi=lo if n*n==d*den*den else F(n+1,den);out.append((lo,hi))
 return out
def sign(a):
 if a==ZERO:return 0
 for bits in(48,96,192,384,768):
  lo=hi=F(0)
  for x,(l,u)in zip(a,bounds(bits),strict=True):
   lo+=x*(l if x>=0 else u);hi+=x*(u if x>=0 else l)
  if lo>0:return 1
  if hi<0:return -1
 raise AssertionError(('unresolved sign',a))
def real(a):return sum(float(x)*math.sqrt(d)for x,d in zip(a,RAD,strict=True))
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def inv(A):
 n=len(A);M=[[F(x)for x in r]+[F(i==j)for j in range(n)]for i,r in enumerate(A)]
 for j in range(n):
  i=next(i for i in range(j,n)if M[i][j]);M[i],M[j]=M[j],M[i];q=M[j][j];M[j]=[x/q for x in M[j]]
  for i in range(n):
   if i!=j:
    q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def pd(A):
 A=[[F(x)for x in r]for r in A];out=[]
 while A:
  q=A[0][0];assert q>0;out.append(q)
  A=[[A[i][j]-A[i][0]*A[0][j]/q for j in range(1,len(A))]for i in range(1,len(A))]
 return out
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(y)for y in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def fdot(a,b):
 z=ZERO
 for x,y in zip(a,b,strict=True):z=add(z,mul(x,y))
 return z
def flinear(a,c):
 z=ZERO
 for x,y in zip(a,c,strict=True):z=add(z,scale(x,y))
 return z

@lru_cache(None)
def contact_options(di,dj,g):
 ai=sqrtq(F(di,2));out=set()
 for c in(F(0),F(1,4),F(1,2)):
  A=F(di,2)*(1+c);disc=F(g*g*di,2)-A*(g*g-c*F(di*dj,2))
  if disc<0:continue
  for branch in(-1,1):
   aj=scale(add(scale(ai,g),scale(sqrtq(disc),branch)),1/A)
   if not sign(aj):continue
   nj=sub(ff(dj),mul(aj,aj));p=sub(ff(g),mul(ai,aj))
   if sign(nj)<=0 or sign(p)>0:continue
   assert mul(p,p)==scale(nj,c*F(di,2))
   out.add(aj)
 return tuple(sorted(out))

@lru_cache(None)
def local_candidates(d,A,J,boundary):
 k=len(d);opts=[(sqrtq(F(d[0],2)),)]+[contact_options(d[0],d[j],A[0][j])for j in range(1,k)]
 accepted=[];trials=0
 for a in product(*opts):
  trials+=1
  if any(mul(a[i],a[i])!=ff(F(d[i],2))for i in boundary):continue
  q=ZERO
  for i in range(k):
   for j in range(k):q=add(q,scale(mul(a[i],a[j]),J[i][j]))
  if q!=ONE:continue
  norms=[sub(ff(d[i]),mul(a[i],a[i]))for i in range(k)]
  assert all(sign(z)>0 for z in norms)
  okay=True
  for i,j in combinations(range(k),2):
   p=sub(ff(A[i][j]),mul(a[i],a[j]));nn=mul(norms[i],norms[j])
   if sign(p)>0 or not any(mul(p,p)==scale(nn,c)for c in(F(0),F(1,4),F(1,2))):okay=False;break
  if okay:accepted.append(a)
 return tuple(accepted),trials

def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']]
 W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']]
 C=[[dot(q,w)for q in Q]for w in W];G=[[dot(w,c)for c in C]for w in W]
 assert Counter(G[i][i]for i in range(210))=={2:42,4:168}
 assert all(G[i][j]<=0 for i,j in combinations(range(210),2))
 walk=json.loads((ROOT/WALK).read_text());assert walk['sources'][RAW]==sha(RAW)and not walk['pending_orbit_ids']
 cusps=[r for r in walk['corner_orbits']if r['type']=='ideal'];assert len(cusps)==18
 rows=[];counts=Counter();unique={};deletions=[];excluded_small=[];systems=set()
 for cusp in cusps:
  assert sum(f['rank']for f in cusp['affine_components'])==20
  for fac in cusp['affine_components']:
   ids=fac['walls'];n=len(ids)-1;marks=list(map(F,fac['positive_kernel']))
   assert min(marks)>0 and all(sum(G[i][j]*marks[k]for k,j in enumerate(ids))==0 for i in ids)
   D=[sum(marks[j]*W[i][d]for j,i in enumerate(ids))for d in range(22)]
   assert D==[F(fac['ray_multiple'])*F(x)for x in cusp['primitive_ray']]
   assert dot(D,[dot(q,D)for q in Q])==0
   counts['all_factors']+=1
   if n<4:
    excluded_small.append({'cusp':cusp['id'],'walls':ids,'rank':n});continue
   counts['rank_at_least_four_factors']+=1
   for omit in range(n+1):
    rr=[i for i in range(n+1)if i!=omit];walls=[ids[i]for i in rr];m=[marks[i]for i in rr]
    B=[[F(G[i][j])for j in walls]for i in walls];piv=pd(B);I=inv(B)
    assert min(x for r in I for x in r)>=0
    assert all(B[i][j]**2/(B[i][i]*B[j][j])in(F(0),F(1,4),F(1,2))for i,j in combinations(range(n),2))
    N=[{j for j in range(n)if i!=j and B[i][j]}for i in range(n)]
    components=[];remaining=set(range(n))
    while remaining:
     c={min(remaining)};todo=list(c)
     for i in todo:
      for j in N[i]-c:c.add(j);todo.append(j)
     remaining-=c;components.append(c)
    assert all(sum(len(N[i]&c)for i in c)==2*(len(c)-1)for c in components)
    counts['deletions']+=1;before=len(rows)
    for size in range(1,4):
     for ss in combinations(range(n),size):
      counts['supports']+=1;S=set(ss)
      if any(not S&c for c in components):counts['uncovered_component']+=1;continue
      boundary=[i for i in ss if N[i]-S]
      assert boundary # Otherwise the unchanged outside support is disconnected.
      if any(B[i][j]**2*4!=B[i][i]*B[j][j]for i in ss for j in N[i]-S):counts['old_label_four_boundary']+=1;continue
      order=[boundary[0]]+[i for i in ss if i!=boundary[0]]
      key=(tuple(B[i][i]for i in order),tuple(tuple(B[i][j]for j in order)for i in order),tuple(tuple(I[i][j]for j in order)for i in order),tuple(j for j,i in enumerate(order)if i in boundary))
      systems.add(key);local,nt=local_candidates(*key);counts['branch_trials_with_repetition']+=nt
      for aa in local:
       a=[ZERO]*n
       for i,z in zip(order,aa,strict=True):a[i]=z
       k=[flinear(a,r)for r in I];ks=[sign(z)for z in k]
       if min(ks)>0:pass
       elif max(ks)<0:a=list(map(neg,a));k=list(map(neg,k))
       else:counts['nonpositive_kernel']+=1;continue
       assert fdot(a,k)==ONE
       # All entries of the full projected Gram are checked, including
       # outside-support pairs; this independently enforces irreducibility.
       H=[[sub(ff(B[i][j]),mul(a[i],a[j]))for j in range(n)]for i in range(n)]
       assert all(sign(H[i][i])>0 for i in range(n))
       edge=[]
       for i,j in combinations(range(n),2):
        assert sign(H[i][j])<=0
        cs=[c for c in(F(0),F(1,4),F(1,2))if mul(H[i][j],H[i][j])==scale(mul(H[i][i],H[j][j]),c)]
        assert len(cs)==1
        if cs[0]:edge.append([walls[i],walls[j],cs[0]])
       assert all(fdot(r,k)==ZERO for r in H)
       e=[flinear(k,[W[i][d]for i in walls])for d in range(22)]
       assert fdot(e,[flinear(e,q)for q in Q])==ONE
       assert flinear(e,[dot(q,D)for q in Q])==ZERO
       pairfac=[flinear(e,c)for c in [C[i]for i in ids]]
       assert pairfac[omit]==scale(flinear(a,m),-1/marks[omit])
       low=0
       for i in range(1,n):
        if sign(sub(scale(k[i],1/m[i]),scale(k[low],1/m[low])))<0:low=i
       L=scale(k[low],1/m[low]);assert all(sign(sub(scale(k[i],1/m[i]),L))>=0 for i in range(n))
       # Exact projective direction key, using normalization by the
       # positive omitted contact instead of division in the number field.
       match=None
       for oldid in unique.get((cusp['id'],tuple(ids)),[]):
        other=rows[oldid]
        if all(mul(x,other['_pairfac'][omit])==mul(y,pairfac[omit])for x,y in zip(pairfac,other['_pairfac'],strict=True)):
         match=oldid;break
       idx=len(rows);entry={'cusp':cusp['id'],'factor':ids,'omitted_wall':ids[omit],'retained_walls':walls,
        'support':[walls[i]for i in range(n)if a[i]!=ZERO],'retained_pairings':a,'positive_kernel':k,
        'unit_normal_constant':e,'null_shift':D,
        'parameter_lower_magnitude':L,'lower_bound_witness_wall':walls[low],
        'projected_edges_squared_cosines':edge,'same_spatial_plane_as_row':match,'_pairfac':pairfac}
       rows.append(entry)
       if match is None:unique.setdefault((cusp['id'],tuple(ids)),[]).append(idx)
    deletions.append({'cusp':cusp['id'],'factor':ids,'omitted_wall':ids[omit],'retained_LDL':piv,'local_directions':len(rows)-before})
 counts['accepted_omission_directions']=len(rows);counts['distinct_spatial_planes']=sum(map(len,unique.values()))
 counts['distinct_local_systems']=local_candidates.cache_info().currsize
 counts['distinct_branch_trials']=sum(local_candidates(*key)[1]for key in systems)
 for r in rows:r.pop('_pairfac')
 return enc({'status':'PASS_EXACT_RANK_AT_LEAST_FOUR_ONE_OMITTED_LOCAL_CLASSIFICATION',
  'sources':{p:sha(p)for p in[RAW,WALK,f'scripts/{TAG}.py']},'field_basis_square_roots':RAD,
  'counts':dict(counts),'unclassified_small_factors':excluded_small,'deletions':deletions,'local_directions':rows,
  'scope':'Complete local direction census for exactly one omitted facet of one active affine factor of rank at least four. Rank-three factors and global Coxeter ridge conditions are not certified by this file.'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['counts'])
