#!/usr/bin/env python3
"""Bounded classical diagonal H13 source and exact finite-corner proposals."""
import argparse,json,time
from pathlib import Path
from fractions import Fraction as F
from math import gcd,isqrt
from itertools import combinations
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h13_a11a2_completion_20261014_diagonal_scout'
@lru_cache(None)
def parts(total,n,maximum):
 if not n:return ((),)if total==0 else ()
 out=[]
 for a in range(min(maximum,isqrt(total)),-1,-1):
  if total-a*a>(n-1)*a*a:continue
  for tail in parts(total-a*a,n-1,a):out.append((a,)+tail)
 return tuple(out)
def run(d,height):
 def ip(x,y):return sum(a*b for a,b in zip(x[1:],y[1:]))-d*x[0]*y[0]
 N=[]
 for i in range(12):
  v=[0]*14;v[i+1]=-1;v[i+2]=1;N.append(tuple(v))
 N.append((0,)+(0,)*12+(-1,))
 pool=[q for q in range(1,2*d+1)if 2*d%q==0];history=[]
 for h,t,q in sorted((F(t*t,q),t,q)for q in pool for t in range(1,isqrt(height*q)+1)if 2*d*t%q==0):
  div=q//gcd(q,2);S=d*t*t+q
  if S%(div*div):continue
  cand=[(t,)+tuple(div*x for x in v)for v in parts(S//(div*div),13,isqrt(S)//div)]
  accepted=[]
  for r in cand:
   if gcd(*r)!=1:continue
   assert ip(r,r)==q
   if all(ip(r,s)<=0 for s in N):accepted.append(r)
  assert all(ip(r,s)<=0 for i,r in enumerate(accepted)for s in accepted[:i]);N+=accepted
  history.append(dict(height=str(h),q=q,t=t,candidates=len(cand),accepted=accepted))
 G=[[ip(r,s)for s in N]for r in N];norm=[G[i][i]for i in range(len(N))]
 # Each connected pi3 component has one norm, but the two components
 # can have different norms (notably the norm-six roots for d=3).
 hits=[];tested=0
 ix=[i for i in range(len(N))if any(i!=j and G[i][j]<0 and 4*G[i][j]**2==norm[i]*norm[j]for j in range(len(N)))]
 for J in combinations(ix,13):
   tested+=1
   if any(G[i][j]and (G[i][j]>0 or 4*G[i][j]**2!=norm[i]*norm[j])for i,j in combinations(J,2)):continue
   adj={i:{j for j in J if i!=j and G[i][j]}for i in J};unseen=set(J);components=[]
   while unseen:
    todo=[min(unseen)];C=set()
    while todo:
     i=todo.pop()
     if i in C:continue
     C.add(i);todo+=list(adj[i]-C)
    unseen-=C;degrees=[len(adj[i])for i in C]
    if max(degrees)>2 or sum(degrees)!=2*(len(C)-1):break
    components.append(sorted(C))
   else:
    if sorted(map(len,components))==[2,11]:hits.append(dict(incident=J,component_norms=[norm[c[0]]for c in components],components=components))
 print('d',d,'roots',len(N),'extra',N[13:],'hits',hits,'subsets',tested,flush=True)
 out=dict(status='BOUNDED_DIAGONAL_ROOT_PREFIX_AND_A11_A2_PROBE',d=d,height=height,normals=N,gram=G,norms=norm,levels=history,hits=hits,subsets=tested)
 (ROOT/f'results/{TAG}_d{d}.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--d',type=int,default=1);ap.add_argument('--height',type=int,default=16);a=ap.parse_args();run(a.d,a.height)
