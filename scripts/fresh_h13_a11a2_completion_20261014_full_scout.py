#!/usr/bin/env python3
"""Bounded full integral reflective-root proposal in A11+A2+<-2>.

The time-zero sector is A11+G2, explicitly including the norm-six mirror.
"""
import argparse,json,time
from pathlib import Path
from functools import lru_cache
from math import isqrt,gcd
from fractions import Fraction as F
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h13_a11a2_completion_20261014_full_scout'
def ip(x,y):return sum(a*b for a,b in zip(x[1:],y[1:]))-2*x[0]*y[0]
@lru_cache(None)
def seq(n,q,s,lo,d,r):
 if q<0 or n*q<s*s or s<n*lo or (s-n*r)%d:return ()
 if n==1:return ((s,),)if s>=lo and s*s==q and(s-r)%d==0 else ()
 out=[];start=max(lo,-isqrt(q));start+=(r-start)%d
 for a in range(start,min(s//n,isqrt(q))+1,d):
  if (q-a*a)*(n-1)<(s-a)**2:continue
  for tail in seq(n-1,q-a*a,s-a,a,d,r):out.append((a,)+tail)
 return tuple(out)
def vectors(n,q,d):
 return tuple(x for r in range(d)if n*r%d==0 for x in seq(n,q,0,-isqrt(q),d,r))
def run(height):
 N=[]
 for i in range(11):
  v=[0]*16;v[1+i]=1;v[2+i]=-1;N.append(tuple(v))
 N += [(0,)+(0,)*12+(0,1,-1),(0,)+(0,)*12+(1,-2,1)]
 history=[];begin=time.monotonic();pool=[2,4,6,8,12,24]
 levels=sorted((F(t*t,q),t,q)for q in pool for t in range(1,isqrt(height*q)+1)if (4*t)%q==0)
 for h,t,q in levels:
  total=2*t*t+q;accepted=[];count=0;d=q//2
  for qx in range(0,total+1,2):
   xs=vectors(12,qx,d);ys=[y for y in vectors(3,total-qx,d)if y[1]>=0]
   for x in xs:
    for y in ys:
     v=(t,)+x+y
     if gcd(*v)!=1:continue
     assert ip(v,v)==q and all(2*(x[i]-x[i+1])%q==0 for i in range(11))and all(2*(y[i]-y[i+1])%q==0 for i in range(2));count+=1
     if all(ip(v,r)<=0 for r in N):accepted.append(v)
  assert all(ip(a,b)<=0 for i,a in enumerate(accepted)for b in accepted[:i])
  N+=accepted;history.append(dict(height=str(h),t=t,q=q,candidates=count,accepted=accepted))
  print(h,t,q,count,len(accepted),len(N),round(time.monotonic()-begin,2),flush=True)
 out=dict(status='BOUNDED_EXACT_FULL_INTEGRAL_REFLECTION_PREFIX',height=height,normals=N,gram=[[ip(a,b)for b in N]for a in N],levels=history,norm_pool=pool,
  coordinates='time,12 integral sumzero A11 coordinates,3 integral sumzero A2 coordinates; form -2time²+Euclidean',
  scope='Complete dominant integral reflective-root enumeration in safe norm pool through saved distance height; zero group A11+G2. No finite-volume or full-reflectivity conclusion.')
 path=ROOT/f'results/{TAG}.json';out['sources']={f'scripts/{TAG}.py':sha256(Path(__file__).read_bytes()).hexdigest()};path.write_text(json.dumps(out,indent=2)+'\n')
 return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--height',type=int,default=10);a=ap.parse_args();run(a.height)
