#!/usr/bin/env python3
"""Exact bounded norm-two subgroup roots; numerical global proposal only."""
import argparse,json,time
from pathlib import Path
from functools import lru_cache
from math import isqrt
from fractions import Fraction as F
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h13_a11a2_completion_20261014_scout'
def ip(x,y):return sum(a*b for a,b in zip(x[1:],y[1:]))-2*x[0]*y[0]
@lru_cache(None)
def sorted_vectors(n,q,s,lo):
 if q<0 or n*q<s*s or s<n*lo:return ()
 if n==1:return ((s,),)if s>=lo and s*s==q else ()
 out=[]
 for a in range(max(lo,-isqrt(q)),min(s//n,isqrt(q))+1):
  if (q-a*a)*(n-1)<(s-a)**2:continue
  for tail in sorted_vectors(n-1,q-a*a,s-a,a):out.append((a,)+tail)
 return tuple(out)
def run(height):
 N=[]
 for start,n in [(1,12),(13,3)]:
  for i in range(n-1):
   v=[0]*16;v[start+i]=1;v[start+i+1]=-1;N.append(tuple(v))
 history=[];begin=time.monotonic()
 for t in range(1,height+1):
  total=2*t*t+2;accepted=[];records=[];count=0
  for q in range(0,total+1,2):
   xs=sorted_vectors(12,q,0,-isqrt(q));ys=sorted_vectors(3,total-q,0,-isqrt(total-q))
   for x in xs:
    for y in ys:
     v=(t,)+x+y;assert ip(v,v)==2;count+=1
     block=next((i for i,r in enumerate(N)if ip(v,r)>0),None)
     if block is None:accepted.append(v)
  assert all(ip(a,b)<=0 for i,a in enumerate(accepted)for b in accepted[:i])
  N+=accepted;history.append(dict(t=t,candidates=count,accepted=accepted))
  print(t,count,len(accepted),len(N),round(time.monotonic()-begin,2),flush=True)
 out=dict(status='BOUNDED_EXACT_NORM_TWO_SUBGROUP_PREFIX',height=height,normals=N,gram=[[ip(a,b)for b in N]for a in N],levels=history,
  coordinates='time,12 integral sumzero A11 coordinates,3 integral sumzero A2 coordinates; form -2time²+Euclidean',
  scope='Complete norm-two dominant root enumeration through saved integer time; no full-lattice reflection completeness, no finite-volume conclusion.')
 path=ROOT/f'results/{TAG}.json';out['sources']={f'scripts/{TAG}.py':sha256(Path(__file__).read_bytes()).hexdigest()};path.write_text(json.dumps(out,indent=2)+'\n')
 return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--height',type=int,default=8);a=ap.parse_args();run(a.height)
