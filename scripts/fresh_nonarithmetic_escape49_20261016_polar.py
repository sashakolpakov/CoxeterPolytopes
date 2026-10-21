#!/usr/bin/env python3
"""Eighty-digit discovery screen on every new48-wall exterior polar."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_joint_regions_20261016_batch as B
D,A,R,C,K,Z,O=B.D,B.A,B.R,B.C,B.K,B.Z,B.O
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_escape49_20261016_polar'
BATCH='results/fresh_nonarithmetic_common_escape_20261016_cap.json';CONE='results/fresh_nonarithmetic_common_escape_20261016_cone.json'
mp.mp.dps=80
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def scalar(raw):
 z=B.iv(tuple(A.B.db(x)for x in raw));v=(z.l+z.h)/2
 return mp.mpf(v.numerator)/v.denominator
def dot(v,w):return sum((a*b for a,b in zip(v,w)),mp.mpf(0))
def mv(Q,v):return[dot(row,v)for row in Q]
def null(rows):
 a=[list(x)for x in rows];piv=[];r=0
 for c in range(10):
  j=max(range(r,len(a)),key=lambda j:abs(a[j][c]))
  if abs(a[j][c])<mp.mpf('1e-38'):continue
  a[r],a[j]=a[j],a[r];q=a[r][c];a[r]=[x/q for x in a[r]]
  for j in range(len(a)):
   if j!=r:
    q=a[j][c];a[j]=[x-q*y for x,y in zip(a[j],a[r])]
  piv.append(c);r+=1
  if r==len(a):break
 assert r==9,(len(rows),r)
 free=next(c for c in range(10)if c not in piv);v=[mp.mpf(0)]*10;v[free]=1
 for i,c in enumerate(piv):v[c]=-a[i][free]
 assert max(abs(dot(row,v))for row in rows)<mp.mpf('1e-33')
 return v
def run():
 b=B.read(BATCH);c=B.read(CONE)
 Q=[[scalar(D.encfield(D.lift(A.B.db(x))))for x in row]for row in b['form']]
 N=[[scalar(x)for x in row]for row in b['normals']];H=[scalar(x)for x in b['center']];E=[scalar(x)for x in c['strict_escape']['point']]
 E=[x/-dot(E,mv(Q,H))for x in E]
 H=[x/mp.sqrt(-dot(H,mv(Q,H)))for x in H]
 normalized=[]
 for n in N:
  q=dot(n,mv(Q,n));normalized.append([x*mp.sqrt(2/q)for x in n])
 N=normalized;cov=[mv(Q,n)for n in N];hcov=mv(Q,H);ecov=mv(Q,E)
 rows=[];cutting_screens=[];stats=Counter();allrays=[];maxerr=mp.mpf(0)
 for i,ids in enumerate(c['complete_census']['hyperideal_rays']):
  v=null([cov[j]for j in ids]);v=[x/-dot(v,hcov)for x in v];q=dot(v,mv(Q,v));assert q>0
  old=[dot(v,x)for x in cov];assert max(old)<mp.mpf('1e-30');maxerr=max(maxerr,max(old));allrays.append(v)
  r=[x*mp.sqrt(2/q)for x in v]
  if dot(r,ecov)<=0:stats['noncutting']+=1;continue
  ps=[dot(r,x)for x in cov];finite=[];bad=[]
  for j,g in enumerate(ps):
   if abs(g)<mp.mpf('1e-30'):finite.append([j,2]);continue
   if g<-2-mp.mpf('1e-30'):continue
   if abs(g+2)<mp.mpf('1e-30'):bad.append([j,'near_parallel']);continue
   m=mp.pi/mp.acos(-g/2);mm=int(mp.nint(m))
   if mm>=2 and abs(m-mm)<mp.mpf('1e-28'):finite.append([j,mm])
   else:bad.append([j,mp.nstr(m,30)])
  stats['cutting']+=1
  cutting_screens.append(dict(ray_index=i,active=ids,escape_pairing=mp.nstr(dot(r,ecov),40),finite=finite,unresolved=bad))
  if not bad:
   row=dict(ray_index=i,active=ids,finite=finite,escape_pairing=mp.nstr(dot(r,ecov),40),point=[mp.nstr(x,65)for x in r]);rows.append(row)
   print('proposal',i,'E',row['escape_pairing'],'finite',finite,flush=True)
  if i%200==0:print('screened',i,flush=True)
 # Coverage is a lower-precision ranking only, not a claim of exact signs.
 Qf=np.array([[float(x)for x in row]for row in Q]);V=np.array([[float(x)for x in v]for v in allrays])
 for row in rows:
  r=np.array([float(x)for x in row['point']]);g=V@Qf@r;row['numerical_cut_rays']=np.flatnonzero(g>1e-7).tolist();row['numerical_cut_count']=len(row['numerical_cut_rays'])
 paths=[f'scripts/{TAG}.py',BATCH,CONE,f'scripts/{B.TAG}.py']
 return dict(status='NUMERICAL_80_DIGIT_COMPLETE_DIRECTION_SCREEN_ONLY',sources={p:sha(p)for p in paths},direction_count=len(allrays),proposals=rows,cutting_screens=cutting_screens,escape_normalization='E/(-E.H)',stats=dict(stats),max_reconstruction_positive_residual=mp.nstr(maxerr,40),scope='Numerical discovery at80decimal digits from160bit outward input enclosures. Exact acceptance and all-label rejection remain unproved; near-parallel contacts unresolved. Coverage scores use double precision only.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved next-polar numerical replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['proposals']))
