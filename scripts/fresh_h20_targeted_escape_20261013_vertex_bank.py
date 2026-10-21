#!/usr/bin/env python3
"""All19-of24 short-wall vertex directions, then bounded height-cap profiles.

Vertex positive-definiteness and projection coefficients use Fraction only.
Cap angle proposals and bank compatibility in this script remain numerical.
"""
import argparse,hashlib,itertools,json,math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from fresh_h20_targeted_escape_20261013_second_volume import extended_setup
from fresh_h20_targeted_escape_20261013 import ev

ROOT=Path(__file__).resolve().parents[1]
VERTICES='data/enumerations/fresh_h20_targeted_escape_20261013_short_vertices.json'
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_vertex_bank.json'
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
VOL='data/enumerations/fresh_h20_targeted_escape_20261013_second_volume.json'
LABELS=[2,3,4,5,6,8,10,12]
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def positive_solve(M,rhs):
 """Sparse exact symmetric elimination: return M^-1 rhs iff M is PD."""
 n=len(M);A={i:{j:F(M[i][j])for j in range(n)if M[i][j]}for i in range(n)}
 b={i:F(rhs[i])for i in range(n)};left=set(range(n));steps=[]
 while left:
  i=min(left,key=lambda j:(len(A[j]),j));d=A[i].get(i,F(0))
  if d<=0:return None
  neighbors={j:x for j,x in A[i].items()if j!=i and j in left}
  steps.append((i,d,b[i],neighbors))
  for j,x in neighbors.items():
   b[j]-=x*b[i]/d
   for k,y in neighbors.items():
    a=A[j].get(k,F(0))-x*y/d
    if a:A[j][k]=a
    else:A[j].pop(k,None)
   A[j].pop(i,None)
  left.remove(i)
 x={}
 for i,d,b,neighbors in reversed(steps):x[i]=(b-sum((a*x[j]for j,a in neighbors.items()),F(0)))/d
 return [x[i]for i in range(n)]

def vertices():
 b=read(BASE);Q=[[F(x[0])for x in row[:20]]for row in b['form'][:20]]
 assert all(F(x[1])==0 for row in b['form']for x in row)
 oldR=[[F(x[0])for x in r[:20]]for r in b['old_normals']];R=oldR[:24]
 h=[F(x[0])for x in b['strict_center'][:20]]
 dot=lambda a,b:sum((x*y for x,y in zip(a,b,strict=True)),F(0))
 mv=lambda M,x:[dot(row,x)for row in M]
 cov=[mv(Q,r)for r in R];longcov=[mv(Q,r)for r in oldR[24:]];G=[[dot(a,c)for c in cov]for a in R];hc=mv(Q,h);bvec=[dot(r,hc)for r in R]
 out=[];count=0
 for ids in itertools.combinations(range(24),19):
  count+=1;M=[[G[i][j]for j in ids]for i in ids];rhs=[bvec[i]for i in ids]
  sol=positive_solve(M,rhs)
  if sol is None:continue
  x=[h[k]-sum((s*R[i][k]for s,i in zip(sol,ids,strict=True)),F(0))for k in range(20)]
  xx=dot(x,mv(Q,x));assert xx<0
  ps=[dot(x,c)for c in cov];assert all(p==0 if i in ids else p<0 for i,p in enumerate(ps))
  assert all(dot(x,c)<0 for c in longcov)
  out.append({'incident_short_walls':list(ids),'projection':[str(y)for y in x],
   'norm':str(xx),'gram_solution':list(map(str,sol))})
 assert count==42504 and len(out)==2046
 d={'status':'EXACT_COMPLETE_19_OF24_SHORT_VERTEX_DIRECTION_CENSUS','source':BASE,'source_sha256':sha(BASE),
  'subsets':count,'positive_definite_subsets':len(out),'vertices':out,
  'scope':'Every19-subset is decided by exact symmetric elimination; each accepted short-wall projection is timelike, annihilates exactlythose19shortwalls and isstrict onother5. This is the short-only finitevertex direction class, not all8346oldH19finitevertices.'}
 (ROOT/VERTICES).write_text(json.dumps(d,indent=2)+'\n');return d

def bank():
 data=read(VERTICES);assert data['source_sha256']==sha(BASE)
 Q,R,h,_,_,_,_=extended_setup();cov=R@Q;norms=np.einsum('ij,jk,ik->i',R,Q,R)
 w=np.array([ev(x)for x in read(VOL)['certificate']['vector']])
 finite=np.cos(np.pi/np.arange(2,121));labels=np.array(LABELS);c2=np.cos(np.pi/labels)**2
 anchors=np.arange(82);z=cov[anchors,-1];counts=Counter();found=[];seen={}
 for index,vertex in enumerate(data['vertices']):
  original=np.array([float(F(x))for x in vertex['projection']]+[0.]);q0=float(F(vertex['norm']))
  scale=np.sqrt(-2/q0);v=original*scale;q0=-2.;a=cov[anchors]@v
  aa=z[:,None]**2-3*norms[anchors,None]*c2
  bb=2*a[:,None]*z[:,None]*np.ones_like(c2)
  cc=a[:,None]**2-q0*norms[anchors,None]*c2
  disc=bb*bb-4*aa*cc;valid=(disc>=-1e-8)&(abs(aa)>1e-10)
  counts['quadratic_profiles']+=aa.size
  for branch in(-1,1):
   js,ms=np.nonzero(valid);tt=(-bb[js,ms]+branch*np.sqrt(np.maximum(0,disc[js,ms])))/(2*aa[js,ms])
   q=q0+3*tt*tt;keep=(q>1e-8)&(a[js]+z[js]*tt<=1e-7)
   js,ms,tt,q=js[keep],ms[keep],tt[keep],q[keep]
   if not len(js):continue
   rr=np.tile(v,(len(js),1));rr[:,-1]=tt;rr*=np.sqrt(2/q)[:,None]
   ratios=-(rr@cov.T)/np.sqrt(2*norms)
   ix=np.clip(np.searchsorted(finite,ratios),1,len(finite)-1)
   err=np.minimum(abs(ratios-finite[ix]),abs(ratios-finite[ix-1]));err=np.where(ratios>=1-1e-7,0,err)
   ok=np.flatnonzero(np.max(err,axis=1)<1e-7)
   for k in ok:
    r=rr[k];key=tuple(np.round(r/np.linalg.norm(r),10))
    if key in seen:counts['duplicate_profiles']+=1;continue
    seen[key]=len(found);found.append({'vertex_index':index,'anchor':int(js[k]),'label':int(labels[ms[k]]),
     'branch':branch,'height_for_normalized_projection':float(tt[k]),'normal':r.tolist(),
     'norm_error':abs(float(r@Q@r)-2),'escape_pairing':float(r@Q@w),
     'pair_ratios':ratios[k].tolist()})
  if index%100==0:print({'vertex':index,'caps':len(found),'cutters':sum(c['escape_pairing']>1e-7 for c in found)},flush=True)
 out={'status':'NUMERICAL_FINITE_VERTEX_DIRECTION_BOUNDED_HEIGHT_BANK',
  'source':VERTICES,'source_sha256':sha(VERTICES),'old_walls':87,'labels':LABELS,'anchor_indices':list(map(int,anchors)),
  'source_vertex_status':data['status'],'vertex_directions':len(data['vertices']),
  'counts':dict(counts),'caps':found,'escape_cutters':[i for i,c in enumerate(found)if c['escape_pairing']>1e-7],
  'scope':'Every direction in the referenced exact vertexcensus; solve finite angleevents against82basewalls for8labels, varyinglastcoordinate; testall87oldwalls numerically. Null quadratic leadingcoefficients skipped. Capentries are numerical, dedup isapproximate; noexactnewbank, all-angle completeness or volumeclaim.'}
 (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print({'done':True,'caps':len(found),'cutters':len(out['escape_cutters'])})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--vertices',action='store_true');p.add_argument('--bank',action='store_true');a=p.parse_args()
 if a.vertices:vertices()
 if a.bank:bank()
