#!/usr/bin/env python3
"""Independent stdlib Fraction reconstruction of the31-wall H22 cusp link.

No producer, project arithmetic, NumPy or SymPy imports. Global orbit arrays
and distant pair gates are outside this separate local reconstruction.
"""
import argparse,hashlib,json
from fractions import Fraction as F
from collections import Counter
from functools import reduce
from math import gcd
from pathlib import Path
P=Path(__file__).resolve().parents[1]
PREFIX='fresh_h22_transported_cusp_completion_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
PROPOSAL='data/enumerations/'+PREFIX+'_local_proposal.json'
SOURCE='results/'+PREFIX+'_verify.json'
OUT='results/'+PREFIX+'_local_independent.json'
SHIFTS=[3,4,3,4,4,4,4,4,5,8,7,9,6]
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),F(0))
def mv(M,x):return[dot(r,x)for r in M]
def rank(M):
 A=[list(map(F,r))for r in M];k=0
 for j in range(len(A[0])):
  i=next((i for i in range(k,len(A))if A[i][j]),None)
  if i is None:continue
  A[i],A[k]=A[k],A[i];q=A[k][j];A[k]=[x/q for x in A[k]]
  for i in range(k+1,len(A)):
   q=A[i][j];A[i]=[x-q*y for x,y in zip(A[i],A[k],strict=True)]
  k+=1
  if k==len(A):break
 return k
def inertia(M):
 A=[list(map(F,r))for r in M];pos=neg=zero=0
 while A:
  n=len(A);i=next((i for i in range(n)if A[i][i]),None)
  if i is None:
   ij=next(((i,j)for i in range(n)for j in range(i)if A[i][j]),None)
   if ij is None:zero+=n;break
   i,j=ij
   for k in range(n):A[i][k]+=A[j][k]
   for k in range(n):A[k][i]+=A[k][j]
  if i:
   A[0],A[i]=A[i],A[0]
   for row in A:row[0],row[i]=row[i],row[0]
  q=A[0][0];assert q;pos+=q>0;neg+=q<0
  A=[[A[i][j]-A[i][0]*A[0][j]/q for j in range(1,n)]for i in range(1,n)]
 return pos,neg,zero
def kind(p,q,n):
 assert p<=0
 if p==0:return'right'
 c=p*p/(q*n)
 if c in(F(1,4),F(1,2),F(3,4)):return{F(1,4):'pi/3',F(1,2):'pi/4',F(3,4):'pi/6'}[c]
 assert c>=1
 return'parallel'if c==1 else'dotted'
def run():
 raw=read(RAW)['coxeter_input'];proposal=read(PROPOSAL);primary=read(SOURCE)
 assert primary['status']=='PASS_EXACT_COMPLETE_TRANSPORTED_H22_CUSP_REPAIR'
 assert all(sha(p)==h for p,h in primary['sources'].items())
 Q=[[-F(x)for x in row]+[F(0)]for row in raw['basis_gram']]+[[F(0)]*22+[F(11)]]
 W=[[F(x)for x in row]+[F(0)]for row in raw['wall_coordinates_as_columns_transposed']]
 h=[sum(int(r[j])for r in W)for j in range(23)];g=reduce(gcd,h);h=[F(x,g)for x in h];hc=mv(Q,h)
 assert inertia(Q)==(22,1,0)and dot(h,hc)==-14
 assert all(dot(w,hc)<0 for w in W)
 delta=list(map(F,proposal['target']['delta']));oldcap=list(map(F,proposal['target']['cap']));incident=proposal['target']['incident'];extra=proposal['target']['extra']
 assert len(incident)==17 and len(set(incident))==17 and extra in incident
 assert dot(delta,mv(Q,delta))==0 and dot(delta,hc)<0
 assert [a+b for a,b in zip(oldcap,W[extra],strict=True)]==delta
 assert dot(oldcap,mv(Q,oldcap))==2 and dot(oldcap,hc)<0
 olddelta=[dot(w,mv(Q,delta))for w in W];assert all(x<=0 for x in olddelta)and [i for i,x in enumerate(olddelta)if x==0]==incident
 bases=[list(map(F,row))for row in proposal['bases']];assert len(bases)==13
 caps=[[x+t*y for x,y in zip(b,delta,strict=True)]for b,t in zip(bases,SHIFTS,strict=True)]
 norms=[dot(r,mv(Q,r))for r in caps];assert norms==list(map(F,proposal['base_norms']))==list(map(F,primary['new_norms']))
 assert all(dot(r,hc)<0 and dot(r,mv(Q,delta))==0 and dot(r,mv(Q,oldcap))==0 for r in caps)
 den=primary['denominator'];assert primary['shifts']==SHIFTS
 assert caps==[[F(x,den)for x in row]+[F(z,den)]for row,z in zip(primary['representative_old_numerators'],primary['representative_height_numerators'],strict=True)]
 oldnorms=[dot(w,mv(Q,w))for w in W];oldcounts=[]
 for r,n in zip(caps,norms,strict=True):
  cc=Counter()
  for i,(w,m)in enumerate(zip(W,oldnorms,strict=True)):
   typ=kind(dot(r,mv(Q,w)),n,m);cc[typ]+=1
   if i not in incident:assert typ=='dotted'
  oldcounts.append(dict(cc))
 assert oldcounts==primary['old210_pair_counts']
 local=[W[i]for i in incident]+[oldcap]+caps;cov=[mv(Q,r)for r in local];G=[[dot(r,c)for c in cov]for r in local]
 assert len(local)==31 and rank(local)==22 and inertia(G)==(21,0,10)
 kinds=Counter(kind(G[i][j],G[i][i],G[j][j])for i in range(31)for j in range(i))
 newkinds=Counter(kind(G[18+i][18+j],norms[i],norms[j])for i in range(13)for j in range(i));assert dict(newkinds)==primary['local_new_pair_counts']
 assert [sum(r[j]for r in local)for j in range(23)]==[66*x for x in delta]
 return {'status':'PASS_INDEPENDENT_STDLIB_FRACTION_TRANSPORTED_H22_LOCAL_LINK','sources':{p:sha(p)for p in(RAW,PROPOSAL,SOURCE)},'ambient_signature':[22,1,0],'center_norm':'-14','field_coordinate':'last coordinate is coefficient ofsqrt11; rationalformdiag(Q,11)','new_norms':list(map(str,norms)),'local_normal_count':31,'normal_rank':22,'gram_inertia':[21,0,10],'positive_weights':[1]*31,'positive_null_relation_multiple':66,'new_old210_pair_counts':oldcounts,'local_pair_counts':dict(kinds),'new_new_pair_counts':dict(newkinds),'scope':'Independent rational reconstruction ofrepresentative13caps,all2730old210pairs,all465localpairs,strictcenter,nullincidence,signature,positivefullsupportrelationandranks. It proves thelocal21dimensionalEuclideancusplinkbounded. Completeorbittransfer/distantpairs dependontheprimaryexactcertificateandsource-review; theyarenotrecomputed here.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(P/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
