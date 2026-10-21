#!/usr/bin/env python3
"""Exact closer placement of all eleven equivariant cusp-repair wall types."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from math import isqrt
import numpy as np
from fresh_h22_new_cusp_completion_20261013_verify import setup,families,PRE,P,read,sha
DEN=26040
SOURCE='results/'+PRE+'_verify.json';OUT='results/'+PRE+'_tight.json'
def upper(q):
 scale=10**6;n=isqrt(q.numerator*scale*scale//q.denominator)+1
 return F(n,scale)
def shifted(D,V,B,A,ids,sizes,Q=None,lower=False):
 for index,(name,base,ix)in enumerate(families(D,V,B,A,ids)):
  heights=np.where(sizes>=4,60,80)if index==0 else np.full(len(ix),[200,100,1,1][index-1],dtype=np.int64)
  if index==0 and lower:
   assert Q is not None
   n=np.einsum('ij,ij->i',A,A@Q)
   heights=np.where(sizes>=4,60,np.where(sizes==3,80,np.where(n==2,20,40)))
  old=base*5+heights[:,None]*(DEN//120)*D[ix]
  yield name,old,heights,ix,index

def run(lower=False):
 prior,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup()
 source=read(SOURCE);assert all(sha(p)==s for p,s in source['sources'].items())
 si=np.array([60,60,80,80,80,80,80,200,100,1,1],dtype=np.int64)
 if lower:si[4:7]=[20,20,40]
 bases=np.array(source['base_old_numerators'],dtype=np.int64)*5
 C=bases+si[:,None]*(DEN//120)*D[0];cov=C@Q
 norms=[F(x)for x in source['new_cap_norms']]
 nn=np.einsum('ij,ij->i',C,cov)+11*si*si*(DEN//10)**2
 assert [F(int(x),DEN**2)for x in nn]==norms
 assert np.all(np.einsum('ij,ij->i',C,cov)<0)and max(C@(Q@h))<0
 local=bases@Q@bases.T;maxima={}
 for name,roots,t,ix,index in shifted(D,V,B,A,ids,sizes,Q,lower):
  assert 22*int(abs(roots).max())*int(abs(cov).max())+11*int(max(si))*int(max(t))*(DEN//10)**2<2**63
  pair=roots@cov.T+11*t[:,None]*si*(DEN//10)**2
  for r,row in zip(roots[ix==0],pair[ix==0]):
   match=np.flatnonzero(np.all(C==r,axis=1));assert len(match)==1 and np.array_equal(row,local[match[0]])
  classes=np.einsum('ij,ij->i',A,A@Q)if index==0 else np.zeros(len(ix),dtype=np.int8)
  results={}
  for cls in np.unique(classes):
   nq=F(int(cls))if index==0 else norms[6+index]
   mx=pair[(ix!=0)&(classes==cls)].max(axis=0)
   for j,p in enumerate(mx):assert p<0 and F(int(p)**2,DEN**4)>norms[j]*nq
   results[str(nq)]=[str(F(int(x),DEN**2))for x in mx]
  maxima[name]=results
 oldnorm=np.diag(W@Q@W.T);pair=W@cov.T;counts=[]
 for j,q in enumerate(norms):
  counter=Counter()
  for value,n in zip(pair[:,j],oldnorm):
   a=F(int(value),DEN);b=a*a/(q*int(n))
   if not a:kind='right'
   elif a<0 and b in(F(1,4),F(1,2),F(3,4)):kind={F(1,4):'pi/3',F(1,2):'pi/4',F(3,4):'pi/6'}[b]
   elif a<0 and b==1:kind='parallel'
   else:assert a<0 and b>1;kind='dotted'
   counter[kind]+=1
  counts.append(dict(counter))
 U=np.load(P/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')['twice_old_projections']
 z=np.load(P/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz');K=z['old_null_projections'];four=z['four_heights']
 R=np.load(P/'data/enumerations/fresh_h22_gap_face_orbit_20261013.npz')['six_old_projections']
 bounds=[]
 for j,q in enumerate(norms):
  cc=cov[j];s=int(si[j]);assert 22*max(int(abs(x).max())for x in(U,K,R))*int(abs(cc).max())<2**63
  mu=int(max(U@cc));assert F(mu,2*DEN)+F(21*s,40)<-upper(q)
  dp=K@cc;ds={}
  for t in(1,4):
   v=int(max(dp[four==t]));assert F(4*v,DEN*t)+F(s,3)<-upper(q);ds[str(t)]=v
  hp=int(h@cc);assert F(hp,DEN)+F(4*s,3)<-upper(2*q)
  rp=R@cc+33*s*(DEN//5);zeros=np.flatnonzero(rp==0).tolist();assert zeros==[0,prior['new_cusp']['partner_index']]
  mr=int(max(rp[rp!=0]));assert mr<0 and F(mr*mr,36*DEN*DEN)>F(2,3)*q
  bounds.append({'twice_timelike_max':mu,'null_maxima':ds,'center_pair':hp,'six_k6_nonincident_max':mr})
 return {'status':'PASS_EXACT_'+('LOWER'if lower else 'TIGHT')+'_COMPLETE_H22_CUSP_REPAIR','source':SOURCE,'source_sha256':sha(SOURCE),'common_denominator':DEN,'twenty_shifts':si.tolist(),'shifts':[str(F(int(x),20))for x in si],'new_cap_old_numerators':C.tolist(),'new_cap_norms':list(map(str,norms)),'positive_cusps':len(D),'completed_signed_cusps':2*len(D),'total_facets':prior['total_facets']+22*len(D),'new_pair_nonlocal_maxima':maxima,'old210_pair_counts':counts,'old_family_bounds':bounds,'scope':'A different placement of the same complete eleven-type cusp-repair orbit; replaces the uniform-shift10 placement. Retains the2,736,428-wall k6 chamber and completes all241920 new cusp links. All finite pairs, actual facets and all m>=3 center labels pass. No finite-volume conclusion without a separate global certificate.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--lower',action='store_true');a=p.parse_args();d=run(a.lower);out=OUT.replace('_tight.json','_lower.json')if a.lower else OUT
 if a.verify:assert d==read(out)
 else:(P/out).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['total_facets'])
