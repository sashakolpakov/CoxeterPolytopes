#!/usr/bin/env python3
"""Bounded numeric global Gram surgery guided by exact affine circuit covers.

Changes every pair in each minimal finite-edge support to one of labels2..6,
excluding its old label. All old parallel pairs are freed as dotted distances.
The ten-wall Schur basis contains no free pair. This is a numerical scout;
any positive survivor requires exact field and full compactness verification.
"""
import argparse,hashlib,json,math,time
from itertools import combinations,product
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_h9_incidence_surgery_20261013.json'
OUT=ROOT/'results/fresh_h9_incidence_surgery_20261013_search.json'
BASES={
 'A14_8_1_C02':[0,2,4,5,6,7,8,10,11,12],
 'A14_10_1_C02':[0,1,2,3,4,7,8,9,11,13],
 'A16_14_1_C01':[1,2,3,4,5,6,8,9,10,13]}
def solve_remaining(H,B,P,inv):
 n=len(H);T=[i for i in range(n)if i not in B];F=np.array(H[np.ix_(T,B)])
 free=[(i,j)for i in range(len(T))for j in range(10)if tuple(sorted((T[i],B[j])))in P]
 assert len(free)<=1
 if not free:options=[[]]
 else:
  i,j=free[0];f=F[i].copy();f[j]=0
  aa=inv[j,j];bb=2*f@inv[:,j];cc=f@inv@f-1
  if abs(aa)<1e-12:
   options=[[(-cc/bb)]]if abs(bb)>1e-12 else []
  else:
   D=bb*bb-4*aa*cc
   options=[]if D<-1e-9 else [[(-bb+sgn*math.sqrt(max(0,D)))/(2*aa)]for sgn in[-1,1]]
 out=[]
 for x in options:
  V=F.copy()
  for (i,j),v in zip(free,x):V[i,j]=v
  S=V@inv@V.T
  if max(abs(np.diag(S)-1))>1e-7:continue
  new=H.copy();ok=True
  for i,j in combinations(range(len(T)),2):
   a,b=T[i],T[j]
   if (min(a,b),max(a,b))in P:
    if S[i,j]>=-1-1e-8:ok=False;break
    new[a,b]=new[b,a]=S[i,j]
   elif abs(S[i,j]-H[a,b])>1e-7:ok=False;break
  for (i,j),v in zip(free,x):
   if v>=-1-1e-8:ok=False;break
   new[T[i],B[j]]=new[B[j],T[i]]=v
  if ok:out.append(new.tolist())
 return out

def search(row,limit_supports=None):
 started=time.time();name=row['id'];H=np.array(row['source_positive_gram'],float);H/=np.sqrt(np.diag(H)[:,None]*np.diag(H)[None,:]);n=len(H)
 B=BASES[name];T=[i for i in range(n)if i not in B];P={tuple(p)for p in row['parallel_pairs']}
 assert not any(i in B and j in B for i,j in P)
 zero=[i for i in T if not any((min(i,j),max(i,j))in P for j in B)]
 pos={a:i for i,a in enumerate(B)}
 supports=row['finite_edge_cover']['minimum_supports'];assert row['finite_edge_cover']['enumeration_complete']
 if limit_supports is not None:supports=supports[:limit_supports]
 fallbacks=[]
 for C in combinations(range(n),10):
  if list(C)==B or any(i in C and j in C for i,j in P):continue
  TT=[i for i in range(n)if i not in C]
  cc=[sum((min(i,j),max(i,j))in P for j in C)for i in TT]
  if max(cc)>1 or cc.count(0)==0:continue
  ee=np.linalg.eigvalsh(H[np.ix_(C,C)])
  if min(abs(ee))<1e-8 or sum(ee<0)!=1:continue
  fallbacks.append((sum(cc),tuple(C)))
 fallbacks.sort();fallbacks=[list(x[1])for x in fallbacks[:32]]
 counts={'assignments':0,'basis_nonlorentz':0,'basis_singular_or_ill_conditioned':0,'singular_psd':0,'singular_one_negative':0,'singular_multiple_negative':0,'fixed_outside_norm_failure':0,'schur_distance_failure':0,'global_gram_survivors':0,'fallback_resolved_nonlorentz':0,'fallback_norm_failure':0,'fallback_metric_failure':0,'fallback_unresolved':0}
 near=[];survivors=[];singular=[]
 for si,support in enumerate(supports):
  choices=[]
  for i,j in support:
   old=round(math.pi/math.acos(-H[i,j]));choices.append([m for m in range(2,7)if m!=old])
  words=list(product(*choices))
  for start in range(0,len(words),1024):
   batch=words[start:start+1024];N=len(batch);A=np.repeat(H[None],N,axis=0)
   for k,(i,j)in enumerate(support):
    vals=-np.cos(np.pi/np.array([w[k]for w in batch]));A[:,i,j]=A[:,j,i]=vals
   counts['assignments']+=N;BB=A[:,B][:,:,B];ev=np.linalg.eigvalsh(BB)
   small=np.min(abs(ev),axis=1)<1e-9
   counts['basis_singular_or_ill_conditioned']+=int(sum(small))
   neg=np.sum(ev<-1e-9,axis=1)
   counts['singular_psd']+=int(sum(small&(neg==0)))
   counts['singular_one_negative']+=int(sum(small&(neg==1)))
   counts['singular_multiple_negative']+=int(sum(small&(neg>1)))
   for ii in np.where(small)[0][:max(0,30-len(singular))]:singular.append({'support_index':si,'labels':batch[ii]})
   good=(np.sum(ev<-1e-9,axis=1)==1)&(~small)
   counts['basis_nonlorentz']+=int(sum(~good&~small))
   inds=np.where(good)[0]
   pending=np.where(small)[0]
   for C in fallbacks:
    if not len(pending):break
    CC=A[pending][:,C][:,:,C];ce=np.linalg.eigvalsh(CC);cs=np.min(abs(ce),axis=1)<1e-9
    cg=(np.sum(ce<-1e-9,axis=1)==1)&(~cs)
    counts['fallback_resolved_nonlorentz']+=int(sum(~cg&~cs))
    ci=np.where(cg)[0]
    if len(ci):
     CI=np.linalg.inv(CC[ci]);ZZ=[i for i in range(n)if i not in C and not any((min(i,j),max(i,j))in P for j in C)]
     CF=A[pending[ci]][:,ZZ][:,:,C]
     cr=np.max(abs(np.einsum('bri,bij,brj->br',CF,CI,CF)-1),axis=1)
     counts['fallback_norm_failure']+=int(sum(cr>=1e-7))
     for z in np.where(cr<1e-7)[0]:
      # A fixed-basis singular branch gets an explicit pending record if
      # its alternate basis has several independent distance variables.
      TT=[i for i in range(n)if i not in C]
      freecount=sum((min(i,j),max(i,j))in P for i in TT for j in C)
      if freecount>1:
       counts['fallback_unresolved']+=1
       if len(singular)<30:singular.append({'support_index':si,'labels':batch[pending[ci[z]]],'reason':'multiple free rows after norm gate'})
      else:
       sols=solve_remaining(A[pending[ci[z]]],C,P,CI[z])
       if not sols:counts['fallback_metric_failure']+=1
       for gg in sols:
        survivors.append({'support_index':si,'changed_pairs':support,'labels':batch[pending[ci[z]]],'normalized_gram':gg});counts['global_gram_survivors']+=1
    pending=pending[cs]
   counts['fallback_unresolved']+=len(pending)
   if not len(inds):continue
   II=np.linalg.inv(BB[inds]);F=A[inds][:,zero][:,:,B]
   res=np.max(abs(np.einsum('bri,bij,brj->br',F,II,F)-1),axis=1)
   for kk in np.argsort(res)[:3]:
    near.append({'residual':float(res[kk]),'support_index':si,'labels':batch[inds[kk]]})
   near=sorted(near,key=lambda r:r['residual'])[:20]
   passed=res<1e-7;counts['fixed_outside_norm_failure']+=int(sum(~passed))
   for kk in np.where(passed)[0]:
    ii=inds[kk];sol=solve_remaining(A[ii],B,P,II[kk])
    if not sol:counts['schur_distance_failure']+=1
    for g in sol:survivors.append({'support_index':si,'changed_pairs':support,'labels':batch[ii],'normalized_gram':g});counts['global_gram_survivors']+=1
  print(name,'support',si+1,'/',len(supports),'assignments',counts['assignments'],'survivors',len(survivors),'near',near[0]['residual']if near else None,flush=True)
 return {'id':name,'basis':B,'fallback_bases':fallbacks,'outside_fixed_norm_walls':zero,'freed_parallel_pairs':row['parallel_pairs'],'tested_support_count':len(supports),'counts':counts,'closest_fixed_norm_proposals':near,'singular_basis_examples':singular,'survivors':survivors,'elapsed_seconds':time.time()-started}
def main():
 p=argparse.ArgumentParser();p.add_argument('--case');p.add_argument('--limit-supports',type=int);a=p.parse_args()
 data=json.loads(SOURCE.read_text());rows=[]
 order=['A16_14_1_C01','A14_10_1_C02','A14_8_1_C02']
 for name in order:
  if a.case and name!=a.case:continue
  rows.append(search(next(r for r in data['cases']if r['id']==name),a.limit_supports))
  out=OUT if not a.case else OUT.with_name(OUT.stem+'_'+name+'.json')
  out.write_text(json.dumps({'status':'BOUNDED_NUMERICAL_GLOBAL_GRAM_SURGERY','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'mutation_alphabet':[2,3,4,5,6],'cases':rows,'scope':'Minimal finite-edge covers, with all old parallel entries freed to strict dotted values. Near-singular chosen bases remain unresolved; failures are numerical, not an exact no-existence theorem.'},indent=2)+'\n')
if __name__=='__main__':main()
