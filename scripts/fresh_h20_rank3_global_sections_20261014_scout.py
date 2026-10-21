#!/usr/bin/env python3
"""Numerical proposal scout; no exact or exhaustive global claim."""
from pathlib import Path
import json,hashlib
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_rank3_global_sections_20261014'
LOCAL='results/fresh_h20_rank3_omitted_sections_20261014.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
mp.mp.dps=70

def load():
 p=json.loads((ROOT/LOCAL).read_text());r=json.loads((ROOT/RAW).read_text())['coxeter_input']
 Q=-np.array(r['basis_gram'],dtype=np.int64);W=np.array(r['wall_coordinates_as_columns_transposed'],dtype=np.int64);G=W@Q@W.T
 return p,Q,W,G

def values(row):
 c=[mp.mpf(str(F(x).numerator))/F(x).denominator for x in row['minimal_polynomial']]
 lo,hi=[mp.mpf(str(F(x).numerator))/F(x).denominator for x in row['root_interval']]
 root=lo if lo==hi else mp.findroot(lambda x:mp.polyval(c,x),(lo,hi),solver='secant')
 def val(cs):return float(mp.polyval([mp.mpf(str(F(x).numerator))/F(x).denominator for x in cs],root))
 return [val(x) for x in row['positive_kernel']],val(row['normal_square']),val(row['parameter_lower_magnitude'])

def run():
 p,Q,W,G=load()
 @lru_cache(None)
 def pd(ids):return np.linalg.eigvalsh(G[np.ix_(ids,ids)].astype(float))[0]>1e-8
 pairs=[];layout_cache={}
 for sid,s in enumerate(p['systems']):
  R=s['retained_walls'];fac=s['factor'];omit=s['omitted_wall']
  for low in R:
   key=(tuple(R),tuple(fac),omit,low)
   if key in layout_cache:
    pairs.append((sid,low,layout_cache[key]));continue
   out=[]
   for i,j in combinations(range(210),2):
    if omit in(i,j):continue
    A=tuple(sorted(set(R+[i,j])));B=tuple(sorted(set([x for x in fac if x!=low]+[i,j])))
    if pd(A) and pd(B):out.append((i,j,False))
   # Include direct old-facet endpoint face when all coefficients vanish.
   for i,j in combinations(range(210),2):
    if omit in(i,j):continue
    A=tuple(sorted(set(R+[i,j])));B=tuple(sorted(set([omit,i,j])))
    if pd(A) and pd(B):out.append((i,j,True))
   layout_cache[key]=out
   pairs.append((sid,low,out))
  if sid%20==0:print('profiles',sid,flush=True)
 pairdict={(sid,low):ps for sid,low,ps in pairs}
 result=[];us=np.linspace(.000001,.999999,101)
 for ix,row in enumerate(p['directions']):
  s=p['systems'][row['system']];R=s['retained_walls'];k,q,L=values(row);k=np.array(k);D=np.array([float(F(x)) for x in row['null_shift']]);a=k@G[R,:];b=-L*(D@Q@W.T)
  pp=a[None,:]+us[:,None]*b[None,:];diag=np.diag(G);n=diag[None,:]-pp*pp/q
  found=[];candidates=[]
  for i,j,direct in pairdict[row['system'],row['lower_bound_witness_wall']]:
   if direct:
    # e0(x1)=-L only if e0-LD is supported on omitted wall.
    fac=s['factor'];marks=np.linalg.svd(G[np.ix_(fac,fac)].astype(float))[2][-1];marks=abs(marks);marks/=marks[fac.index(s['omitted_wall'])]
    dscale=np.linalg.norm(D)/np.linalg.norm(marks@W[fac,:]);marks*=dscale
    if max(abs(k-L*np.array([marks[fac.index(x)]for x in R])))>1e-8:continue
   if min(n[:,i].min(),n[:,j].min())<=1e-10:continue
   pair=G[i,j]-pp[:,i]*pp[:,j]/q;cs=pair*pair/(n[:,i]*n[:,j])
   mode=None;gap=None
   if pair.min()>1e-11:mode='obtuse'
   elif pair.max()< -1e-11:
    for l,h in zip([0,.25,.5,.75],[.25,.5,.75,1]):
     if cs.min()>l+1e-11 and cs.max()<h-1e-11:mode='gap';gap=[l,h];break
   if mode:
    found={'ridge':[i,j],'mode':mode,'gap':gap,'direct_omitted_endpoint':direct,'sample_cosine_range':[float(cs.min()),float(cs.max())],'sample_pair_range':[float(pair.min()),float(pair.max())]};break
   if np.max(cs)<1-1e-7 and np.min(pair)<-1e-6 and np.max(cs)-np.min(cs)>1e-6:
    candidates.append({'ridge':[i,j],'direct_omitted_endpoint':direct,'sample_cosine_range':[float(cs.min()),float(cs.max())],'sample_pair_range':[float(pair.min()),float(pair.max())]})
  out={'row':ix,'cusp':s['cusp'],'factor':s['factor'],'omit':s['omitted_wall'],'target':s['target'],'degree':len(row['minimal_polynomial'])-1,'q':q,'L':L,'witness':found,'fallback_pairs':candidates[:20]}
  result.append(out);print(ix,s['target'],found or ('NO_UNIFORM',len(candidates)),flush=True)
 out={'status':'NUMERICAL_ACTUAL_RIDGE_INTERPOLATION_PROPOSALS','sources':{x:hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in[LOCAL,RAW,f'scripts/{TAG}_scout.py']},'rows':result,'uniform_proposals':sum(bool(r['witness']) for r in result),'scope':'Numerical proposals over101 interior samples only. Exact certificates required.'}
 (ROOT/'results'/f'{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n')
 print(out['uniform_proposals'],'of',len(result),flush=True)
if __name__=='__main__':run()
