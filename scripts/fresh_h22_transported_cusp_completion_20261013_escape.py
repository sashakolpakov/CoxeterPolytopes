#!/usr/bin/env python3
"""Exact strict escape for all11,688,030 transported-repair H22 facets.

Only proposal generation is numerical. Saved replay uses integer vectors,
checked native integer contractions and exact integer square comparisons.
"""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from fresh_h22_incidence_binary_20261013_escape import check as old_check
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_transported_cusp_completion_20261013'
SOURCE='results/'+PRE+'_verify.json';FRAME='data/enumerations/'+PRE+'_frames.npz';RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
SCOUT='data/enumerations/'+PRE+'_volume_scout.json';OUT='results/'+PRE+'_escape.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def check(nums,den):
 assert len(nums)==23 and all(type(x)is int for x in nums)and type(den)is int and den>0
 prior=read(SOURCE);assert prior['status']=='PASS_EXACT_COMPLETE_TRANSPORTED_H22_CUSP_REPAIR'and all(sha(p)==h for p,h in prior['sources'].items())
 inherited=old_check(nums,den);assert inherited['facets']==5398110
 raw=read(RAW)['coxeter_input'];Q=[[-int(x)for x in row]for row in raw['basis_gram']]
 cov=[sum(a*b for a,b in zip(row,nums[:22],strict=True))for row in Q];assert max(map(abs,cov))<2**63
 square=F(sum(a*b for a,b in zip(nums[:22],cov,strict=True))+nums[22]*nums[22],den*den);assert square>0 and square==F(inherited['norm'])
 z=np.load(P/FRAME);bases=z['old_base_numerators'];E=z['six_delta_old'];height=z['height_numerators'];capden=int(z['denominator']);assert bases.shape==(241920,9,22)and E.shape==(241920,22)and capden==32760
 assert bases.dtype==E.dtype==height.dtype==np.int64
 specs=[(i,1)for i in range(5)]+[(i,s)for i in range(5,9)for s in(1,-1)];shifts=prior['shifts'];assert len(shifts)==13 and all(type(t)is int for t in shifts)
 cc=np.array(cov,dtype=np.int64);records=[]
 for k,((j,s),t)in enumerate(zip(specs,shifts,strict=True)):
  arithmetic_bound=int(abs(bases[:,j]).max())+abs(t)*(capden//6)*int(abs(E).max());assert arithmetic_bound<2**63
  rows=s*bases[:,j]+t*(capden//6)*E;zz=s*int(height[j])+2*t*capden
  assert 22*int(abs(rows).max())*max(map(abs,cov))<2**63
  values=rows@cc;idx=int(np.argmax(values));a=int(values[idx]);b=abs(zz*nums[22]);margin=a*a-11*b*b;assert a<0 and margin>0
  records.append({'profile':k,'maximum_old_integer_pairing':str(a),'maximizer':idx,'height_numerator':zz,'absolute_sqrt11_numerator':str(b),'positive_square_margin':str(margin),'denominator':capden*den,'covered_signed_walls':2*len(rows)})
 assert sum(r['covered_signed_walls']for r in records)==prior['new_caps']==6289920
 return {'status':'PASS_EXACT_TRANSPORTED_H22_ALL_FACETS_STRICT_ESCAPE','sources':{p:sha(p)for p in(SOURCE,FRAME,RAW)},'facets':prior['facets'],'numerator':nums,'denominator':den,'norm':str(square),'center_pairing':inherited['inherited_bank_escape_replay']['center_pairing'],'inherited_5398110_replay':inherited,'new_family_integer_maxima':records,'scope':'An exact rational future spacelike point is strict on every retained5,398,110 wall and all6,289,920 signed new walls. Integer maxima and A<0,A^2>11B^2 prove each full-family sign. The11,688,030-facet chamber therefore has infinite volume despite483,840 specified bounded cusp links. This is not an unrestricted construction obstruction.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 if a.verify:
  d=read(OUT);assert check(d['numerator'],d['denominator'])==d
 else:
  scout=read(SCOUT);x=np.array(scout['runs'][0]['point'],dtype=float);h=np.array(read('results/fresh_h22_incidence_core_20261013.json')['center'],dtype=float);x*=.999;x[:22]+=.001*h/np.sqrt(14);den=10**6;nums=np.rint(x*den).astype(np.int64).tolist();d=check(nums,den);(P/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['norm'])
