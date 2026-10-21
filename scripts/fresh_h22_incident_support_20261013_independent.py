#!/usr/bin/env python3
"""Independent Fraction coefficient replay of the orthogonal-support identity.

No producer arithmetic, optimizer, matrix solve or symbolic package is used.
The only irrationality is the positive square root of eleven.
"""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1]
SOURCE='results/fresh_h22_incident_support_20261013_orthogonal.json'
OUT='results/fresh_h22_incident_support_20261013_independent.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]+11*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def rat(x):return(F(x),F(0))
def run():
 d=read(SOURCE);assert d['status']=='PASS_EXACT_ORTHOGONAL_INCIDENT_ACUTE_NORMALS_MISS_FIXED_ESCAPE'
 assert all(sha(p)==h for p,h in d['sources'].items())
 raw=read('data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json')['coxeter_input']
 W=raw['wall_coordinates_as_columns_transposed'];Q=[[-F(v)for v in r]+[F(0)]for r in raw['basis_gram']]+[[F(0)]*22+[F(11)]]
 local=read('results/fresh_h22_lower_new_cusp_20261013.json')['representatives'][0]
 delta=list(map(F,local['null_ray']));free=set(local['incident_old_wall_indices'])
 assert sorted(free)==d['free_incident_old_coefficients']and len(free)==17
 def dot(x,y):return sum((a*b for a,b in zip(x,y,strict=True)),F(0))
 dc=[dot(row,delta)for row in Q];assert dot(delta,dc)==0
 assert {i for i,w in enumerate(W)if dot(list(map(F,w))+[F(0)],dc)==0}==free
 ex=read('results/fresh_h22_new_cusp_completion_20261013_lower_escape.json')
 target=[rat(F(x,ex['denominator']))for x in ex['numerator'][:-1]]+[(F(0),F(ex['numerator'][-1],11*ex['denominator']))]
 with np.load(P/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz')as z:
  nulls=z['old_null_projections'];heights=z['four_heights']
 total=[rat(0)for _ in range(23)];positive=[]
 lo=F(3316624790,10**9);hi=F(3316624791,10**9);assert lo*lo<11<hi*hi
 for term in d['terms']:
  name,i,s=term['source'];a,b=map(F,term['coefficient']);co=(a,b)
  if name=='old':r=list(map(rat,W[i]))+[rat(0)]
  elif name=='delta':r=list(map(rat,delta))
  else:
   assert name=='old_null'and i in(95611,95665)and s==1
   old=[F(int(x))for x in nulls[i]];assert dot(old,[dot(row[:22],old)for row in Q[:22]])==0
   r=list(map(rat,old))+[(F(0),F(int(heights[i]),44))]
  if name!='delta'and not(name=='old'and i in free):
   lower=a+b*(lo if b>=0 else hi);assert lower>0
   positive.append({'wall':term['source'],'coefficient_lower_bound':str(lower)})
  total=[add(x,mul(co,y))for x,y in zip(total,r,strict=True)]
 assert total==target and len(positive)==5
 assert {tuple(t['wall'])for t in positive}=={('old',6,1),('old',8,1),('old',38,1),('old_null',95611,1),('old_null',95665,1)}
 return {'status':'PASS_INDEPENDENT_EXACT_ORTHOGONAL_CUSP_SUPPORT_IDENTITY','source':SOURCE,'source_sha256':sha(SOURCE),'identity_coordinates':23,'free_old_wall_count':17,'positive_terms':positive,'scope':'Independent coordinate identity and sign proof. It excludes cutting the fixed escape using normals perpendicular to delta4 and the17 old incident walls, with no center-label or norm restriction. It does not exclude other incident profiles or other cusps.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
