#!/usr/bin/env python3
"""Numerical LP proposals for exact all-real center/facet section gates."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json, sys
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import fresh_h20_global_fraction_verify_20261013 as S
Q,W,G,h=S.source();Q=np.array(Q,dtype=np.int64);W=np.array(W,dtype=np.int64);G=np.array(G,dtype=np.int64);h=np.array(h,dtype=np.int64);C=W@Q;hc=h@Q
out=[]
for i in [0,next(i for i in range(210)if G[i,i]==4)]:
 lp=linprog(C[i],A_ub=C,b_ub=np.zeros(210),A_eq=hc[None],b_eq=[-1],bounds=[(None,None)]*22,method='highs');assert lp.success
 x=[F(float(v)).limit_denominator(1000000)for v in lp.x]
 weights=[F(float(-v)).limit_denominator(1000000)for v in lp.ineqlin.marginals];T=F(float(-lp.fun)).limit_denominator(1000000)
 dot=lambda a,b:sum(x*y for x,y in zip(a,b,strict=True))
 assert dot(hc,x)==-1 and dot(C[i],x)==-T and all(dot(c,x)<=0 for c in C)
 assert min(weights)>=0 and all(T*h[j]-W[i,j]==sum(weights[k]*W[k,j]for k in range(210))for j in range(22))
 cand=[j for j in range(210)if G[j,j]==2 and G[i,j]==0 and dot(C[j],x)==0]
 j,k=next((j,k)for j,k in combinations(cand,2)if G[j,k]==0)
 out.append({'wall':i,'T':str(T),'maximum_point':list(map(str,x)),'dual_weights':[[k,str(v)]for k,v in enumerate(weights)if v],'ridge_pair':[j,k]})
 print('exact proposal',i,T,[j,k],flush=True)
p=ROOT/'results/fresh_h20_borcherds_sections_20261013_facet_recipes.json';p.write_text(json.dumps(out,indent=2)+'\n')
