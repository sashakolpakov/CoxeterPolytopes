#!/usr/bin/env python3
"""Target nonpolar Coxeter caps after releasing wall 4 of the H8 lift.

Each proposal is specified by exact pairings against ten independent seed
walls. A floating norm filter only proposes survivors; all accepted caps
are checked over Q(phi). This is a bounded pairing dictionary, not a full
reflection enumeration or a compactness theorem.
"""
import argparse
from itertools import combinations, product
import json
from pathlib import Path
import time
import numpy as np
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,pair,matrix,
    inverse,rref,
)

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h8_lift_h9_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h8_release_h9_caps_20261013.json'
real=lambda x:float(x[0])+(1+5**.5)/2*float(x[1])
encoded=lambda rows:[[enc(x)for x in row]for row in rows]

def run(max_support=5,output=OUTPUT,norms=(2,)):
    start=time.monotonic();raw=json.loads(SOURCE.read_text())['seeds'][0]
    assert raw['new_label']==3
    Q=matrix(raw['host_gram']);allroots=matrix(raw['wall_normals']);ids=[i for i in range(11)if i!=4]
    roots=tuple(allroots[i]for i in ids);G=tuple(tuple(pair(r,Q,s)for s in roots)for r in roots);I=inverse(G)
    assert len(rref(roots)[1])==10
    In=np.array([[real(x)for x in row]for row in I])
    magnitudes=(O,P,scale(O,2),add(O,P),scale(P,2))
    vals=np.array([real(x)for x in magnitudes]);caps=[];seen=set();counts=[]
    for k in range(1,max_support+1):
        choices=np.array(list(product(range(len(vals)),repeat=k)),dtype=int);p=vals[choices]
        proposed=0;tested=0
        for support in combinations(range(10),k):
            H=In[np.ix_(support,support)];squares=np.einsum('bi,ij,bj->b',p,H,p)
            close=np.zeros(len(squares),dtype=bool)
            for target in norms:close|=abs(squares-target)<1e-7
            for ci in np.flatnonzero(close):
                proposed+=1;v=[Z]*10
                for j,m in zip(support,choices[ci]):v[j]=neg(magnitudes[m])
                v=tuple(v);coeffs=mv(I,v);q=dot(v,coeffs)
                if q not in tuple(scale(O,t)for t in norms):continue
                r=tuple(total(mul(a,w[j])for a,w in zip(coeffs,roots))for j in range(10))
                assert pair(r,Q,r)==q
                assert tuple(pair(r,Q,w)for w in roots)==v
                assert all(sign(x)<=0 for x in v)
                for x in v:
                    ratio=div(scale(mul(x,x),2),q)
                    assert sign(sub(ratio,scale(O,4)))>=0 or ratio in (Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P))
                if r in seen:continue
                seen.add(r);tested+=1
                caps.append({'normal':[enc(x)for x in r],'norm':enc(q),'seed_pairings':[enc(x)for x in v],
                             'seed_coefficients':[enc(x)for x in coeffs],
                             'cut_hyperideal_seed_vertices':[ids[j]for j in range(10)if sign(I[j][j])>0 and sign(coeffs[j])<0],
                             'nonzero_pairings':k})
        row={'support_size':k,'patterns':len(list(combinations(range(10),k)))*len(choices),'numeric_norm_proposals':proposed,'exact_caps':tested}
        counts.append(row);print(row,flush=True)
    result={'schema':'released-H8-nonpolar-H9-caps-v1','source':str(SOURCE.relative_to(ROOT)),
            'seed_new_label':3,'deleted_original_wall':4,'retained_original_walls':ids,
            'host_gram':encoded(Q),'seed_normals':encoded(roots),'seed_gram':encoded(G),
            'seed_gram_inverse':encoded(I),'norm_pool':[enc(scale(O,t))for t in norms],
            'pairing_magnitudes':[enc(x)for x in magnitudes],
            'maximum_pairing_support':max_support,'counts':counts,'caps':caps,
            'scope':'Bounded pairing dictionary with exact accepted cap norm/angle/realization checks; numeric proposal filter is not an exhaustion proof.',
            'seconds':round(time.monotonic()-start,3)}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print({'caps':len(caps),'seconds':result['seconds']},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--max-support',type=int,default=5);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--norms',type=int,nargs='+',default=[2]);a=p.parse_args();run(a.max_support,a.output,tuple(a.norms))
