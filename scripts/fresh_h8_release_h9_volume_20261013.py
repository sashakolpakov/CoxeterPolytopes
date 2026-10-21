#!/usr/bin/env python3
"""Exact replay of LP-proposed escapes in a released H8 seed and cap bank.

Coordinates y_i=-(x,r_i) are nonnegative. The section sum(y)=1 is a bounded
simplex cut by the cap inequalities, so a complete future-timelike vertex
census would prove compactness. This script only certifies found escapes.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,pair,matrix,
    inverse,rref,
)

ROOT=Path(__file__).resolve().parents[1]
DEFAULT=ROOT/'data/enumerations/fresh_h8_release_h9_caps_20261013.json'
real=lambda x:float(x[0])+(1+5**.5)/2*float(x[1])
encvec=lambda v:[enc(x)for x in v]

def gate(path,selected=None,tries=100):
    d=json.loads(path.read_text());G=matrix(d['seed_gram']);I=matrix(d['seed_gram_inverse'])
    A=[tuple(map(parse,row['seed_coefficients']))for row in d['caps']]
    if selected is None:selected=list(range(len(A)))
    A=[A[i]for i in selected]
    gn=np.array([[real(x)for x in row]for row in G]);qn=np.array([[real(x)for x in row]for row in I]);an=np.array([[real(x)for x in row]for row in A])
    _,vec=np.linalg.eigh(gn);v=vec[:,0]
    if v[0]<0:v=-v
    assert min(v)>0
    coefficients=tuple(scale(O,int(round(x*1000000)))for x in v)
    assert all(sign(x)>0 for x in coefficients)
    ycenter=tuple(neg(x)for x in mv(G,coefficients))
    assert all(sign(x)>0 for x in ycenter) and sign(pair(ycenter,I,ycenter))<0
    ycenter=tuple(div(x,total(ycenter))for x in ycenter)
    assert all(sign(dot(a,ycenter))>0 for a in A)
    rows=[tuple(O if i==j else Z for j in range(10))for i in range(10)]+A
    nums=np.array([[real(x)for x in row]for row in rows]);rng=np.random.default_rng(20261013)
    for trial in range(tries):
        p=linprog(rng.normal(size=10),A_ub=-an,b_ub=np.zeros(len(A)),A_eq=np.ones((1,10)),b_eq=[1],bounds=(0,None),method='highs')
        assert p.success,p.message
        if p.x@qn@p.x<=1e-9:continue
        ids=np.flatnonzero(abs(nums@p.x)<1e-8)
        active=[];system=[(O,)*10]
        for i in ids:
            if len(rref([*system,rows[i]])[1])>len(system):system.append(rows[i]);active.append(int(i))
            if len(system)==10:break
        if len(system)!=10:continue
        y=mv(inverse(system),(O,*([Z]*9)));q=pair(y,I,y)
        values=tuple(dot(a,y)for a in rows)
        if sign(q)<=0 or any(sign(x)<0 for x in values):continue
        assert sign(pair(y,I,ycenter))<0
        epsilon=F(1)
        while True:
            strict=tuple(add(x,scale(c,epsilon))for x,c in zip(y,ycenter));sq=pair(strict,I,strict)
            if sign(sq)>0:break
            epsilon/=2
        assert all(sign(dot(a,strict))>0 for a in rows) and sign(pair(strict,I,ycenter))<0
        return {'status':'EXACT_COMMON_SPACELIKE_ESCAPE','source':str(path.relative_to(ROOT)),
                'selected_cap_indices':selected,'seed_facet_count':10,'cap_count':len(A),
                'y_coordinates':encvec(y),'norm':enc(q),'constraint_values':encvec(values),
                'active_constraints':active,'strict_center_y':encvec(ycenter),
                'positive_center_seed_coefficients':encvec(coefficients),
                'strict_escape_y':encvec(strict),'strict_escape_norm':enc(sq),
                'strict_escape_center_coefficient':str(epsilon),
                'scope':'Exact escape satisfies every selected inequality, hence every subsystem. It rules out this finite cap dictionary, not caps outside it.'}
    return {'status':'UNDECIDED_NO_LP_ESCAPE','source':str(path.relative_to(ROOT)),
            'selected_cap_indices':selected,'tries':tries}

def main():
    p=argparse.ArgumentParser();p.add_argument('source',type=Path,nargs='?',default=DEFAULT)
    p.add_argument('--output',type=Path,default=ROOT/'data/enumerations/fresh_h8_release_h9_volume_20261013.json');a=p.parse_args()
    out=gate(a.source.resolve());a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'cap_count':out.get('cap_count'),'norm':out.get('norm'),'output':str(a.output)},indent=2))

if __name__=='__main__':main()
