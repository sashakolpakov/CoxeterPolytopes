#!/usr/bin/env python3
"""Exact saved strict-escape gate for the coherent94-wall package."""
import argparse,json
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from fresh_h20_targeted_escape_20261013_volume import check,read,sha,ROOT,E,parse,vec,Extension,B
from fresh_h20_edge_family_completion_20261013_volume import setup as prior_setup
from fresh_h20_targeted_escape_20261013 import ev
CAPS='data/enumerations/fresh_h20_edge_family_completion_20261013_four_exact.json'
OUT='data/enumerations/fresh_h20_edge_family_completion_20261013_four_volume.json'

def setup():
    Qn,Rn,hn,Q,h,hc,roots=prior_setup();data=read(CAPS)
    for c in data['new_caps']:
        d=ev(c['field']['radicand']);r=np.array([ev(x[0])+np.sqrt(d)*ev(x[1])for x in c['normal']]);Rn=np.vstack([Rn,r])
        f=Extension(B,parse(c['field']['radicand']));roots.append((f,tuple(map(f.parse,c['normal']))))
    assert len(roots)==94
    return Qn,Rn,hn,Q,h,hc,roots

def run(attempts=200):
    Qn,Rn,hn,Q,h,hc,roots=setup();ns=np.einsum('ij,jk,ik->i',Rn,Qn,Rn);A=Rn@Qn/np.sqrt(ns[:,None]);center=hn/np.sqrt(-hn@Qn@hn)
    rng=np.random.default_rng(20262094);trials=[];cert=None
    for i in range(attempts):
        lp=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(center@Qn)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
        if not lp.success:continue
        q=float(lp.x@Qn@lp.x);ids=np.flatnonzero(abs(A@lp.x)<1e-7).tolist();trials.append({'attempt':i,'numeric_norm':q,'active_indices':ids})
        if q<=1e-7:continue
        scale=max(abs(lp.x));base=tuple(E(F(round(x/scale*10**9),10**9))for x in lp.x)
        for power in range(12,61):
            w=tuple(x+E(F(1,2**power))*y for x,y in zip(base,h,strict=True))
            try:cert=check(w,Q,hc,roots)
            except AssertionError:continue
            cert.update({'source_attempt':i,'perturbation_power':power,'numeric_active_indices':ids});break
        if cert:break
    return {'status':'EXACT_STRICT_ESCAPE_FOR_94_WALL_H20'if cert else'BOUNDED_94_WALL_VOLUME_UNDECIDED','source':CAPS,'source_sha256':sha(CAPS),
            'attempt_limit':attempts,'trials':trials,'certificate':cert,'scope':'Numerical LP proposals, exact rational/Qsqrt6 strict future spacelike gate against all94 quadratic-field walls. Saved replay needs no LP.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    if a.verify:
        out=read(OUT);assert out['source_sha256']==sha(CAPS);_,_,_,Q,h,hc,roots=setup();c=out['certificate'];assert c
        actual=check(vec(c['vector']),Q,hc,roots);assert all(c[k]==v for k,v in actual.items())
    else:out=run();(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],None if out['certificate']is None else out['certificate']['norm'])
