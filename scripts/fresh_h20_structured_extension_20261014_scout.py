#!/usr/bin/env python3
"""Continuous numerical scout of four paired PD18 complements; no exact gate."""
import json
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from itertools import combinations, product
import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_structured_extension_20261014'
F=lambda x:float(Fraction(x))

def intervals(a,b,c):
    eps=1e-8
    a,b,c=[0. if abs(x)<eps else float(x) for x in (a,b,c)]
    if not a:
        if not b:return [(0.,np.inf)] if c<=0 else []
        z=-c/b
        return [(0.,z)] if b>0 and z>0 else ([(max(0.,z),np.inf)] if b<0 else [])
    roots=np.roots([a,b,c]); rr=sorted(float(z.real) for z in roots if abs(z.imag)<1e-6 and z.real>0)
    bounds=[0.]+rr+[np.inf];out=[]
    for l,u in zip(bounds,bounds[1:]):
        p=(l+u)/2 if np.isfinite(u) else max(1.,2*l)
        if (a*p+b)*p+c<=eps:out.append((l,u))
    return out

def intersect(A,B):
    return [(max(a,c),min(b,d)) for a,b in A for c,d in B if max(a,c)<=min(b,d)]

def run():
    path=f'results/fresh_h20_coupled_lift_20261014.json';d=json.loads((ROOT/path).read_text())
    Q=np.array([[F(x) for x in row]for row in d['form']])
    N=np.array([[F(a)+np.sqrt(F(r['d']))*F(b)for a,b in zip(r['a'],r['b'])]for r in d['normals'][:60]])
    H=np.array(list(map(F,d['strict_center'])))
    old=json.loads((ROOT/'results/fresh_dimension_construction_routes_20261014.json').read_text())
    v=np.array([F(x['a'])for x in old['strict_escape']['boundary_point']]);E=v/5
    frames=[];ranges=[]
    for i,rec in enumerate(d['base_recipes']):
        h=np.array(list(map(F,rec['h0'])));w=np.array(list(map(F,rec['w'])))
        W=w/np.sqrt(w@Q@w);T=h-(h@Q@E)*E;T/=np.sqrt(-T@Q@T)
        if T@Q@H>0:T=-T
        for branch in (1,-1):
            EE=branch*E;coeff=np.array([N@Q@(EE+T),2*abs(N@Q@W),N@Q@(EE-T)]).T
            ints=[(0.,np.inf)]
            for aa,bb,cc in coeff:ints=intersect(ints,intervals(aa,bb,cc))
            ints=intersect(ints,intervals(H@Q@(EE+T),0,H@Q@(EE-T)))
            ranges.append(dict(support=i,branch=branch,t_intervals=ints))
            for l,u in ints:
                assert l>0
                # An infinite upper limit means an ideal asymptote. The finite
                # scout truncation below is explicit and never an exclusion.
                frames.append((i,branch,EE,T,W,(np.log(l),np.log(u) if np.isfinite(u) else np.log(l)+15)))
    def cap(frame,s,sgn):
        _,_,E,T,W,_=frame
        return np.cosh(s)*E+np.sinh(s)*T+sgn*W
    maximum=[]
    for fi,fj in combinations(range(len(frames)),2):
        x,y=frames[fi],frames[fj]
        if x[0]==y[0]:continue
        for si,sj in product((1,-1),repeat=2):
            fun=lambda z: -cap(x,z[0],si)@Q@cap(y,z[1],sj)/2
            z=differential_evolution(fun,[x[-1],y[-1]],seed=20261014,tol=1e-10,popsize=10,maxiter=120,polish=True)
            maximum.append(dict(supports=[x[0],y[0]],branches=[x[1],y[1]],signs=[si,sj],
                                max_normalized_pair=-float(z.fun),s=z.x.tolist(),
                                at_acute_boundary=[bool(abs(z.x[k]-[x,y][k][-1][0])<1e-5) for k in range(2)]))
    out=dict(status='NUMERICAL_CONTINUOUS_ACUTE_PAIRED_COMPLEMENT_SCOUT',
             source_sha256=sha256((ROOT/path).read_bytes()).hexdigest(),ranges=ranges,
             pair_maxima=maximum,
             scope='All real branches are parameterized. Acute ranges are numerical; unbounded log-parameter intervals are truncated at lower+15 for pair optimization. No Coxeter acceptance or all-real exclusion is claimed.')
    (ROOT/f'results/{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n')
    print('ranges',ranges)
    print('pair maxima',[(r['supports'],r['signs'],round(r['max_normalized_pair'],8))for r in maximum])
    return out

if __name__=='__main__':run()
