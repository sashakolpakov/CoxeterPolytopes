#!/usr/bin/env python3
"""Standalone integer-only replay of 132 explicit infinite-volume prefixes.

No search imports, no floating point, and no reliance on root-list exhaustion.
Each verdict follows from an explicit feasible spacelike point. This does NOT
exclude later mirrors or prove nonreflectivity of any complete host.
"""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]

def verify_all():
    data=json.loads((ROOT/'data/enumerations/quadratic_diagonal_h9_2026_10_13.json').read_text())['records']
    extra=json.loads((ROOT/'data/enumerations/quadratic_diagonal_h9_geometry_2026_10_13.json').read_text())
    audits={(r['d'],tuple(r['alpha'])):r for r in extra}
    count=0
    for row in data:
        d=row['d']; half=d%4==1; alpha=tuple(row['alpha']); n=row['dimension']; H=row['height_bound']
        def mul(x,y):
            a,b=x;c,e=y
            return (a*c+b*e*((d-1)//4 if half else d),a*e+b*c+(b*e if half else 0))
        def add(x,y):return (x[0]+y[0],x[1]+y[1])
        def neg(x):return (-x[0],-x[1])
        def times(x,k):return (x[0]*k,x[1]*k)
        def sign(x,conj=False):
            # Compare rational and radical parts after clearing the denominator.
            A=(2*x[0]+x[1]) if half else x[0]
            B=x[1]*(-1 if conj else 1)
            if B==0:return (A>0)-(A<0)
            if A>=0 and B>=0:return 1
            if A<=0 and B<=0:return -1
            if A>0:return 1 if A*A>d*B*B else -1
            return 1 if A*A<d*B*B else -1
        assert sign(alpha)>0 and sign(alpha,True)<0
        def dot(y,z):
            v=(0,0)
            for a,b in zip(y,z):v=add(v,mul(a,b))
            return v
        walls=row['walls']
        for i,r in enumerate(walls):
            q=r['norm'];x=r['x'];ys=r['spatial']
            assert q in (1,2) and len(ys)==n and sign(x)>0
            assert sign(ys[-1])>=0 and all(sign(add(a,neg(b)))>=0 for a,b in zip(ys,ys[1:]))
            ax2=mul(alpha,mul(x,x))
            assert add(dot(ys,ys),neg(ax2))==(q,0)
            assert sign(add((q*H,0),neg(ax2)))>=0
            # 2/q is an integer, so all reflection coefficients are integral.
            for s in walls[:i]:assert sign(add(dot(ys,s['spatial']),neg(mul(alpha,mul(x,s['x'])))))<=0
        if row['spacelike_witnesses']:
            w=row['spacelike_witnesses'][0];k=w['k']; num=w['t_numerator'];den=w['t_denominator']
            assert sign(num)>=0 and sign(den)>0 and 1<=k<=n
            ys=[num]*k+[(0,0)]*(n-k)
        else:
            w=audits[(d,alpha)]['exact_feasible_spacelike_point']
            den=(w['denominator'],0);ys=[(a,0) for a in w['spatial_numerators']]
            assert sign(den)>0 and len(ys)==n
            assert sign(ys[-1])>=0 and all(sign(add(a,neg(b)))>=0 for a,b in zip(ys,ys[1:]))
        assert sign(add(dot(ys,ys),neg(mul(alpha,mul(den,den)))))>0
        for r in walls:
            assert sign(add(dot(ys,r['spatial']),neg(mul(alpha,mul(r['x'],den)))))<=0
        count+=1
    assert len(audits)==sum(not r['spacelike_witnesses'] for r in data)
    return count

if __name__=='__main__': print(f'PASS: {verify_all()} exact infinite-volume prefixes; complete hosts remain undecided.')
