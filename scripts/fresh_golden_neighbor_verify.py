#!/usr/bin/env python3
"""Standalone Fraction replay of a golden determinant-neighbor escape."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json

Z=(F(0),F(0))
def parse(x): return F(x['a']),F(x['b'])
def add(x,y): return x[0]+y[0],x[1]+y[1]
def scale(x,k): return x[0]*k,x[1]*k
def mul(x,y): return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def inv(x):
    n=x[0]*x[0]+x[0]*x[1]-x[1]*x[1]
    return (x[0]+x[1])/n,-x[1]/n
def sign(x):
    a=2*x[0]+x[1];b=x[1]
    if not b:return (a>0)-(a<0)
    if not a or (a>0)==(b>0):return (b>0)-(b<0)
    s=a*a-5*b*b
    return ((a>0)-(a<0))*((s>0)-(s<0))
EDGES=[(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(3,8)]
def pair(x,y,alpha):
    v=scale(mul(alpha,mul(x[0],y[0])),-1)
    for i in range(1,9):v=add(v,scale(mul(x[i],y[i]),2))
    v=add(v,mul(x[9],y[9]))
    for i,j in EDGES:v=add(v,scale(add(mul(x[i],y[j]),mul(x[j],y[i])),-1))
    return v
def verify(path):
    d=json.loads(path.read_text());alpha=parse(d['alpha'])
    walls=[tuple(map(parse,r)) for r in d['all_walls']]
    inside=tuple(map(parse,d['strict_interior']))
    assert sign(pair(inside,inside,alpha))<0
    norms=set(map(parse,d['norm_pool'])); previous=Z;bound=parse(d['completed_height'])
    for j,r in enumerate(walls):
        assert all(a.denominator==b.denominator==1 for a,b in r)
        q=pair(r,r,alpha);assert q in norms and sign(q)>0
        for i in range(10):
            e=[Z]*10;e[i]=(F(1),F(0))
            c=scale(mul(pair(r,e,alpha),inv(q)),2)
            assert all(v.denominator==1 for v in c)
        h=mul(mul(alpha,mul(r[0],r[0])),inv(q))
        assert sign(add(h,scale(previous,-1)))>=0 and sign(add(bound,scale(h,-1)))>=0
        previous=h
        assert sign(pair(inside,r,alpha))<0
        assert all(sign(pair(r,s,alpha))<=0 for s in walls[:j])
    w=d['hyperideal_witness'];r=tuple(map(parse,w['vector']))
    q=pair(r,r,alpha)
    assert q==parse(w['norm']) and sign(q)>0
    values=[pair(r,s,alpha) for s in walls]
    assert values==list(map(parse,w['wall_pairings'])) and all(sign(v)<=0 for v in values)
    print(f'PASS {path}: {len(walls)} integral acute mirrors, exact spacelike escape {q}.')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('artifact',type=Path)
    verify(p.parse_args().artifact)
