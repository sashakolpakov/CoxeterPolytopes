#!/usr/bin/env python3
"""Standalone Fraction arithmetic replay of a mixed-norm H9 escape.

Uses neither project field/matrix code nor the root search. Checks exact root
norms, reflection integrality, chronological height, strict timelike interior,
and a feasible spacelike point. Does not certify height-search exhaustion.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

Z=(F(0),F(0)); PHI=(F(0),F(1))
def parse(x):return F(x['a']),F(x['b'])
def plus(x,y):return x[0]+y[0],x[1]+y[1]
def times(x,k):return x[0]*k,x[1]*k
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def sign(x):
    a=2*x[0]+x[1]; b=x[1]
    if not b:return (a>0)-(a<0)
    if not a or (a>0)==(b>0):return (b>0)-(b<0)
    c=a*a-5*b*b
    return ((a>0)-(a<0))*((c>0)-(c<0))

EDGES=[(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(3,8)]
def pairing(x,y,a1):
    total=times(mul(PHI,mul(x[0],y[0])),-2)
    for i in range(1,9):total=plus(total,times(mul(x[i],y[i]),2))
    total=plus(total,times(mul(x[9],y[9]),a1))
    for i,j in EDGES:total=plus(total,times(plus(mul(x[i],y[j]),mul(x[j],y[i])),-1))
    return total

def verify(path):
    proof=json.loads(path.read_text())
    source=json.loads(Path(proof['source']).read_text())
    assert source['dimension']==9
    a1=source.get('a1_norm',2)
    walls=[tuple(map(parse,r)) for r in source['all_walls']]
    witness=tuple(map(parse,proof['vector']))
    interior=tuple(map(parse,source['strict_interior']))
    H=parse(source['height_bound'])
    assert a1 in (1,2) and len(walls)==source['wall_count']==proof['wall_count']
    assert len(witness)==len(interior)==10 and witness[0][0]>0 and witness[0][1]==0
    assert sign(pairing(interior,interior,a1))<0 and sign(interior[0])>0
    previous=Z
    norms={1,2} if a1==1 else {2,4}
    for r in walls:
        assert len(r)==10 and all(a.denominator==b.denominator==1 for a,b in r)
        q=pairing(r,r,a1)
        assert q[1]==0 and q[0] in norms
        for i in range(10):
            e=[Z]*10;e[i]=(F(1),F(0))
            c=times(pairing(r,e,a1),F(2)/q[0])
            assert c[0].denominator==c[1].denominator==1
        height=times(mul(PHI,mul(r[0],r[0])),F(2)/q[0])
        assert sign(plus(height,times(previous,-1)))>=0
        assert sign(plus(H,times(height,-1)))>=0
        previous=height
        assert sign(pairing(interior,r,a1))<0
    pairs=[pairing(witness,r,a1) for r in walls]
    assert all(sign(p)<=0 for p in pairs)
    assert pairs==list(map(parse,proof['wall_pairings']))
    norm=pairing(witness,witness,a1)
    assert norm==parse(proof['norm']) and sign(norm)>0
    print(f'PASS: {len(walls)} integral mirrors, strict interior, exact feasible spacelike point of norm {norm}.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('witness',type=Path)
    verify(p.parse_args().witness)
