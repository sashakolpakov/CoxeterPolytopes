"""Exact thin-box enumeration in the cubic-seven integer ring."""
from fractions import Fraction as F
from functools import lru_cache
from coxeter_polytopes.exact.cubic7 import Cubic7 as K,THETA7 as T

DUAL=((K(1,-F(1,7),-F(3,7))),K(-F(1,7),F(2,7),F(1,7)),K(-F(3,7),F(1,7),F(2,7)))
BOUNDS=((F(1,6),F(1),F(1,7)),(F(4,9),F(1,4),F(1,5)),(F(1,5),F(4,9),F(1,4)))
INTERVALS=((F(1246979603717,10**12),F(1246979603718,10**12)),
           (F(-445041867913,10**12),F(-445041867912,10**12)),
           (F(-1801937735805,10**12),F(-1801937735804,10**12)))
for e,(lo,hi) in enumerate(INTERVALS):
    assert (T-lo).sign(e)>0 and (hi-T).sign(e)>0
for d,bs in zip(DUAL,BOUNDS):
    assert all((K(b)-d).sign(e)>0 and (K(b)+d).sign(e)>0 for e,b in enumerate(bs))

def ceiling(x):return -((-x.numerator)//x.denominator)
def interval_scale(a,lo,hi):return (a*lo,a*hi) if a>=0 else (a*hi,a*lo)

@lru_cache(None)
def thin_alphabet(caps):
    if any(c.sign(e)<0 for e,c in enumerate(caps)):return ()
    ms=[]
    for e,c in enumerate(caps):
        m=0 if not c else 1
        while (K(m*m)-c).sign(e)<0:m*=2
        ms.append(m)
    bounds=[ceiling(sum(b*m for b,m in zip(row,ms))) for row in BOUNDS]
    out=[]
    for b in range(-bounds[1],bounds[1]+1):
        for c in range(-bounds[2],bounds[2]+1):
            low,high=-bounds[0],bounds[0]
            for e,(lo,hi) in enumerate(INTERVALS):
                blo,bhi=interval_scale(b,lo,hi)
                clo,chi=interval_scale(c,min(lo*lo,hi*hi),max(lo*lo,hi*hi))
                low=max(low,ceiling(-ms[e]-bhi-chi))
                bound=ms[e]-blo-clo
                high=min(high,bound.numerator//bound.denominator)
                if low>high:break
            for a in range(low,high+1):
                x=K(a,b,c)
                if all((cap-x*x).sign(e)>=0 for e,cap in enumerate(caps)):out.append(x)
    return tuple(sorted(out))
