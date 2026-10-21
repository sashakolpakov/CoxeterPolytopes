"""Exact three-embedding integer alphabet in the cubic field of discriminant81."""
from fractions import Fraction as F
from functools import lru_cache
from coxeter_polytopes.exact.cubic9 import Cubic9 as K,THETA9 as T

DUAL=(K(F(11,9),F(2,9),-F(4,9)),K(F(2,9),F(2,9),-F(1,9)),K(-F(4,9),-F(1,9),F(2,9)))
INTERVALS=((F(1879385241571,10**12),F(1879385241572,10**12)),
           (F(-347296355334,10**12),F(-347296355333,10**12)),
           (F(-1532088886238,10**12),F(-1532088886237,10**12)))
assert [[(d*x).trace() for x in (K(1),T,T*T)] for d in DUAL]==[[1,0,0],[0,1,0],[0,0,1]]
for e,(lo,hi) in enumerate(INTERVALS):assert (T-lo).sign(e)>0 and (hi-T).sign(e)>0
BOUNDS=[]
for d in DUAL:
    row=[]
    for e in range(3):
        b=F(1,12)
        while (b-d).sign(e)<=0 or (b+d).sign(e)<=0:b+=F(1,12)
        row.append(b)
    BOUNDS.append(tuple(row))
BOUNDS=tuple(BOUNDS)

def ceiling(x):return -((-x.numerator)//x.denominator)
def interval_scale(a,lo,hi):return (a*lo,a*hi) if a>=0 else (a*hi,a*lo)

def coefficient_bounds(caps):
    if any(c.sign(e)<0 for e,c in enumerate(caps)):return None
    ms=[]
    for e,c in enumerate(caps):
        m=0 if not c else 1
        while (K(m*m)-c).sign(e)<0:m*=2
        ms.append(m)
    return tuple(ms),tuple(ceiling(sum(b*m for b,m in zip(row,ms))) for row in BOUNDS)

@lru_cache(None)
def alphabet(caps):
    bounds=coefficient_bounds(caps)
    if bounds is None:return ()
    ms,bs=bounds;out=[]
    for b in range(-bs[1],bs[1]+1):
        for c in range(-bs[2],bs[2]+1):
            low,high=-bs[0],bs[0]
            for e,(lo,hi) in enumerate(INTERVALS):
                blo,bhi=interval_scale(b,lo,hi)
                clo,chi=interval_scale(c,min(lo*lo,hi*hi),max(lo*lo,hi*hi))
                low=max(low,ceiling(-ms[e]-bhi-chi))
                bound=ms[e]-blo-clo;high=min(high,bound.numerator//bound.denominator)
                if low>high:break
            for a in range(low,high+1):
                x=K(a,b,c)
                if all((cap-x*x).sign(e)>=0 for e,cap in enumerate(caps)):out.append(x)
    return tuple(sorted(out))
