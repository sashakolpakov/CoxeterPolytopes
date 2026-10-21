#!/usr/bin/env python3
"""Exact scalar probe only; it supplies no Coxeter contact classification.

For each displayed prefix, bound every isolated-foot parameter 0<=u<=1
by covering that interval with 1000 rational closed cells. Square roots
are enclosed using integer square roots. Independent worst-cell errors
are added; this intentionally discards all geometric compatibility.
Termination is the stated finite loop. These scalar outputs are not a
dimension theorem or a claim that every actual cap has the tested form.
"""
from fractions import Fraction as F
from math import isqrt
import json

SCALE = 10**15


class I:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def box(x):
        return x if isinstance(x, I) else I(x)

    def __add__(self, other):
        o = I.box(other)
        return I(self.lo+o.lo, self.hi+o.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -I.box(other)

    def __rsub__(self, other):
        return I.box(other) + -self

    def __mul__(self, other):
        o = I.box(other)
        vals = [a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return I(min(vals),max(vals))

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = I.box(other)
        assert o.lo*o.hi > 0
        return self * I(1/o.hi,1/o.lo)

    def __rtruediv__(self, other):
        return I.box(other)/self

    def sqrt(self):
        assert self.lo >= 0
        a = isqrt(self.lo.numerator*SCALE*SCALE//self.lo.denominator)
        b = isqrt(self.hi.numerator*SCALE*SCALE//self.hi.denominator)
        return I(F(a,SCALE),F(b+1,SCALE))


def square(x):
    assert x.lo >= 0
    return I(x.lo*x.lo,x.hi*x.hi)


def root(n):
    return I(n).sqrt()


def outward(x):
    unit = 10**9
    lo = x.lo.numerator*unit//x.lo.denominator
    hi = -((-x.hi.numerator*unit)//x.hi.denominator)
    return dict(lo=str(F(lo,unit)),hi=str(F(hi,unit)))


def foot_error(beta, a, loss):
    best = None
    for j in range(1000):
        lo,hi = F(j,1000),F(j+1,1000)
        at = min(max(F(1,2),lo),hi)
        radical = I(at*(1-at)).sqrt()
        # 2u-u^2 increases on [0,1]; u(1-u) peaks at the clamped half.
        bound = 2*beta*loss - 2*beta*a*hi - 2*beta*radical
        bound += square(beta)*(2*lo-lo*lo)
        if best is None or bound.lo < best:
            best = bound.lo
    return best


def prefix(k, loss, final_loss=None):
    last = loss if final_loss is None else final_loss
    floors = {i:root(i)-I(last if i==k else loss)/root(2)
              for i in range(5,k+1)}
    slopes = [floors[5]/5]*5 + [floors[i]-floors[i-1]
                               for i in range(6,k+1)]
    assert all(w.lo>0 for w in slopes)
    assert all(slopes[i].lo>=slopes[i+1].hi for i in range(len(slopes)-1)
               if i>=4)
    energy = sum((square(w) for w in slopes),I(0))
    betas = {i:(slopes[i-1]-(slopes[i] if i<k else I(0)))/root(2)
             for i in range(5,k+1)}
    errors=[]
    for i,beta in betas.items():
        assert beta.lo>0
        q=root(2*i)
        a=1/(q+(square(q)-1).sqrt())
        error=foot_error(beta,a,last if i==k else loss)
        errors.append(dict(index=i,beta=outward(beta),error_lower=str(error)))
    total=energy+sum((F(x['error_lower']) for x in errors),F(0))
    return dict(rank=k,loss=str(loss),last_loss=str(last),
                reference_energy=outward(energy),
                all_isolated_foot_effective_energy=outward(total),errors=errors)


def main():
    cases=[(15,F(12,25),F(11,20)),(17,F(1,2),F(11,20)),
           (17,F(12,25),F(11,20)),(17,F(23,50),F(11,20))]
    rows=[prefix(*x) for x in cases]
    print(json.dumps(dict(status='EXACT_SCALAR_PROBE_NO_GEOMETRIC_CLASSIFICATION',
          cells_per_foot=1000,root_scale=SCALE,rows=rows),sort_keys=True))


if __name__=='__main__':
    main()
