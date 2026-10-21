#!/usr/bin/env python3
"""Exact, bounded norm-1/2 searches in new quadratic diagonal H9 hosts.

This is NOT an exhaustive search of reflective lattices. It exhausts the
norm-1 and norm-2 vectors in each stated height bound, modulo the B_n
stabilizer, then greedily selects acute walls in height order. Other norm
classes are not searched. A feasible spacelike point certifies infinite
volume only of the displayed prefix, never nonreflectivity of the host.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
from functools import cmp_to_key, lru_cache
from itertools import product
from math import isqrt, sqrt
from pathlib import Path
import json
import time

Z = (0, 0)
O = (1, 0)
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def sub(x,y): return add(x,neg(y))
def scale(x,k): return (x[0]*k,x[1]*k)

@dataclass(frozen=True)
class Field:
    d: int
    @property
    def half(self): return self.d % 4 == 1
    def mul(self,x,y):
        a,b=x; c,e=y
        return (a*c+b*e*((self.d-1)//4 if self.half else self.d),a*e+b*c+(b*e if self.half else 0))
    def sq(self,x): return self.mul(x,x)
    def sign(self,x,place=0):
        a,b=x
        if self.half: a,b=2*a+b,b
        if place: b=-b
        if b==0: return (a>0)-(a<0)
        if a==0 or (a>0)==(b>0): return (b>0)-(b<0)
        c=a*a-self.d*b*b
        return ((c>0)-(c<0))*((a>0)-(a<0))
    def val(self,x,place=0):
        w=(1+(-1 if place else 1)*sqrt(self.d))/2 if self.half else (-1 if place else 1)*sqrt(self.d)
        return x[0]+x[1]*w
    def bounded(self,U,V,predicate):
        # x0-x1=b*sqrt(d) in the half basis, else 2*b*sqrt(d).
        gap=isqrt(self.d)*(1 if self.half else 2)
        B=(U+V+gap-1)//gap
        A=V+B*(isqrt(self.d)+1)
        return [x for a in range(-A,A+1) for b in range(-B,B+1)
                if predicate(x:=(a,b))]
    def controls(self,alpha,q,H):
        U=1
        while self.sign(sub(scale(alpha,U*U),(q*H,0)))<0: U*=2
        V=1
        while self.sign(sub(scale(neg(alpha),V*V),(q,0)),1)<0: V*=2
        return self.bounded(U,V,lambda x: self.sign(x)>0 and self.sign(sub((q*H,0),self.mul(alpha,self.sq(x))))>=0 and self.sign(add((q,0),self.mul(alpha,self.sq(x))),1)>=0)
    def spatial(self,target):
        bounds=[]
        for p in (0,1):
            U=1
            while self.sign(sub((U*U,0),target),p)<0: U*=2
            bounds.append(U)
        return sorted(self.bounded(*bounds,lambda x:self.sign(x)>=0 and all(self.sign(sub(target,self.sq(x)),p)>=0 for p in (0,1))),key=cmp_to_key(lambda x,y:self.sign(sub(x,y))),reverse=True)


def search(d,alpha,n,H):
    started=time.monotonic(); K=Field(d)
    assert K.sign(alpha)>0 and K.sign(alpha,1)<0
    roots=[]; controls=0; raw=0
    for q in (1,2):
        for x in K.controls(alpha,q,H):
            controls+=1
            target=add((q,0),K.mul(alpha,K.sq(x)))
            alphabet=K.spatial(target); squares=[K.sq(y) for y in alphabet]
            @lru_cache(None)
            def partitions(start,slots,a,b):
                rem=(a,b)
                if rem==Z: return ((Z,)*slots,)
                if slots==0: return ()
                if any(K.sign(rem,p)<0 for p in (0,1)): return ()
                out=[]
                for j in range(start,len(alphabet)):
                    y=alphabet[j]
                    if y==Z: continue
                    sq=squares[j]
                    if K.sign(sub(scale(sq,slots),rem))<0: break
                    rest=sub(rem,sq)
                    if any(K.sign(rest,p)<0 for p in (0,1)): continue
                    for tail in partitions(j,slots-1,*rest): out.append((y,)+tail)
                return tuple(out)
            ys=partitions(0,n,*target)
            raw+=len(ys)
            roots.extend((q,x,y) for y in ys)
    def height_cmp(r,s): return K.sign(sub(scale(K.sq(r[1]),s[0]),scale(K.sq(s[1]),r[0])))
    roots.sort(key=cmp_to_key(height_cmp))
    def pair(r,s):
        val=neg(K.mul(alpha,K.mul(r[1],s[1])))
        for a,b in zip(r[2],s[2]): val=add(val,K.mul(a,b))
        return val
    walls=[]
    for r in roots:
        if all(K.sign(pair(r,s))<=0 for s in walls): walls.append(r)
    # Verify all root equations and integral reflection coefficients.
    for q,x,y in walls:
        assert pair((q,x,y),(q,x,y))==(q,0)
        coeff=[scale(K.mul(alpha,x),-2),*(scale(t,2) for t in y)]
        assert all(a%q==0 and b%q==0 for a,b in coeff)
    # Exact feasible test points: y=t*(1^k,0^(n-k)), x0=1.
    # t=min(alpha*x/sum(y[:k])) represented by a numerator/denominator pair.
    witnesses=[]
    for k in range(1,n+1):
        limits=[]
        for q,x,y in walls:
            denom=Z
            for a in y[:k]: denom=add(denom,a)
            if denom!=Z: limits.append((K.mul(alpha,x),denom))
        if not limits:
            # Choose t integer with k*t^2>alpha.
            t=1
            while K.sign(sub((k*t*t,0),alpha))<=0:t*=2
            num,den=(t,0),O
        else:
            num,den=min(limits,key=cmp_to_key(lambda a,b:K.sign(sub(K.mul(a[0],b[1]),K.mul(b[0],a[1])))))
        norm_numerator=sub(scale(K.sq(num),k),K.mul(alpha,K.sq(den)))
        if K.sign(norm_numerator)>0:
            for q,x,y in walls:
                sy=Z
                for a in y[:k]: sy=add(sy,a)
                assert K.sign(sub(K.mul(sy,num),K.mul(K.mul(alpha,x),den)))<=0
            witnesses.append({'k':k,'t_numerator':num,'t_denominator':den,'norm_numerator':norm_numerator})
    return {'d':d,'integral_basis':'(1+sqrt(d))/2' if K.half else 'sqrt(d)', 'alpha':alpha,'dimension':n,'height_bound':H,'norms':[1,2], 'control_count':controls,'dominant_root_count':raw,'positive_height_wall_count':len(walls),'wall_count':len(walls)+n,'walls':[{'norm':q,'x':x,'spatial':y} for q,x,y in walls], 'spacelike_witnesses':witnesses,'prefix_finite_volume':'false' if witnesses else 'undecided','host_reflectivity':'undecided','elapsed_seconds':round(time.monotonic()-started,4)}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--height',type=int,default=300)
    parser.add_argument('--dimension',type=int,default=9)
    parser.add_argument('--fields',type=int,nargs='+',default=[2,3,5,6,7,10,11,13,17,19,21,29,33])
    parser.add_argument('--output',type=Path,default=Path('data/enumerations/quadratic_diagonal_h9_2026_10_13.json'))
    args=parser.parse_args(); rows=[]
    for d in args.fields:
        K=Field(d)
        for a,b in product(range(-2,4),[1,2]):
            alpha=(a,b)
            if K.sign(alpha)<=0 or K.sign(alpha,1)>=0:continue
            row=search(d,alpha,args.dimension,args.height); rows.append(row)
            print(json.dumps({k:v for k,v in row.items() if k not in ('walls','spacelike_witnesses')}),flush=True)
            args.output.write_text(json.dumps({'scope':__doc__,'records':rows},indent=2)+'\n')
if __name__=='__main__': main()
