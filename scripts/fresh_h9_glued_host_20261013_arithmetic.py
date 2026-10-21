#!/usr/bin/env python3
"""Exact arithmetic in Q(T), T=2cos(pi/8), T^4-4T^2+2=0."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import ceil,isqrt,sqrt
from coxeter_polytopes.exact.silver import Silver as S

@dataclass(frozen=True)
class K:
    a:F=F(0);b:F=F(0);c:F=F(0);d:F=F(0)
    def __post_init__(self):
        for key in('a','b','c','d'):object.__setattr__(self,key,F(getattr(self,key)))
    @property
    def coeffs(self):return self.a,self.b,self.c,self.d
    @staticmethod
    def coerce(x):return x if isinstance(x,K)else K(x)
    def __add__(self,y):
        y=K.coerce(y);return K(*(a+b for a,b in zip(self.coeffs,y.coeffs)))
    __radd__=__add__
    def __neg__(self):return K(*(-x for x in self.coeffs))
    def __sub__(self,y):return self+-K.coerce(y)
    def __rsub__(self,y):return K.coerce(y)+-self
    def __mul__(self,y):
        y=K.coerce(y);v=[F(0)]*7
        for i,a in enumerate(self.coeffs):
            for j,b in enumerate(y.coeffs):v[i+j]+=a*b
        for j in(6,5,4):v[j-2]+=4*v[j];v[j-4]-=2*v[j]
        return K(*v[:4])
    __rmul__=__mul__
    def tower(self):return S(self.a+2*self.c,self.c),S(self.b+2*self.d,self.d)
    def inverse(self):
        u,v=self.tower();den=u*u-S(2,1)*v*v
        u,v=u/den,-v/den
        return K(u.a-2*u.b,v.a-2*v.b,u.b,v.b)
    def __truediv__(self,y):return self*K.coerce(y).inverse()
    def __rtruediv__(self,y):return K.coerce(y)*self.inverse()
    def __pow__(self,n):
        if n<0:return self.inverse()**-n
        out=K(1);x=self
        while n:
            if n%2:out*=x
            n//=2;x*=x
        return out
    def __bool__(self):return any(self.coeffs)
    def __eq__(self,y):
        if isinstance(y,(int,F)):return self.coeffs==(F(y),F(0),F(0),F(0))
        return isinstance(y,K)and self.coeffs==y.coeffs
    def __hash__(self):return hash(self.a)if not any((self.b,self.c,self.d))else hash(self.coeffs)
    def sign(self,embedding=0):
        assert embedding in range(4)
        u,v=self.tower();e=1 if embedding in(1,2)else 0
        su=u.sign(e);sv=v.sign(e)*(1 if embedding in(0,1)else-1)
        if not su:return sv
        if not sv or su==sv:return su
        return su*(u*u-S(2,1)*v*v).sign(e)
    def __lt__(self,y):return(self-y).sign()<0
    def __le__(self,y):return(self-y).sign()<=0
    def __gt__(self,y):return(self-y).sign()>0
    def __ge__(self,y):return(self-y).sign()>=0
    @property
    def is_integral(self):return all(x.denominator==1 for x in self.coeffs)
    def trace(self):return 4*self.a+8*self.c
    def norm(self):
        u,v=self.tower();return(u*u-S(2,1)*v*v).norm()
    def real(self,e=0):
        t=(1 if e in(0,1)else-1)*sqrt(2+(1 if e in(0,3)else-1)*sqrt(2))
        return sum(float(x)*t**i for i,x in enumerate(self.coeffs))
    def to_json(self):return list(map(str,self.coeffs))
    @classmethod
    def from_json(cls,x):return cls(*map(F,x))
    def __str__(self):return ','.join(map(str,self.coeffs))

Z=K(0);O=K(1);T=K(0,1);SQRT2=T*T-2;ALPHA=T-1
PI=2+T;NORMS=(O,PI,T*T,2/PI,K(2))

def inverse(g):
    n=len(g);a=[list(row)+[K(i==j)for j in range(n)]for i,row in enumerate(g)]
    for i in range(n):
        p=next(j for j in range(i,n)if a[j][i]);a[i],a[p]=a[p],a[i]
        q=a[i][i];a[i]=[x/q for x in a[i]]
        for j in range(n):
            if j!=i:
                q=a[j][i];a[j]=[x-q*y for x,y in zip(a[j],a[i])]
    return tuple(tuple(row[n:])for row in a)

def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),Z)
def mv(g,v):return tuple(dot(row,v)for row in g)
def pair(a,g,b):return dot(a,mv(g,b))
def enc(v):return[x.to_json()for x in v]
def menc(g):return[enc(row)for row in g]
def le_all(a,b):return all((b-a).sign(e)>=0 for e in range(4))
def ray_key(v):
    first=next(x for x in v if x);return first.sign(),tuple(x/first for x in v)

@lru_cache(None)
def trace_ball(bound):
    """Trace(x^2)=4(a+2c)^2+8(b+3d)^2+8c^2+8d^2."""
    result=[]
    for A in range(-isqrt(bound//4),isqrt(bound//4)+1):
        n=bound-4*A*A
        for B in range(-isqrt(n//8),isqrt(n//8)+1):
            m=n-8*B*B
            for C in range(-isqrt(m//8),isqrt(m//8)+1):
                k=m-8*C*C
                for D in range(-isqrt(k//8),isqrt(k//8)+1):
                    x=K(A-2*C,B-3*D,C,D)
                    if x.sign()>=0:result.append((x,x*x))
    return tuple(result)

def ceil_embedding(x,e):
    n=ceil(x.real(e))
    while(x-n).sign(e)>0:n+=1
    while(x-(n-1)).sign(e)<=0:n-=1
    return n

@lru_cache(None)
def alphabet(caps):
    if any(c.sign(e)<0 for e,c in enumerate(caps)):return()
    bound=ceil(caps[0].trace())if all(c==caps[0]for c in caps)else sum(ceil_embedding(c,e)for e,c in enumerate(caps))
    return tuple(sorted(x for x,sq in trace_ball(bound)if all((c-sq).sign(e)>=0 for e,c in enumerate(caps))))

F4=tuple(tuple(K(2)if i==j else -SQRT2 if{i,j}=={1,2}else -O if abs(i-j)==1 else Z for j in range(4))for i in range(4))
FI=inverse(F4)
SCHUR=tuple(inverse(tuple(row[:k]for row in F4[:k]))for k in range(1,5))
ALLOWED={Z,O,K(2),K(3),2+SQRT2,2+T}
