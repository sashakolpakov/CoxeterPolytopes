"""Exact inherited radical-tower arithmetic for the coupled-role construction."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
RAD=[1,2,79,158,11,22,869,1738,3,6,237,474,33,66,2607,5214]
def parse(x,k):return tuple(map(F,x))if k==0 else tuple(parse(y,k-1)for y in x)
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(tuple,list)):return[enc(y)for y in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
class Field:
    def __init__(self,old,bridge,edge):
        assert old['field']['base_square_root_basis']==RAD
        self.z=[B.Z]
        for _ in range(3):self.z.append((self.z[-1],self.z[-1]))
        self.rad=[parse(old['field']['discriminant'],0),parse(bridge['field']['new_positive_sqrt_radical'],1),parse(edge['field']['new_positive_sqrt_radical'],2)]
        self.Q=[[F(q)for q in row]for row in old['form']]
        self.rows=[[(j,q)for j,q in enumerate(row)if q]for row in self.Q]
        self.maxbits=0;self.cache={}
    def lift(self,x,k,target=3):
        for j in range(k,target):x=(x,self.z[j])
        return x
    def c(self,q,k=3):return self.lift((F(q),)+(F(0),)*15,0,k)
    def add(self,x,y,k=3):return B.add(x,y)if k==0 else tuple(self.add(a,b,k-1)for a,b in zip(x,y))
    def scale(self,x,q,k=3):return B.scale(x,q)if k==0 else tuple(self.scale(a,q,k-1)for a in x)
    def mul(self,x,y,k=3):
        if x==self.z[k]or y==self.z[k]:return self.z[k]
        if k==0:return B.mul(x,y)
        a,b=x;c,d=y;z=self.z[k-1]
        if b==z:return self.mul(a,c,k-1),self.mul(a,d,k-1)
        if d==z:return self.mul(a,c,k-1),self.mul(b,c,k-1)
        return(self.add(self.mul(a,c,k-1),self.mul(self.rad[k-1],self.mul(b,d,k-1),k-1),k-1),self.add(self.mul(a,d,k-1),self.mul(b,c,k-1),k-1))
    def total(self,xs,k=3):
        r=self.z[k]
        for x in xs:r=self.add(r,x,k)
        return r
    def vlift(self,v,k,target=3):return tuple(self.lift(x,k,target)for x in v)
    def lin(self,vs,cs,k=3):return tuple(self.total((self.mul(c,v[j],k)for c,v in zip(cs,vs)),k)for j in range(21))
    def dot(self,v,w,k=3):
        key=k,w
        if key not in self.cache:self.cache[key]=tuple(self.total((self.scale(w[j],q,k)for j,q in row),k)for row in self.rows)
        return self.total((self.mul(x,y,k)for x,y in zip(v,self.cache[key])),k)
    @lru_cache(None)
    def rb(self,bits):
        den=1<<bits;out=[]
        for d in RAD:
            a=isqrt(d*den*den);out.append((F(a,den),F(a+int(a*a!=d*den*den),den)))
        return out
    def bounds(self,x,k=3,bits=256):
        if k==0:
            lo=hi=F(0)
            for q,(l,h)in zip(x,self.rb(bits)):
                lo+=q*(l if q>=0 else h);hi+=q*(h if q>=0 else l)
            return lo,hi
        l,h=self.bounds(x[0],k-1,bits)
        if x[1]==self.z[k-1]:return l,h
        a,b=self.bounds(x[1],k-1,bits);c,d=self.sqrtbounds(self.rad[k-1],k-1,bits);vv=[a*c,a*d,b*c,b*d]
        return l+min(vv),h+max(vv)
    @lru_cache(None)
    def sqrtbounds(self,x,k,bits):
        l,h=self.bounds(x,k,bits);assert l>0
        den=1<<bits;a=isqrt(l.numerator*den*den//l.denominator);b=isqrt(h.numerator*den*den//h.denominator)+1
        assert F(a*a,den*den)<=l and F(b*b,den*den)>h
        return F(a,den),F(b,den)
    def sign(self,x,k=3):
        if x==self.z[k]:return 0
        for bits in[192,384,768,1536,3072,6144]:
            try:l,h=self.bounds(x,k,bits)
            except AssertionError:continue
            if l>0 or h<0:self.maxbits=max(self.maxbits,bits);return 1 if l>0 else -1
        raise AssertionError('unresolved exact scalar sign')
    def rational(self,x,k=1):
        while k:assert x[1]==self.z[k-1];x=x[0];k-=1
        assert not any(x[1:]);return x[0]
    def normalized_bounds(self,coeff,pairs,squares,k=3,bits=256):
        lo=hi=F(0)
        for c,g,q in zip(coeff,pairs,squares):
            if c==B.Z or g==self.z[k]:continue
            a,b=self.bounds(self.mul(self.lift(c,0,k),g,k),k,bits)
            l,h=self.sqrtbounds(q,1,bits);vv=[a/l,a/h,b/l,b/h]
            lo+=min(vv);hi+=max(vv)
        return lo,hi
    def normalized_sign(self,coeff,pairs,squares,k=3):
        if all(c==B.Z or g==self.z[k]for c,g in zip(coeff,pairs)):return 0,(F(0),F(0))
        for bits in[256,512,1024,2048,4096,8192]:
            try:l,h=self.normalized_bounds(coeff,pairs,squares,k,bits)
            except AssertionError:continue
            if l>0 or h<0:self.maxbits=max(self.maxbits,bits);return(1 if l>0 else -1),(l,h)
        raise AssertionError('unresolved normalized-basis sign')
    def normalized_dotted(self,coeff,pairs,squares,oldnorm):
        for bits in[256,512,1024,2048,4096,8192]:
            try:l,h=self.normalized_bounds(coeff,pairs,squares,3,bits)
            except AssertionError:continue
            qlo,qhi=self.bounds(oldnorm,3,bits)
            if h<0 and h*h>2*qhi:self.maxbits=max(self.maxbits,bits);return(l,h)
        raise AssertionError('normalized cap is not certified dotted')
def pd(M):
    M=[list(row)for row in M]
    for k in range(len(M)):
        p=M[k][k]
        if p<=0:return False
        for i in range(k+1,len(M)):
            for j in range(k+1,len(M)):M[i][j]-=M[i][k]*M[k][j]/p
    return True
