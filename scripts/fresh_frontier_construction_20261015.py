#!/usr/bin/env python3
"""Complete bounded Vinberg pilot on a new Q(sqrt2), discriminant-31 host.

Self-contained Fraction arithmetic; all search acceptance is exact.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import cmp_to_key
from hashlib import sha256
from math import isqrt, lcm, ceil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_frontier_construction_20261015'
HEIGHT = F(20)

def add(x,y): return x[0]+y[0], x[1]+y[1]
def sub(x,y): return x[0]-y[0], x[1]-y[1]
def mul(x,y): return x[0]*y[0]+2*x[1]*y[1], x[0]*y[1]+x[1]*y[0]
def scale(x,k): return x[0]*k,x[1]*k
def sign(x,conjugate=False):
    a,b=x
    if conjugate: b=-b
    if not b: return (a>0)-(a<0)
    if not a: return (b>0)-(b<0)
    if (a>0)==(b>0): return (a>0)-(a<0)
    z=a*a-2*b*b
    assert z
    return ((a>0)-(a<0))*((z>0)-(z<0))

class K:
    def __init__(self,a=0,b=0):
        if isinstance(a,K): self.a,self.b=a.a,a.b
        else: self.a,self.b=F(a),F(b)
    def __add__(self,y):
        y=K(y);return K(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,y):return self+-K(y)
    def __rsub__(self,y):return K(y)+-self
    def __mul__(self,y):
        y=K(y);return K(*mul((self.a,self.b),(y.a,y.b)))
    __rmul__=__mul__
    def conjugate(self):return K(self.a,-self.b)
    def __truediv__(self,y):
        y=K(y);d=y.a*y.a-2*y.b*y.b;assert d
        return self*K(y.a/d,-y.b/d)
    def __rtruediv__(self,y):return K(y)/self
    def __eq__(self,y):
        y=K(y);return self.a==y.a and self.b==y.b
    def __bool__(self):return bool(self.a or self.b)
    def __hash__(self):return hash((self.a,self.b))
    def __lt__(self,y):return (self-y).sign()<0
    def sign(self,conjugate=False):return sign((self.a,self.b),conjugate)
    def integral(self):return self.a.denominator==self.b.denominator==1
    def ri(self):assert self.integral();return int(self.a),int(self.b)
    def enc(self):return dict(a=str(self.a),b=str(self.b))

S=K(0,1);UNIT=1+S;DELTA=1+4*S;BETA=9+5*S
C=[[2,0,-1,0,0,0,0,0], [0,2,0,-1,0,0,0,0],
   [-1,0,2,-1,0,0,0,0], [0,-1,-1,2,-1,0,0,0],
   [0,0,0,-1,2,-1,0,0], [0,0,0,0,-1,2,-1,0],
   [0,0,0,0,0,-1,2,-1], [0,0,0,0,0,0,-1,2]]

def inverse(M):
    n=len(M);a=[[F(x)for x in r]+[F(i==j)for j in range(n)]for i,r in enumerate(M)]
    for i in range(n):
        j=next(j for j in range(i,n)if a[j][i]);a[i],a[j]=a[j],a[i]
        q=a[i][i];a[i]=[x/q for x in a[i]]
        for j in range(n):
            if j!=i:
                q=a[j][i];a[j]=[x-q*y for x,y in zip(a[j],a[i])]
    return [r[n:]for r in a]
def dot(x,y):return sum((a*b for a,b in zip(x,y)),K())
def mv(M,v):return [dot(row,v)for row in M]
def pair(x,y,G):return dot(x,mv(G,y))
def dec(x):return K(x['a'],x['b'])
def enc(x):
    if isinstance(x,K):return x.enc()
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def sqrt_ceil(x):
    x=F(x);n=isqrt(x.numerator//x.denominator)
    while n*n<x:n+=1
    return n
def rational_ldl(M):
    M=[[F(x)for x in r]for r in M];p=[]
    for k in range(len(M)):
        q=M[k][k];assert q>0;p.append(q)
        for i in range(k+1,len(M)):
            for j in range(k+1,len(M)):M[i][j]-=M[i][k]*M[k][j]/q
    return p

def zero_shell(J):
    forms=[inverse([r[:k]for r in J[:k]])for k in range(1,9)]
    bounds=[sqrt_ceil(2*J[i][i])for i in range(8)]
    roots=[];nodes=[0]*9
    def visit(v):
        k=len(v);nodes[k]+=1
        if k:
            q=sum(v[i]*forms[k-1][i][j]*v[j]for i in range(k)for j in range(k))
            if q>2:return
        if k==8:
            if q==2:roots.append(v)
            return
        for z in range(-bounds[k],bounds[k]+1):visit(v+[z])
    visit([])
    assert len(roots)==240
    assert all(sum(v[i]*C[i][j]*v[j]for i in range(8)for j in range(8))==2 for v in roots)
    return dict(integer_E8_roots=roots,nodes=nodes,total_roots_E8_A1=242)

def run(height_bound=HEIGHT):
    T=F(height_bound);assert T>0
    J=inverse(C);piv=rational_ldl(C)
    assert all(x.denominator==1 and x>=0 for r in J for x in r)
    determinant=F(1)
    for x in piv:determinant*=x
    assert determinant==1
    assert UNIT*UNIT.conjugate()==-1 and BETA==UNIT*DELTA
    assert DELTA*DELTA.conjugate()==-31 and BETA*BETA.conjugate()==31
    assert all(31%p for p in (2,3,5))
    assert BETA.sign()>0 and BETA.sign(True)>0
    assert (DELTA-6).sign()>0 and (-DELTA.conjugate()-4).sign()>0
    norms=[K(2),2*BETA]
    assert all((33-q).sign()>0 and (4-q.conjugate()).sign()>0 for q in norms)
    G=[[K()for _ in range(10)]for _ in range(10)]
    G[0][0]=-2*S;G[0][9]=G[9][0]=K(-1);G[9][9]=K(2)
    for i in range(8):
        for j in range(8):G[i+1][j+1]=K(C[i][j])
    H=[K(2)]+[K()]*8+[K(1)]
    assert pair(H,H,G)==-2*DELTA and mv(G,H)==[-DELTA]+[K()]*9
    simple=[[K(-int(i==j))for i in range(10)]for j in range(1,9)]
    simple.append([K(int(i==9))for i in range(10)])
    zero=zero_shell(J)
    # x²<11T; |x'|<sqrt2<2. For every positive coordinate z,
    # z²<66(1+T), |z'|<sqrt8<3. Since 2sqrt2>2,
    # |b|=|z-z'|/(2sqrt2)<(Z+3)/2, |a|<(Z+3)/2.
    X=sqrt_ceil(11*T);Z=sqrt_ceil(66*(1+T))
    xb=xa=ceil(F(X+2,2));zb=za=ceil(F(Z+3,2))
    controls=[]
    for q in norms:
        for a in range(-xa,xa+1):
            for b in range(-xb,xb+1):
                x=K(a,b);target=2*q+DELTA*x*x
                if x.sign()>0 and (T-DELTA*x*x/(2*q)).sign()>=0 and target.sign(True)>=0:
                    controls.append((DELTA*x*x/(2*q),q,x))
    controls.sort(key=lambda row:row[0])
    assert all(a[0]!=b[0]for a,b in zip(controls,controls[1:]))
    order=list(range(7,-1,-1))
    Finv=[[F(C[i][j],2)for j in order]+[F(0)]for i in order]+[[F(0)]*8+[F(1)]]
    forms=[]
    for n in range(1,10):
        M=inverse([r[:n]for r in Finv[:n]])
        den=lcm(*(z.denominator for r in M for z in r))
        forms.append((den,[[int(z*den)for z in r]for r in M]))
    previous=[];levels=[];total=Counter()
    for index,(height,q,x)in enumerate(controls):
        target=2*q+DELTA*x*x;ti=target.ri();xi=x.ri();div=q/2
        alphabet=[]
        for a in range(-za,za+1):
            for b in range(-zb,zb+1):
                z=K(a,b);rem=target-z*z
                if z.sign()>=0 and rem.sign()>=0 and rem.sign(True)>=0 and (z/div).integral():alphabet.append((a,b))
        alphabet.sort(key=cmp_to_key(lambda a,b:sign(sub(a,b))))
        co=[];bounds=[]
        for r in previous:
            c=[2*r[i+1]for i in order]+[r[0]-2*r[9]]
            assert all(z.sign()>=0 for z in c)
            co.append([z.ri()for z in c]);bounds.append((DELTA*x*r[0]).ri())
        visited=[0]*10;schur=[0]*10;walls=[0]*10;roots=[];parity=[0]
        def recurse(prefix,minimum,partial):
            depth=len(prefix);visited[depth]+=1
            if depth:
                num,den=minimum;rem=sub(scale(ti,den),num)
                if sign(rem)<0 or sign(rem,True)<0:schur[depth]+=1;return
                z=prefix[-1];partial=[add(p,mul(z,c[depth-1]))for p,c in zip(partial,co)]
                if any(sign(sub(bound,p))<0 for bound,p in zip(bounds,partial)):walls[depth]+=1;return
            if depth==9:
                num,den=minimum
                if num!=scale(ti,den):return
                l=K(*prefix[-1]);lam=[K(*z)for z in reversed(prefix[:8])]
                y=mv(J,lam);t=(x-l)/2;r=[x,*y,t]
                assert all(z.integral()for z in r)and pair(r,r,G)==q
                assert all((2*z/q).integral()for z in mv(G,r))
                assert all(pair(r,n,G).sign()<=0 for n in simple+previous)
                roots.append(r);return
            den,M=forms[depth];base=(0,0);linear=(0,0)
            for i in range(depth):
                linear=add(linear,scale(prefix[i],2*M[i][depth]))
                for j in range(depth):base=add(base,scale(mul(prefix[i],prefix[j]),M[i][j]))
            for z in alphabet:
                if depth==8 and any((z[i]-xi[i])%2 for i in range(2)):
                    parity[0]+=1;continue
                num=add(add(base,mul(linear,z)),scale(mul(z,z),M[depth][depth]))
                recurse(prefix+[z],(num,den),partial)
        recurse([],None,[(0,0)]*len(previous))
        assert all(pair(a,b,G).sign()<=0 for i,a in enumerate(roots)for b in roots[i+1:])
        previous+=roots;total['nodes']+=sum(visited);total['levels']+=1
        levels.append(dict(height=height,norm=q,control=x,target=target,alphabet=[K(*z)for z in alphabet],
            nodes=visited,schur_prunes=schur,wall_prunes=walls,parity_prunes=parity[0],roots=roots))
        print('level',index+1,'of',len(controls),'roots',len(roots),'walls',9+len(previous),'nodes',sum(visited),flush=True)
    N=simple+previous;gram=[[pair(r,s,G)for s in N]for r in N]
    angles={K():2,K(1):3,K(2):4,K(3):6,2+S:8,K(4):'parallel'}
    relations=[];counts=Counter()
    for i in range(len(N)):
        for j in range(i):
            assert gram[i][j].sign()<=0
            v=4*gram[i][j]*gram[i][j]/(gram[i][i]*gram[j][j])
            if (v-4).sign()>0:kind='dotted'
            else:
                assert v in angles,(i,j,v.enc())
                m=angles[v];kind='parallel'if m=='parallel'else'right'if m==2 else f'pi/{m}'
            counts[kind]+=1;relations.append(dict(pair=[j,i],cartan_cycle=v,kind=kind))
    direction=[K()]+mv(J,[K(1)]*8)+[K(-1)]
    epsilon=F(1)
    while True:
        center=[a+epsilon*b for a,b in zip(H,direction)]
        if pair(center,center,G).sign()<0 and all(pair(center,n,G).sign()<0 for n in N):break
        epsilon/=2
    facets=[]
    for i,n in enumerate(N):
        point=[a-pair(center,n,G)/gram[i][i]*b for a,b in zip(center,n)]
        assert pair(point,point,G).sign()<0 and pair(point,H,G).sign()<0
        assert all(pair(point,r,G).sign()==(0 if j==i else -1)for j,r in enumerate(N))
        facets.append(point)
    return enc(dict(status='PASS_COMPLETE_BOUNDED_DISCRIMINANT31_PREFIX',
        sources={f'scripts/{TAG}.py':sha(f'scripts/{TAG}.py')},height_bound=T,
        field='Q(sqrt2); every value encoded a+b sqrt2',form=G,basepoint=H,determinant=-DELTA,
        norm_classes=norms,totally_positive_prime=BETA,zero_simple_walls=simple,zero_shell=zero,
        inverse_E8=J,E8_ldl=piv,positive_coordinate_order=[f'lambda{i+1}'for i in order]+['l'],
        uniform_bounds=dict(control_physical=X,control_conjugate=2,control_coefficient_box=[xa,xb],
            coordinate_physical=Z,coordinate_conjugate=3,coordinate_coefficient_box=[za,zb]),
        levels=levels,roots=N,gram=gram,strict_center=center,center_epsilon=epsilon,
        actual_facet_points=facets,pair_relations=relations,pair_counts=dict(counts),
        counts=dict(total,positive_walls=len(previous),walls=len(N)),
        completeness=dict(norms='PID, evenness, primitive Bezout and adjugate imply q/2 divides determinant; norm31 prime and total positivity leave 1,beta modulo unit squares.',
            zero='x=0 norm2: u²+2v²=2 in even E8+A1 forces v=0; complete 242 ordinary roots. Long roots at x0 would be divisible by beta.',
            controls='x²<11T and |x_conjugate|<2; explicit coefficient boxes then exact both-place tests.',
            coordinates='2lambda.J.lambda+l²=2q+delta*x²; inverse positive form diagonal1 implies each coordinate bounded by target at both places.',
            congruences='All lambda and l divisible by q/2; l=x mod2O; explicit reflection integrality checked.',
            pruning='Rational prefix Schur minimum at both embeddings; all remaining prior-wall coefficients nonnegative at physical embedding.'),
        scope='Complete stated height prefix for one new lattice, not a reflectivity or finite-volume assertion.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved discriminant31 prefix replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['counts'])
