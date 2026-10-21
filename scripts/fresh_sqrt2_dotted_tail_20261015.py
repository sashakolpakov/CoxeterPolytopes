#!/usr/bin/env python3
"""Exact rank and integrality gates for two C5 opposite-pair attachment families."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt2_dotted_tail_20261015'

def sg(x):return (x>0)-(x<0)
class K:
    def __init__(self,a=0,b=0):
        if isinstance(a,K):self.a,self.b=a.a,a.b
        else:self.a,self.b=F(a),F(b)
    def __add__(self,o):o=K(o);return K(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,o):return self+-K(o)
    def __rsub__(self,o):return K(o)+-self
    def __mul__(self,o):o=K(o);return K(self.a*o.a+2*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=K(o);n=o.a*o.a-2*o.b*o.b;assert n
        return self*K(o.a/n,-o.b/n)
    def __rtruediv__(self,o):return K(o)/self
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,o):o=K(o);return self.a==o.a and self.b==o.b
    def sign(self):
        if not self.b:return sg(self.a)
        if not self.a or sg(self.a)==sg(self.b):return sg(self.b)
        return sg(self.a)*sg(self.a*self.a-2*self.b*self.b)
    def conj(self):return K(self.a,-self.b)
    def enc(self):return [str(self.a),str(self.b)]
    def integral(self):return self.a.denominator==self.b.denominator==1

def rational_sqrt(q):
    if q<0:return None
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    return F(a,b) if a*a==q.numerator and b*b==q.denominator else None

def field_sqrt(z):
    """Return its positive physical square root if in K, otherwise None."""
    assert z.sign()>0
    if not z.b:
        a=rational_sqrt(z.a)
        if a is not None:return K(a)
        b=rational_sqrt(z.a/2)
        return None if b is None else K(0,b)
    norm=rational_sqrt(z.a*z.a-2*z.b*z.b)
    if norm is None:return None
    for v in {norm,-norm}:
        a=rational_sqrt((z.a+v)/2)
        if a:
            root=K(a,z.b/(2*a));assert root*root==z
            return root if root.sign()>0 else -root
    return None

def sub(M,ix):return [[M[i][j] for j in ix] for i in ix]
def inertia(M):
    A=[[K(x) for x in row] for row in M];pos=neg=zero=0
    while A:
        n=len(A);pivot=next((i for i in range(n) if A[i][i]),None)
        if pivot is not None:
            ix=[pivot]+[i for i in range(n) if i!=pivot];A=sub(A,ix)
            d=A[0][0];pos+=d.sign()>0;neg+=d.sign()<0
            A=[[A[i][j]-A[i][0]*A[0][j]/d for j in range(1,n)] for i in range(1,n)]
        else:
            pair=next(((i,j) for i in range(n) for j in range(i+1,n) if A[i][j]),None)
            if pair is None:zero+=n;break
            ix=list(pair)+[i for i in range(n) if i not in pair];A=sub(A,ix);b=A[0][1]
            pos+=1;neg+=1
            A=[[A[i][j]-(A[i][0]*A[1][j]+A[i][1]*A[0][j])/b for j in range(2,n)] for i in range(2,n)]
    return [pos,neg,zero]

def literal(ell=0,double=False,contact=None):
    """Norm-two matrix, with the single-family y rescaled by sqrt(D)."""
    n=7+ell;M=[[K(2*int(i==j)) for j in range(n)] for i in range(n)]
    edges=[(2,3,K(0,1)),(2,4,K(1)),(4,5,K(1)),(5,6,K(1)),(6,3,K(1)),(0,5,K(1))]
    if double:edges.extend([(1,5,K(1)),(0,1,contact)])
    else:
        D=contact;M[1][1]=2*D;edges.append((0,1,D))
    for k in range(ell):edges.append((0 if k==0 else 6+k,7+k,K(1)))
    for i,j,c in edges:M[i][j]=M[j][i]=-c
    return M

def verify_local(ell,double,contact):
    M=literal(ell,double,contact);n=len(M);W=[4,5,6]+list(range(7,n))
    corners=[]
    for pair in [(0,2),(2,1),(1,3),(3,0)]:
        ix=W+list(pair);assert inertia(sub(M,ix))==[len(ix),0,0]
        corners.append(dict(marked=list(pair),support=ix))
    assert inertia(M)==[n-2,1,1]
    assert inertia([[x.conj() for x in row] for row in M])==[n-1,0,1]
    for pair in [(0,1),(2,3)]:assert inertia(sub(M,W+list(pair)))[1]==1
    tested=0
    for size in range(1,n+1):
        for ix in combinations(range(n),size):
            iv=inertia(sub(M,ix));assert not(iv[1]==0 and iv[2]>0),(ell,double,ix)
            tested+=1
    return dict(ell=ell,double=double,contact=contact.enc(),gram=[[x.enc() for x in row] for row in M],
                corners=corners,physical_inertia=[n-2,1,1],conjugate_inertia=[n-1,0,1],no_affine_subset_tests=tested)

def run():
    p,q,x=s.symbols('p q x');r2=s.sqrt(2);M=2*s.eye(7);M[0,0]=2*p;M[1,1]=2*q
    for i,j,c in [(2,3,r2),(2,4,1),(4,5,1),(5,6,1),(6,3,1),(0,5,1),(0,1,x)]:M[i,j]=M[j,i]=-c
    d1=s.expand(M.det());assert s.simplify(d1-2*(1+r2)*(x*x-q*(4*p-1+r2)))==0
    M[1,5]=M[5,1]=-1;t=(1-r2)/2
    d2=s.expand(M.det());f=x*x+2*t*x-4*p*q+2*t*(p+q)
    assert s.simplify(d2-2*(1+r2)*f)==0
    tail=[]
    for ell in range(5):
        pg=F(ell+2,2*(ell+1))
        if ell:
            A=2*s.eye(ell)
            for i in range(ell-1):A[i,i+1]=A[i+1,i]=-1
            assert s.Rational(pg)==1-A.inv()[0,0]/2
        tail.append(dict(length=ell,reduced_unit_diagonal=str(pg)))
    # Affine B3 fork forced by any pendant contact at b or d.
    fork=[[K(1),K(0,F(-1,2)),K(F(-1,2)),K(F(-1,2))],
          [K(0,F(-1,2)),K(1),K(0),K(0)],
          [K(F(-1,2)),K(0),K(1),K(0)],
          [K(F(-1,2)),K(0),K(0),K(1)]]
    kernel=[K(1),K(0,F(1,2)),K(F(1,2)),K(F(1,2))]
    assert all(sum((fork[i][j]*kernel[j] for j in range(4)),K())==0 for i in range(4))
    assert inertia(fork)==[3,0,1]
    # The a-tail corners are D5,E6,E7,E8, followed by an affine E8.
    corner_tails=[]
    for ell in range(5):
        M0=literal(ell,False,K(3,1));ix=[4,5,6]+list(range(7,7+ell))+[0,2]
        iv=inertia(sub(M0,ix));assert iv==([5+ell,0,0] if ell<4 else [8,0,1])
        corner_tails.append(dict(length=ell,inertia=iv))
    singles=[]
    for ell in range(4):
        R=F(ell+3,ell+1);D=K(R,1);special=K(F(1,2),F(-1,4))*D
        assert not special.integral()
        singles.append(dict(length=ell,contact_square=D.enc(),integral=D.integral(),
                            I2_8_y_attachment_contact_square=special.enc()))
    assert [r['length'] for r in singles if r['integral']]==[0,1]
    doubles=[];tt=K(F(1,2),F(-1,2))
    for ell in range(4):
        for m in range(4):
            pp=F(ell+2,2*(ell+1));qq=F(m+2,2*(m+1));delta=(2*pp-tt)*(2*qq-tt)
            root=field_sqrt(delta);contact=None if root is None else -tt+root
            if contact is not None:assert contact*contact+2*tt*contact-4*pp*qq+2*tt*(pp+qq)==0
            doubles.append(dict(tails=[ell,m],discriminant=delta.enc(),field_sqrt=None if root is None else root.enc(),
                                contact=None if contact is None else contact.enc(),integral=contact is not None and contact.integral()))
    assert [r['tails'] for r in doubles if r['integral']]==[[0,0]]
    allowed=[verify_local(0,False,K(3,1)),verify_local(1,False,K(2,1)),verify_local(0,True,K(1,1))]
    deps=['fresh_sqrt2_lanner_catalogue_20261015','fresh_arithmetic_sqrt2_quad_rank_20261015']
    files=[f'scripts/{TAG}.py']+[p for tag in deps for p in (f'scripts/{tag}.py',f'results/{tag}.json',f'results/{tag}.md')]
    digest=lambda p:sha256((ROOT/p).read_bytes()).hexdigest()
    for p in files:
        if p.endswith('.json'):
            for f,h in json.loads((ROOT/p).read_text()).get('sources',{}).items():assert digest(f)==h
    return dict(status='PASS_C5_OPPOSITE_PAIR_PENDANT_TREE_RANK_AND_INTEGRALITY_CLASSIFICATION',
                sources={p:digest(p) for p in files},single_determinant=str(s.factor(d1)),double_determinant=str(s.factor(d2)),
                tail_diagonals=tail,corner_tails=corner_tails,affine_fork_inertia=[3,0,1],
                single_family_cases=singles,double_family_cases=doubles,allowed_local_stars=allowed,
                scope='Complete under the fixed C5 core and single/double middle attachment, with arbitrary pendant white trees attached once at any core node and orthogonal white components. Shared/external components with multiple core attachments, extra edges and different core labels are outside scope. No global dimension bound or finite-volume example.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    data=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'])
