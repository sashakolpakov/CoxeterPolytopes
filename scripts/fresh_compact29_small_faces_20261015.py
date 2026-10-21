#!/usr/bin/env python3
"""Exact local small-face completions, not compact H29 constructions.

All arithmetic is Fraction in K=Q(sqrt2,sqrt5), with one explicit positive
quadratic extension per quadrilateral. No numeric eigenvalue is used.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact29_small_faces_20261015'


def sg(x): return (x>0)-(x<0)


def q2sign(a,b):
    if not b: return sg(a)
    if not a or sg(a)==sg(b): return sg(b)
    return sg(a)*sg(a*a-2*b*b)


class K:
    """Basis 1,sqrt2,sqrt5,sqrt10."""
    def __init__(self,x=0):
        self.v=x.v if isinstance(x,K) else tuple(map(F,x)) if isinstance(x,(tuple,list)) else (F(x),F(0),F(0),F(0))
    def __add__(self,y):return K(tuple(a+b for a,b in zip(self.v,K(y).v)))
    __radd__=__add__
    def __neg__(self):return K(tuple(-a for a in self.v))
    def __sub__(self,y):return self+-K(y)
    def __rsub__(self,y):return K(y)+-self
    def __mul__(self,y):
        y=K(y);out=[F(0)]*4
        for i,a in enumerate(self.v):
            for j,b in enumerate(y.v):
                mask=i&j;out[i^j]+=a*b*(2 if mask&1 else 1)*(5 if mask&2 else 1)
        return K(out)
    __rmul__=__mul__
    def inv(self):
        a,b,c,d=self.v; conj=K((a,b,-c,-d)); u=self*conj
        assert u.v[2:]==(0,0)
        h=u.v[0]*u.v[0]-2*u.v[1]*u.v[1];assert h
        return conj*K((u.v[0]/h,-u.v[1]/h,0,0))
    def __truediv__(self,y):return self*K(y).inv()
    def __rtruediv__(self,y):return K(y)*self.inv()
    def __bool__(self):return any(self.v)
    def __eq__(self,y):return self.v==K(y).v
    def sign(self):
        a,b,c,d=self.v; su=q2sign(a,b);sv=q2sign(c,d)
        if not sv:return su
        if not su or su==sv:return sv
        diff=K((a,b,0,0))*K((a,b,0,0))-5*K((c,d,0,0))*K((c,d,0,0))
        return su*q2sign(diff.v[0],diff.v[1])
    def enc(self):return [str(x) for x in self.v]


class E:
    D=K(1)
    def __init__(self,a=0,b=0):
        if isinstance(a,E):self.a,self.b=a.a,a.b
        else:self.a,self.b=K(a),K(b)
    def __add__(self,y):y=E(y);return E(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return E(-self.a,-self.b)
    def __sub__(self,y):return self+-E(y)
    def __rsub__(self,y):return E(y)+-self
    def __mul__(self,y):y=E(y);return E(self.a*y.a+self.b*y.b*self.D,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def inv(self):
        den=self.a*self.a-self.D*self.b*self.b;assert den
        return E(self.a/den,-self.b/den)
    def __truediv__(self,y):return self*E(y).inv()
    def __rtruediv__(self,y):return E(y)*self.inv()
    def __bool__(self):return bool(self.a) or bool(self.b)
    def __eq__(self,y):y=E(y);return self.a==y.a and self.b==y.b
    def sign(self):
        sa,sb=self.a.sign(),self.b.sign()
        if not sb:return sa
        if not sa or sa==sb:return sb
        return sa*(self.a*self.a-self.D*self.b*self.b).sign()
    def enc(self):return [self.a.enc(),self.b.enc()]


def inverse(M):
    n=len(M);cls=type(M[0][0]); A=[list(row)+[cls(i==j) for j in range(n)] for i,row in enumerate(M)]
    for j in range(n):
        p=next(i for i in range(j,n) if A[i][j]);A[j],A[p]=A[p],A[j]
        z=A[j][j];A[j]=[x/z for x in A[j]]
        for i in range(n):
            if i!=j:
                z=A[i][j]
                if z:A[i]=[a-z*b for a,b in zip(A[i],A[j])]
    return [row[n:] for row in A]


def det(M):
    A=[list(r) for r in M];cls=type(M[0][0]);d=cls(1);n=len(A)
    for j in range(n):
        p=next((i for i in range(j,n) if A[i][j]),None)
        if p is None:return cls(0)
        if p!=j:A[p],A[j]=A[j],A[p];d=-d
        z=A[j][j];d=d*z
        for i in range(j+1,n):
            t=A[i][j]/z
            if t:
                for k in range(j+1,n):A[i][k]=A[i][k]-t*A[j][k]
    return d


def dot(x,y,M):return sum((x[i]*M[i][j]*y[j] for i in range(len(x)) for j in range(len(y))),type(M[0][0])(0))
def mv(M,x):return [sum((a*b for a,b in zip(row,x)),type(row[0])(0)) for row in M]
def sub(M,ix):return [[M[i][j] for j in ix] for i in ix]


def psd_status(M):
    """1=positive definite,0=PSD singular,-1=not PSD, by exact Schur pivots."""
    A=[list(r) for r in M];singular=False
    for i in range(len(A)):
        s=A[i][i].sign()
        if s<0:return -1
        if s==0:
            if any(A[i][j] for j in range(i+1,len(A))):return -1
            singular=True;continue
        for j in range(i+1,len(A)):
            for k in range(j,len(A)):
                A[k][j]=A[j][k]=A[j][k]-A[j][i]*A[i][k]/A[i][i]
    return 0 if singular else 1


@lru_cache(None)
def cached_psd(entries,n):return psd_status([[K(entries[i*n+j]) for j in range(n)] for i in range(n)])


def no_parabolics(G,dotted=None):
    checked=skipped=0
    for size in range(1,len(G)+1):
        for S in combinations(range(len(G)),size):
            if dotted is not None and set(dotted)<=set(S):
                skipped+=1;continue # Its 2x2 principal minor is negative.
            M=sub(G,S)
            if isinstance(M[0][0],E):
                assert all(not x.b for row in M for x in row)
                M=[[x.a for x in row] for row in M]
            key=tuple(x.v for row in M for x in row)
            status=cached_psd(key,size)
            assert status!=0,(size,S)
            checked+=1
    return dict(exact_principal_tests=checked,subsets_excluded_by_dotted_pair=skipped)


COS={2:K(),3:K(1),4:K((0,1,0,0)),5:K((F(1,2),0,F(1,2),0))}


def diagram(name):
    kind=name[0];n=int(name[1:]);e={(i,i+1):3 for i in range(n-1)}
    if kind=='B':e[0,1]=4
    if kind=='F':e[1,2]=4
    if kind=='H':e[0,1]=5
    if kind=='D':e={(0,1):3,(1,2):3,(1,3):3}
    if kind=='E':e={(i,i+1):3 for i in range(n-2)};e[2,n-1]=3
    return n,e


def base_gram(n,e,size):
    M=[[K(2*(i==j)) for j in range(size)] for i in range(size)]
    for (i,j),m in e.items():M[i][j]=M[j][i]=-COS[m]
    return M


def schur(M,W,B):
    J=inverse(sub(M,W))
    return [[M[i][j]-sum((M[i][W[k]]*J[k][l]*M[W[l]][j] for k in range(len(W)) for l in range(len(W))),type(M[0][0])(0))
             for j in B] for i in B]


def triangle(recipe):
    n,e=diagram(recipe['type']);a,b=recipe['marked'];c=n;W=[i for i in range(n) if i not in [a,b]]
    for i,m in recipe['contacts']:e[i,c]=m
    G=base_gram(n,e,n+1);B=[a,b,c]
    assert psd_status(sub(G,list(range(n))))==1
    assert all(psd_status(sub(G,W+[i,j]))==1 for i,j in combinations(B,2))
    R=schur(G,W,B);dr=det(R);assert dr.sign()<0
    I=inverse(R);rays=[[-I[i][j] for i in range(3)] for j in range(3)]
    dots=[[dot(x,y,R) for y in rays] for x in rays]
    assert all(x.sign()<0 for row in dots for x in row)
    for j,x in enumerate(rays):assert mv(R,x)==[K(-int(i==j)) for i in range(3)]
    return dict(residual_Gram=R,residual_determinant=dr,face_vertices=rays,vertex_pairings=dots,
                parabolic_check=no_parabolics(G),ambient_signature=[n,1,0],corner_count=3,
                center='Sum of the three displayed future vertex rays, then a sufficiently small negative white-dual perturbation.')


def quad(recipe):
    n,e=diagram(recipe['type']);a,b=recipe['marked'];c,d=n,n+1;W=[i for i in range(n) if i not in [a,b]]
    for i,m in recipe['c_contacts']:e[i,c]=m
    for i,m in recipe['d_contacts']:e[i,d]=m
    e[a,c]=recipe['ac_label'];e[c,d]=recipe['cd_label']
    G0=base_gram(n,e,n+2);B=[a,b,c,d]
    def residual(y):
        M=[row[:] for row in G0];M[b][d]=M[d][b]=-K(y)
        return schur(M,W,B)
    p0,p1,pm=[det(residual(y)) for y in [0,1,-1]]
    A=(p1+pm)/2-p0;BB=(p1-pm)/2;C=p0;D=BB*BB-4*A*C
    assert A.sign()>0 and C.sign()<0 and D.sign()>0
    assert (4*A+2*BB+C).sign()<0 # positive root strictly greater than2
    E.D=D;y=E(-BB/(2*A),1/(2*A))
    assert (y-2).sign()>0 and E(A)*y*y+E(BB)*y+E(C)==0
    G=[[E(x) for x in row] for row in G0];G[b][d]=G[d][b]=-y
    R=schur(G,W,B);assert not det(R)
    A3=sub(R,[0,1,2]);assert det(A3).sign()<0
    assert all(psd_status(sub(G,W+[B[i],B[j]]))==1 for i,j in [(0,1),(1,2),(2,3),(3,0)])
    assert det(sub(R,[0,2])).sign()<0 and det(sub(R,[1,3])).sign()<0
    rays=[];pairings=[]
    for i,j in [(0,1),(1,2),(2,3),(3,0)]:
        rows=[A3[k] if k<3 else [R[t][3] for t in range(3)] for k in [i,j]]
        k=next(k for k in range(3) if k not in [i,j]);rows.append(A3[k])
        x=mv(inverse(rows),[E(0),E(0),E(-1)])
        ps=mv(A3,x)+[sum((R[t][3]*x[t] for t in range(3)),E())]
        assert all(not ps[k] if k in [i,j] else ps[k].sign()<0 for k in range(4))
        rays.append(x);pairings.append(ps)
    dots=[[dot(x,z,A3) for z in rays] for x in rays]
    assert all(x.sign()<0 for row in dots for x in row)
    return dict(dotted_pair=[b,d],dotted_magnitude=y,positive_radicand=D,quadratic=[A,BB,C],
                residual_Gram=R,residual_basis=[0,1,2],residual_basis_determinant=det(A3),
                face_vertices=rays,face_vertex_wall_pairings=pairings,vertex_pairings=dots,
                parabolic_check=no_parabolics(G,(b,d)),ambient_signature=[n,1,1],corner_count=4,
                center='Sum of the four displayed future vertex rays, then a sufficiently small negative white-dual perturbation.')


def enc(x):
    if isinstance(x,(K,E)):return x.enc()
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def run():
    path=ROOT/f'results/{TAG}_recipes.json';recipes=json.loads(path.read_text())['recipes'];out=[]
    expected={(name,p) for name in ['A4','B4','D4','F4','H4','E7','E8'] for p in combinations(range(int(name[1:])),2)}
    assert {(r['type'],tuple(r['marked'])) for r in recipes}==expected and len(recipes)==79
    for idx,r in enumerate(recipes):
        proof=triangle(r) if r['kind']=='triangle' else quad(r)
        out.append(enc(dict(recipe=r,proof=proof)))
        print(idx+1,r['type'],r['marked'],r['kind'],flush=True)
    counts={name:dict(Counter(r['kind'] for r in recipes if r['type']==name)) for name in sorted({r['type'] for r in recipes})}
    assert sum(r['kind']=='triangle' for r in recipes)==56
    # The disjoint orthogonal product with spectator roots preserves every
    # finite corner, the unique negative direction, and absence of parabolics.
    # E7+E7+E7+E8 is an actual typed forest in the frozen extremal list.
    old=ROOT/'results/fresh_compact_dimension_bound_20261015.json';data=json.loads(old.read_text())
    assert any(['E7','E7','E7','E8'] in row['typed_forests'] for row in data['extremal_rank_partitions'])
    sources=[f'scripts/{TAG}.py',f'results/{TAG}_recipes.json','results/fresh_compact_dimension_bound_20261015.json',
             'scripts/fresh_compact_dimension_bound_20261015.py']
    return dict(status='PASS_EXACT_LOCAL_SMALL_FACE_BARRIER_79_MARKED_PAIRS',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                base_field='K=Q(sqrt2,sqrt5), basis1,sqrt2,sqrt5,sqrt10; E entries encode a+b*sqrt(D), positive square root, D saved per quadrilateral.',
                cases=out,type_counts=counts,marked_pair_count=79,triangle_count=56,quadrilateral_count=23,
                exact_principal_tests=sum(r['proof']['parabolic_check']['exact_principal_tests'] for r in out),
                subsets_excluded_by_dotted_pair=sum(r['proof']['parabolic_check']['subsets_excluded_by_dotted_pair'] for r in out),
                extremal_witness_forest=['E7','E7','E7','E8'],rank=29,original_Vinberg_weight=87,
                conclusion='Every positive-weight marked pair of E7+E7+E7+E8 individually admits a realizable local compact triangle/quadrilateral face star with finite adjacent corners and no parabolic subdiagram. No individual marked angle of this forest can be excluded by these local conditions alone.',
                scope='Separate realizable local stars, not simultaneous completion of the forest, not compact ambient chambers, not a dimension bound or new finite-volume example. No negative classification is inferred from the numerical scouts.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(path.read_text())==result
    else:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'],result['exact_principal_tests'],flush=True)
