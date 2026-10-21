#!/usr/bin/env python3
"""Independent SymPy-field and complete residual-cone replay of local stars.

No producer arithmetic imports. No-affine checks use recursive principal
minors, rather than the producer's Schur-pivot implementation.
"""
import argparse,json,math
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact29_small_faces_20261015_review'
PRIMARY='results/fresh_compact29_small_faces_20261015.json'
RECIPES='results/fresh_compact29_small_faces_20261015_recipes.json'
K=s.QQ.algebraic_field(s.sqrt(2),s.sqrt(5));Z,O=K.zero,K.one
BASE=[s.S.One,s.sqrt(2),s.sqrt(5),s.sqrt(10)]

class I:
    def __init__(self,l=0,h=None):self.l=F(l);self.h=F(l if h is None else h)
    def __add__(self,o):
        o=o if isinstance(o,I)else I(o);return I(self.l+o.l,self.h+o.h)
    __radd__=__add__
    def __neg__(self):return I(-self.h,-self.l)
    def __sub__(self,o):return self+-o if isinstance(o,I)else self+I(-o)
    def __mul__(self,o):
        o=o if isinstance(o,I)else I(o);v=[a*b for a in(self.l,self.h)for b in(o.l,o.h)]
        return I(min(v),max(v))
    __rmul__=__mul__
    def sqrt(self):
        assert self.l>=0;scale=1<<224
        lo=math.isqrt(self.l.numerator*scale*scale//self.l.denominator)
        hi=math.isqrt(self.h.numerator*scale*scale//self.h.denominator)
        return I(F(lo,scale),F(hi+1,scale))
    def pow(self,n):
        assert n>=0;r=I(1)
        for _ in range(n):r=r*self
        return r
def expr_interval(e):
    if e.is_Rational:return I(F(int(e.p),int(e.q)))
    if e.is_Add:return sum((expr_interval(x)for x in e.args),I())
    if e.is_Mul:
        v=I(1)
        for x in e.args:v=v*expr_interval(x)
        return v
    if e.is_Pow:
        a,b=e.args
        if b==s.Rational(1,2):return expr_interval(a).sqrt()
        if b.is_Integer and b>=0:return expr_interval(a).pow(int(b))
    raise ValueError(e)

class Field:
    def __init__(self,T):
        self.T=T;self.z=T.zero;self.o=T.one;self.theta=expr_interval(T.ext.as_expr())
        b2,b5=T.from_sympy(s.sqrt(2)),T.from_sympy(s.sqrt(5))
        self.basis=[self.o,b2,b5,b2*b5]
    @lru_cache(None)
    def sign(self,x):
        if not x:return 0
        a=I()
        for c in x.to_list():a=a*self.theta+I(F(str(c)))
        assert a.l>0 or a.h<0,('unresolved224bit sign',x)
        return 1 if a.l>0 else-1
    @lru_cache(None)
    def decode_tuple(self,x):return sum((self.T.convert(F(a))*b for a,b in zip(x,self.basis)),self.z)
    def decode_base(self,x):return self.decode_tuple(tuple(x))
    def dm(self,M):return DomainMatrix(M,(len(M),len(M[0])),self.T)
    def det(self,M):return self.dm(M).det()
    def inv(self,M):return self.dm(M).inv().to_list()
    def dot(self,v,w):return sum((a*b for a,b in zip(v,w)),self.z)
    def mv(self,M,v):return[self.dot(row,v)for row in M]
    def pair(self,v,w,G):return self.dot(v,self.mv(G,w))
    def pd(self,M):
        # Sylvester's leading-principal-minor criterion, reversed ordering.
        A=[list(reversed(row))for row in reversed(M)]
        return all(self.sign(self.det([r[:j]for r in A[:j]]))>0 for j in range(1,len(A)+1))
    def schur(self,G,W,B):
        J=self.inv(sub(G,W))
        return [[G[i][j]-self.dot([G[i][k]for k in W],self.mv(J,[G[k][j]for k in W]))for j in B]for i in B]

F0=Field(K)
def sub(M,S):return[[M[i][j]for j in S]for i in S]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def expression(x):return K.to_sympy(x)
def enc(x):return list(map(str,x.to_list()))

def base_matrix(recipe):
    name=recipe['type'];n=int(name[1:]);kind=name[0]
    edges={(i,i+1):3 for i in range(n-1)}
    if kind in('B','H'):edges[0,1]=4 if kind=='B' else 5
    if kind=='F':edges[1,2]=4
    if kind=='D':edges={(1,j):3 for j in(0,2,3)}
    if kind=='E':edges={(i,i+1):3 for i in range(n-2)};edges[2,n-1]=3
    a,b=recipe['marked'];c=n;quad=recipe['kind']=='quadrilateral';size=n+1+quad
    if not quad:
        for i,m in recipe['contacts']:edges[i,c]=m
    else:
        d=n+1
        for i,m in recipe['c_contacts']:edges[i,c]=m
        for i,m in recipe['d_contacts']:edges[i,d]=m
        edges[a,c]=recipe['ac_label'];edges[c,d]=recipe['cd_label']
    cval={2:Z,3:O,4:K.from_sympy(s.sqrt(2)),5:K.from_sympy((1+s.sqrt(5))/2)}
    G=[[K.convert(2*(i==j))for j in range(size)]for i in range(size)]
    for(i,j),m in edges.items():G[i][j]=G[j][i]=-cval[m]
    return n,G

@lru_cache(None)
def determinant(entries,n):return F0.det([[entries[i*n+j]for j in range(n)]for i in range(n)])
def no_affine(G,dotted):
    status={0:1};covered=skipped=evaluated=0;minimal_bad=0
    n=len(G);dm=0 if dotted is None else sum(1<<i for i in dotted)
    for size in range(1,n+1):
        for ss in combinations(range(n),size):
            mask=sum(1<<i for i in ss)
            if dm and mask&dm==dm:status[mask]=-1;skipped+=1;continue
            covered+=1
            if any(status[mask^(1<<i)]<0 for i in ss):status[mask]=-1;continue
            ix=tuple(reversed(ss));M=sub(G,ix);v=determinant(tuple(x for row in M for x in row),size)
            sg=F0.sign(v);evaluated+=1
            # Every proper principal submatrix is PSD. A zero determinant
            # here is therefore exactly a forbidden PSD singular submatrix.
            assert sg!=0,('affine principal submatrix',ss)
            status[mask]=sg;minimal_bad+=sg<0
    return dict(covered=covered,dotted_shortcuts=skipped,determinants_evaluated=evaluated,minimal_indefinite=minimal_bad)

def cross(a,b):return[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def canonical(v):
    t=next(x for x in v if x);return tuple(x/t for x in v)

def case_check(recipe,saved):
    n,G0=base_matrix(recipe);a,b=recipe['marked'];W=[i for i in range(n)if i not in(a,b)]
    assert F0.pd(sub(G0,list(range(n))))
    quad=recipe['kind']=='quadrilateral';B=[a,b,n]+([n+1]if quad else[])
    R0=F0.schur(G0,W,B)
    if quad:
        # Direct determinant cofactors, independently of interpolation at0,+/-1.
        A=-F0.det(sub(R0,[0,2]))
        minor=[[R0[i][j]for j in range(4)if j!=3]for i in range(4)if i!=1]
        BB=-2*F0.det(minor);C=F0.det(R0);disc=BB*BB-4*A*C
        assert [A,BB,C]==list(map(F0.decode_base,saved['quadratic']))
        assert disc==F0.decode_base(saved['positive_radicand'])and F0.sign(disc)>0
        assert F0.sign(A)>0 and F0.sign(C)<0 and F0.sign(4*A+2*BB+C)<0
        de=expression(disc);T=s.QQ.algebraic_field(s.sqrt(2),s.sqrt(5),s.sqrt(de));F1=Field(T)
        @lru_cache(None)
        def cv(x):
            assert K.ext.as_expr()==s.sqrt(2)+s.sqrt(5)
            out=T.zero;theta=F1.basis[1]+F1.basis[2]
            for c in x.to_list():out=out*theta+T.convert(F(str(c)))
            return out
        root=T.from_sympy(s.sqrt(de));y=(-cv(BB)+root)/(2*cv(A))
        assert F1.sign(y-2)>0
        decode=lambda x:F1.decode_base(x[0])+F1.decode_base(x[1])*root
        assert y==decode(saved['dotted_magnitude'])
        G=[[cv(x)for x in row]for row in G0];G[b][n+1]=G[n+1][b]=-y
        R=F1.schur(G,W,B);assert F1.dm(R).rank()==3
        dotted=(b,n+1)
    else:
        F1=F0;G=G0;R=R0;decode=F0.decode_base;dotted=None
        assert F1.det(R)==decode(saved['residual_determinant'])and F1.sign(F1.det(R))<0
    assert R==[[decode(x)for x in row]for row in saved['residual_Gram']]
    J=sub(R,[0,1,2]);assert F1.pd(sub(J,[0,1]))and F1.sign(F1.det(J))<0
    if quad:assert F1.det(J)==decode(saved['residual_basis_determinant'])
    cyc=[(i,(i+1)%len(B))for i in range(len(B))]
    for i,j in cyc:assert F1.pd(sub(G,W+[B[i],B[j]]))
    forms=J+([[R[i][3]for i in range(3)]]if quad else[])
    rays=[];contacts=[];rejected=[]
    for i,j in combinations(range(len(B)),2):
        x=cross(forms[i],forms[j]);assert any(x)
        found=[]
        for v in[x,[-t for t in x]]:
            vals=F1.mv(forms,v)
            if all(F1.sign(t)<=0 for t in vals):found.append((v,vals))
        if found:
            assert len(found)==1;v,vals=found[0]
            assert [k for k,t in enumerate(vals)if not t]==[i,j]
            assert F1.sign(F1.pair(v,v,J))<0
            rays.append(v);contacts.append([i,j])
        else:rejected.append([i,j])
    assert len(rays)==len(B)and {tuple(x)for x in contacts}=={tuple(sorted(x))for x in cyc}
    literal=[[decode(x)for x in v]for v in saved['face_vertices']]
    assert {canonical(v)for v in rays}=={canonical(v)for v in literal}
    dots=[[F1.pair(v,w,J)for w in literal]for v in literal]
    assert dots==[[decode(x)for x in row]for row in saved['vertex_pairings']]
    assert all(F1.sign(t)<0 for row in dots for t in row)
    if quad:assert [F1.mv(forms,v)for v in literal]==[[decode(x)for x in row]for row in saved['face_vertex_wall_pairings']]
    H=[sum((v[i]for v in literal),F1.z)for i in range(3)]
    assert F1.sign(F1.pair(H,H,J))<0 and all(F1.sign(t)<0 for t in F1.mv(forms,H))
    na=no_affine(G0,dotted)
    assert na['covered']==saved['parabolic_check']['exact_principal_tests']
    assert na['dotted_shortcuts']==saved['parabolic_check']['subsets_excluded_by_dotted_pair']
    assert saved['ambient_signature']==[n,1,int(quad)]
    return dict(type=recipe['type'],marked=recipe['marked'],kind=recipe['kind'],ambient_signature=[n,1,int(quad)],
        exact_field_degree=F1.T.ext.minpoly.degree(),pair_intersections=len(list(combinations(B,2))),
        feasible_timelike_rays=len(rays),adjacent_pairs=contacts,rejected_opposite_pairs=rejected,
        no_affine=na)

def run():
    p=read(PRIMARY);rs=read(RECIPES)['recipes']
    for f,h in p['sources'].items():assert sha(f)==h,f
    target={(k,tuple(v))for k in ['A4','B4','D4','F4','H4','E7','E8']for v in combinations(range(int(k[1:])),2)}
    assert len(rs)==79 and {(r['type'],tuple(r['marked']))for r in rs}==target
    out=[]
    for i,(r,row)in enumerate(zip(rs,p['cases'])):
        assert r==row['recipe'];out.append(case_check(r,row['proof']))
        print(i+1,r['type'],r['marked'],'PASS',flush=True)
    assert sum(x['kind']=='triangle'for x in out)==56
    assert sum(x['no_affine']['covered']for x in out)==25457
    assert sum(x['no_affine']['dotted_shortcuts']for x in out)==4864
    targetdata=read('results/fresh_compact_dimension_bound_20261015.json')
    assert any(['E7','E7','E7','E8']in r['typed_forests']for r in targetdata['extremal_rank_partitions'])
    assert p['original_Vinberg_weight']==3*21+24==87
    paths=[PRIMARY,RECIPES,'scripts/fresh_compact29_small_faces_20261015.py',
        'results/fresh_compact_dimension_bound_20261015.json',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT79_LOCAL_SMALL_FACE_STARS',sources={f:sha(f)for f in paths},
        arithmetic='Independent SymPy number fields,224-bit rational outward signs, determinant cofactors, full pair-intersection residual cone and recursive principal-minor PSD tests; no producer imports.',
        cases=out,marked_pairs=79,E7_E8_pairs=49,triangles=56,quadrilaterals=23,
        principal_subsets_covered=25457,dotted_pair_shortcuts=4864,
        actual_determinant_evaluations=sum(x['no_affine']['determinants_evaluated']for x in out),
        conclusion='Every individual marked pair in the stated finite component types has an actual local compact triangle/quadrilateral face star. Orthogonal spectators give each individual E7^3+E8 marked angle a rank29 local realization.',
        scope='Separate local stars only. No simultaneous compact chamber, finite-volume construction, global star compatibility, or dimension improvement.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text());print('PASS saved independent79 local-star review')
    else:p.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
