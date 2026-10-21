#!/usr/bin/env python3
"""Independent Qsqrt5 radical gauge and complete principal-core replay.

No producer imports. Pair coefficients use 1,sqrt5, not the golden-ratio basis.
Edges are reconstructed by exact quadratic square roots of the required square,
not by the producer's radical-mask exponent formula. Cycle constraints use DSU.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_sqrt5_quad_filters_20261015';TAG=BASE+'_review'
Z=(F(0),F(0));ONE=(F(1),F(0));TWO=(F(2),F(0))
RADS=((F(2),F(0)),(F(3),F(0)),(F(5,2),F(1,2)))
BITS={3:0,4:1,5:0,6:2,10:4}
SQUARES={3:ONE,4:TWO,5:(F(3,2),F(1,2)),6:(F(3),F(0)),10:RADS[2]}
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def mul(x,y):return x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def div(x,y):
    d=y[0]*y[0]-5*y[1]*y[1];assert d
    z=mul(x,(y[0],-y[1]));return z[0]/d,z[1]/d

def sign(x):
    a,b=x
    if not b:return (a>0)-(a<0)
    if not a or (a>0)==(b>0):return (b>0)-(b<0)
    v=a*a-5*b*b;return ((a>0)-(a<0))*((v>0)-(v<0))
def conj(x):return x[0],-x[1]
def unpack(x):a,b=map(F,x);return a+b/2,b/2
def enc(x):return [str(x[0]),str(x[1])]
def rroot(x):
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return F(a,b) if a*a==x.numerator and b*b==x.denominator else None

def square_root(x):
    assert sign(x)>0
    if not x[1]:
        a=rroot(x[0])
        if a is not None:return a,F(0)
        b=rroot(x[0]/5)
        if b is not None:return F(0),b
    norm=rroot(x[0]*x[0]-5*x[1]*x[1]);assert norm is not None
    for s in (1,-1):
        a=rroot((x[0]+s*norm)/2)
        if a:
            y=(a,x[1]/(2*a))
            if mul(y,y)==x:return y if sign(y)>0 else neg(y)
    raise AssertionError(('not a quadratic square',x))

def coloring(row):
    n=row['nodes'];parent=list(range(n));offset=[0]*n
    def find(i):
        if parent[i]==i:return i,0
        r,k=find(parent[i]);return r,k^offset[i]
    for (i,j),m in zip(combinations(range(n),2),row['labels']):
        if m==2:continue
        a,ca=find(i);b,cb=find(j);delta=ca^cb^BITS[m]
        if a==b:
            if delta:return None
        else:parent[a]=b;offset[a]=delta
    roots={};colors=[]
    for i in range(n):
        r,c=find(i)
        if r not in roots:roots[r]=c
        colors.append(c^roots[r])
    return colors

def matrix(row,colors):
    scale=[]
    for c in colors:
        value=ONE
        for bit,r in enumerate(RADS):
            if c>>bit&1:value=mul(value,r)
        assert sign(value)>0 and sign(conj(value))>0
        scale.append(value)
    n=row['nodes'];M=[[mul(TWO,scale[i])if i==j else Z for j in range(n)]for i in range(n)]
    for (i,j),m in zip(combinations(range(n),2),row['labels']):
        if m==2:continue
        value=square_root(mul(SQUARES[m],mul(scale[i],scale[j])))
        M[i][j]=M[j][i]=neg(value)
    return M

def inverse(A):
    n=len(A);B=[row[:]+[ONE if i==j else Z for j in range(n)]for i,row in enumerate(A)]
    for j in range(n):
        k=next(i for i in range(j,n)if B[i][j]!=Z)
        B[k],B[j]=B[j],B[k];d=B[j][j];B[j]=[div(x,d)for x in B[j]]
        for i in range(n):
            if i!=j:
                a=B[i][j];B[i]=[sub(x,mul(a,y))for x,y in zip(B[i],B[j])]
    return [r[n:]for r in B]

def inertia(M):
    A=[r[:]for r in M];counts=[0,0,0]
    while A:
        n=len(A);k=next((i for i in range(n)if A[i][i]!=Z),None)
        if k is None:
            edge=next(((i,j)for i,j in combinations(range(n),2)if A[i][j]!=Z),None)
            if edge is None:counts[2]+=n;break
            i,j=edge
            # Congruence replacing e_i by e_i+e_j makes a nonzero diagonal.
            for h in range(n):A[i][h]=add(A[i][h],A[j][h])
            for h in range(n):A[h][i]=add(A[h][i],A[h][j])
            k=i
        d=A[k][k];counts[0 if sign(d)>0 else 1]+=1
        I=[i for i in range(n)if i!=k]
        A=[[sub(A[i][j],div(mul(A[i][k],A[k][j]),d))for j in I]for i in I]
    return counts

def residual(M,C):
    J=inverse([[M[i][j]for j in C]for i in C]);O=[i for i in range(len(M))if i not in C]
    S=[]
    for i in O:
        row=[]
        for j in O:
            v=M[i][j]
            for a,u in enumerate(C):
                for b,w in enumerate(C):v=sub(v,mul(M[i][u],mul(J[a][b],M[w][j])))
            row.append(v)
        S.append(row)
    return S

def run():
    primary=json.loads((ROOT/'results'/f'{BASE}.json').read_text())
    union=json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    counts=Counter();rows=[];radical_edges=0
    for j,cl in enumerate(union['marked_isomorphism_classes']):
        row=union['labelled_union_supports'][cl['representative']];n=row['nodes'];saved=primary['records'][j]
        assert saved['class_index']==j and saved['union_index']==cl['representative']
        if 0 in row['labels']:status='DOTTED_PARAMETER_UNFILTERED';rec={}
        else:
            colors=coloring(row)
            assert (colors is not None)==saved['cyclic_field']['passes']
            if colors is None:status='REJECT_CYCLE_FIELD';rec={}
            else:
                assert colors==saved['cyclic_field']['colors']
                M=matrix(row,colors);radical_edges+=sum(m!=2 for m in row['labels'])
                C=row['white']+[0,1];corner=[[M[i][k]for k in C]for i in C]
                assert inertia(corner)==[n-2,0,0]
                assert inertia([[conj(x)for x in r]for r in corner])==[n-2,0,0]
                S=residual(M,C);assert len(S)==2
                for key,conjugated in [('physical',False),('conjugate',True)]:
                    A=[[conj(x)if conjugated else x for x in r]for r in S]
                    d=sub(mul(A[0][0],A[1][1]),mul(A[0][1],A[1][0]))
                    assert A==[[unpack(x)for x in r]for r in saved[key]['schur']]
                    assert d==unpack(saved[key]['determinant'])
                    assert [sign(A[0][0]),sign(A[1][1]),sign(d)]==saved[key]['signs']
                phys=inertia(M);other=inertia([[conj(x)for x in r]for r in M])
                if phys[1]>=2:status='REJECT_TWO_NEGATIVE_DIRECTIONS'
                elif other[1]:status='REJECT_CONJUGATE_NOT_PSD'
                else:status='PASS_PRINCIPAL_UNION_NECESSARY_TESTS'
                rec=dict(physical_inertia=phys,conjugate_inertia=other)
                if status=='PASS_PRINCIPAL_UNION_NECESSARY_TESTS':
                    assert M==[[unpack(x)for x in r]for r in saved['scaled_gram']]
                    assert saved['singular']==bool(phys[2])
        assert status==saved['status'],(j,status,saved['status']);counts[status]+=1
        rows.append(dict(class_index=j,status=status,**rec))
    assert dict(counts)==primary['counts']
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',
           'scripts/fresh_sqrt5_quad_exhaustion_20261015.py','results/fresh_sqrt5_quad_exhaustion_20261015.json',
           'results/fresh_sqrt5_lanner_catalogue_20261015_review.json']
    return dict(status='PASS_INDEPENDENT_QSQRT5_PRINCIPAL_CORE_GAUGE_REVIEW',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                counts=dict(counts),reconstructed_nonright_entries=radical_edges,records=rows,
                arithmetic='Standalone Fraction basis1,sqrt5; DSU parity; exact quadratic square roots; full Gaussian congruence and inverse Schur',
                scope='All1180 principal core filters; no independent union enumeration here, no external white/rank completion and no new dimension bound.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],out['reconstructed_nonright_entries'])
