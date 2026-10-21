#!/usr/bin/env python3
"""Independent pair-arithmetic check of the affine-A height rectangles."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='review_simple_cusp_affine_a_quadrilateral_20261017'
SOURCE='fresh_simple_cusp_affine_a_quadrilateral_20261017'


def add(x,y):return x[0]+y[0], x[1]+y[1]
def neg(x):return -x[0],-x[1]
def mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def scale(x,a):return a*x[0],a*x[1]
def nonnegative(x):
    a,b=x
    if a>=0 and b>=0:return True
    if a<0 and b<=0:return False
    return a*a>=2*b*b if a>=0 else 2*b*b>=a*a
def pair(x):return tuple(map(F,x))


def run():
    doc=json.loads((ROOT/f'results/{SOURCE}.json').read_text())
    assert doc['status']=='PASS_ALL_RANK_AFFINE_A_QUADRILATERAL_PAYMENT'
    zero=(F(),F()); one=(F(1),F())
    records=[]
    for c,d in product(doc['profiles'],repeat=2):
        for upper_u,upper_v in product((0,1),repeat=2):
            u=pair(c['umax']) if upper_u else zero
            v=pair(d['umax']) if upper_v else zero
            L=[add(mul(u,pair(c['moment'][i])),mul(v,pair(d['moment'][i]))) for i in range(2)]
            quadratic=add(mul(pair(c['c0']),mul(u,u)),mul(pair(d['c0']),mul(v,v)))
            difference=add(L[0],neg(L[1]))
            value=add(one,neg(add(scale(add(L[0],L[1]),F(1,2)),
                       add(scale(quadratic,F(1,2)),scale(mul(difference,difference),F(1,16))))))
            assert nonnegative(value)
            records.append(dict(c=c['name'],d=d['name'],upper_u=upper_u,upper_v=upper_v,
                                F=list(map(str,value))))
    assert records==doc['rectangle_checks']
    assert len(records)==256
    zeros=[r for r in records if pair(r['F'])==zero]
    assert zeros==doc['zero_checks'] and len(zeros)==3
    # Independently read each saved literal affine Gram and test its positive kernel.
    kernels=[]
    for w in doc['affine_exclusion_witnesses']:
        k=list(map(pair,w['positive_kernel'])); n=len(k)
        assert all(nonnegative(x) and x!=zero for x in k)
        products=[scale(x,F(2)) for x in k]
        for i,j,bond in w['edges']:
            products[i]=add(products[i],neg(mul(pair(bond),k[j])))
            products[j]=add(products[j],neg(mul(pair(bond),k[i])))
        assert all(x==zero for x in products)
        kernels.append(dict(name=w['name'],order=n,kernel_checked=True))
    paths=[f'scripts/{TAG}.py',f'scripts/{SOURCE}.py',f'results/{SOURCE}.json',f'results/{SOURCE}.md']
    return dict(status='PASS_INDEPENDENT_AFFINE_A_RECTANGLE_AND_KERNEL_REVIEW',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        rectangle_checks=len(records),zeros=zeros,positive_kernels=kernels,
        scope='Independent Fraction-pair arithmetic for256 saved rectangle corners and four positive kernels. The accompanying report separately reviews classification, height exclusions, Schur elimination and all-rank concavity/monotonicity. This does not duplicate the formal determinant expansion, prove the inherited global affine exclusion or assert a new ambient dimension bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    r=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert r==json.loads(out.read_text())
    else:out.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'],r['rectangle_checks'],len(r['positive_kernels']))
