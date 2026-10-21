#!/usr/bin/env python3
"""Separate Fraction witness replay; completeness is source-reviewed."""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import argparse,json

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_integer_inner_gap_20261014_review'
SOURCE='results/fresh_h20_integer_inner_gap_20261014.json'
def invert(A):
    n=len(A);B=[row[:]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        pivot=B[k][k];assert pivot>0
        B[k]=[x/pivot for x in B[k]]
        for j in range(n):
            if j!=k:
                p=B[j][k];B[j]=[x-p*y for x,y in zip(B[j],B[k])]
    return [r[n:] for r in B]
@lru_cache(None)
def matrix(n):
    A=[[F(i==j)-F(abs(i-j)==1,2) for j in range(n)] for i in range(n)]
    return A,invert(A)
def dot(A,x,y):return sum(x[i]*A[i][j]*y[j] for i in range(len(A)) for j in range(len(A)))
@lru_cache(None)
def reconstruct(gaps):
    n=len(gaps)-1;A,B=matrix(n)
    assert all(type(g) is int and g>=0 for g in gaps) and sum(gaps)==n+1
    p=[F(1-g) for g in gaps[1:]]
    c=[sum(x*y for x,y in zip(r,p)) for r in B]
    z=[sum(r) for r in B]
    assert sum(p)==gaps[0]-1
    return dot(A,c,c),dot(A,c,z),dot(A,z,z)
def run():
    d=json.loads((ROOT/SOURCE).read_text());paths={SOURCE}
    for p,h in d['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p;paths.add(p)
    seen=set();allowed=[]
    for row in d['records']:
        key=tuple(map(tuple,row['gaps']));assert key not in seen;seen.add(key)
        terms=[reconstruct(g) for g in key]
        q,p,zq=[sum(t[i] for t in terms) for i in range(3)]
        assert zq==196 and 28*q==row['norm28'] and 28*p==row['dot_vertex28']
        power=q-F(3,10)*p+F(9,400)*zq-4
        assert power*11200==row['power_numerator'] and power<=F(1,100)
        cosine=(q-8)/8
        compatible=cosine in [0,F(1,2)] or cosine>=1
        assert compatible==row['central_compatible']
        if compatible:
            assert power==F(1,100) and cosine==0
            allowed.append(dict(gaps=row['gaps'],center_square=str(q),
                                center_dot_vertex=str(p),power=str(power)))
    assert len(seen)==651 and len(allowed)==2
    x2=F(9,400)*196;D=F(d['escape']['D']);qE=x2-2*D
    assert qE==F(1,200) and D==F(881,400) and -100-D<0
    assert F(1,100)>qE and x2-4==F(41,100)>qE
    assert 196*F(17,20)**2-121==F(2061,100)>qE
    assert F(17,20)==F(d['escape']['minimum_seed_gap'])
    paths.add(str(Path(__file__).resolve().relative_to(ROOT)))
    return dict(status='PASS_SOURCE_REVIEW_AND_SEPARATE_FRACTION_INNER_GAP_WITNESSES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)},
        saved_recipes_reconstructed=651,central_compatible_witnesses=allowed,
        escape_norm='1/200',seed_gap='17/20',integer_cap_power_lower_bound='1/100',
        central_power='41/100',endpoint_radius_at_most=11,endpoint_power='2061/100',
        scope='Source/completeness review plus fresh primary saved replay; independent Fraction '
              'Gaussian inverses reconstruct all651 archived small-power recipes and both '
              'attaining witnesses. No second enumeration of omitted factor compositions.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    d=run();p=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['saved_recipes_reconstructed'])
