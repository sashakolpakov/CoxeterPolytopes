#!/usr/bin/env python3
"""Independent count/scalar/overflow review; does not repeat the C++ enumeration."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='review_cusp_minimum_factor_rank9_dimension85_20261017'
PRIMARY='fresh_cusp_minimum_factor_rank9_dimension85_20261017'


def run():
    data=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert data['status']=='PASS_EXACT_MINIMUM_CUSP_FACTOR_RANK9_DIMENSION85'
    for p,h in data['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    # A partition into exactly k ranks>=9, after subtracting9 from each,
    # is a partition of d-9k into at most k parts. Conjugate Ferrers diagrams
    # identify this with partitions whose largest part is at most k.
    N=162;small=[[0]*(N+1)for _ in range(N//9+1)];small[0][0]=1
    for k in range(1,len(small)):
        for e in range(N+1):small[k][e]=small[k-1][e]+(small[k][e-k]if e>=k else 0)
    counts=[sum(small[k][d-9*k]for k in range(1,d//9+1))for d in range(N+1)]
    rows=data['enumeration_rows'];assert len(rows)==154
    for n,row in zip(range(10,164),rows):
        assert row['ambient_dimension']==n and row['partitions']==counts[n-1]
        if n>=86:assert int(row['minimum_integer_scaled_gap'])>0
        else:assert int(row['minimum_integer_scaled_gap'])==0
    assert sum(counts)==data['enumerated_nonempty_partitions']==50378399
    assert sum(counts[85:])==data['excluded_range_partitions']==50330195
    # All operations are bounded even if differences are estimated in absolute value.
    V=2**81;G=81*V;max_twice_cost=162
    bounds=dict(link_vertices=V,weighted_coefficient_sum=G,window_addition=2*V,
        preassert_child_product=163*V,absolute_diff=(5*81+324)*V,
        absolute_gap=(5*(5*81+324)+2)*V)
    assert all(v<2**94 for v in bounds.values()) and 2**94<2**127
    assert max(counts)<2**64
    assert all((r+1)**2<=2**r for r in range(9,163))
    margins=[]
    for n in range(86,164):
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        a=beta-7*n;assert a>0;margins.append([n,str(a)])
    assert margins==data['finite_vertex_margins']
    # A separate direct coefficient multiplication checks the stated lower-neighbor obstruction.
    poly={0:1}
    for _ in range(6):
        nxt={}
        for i,a in poly.items():
            for j in range(15):nxt[i+j]=nxt.get(i+j,0)+a
        poly=nxt
    vert=sum(poly.values());g=sum(a*min(j,84-j)for j,a in poly.items())
    obstruction=F(5*g-84*vert,2*vert)-42
    assert obstruction==F(-70391,253125)
    assert str(obstruction)==data['dimension85_coarse_comparison_obstruction']['margin']
    paths=[Path(__file__),ROOT/f'scripts/{PRIMARY}.cpp',ROOT/f'scripts/{PRIMARY}.py',ROOT/f'results/{PRIMARY}.json',
        ROOT/'results/fresh_cusp_link_shape_curvature_20261017.md',ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    return dict(status='PASS_INDEPENDENT_RANK9_CUSP_DIMENSION85_COUNT_SCALAR_REVIEW',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        count_method='Exactly k ranks>=9 -> Ferrers conjugation of partitions of d-9k with at most k parts; independent largest-part dynamic program.',
        per_link_dimension_counts=counts,nonempty_partitions=sum(counts),excluded_range_partitions=sum(counts[85:]),
        maximum_per_dimension_count=max(counts),integer_bounds=bounds,finite_vertex_margins=margins,
        lower_neighbor_obstruction=str(obstruction),producer_enumeration_repeated=False,
        conclusion_dimension_bound=85,hypothesis='Every irreducible affine cusp factor has rank at least9.',
        scope='Independent mathematical/source review plus exact partition counts, integer bounds and scalar composition checks. Inherits the producer full integer-convolution enumeration and frozen geometric distance-seven/link-curvature/unrestricted163 theorems. No external refereeing, formalization or novelty claim.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['nonempty_partitions'],d['excluded_range_partitions'],flush=True)
