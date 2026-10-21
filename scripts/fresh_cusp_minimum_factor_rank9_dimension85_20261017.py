#!/usr/bin/env python3
"""Exact finite rank-partition comparison, inheriting the unrestricted163 proof."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_cusp_minimum_factor_rank9_dimension85_20261017'
INPUTS=['scripts/'+TAG+'.cpp',
 'results/fresh_cusp_factor_payment_framework_20261017.json',
 'results/fresh_cusp_factor_payment_framework_20261017.md',
 'results/fresh_cusp_link_shape_curvature_20261017.json',
 'results/fresh_cusp_link_shape_curvature_20261017.md',
 'results/fresh_nikulin_two_budget_frontier_20261015_dimension163.json',
 'results/fresh_nikulin_two_budget_frontier_20261015_dimension163.md']

def run():
    inherited=json.loads((ROOT/INPUTS[-2]).read_text())
    assert inherited['combined_candidate_dimension_at_most']==163
    framework=json.loads((ROOT/INPUTS[1]).read_text())
    assert framework['hypothesis']=='Every irreducible affine cusp factor has rank at least9.'
    shape=json.loads((ROOT/INPUTS[3]).read_text())
    assert shape['status']=='PASS_EXACT_CUSP_LINK_SHAPE_CURVATURE'
    # Integer-safety proof: V^2 <= 2^sum(r) <= 2^162, so V<=2^81.
    # h_j<=V; G<=81V; all C++ additions/comparisons have magnitude<2^94.
    assert all((r+1)**2<=2**r for r in range(9,163))
    with tempfile.TemporaryDirectory(prefix='coxeter_rank9_exact_') as tmp:
        exe=Path(tmp)/'enumerate'
        subprocess.run(['clang++','-O3','-std=c++17',str(ROOT/INPUTS[0]),'-o',str(exe)],check=True,capture_output=True,text=True)
        raw=subprocess.run([str(exe)],check=True,capture_output=True,text=True)
        rows=json.loads(raw.stdout)
    # Independently count integer partitions by the ordinary coin recurrence.
    counts=[0]*163;counts[0]=1
    for r in range(9,163):
        for d in range(r,163):counts[d]+=counts[d-r]
    assert [r['ambient_dimension'] for r in rows]==list(range(10,164))
    for row in rows:
        d=row['ambient_dimension']-1
        assert row['partitions']==counts[d]
        if d>=85:assert int(row['minimum_integer_scaled_gap'])>0
    tested=sum(r['partitions'] for r in rows if r['ambient_dimension']>=86)
    assert sum(counts[9:])==50378399
    finite=[]
    for n in range(86,164):
        beta=F(n*(n-6),8) if n%2==0 else F((n-1)*(n-5),8)
        margin=beta-7*n;assert margin>0
        finite.append([n,str(margin)])
    # Exact failure of this coarser pointwise comparison one dimension below.
    h=[1]
    for r in [14]*6:
        nxt=[0]*(len(h)+r)
        for i,c in enumerate(h):
            for j in range(r+1):nxt[i+j]+=c
        h=nxt
    V=sum(h);G=sum(c*min(i,84-i)for i,c in enumerate(h))
    obstruction=F(5*G-84*V,2*V)-42
    assert obstruction<0
    paths=[Path(__file__)]+[ROOT/p for p in INPUTS]
    return dict(status='PASS_EXACT_MINIMUM_CUSP_FACTOR_RANK9_DIMENSION85',
      sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
      dimension_bound=85,first_excluded_dimension=86,inherited_unrestricted_bound=163,
      minimum_irreducible_cusp_factor_rank=9,maximum_cusp_factor_count=None,
      arithmetic_hypothesis=False,
      scope='Every finite-sided finite-volume real hyperbolic Coxeter polytope whose irreducible affine cusp factors all have rank at least9. Compact polytopes satisfy the condition vacuously. No conclusion for a polytope with a smaller cusp factor.',
      comparison='For every rank partition of n-1 with parts>=9 and86<=n<=163, kappa(L)/V - sum_i min(7,r_i/2) > 1/5.',
      integer_test='5*(5*G-(d+sum_i min(14,r_i))*V)-2*V > 0',
      enumeration_rows=rows,enumerated_nonempty_partitions=sum(counts[9:]),
      excluded_range_partitions=tested,independent_partition_count_method='unbounded coin recurrence, ranks9through162',
      integer_safety='V<=2^81; every coefficient<=V;G<=81V;all intermediate integer magnitudes<2^94, inside signed128.',
      finite_vertex_margins=finite,
      dimension85_coarse_comparison_obstruction=dict(ranks=[14]*6,margin=str(obstruction),scope='Local curvature/payment model only; no actual ambient chamber asserted.'),
      formalization=False,external_refereeing=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();path=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(path.read_text())
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['enumerated_nonempty_partitions'],'partitions;',d['excluded_range_partitions'],'in excluded range',flush=True)
