#!/usr/bin/env python3
"""Exact product-link curvature accounting and a scoped coarse-payment obstruction."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_cusp_link_shape_curvature_20261017'
INPUTS=['results/fresh_cusp_factor_payment_framework_20261017.json',
        'results/fresh_cusp_factor_payment_framework_20261017.md',
        'results/fresh_simple_cusp_distance_bound_20261017.json',
        'results/fresh_finite_volume_simplex_truncation_20261017.md']

def product_link(ranks):
    h=[1]
    for r in ranks:
        assert isinstance(r,int) and r>=1
        nxt=[0]*(len(h)+r)
        for i,c in enumerate(h):
            for j in range(r+1):nxt[i+j]+=c
        h=nxt
    d=sum(ranks);V=sum(h);E=F(d*V,2)
    assert h==h[::-1] and h[0]==h[-1]==1
    assert E.denominator==1
    assert sum(i*c for i,c in enumerate(h))==E
    g=sum(c*min(i,d-i)for i,c in enumerate(h))
    k=F(5*g,2)-E
    absdev=sum(F(c*abs(2*i-d),2*V)for i,c in enumerate(h))
    assert k/V==F(3*d,4)-F(5,2)*absdev
    var=sum(F(r*(r+2),12)for r in ranks)
    directvar=sum(F(c*(2*i-d)**2,4*V)for i,c in enumerate(h))
    assert var==directvar and absdev**2<=var
    assert E-g==sum(c*(2*i-d)for i,c in enumerate(h)if 2*i>d)
    assert E-g>=d
    return dict(ranks=ranks,link_dimension=d,ambient_dimension=d+1,
        h_coefficients=h,vertices=V,edges=int(E),minimum_uncut_edge_bound=g,
        curvature=str(k),curvature_per_link_vertex=str(k/V),
        mean_absolute_deviation=str(absdev),variance=str(var),
        strict_global_extremum_slack=str(F(5,2)*(E-g)))

def run():
    framework=json.loads((ROOT/INPUTS[0]).read_text())
    for p,h in framework['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    rows=[]
    for d in range(2,181):
        a=product_link([d]);n=d+1
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert F(a['curvature'])==beta
        rows.append(dict(dimension=n,curvature=str(beta)))
    cases=[[1]*d for d in (3,8,16,25,50,100)]+[[6]*q for q in (5,7,10,20)]
    cases += [[30]*5+[13],[30]*5+[12],[30]*5+[10],[9,9],[17,17],[30,30],[30,30,30]]
    records=[product_link(x)for x in cases]
    blocked=next(x for x in records if x['ranks']==[30]*5+[13])
    cost=sum(min(F(r,2),F(15))for r in blocked['ranks'])
    margin=F(blocked['curvature_per_link_vertex'])-cost
    assert cost==F(163,2) and margin<0
    # This can be a literal product of affine A factors: on order31 cycles
    # every pair has distance<=15, and likewise on order14.
    paircounts=[]
    for r in (30,13):
        m=r+1;count=sum(min(j-i,m-(j-i))<=15 for i in range(m)for j in range(i+1,m))
        assert count==m*r//2
        paircounts.append(dict(factor_rank=r,cycle_order=m,paid_pairs=count,
                               cost_per_link_vertex=str(F(count,m))))
    paths=[Path(__file__)]+[ROOT/p for p in INPUTS]
    return dict(status='PASS_EXACT_CUSP_LINK_SHAPE_CURVATURE',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        global_inequality='D5 >= beta_n*v_f + sum_ideal_v ((5/2)*sum_j h_vj*min(j,n-1-j) - (n-1)*V_v/2) + 5*(n-1)',
        h_polynomial='product_i(1+t+...+t^r_i)',
        simplex_identity_checks=rows,product_link_checks=records,
        obstruction=dict(factor_ranks=blocked['ranks'],ambient_dimension=164,
            curvature_per_link_vertex=blocked['curvature_per_link_vertex'],
            distance15_cost_per_link_vertex=str(cost),margin=str(margin),
            cycle_pair_counts=paircounts,
            scope='Local product-link curvature/payment comparison only. The link has a standard Euclidean affine-product realization; no finite-volume ambientH164 chamber is asserted. The strict global-extremum slack is positive but nonuniform relative to arbitrary cusp counts and does not validate a per-vertex domination with this negative coefficient.'),
        scope='Exact stronger global combinatorial inequality and local obstruction to one coarse per-cusp payment comparison. No new dimension theorem or optimizer.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/'results'/f'{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['simplex_identity_checks']),'simplex identities;',len(d['product_link_checks']),'product links;',d['obstruction']['margin'],flush=True)
