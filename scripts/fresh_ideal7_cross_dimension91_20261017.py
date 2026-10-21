#!/usr/bin/env python3
"""Exact cap15 curvature transfer with finite cross-component half-payments."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb,prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ideal7_cross_dimension91_20261017'
BASE='fresh_distance8_ideal_dimension97_20261017'
LOCAL='fresh_ideal7_finite_cross_payments_20261017'
MIXED='fresh_ideal_distance7_mixed_replacement_20261017'

def digest(path):
    h=sha256()
    with path.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def profile(ranks):
    d=sum(ranks);h=[1]
    for r in ranks:
        out=[0]*(len(h)+r)
        for i,a in enumerate(h):
            for j in range(r+1):out[i+j]+=a
        h=out
    V=prod(r+1 for r in ranks);assert sum(h)==V and h==h[::-1]
    G=sum(a*min(j,d-j)for j,a in enumerate(h))
    curvature=F(5*G-d*V,2*V);cost=sum((min(F(15,2),F(r,2))for r in ranks),F())
    moment=0;s=len(ranks)
    for k in range(s+1):
        for chosen in combinations(ranks,k):
            shift=sum(r+1 for r in chosen);m=(d-1)//2-shift
            if m>=0:moment+=(-1)**k*((d-2*shift)*comb(m+s,s)-2*s*comb(m+s,s+1))
    assert curvature==F(3*d,4)-F(5*moment,2*V)
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=V,h_polynomial=h,G=G,
        normalized_curvature=str(curvature),normalized_cost=str(cost),margin=str(curvature-cost))

def scalars():
    rows=[];minima=[]
    for d in range(91,97):
        shapes=set()
        for copies in range(d//15+1):
            for low in range(15):
                high=d-15*copies-low
                if high==0 or high>=16:shapes.add(tuple(sorted([15]*copies+([low]if low else[])+([high]if high else[]))))
        group=[profile(rs)for rs in sorted(shapes)];assert all(F(r['margin'])>0 for r in group)
        rows+=group;n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>9*n
        smallest=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=smallest['margin'],minimizing_ranks=smallest['ranks'],finite_slack=str(beta-9*n)))
    assert min(F(r['margin'])for r in rows)==F(10121,4194304)
    boundary=profile((15,)*6);assert F(boundary['margin'])==F(-428955,2097152)
    return dict(reduced_profiles=rows,dimension_minima=minima,boundary_at_n91=boundary)

def run(scalar_only=False):
    out=scalars()
    if scalar_only:
        return dict(status='PASS_CONDITIONAL_CAP15_SCALARS',profile_count=len(out['reduced_profiles']),dimension_minima=out['dimension_minima'],boundary=out['boundary_at_n91']['margin'])
    inputs={}
    for stem in(BASE,BASE+'_review',LOCAL,MIXED):
        d=json.loads((ROOT/f'results/{stem}.json').read_text());assert d['status'].startswith('PASS')
        for rel,h in d['sources'].items():assert digest(ROOT/rel)==h,(stem,rel)
        inputs[stem]=d
    assert inputs[BASE]['dimension_bound']==97 and inputs[BASE+'_review']['dimension_bound']==97
    assert inputs[LOCAL]['status']=='PASS_IDEAL7_SELECTIVE17_AND_FINITE_CROSS_COMPONENT_PAYMENTS'
    assert{t:r['cost']for t,r in inputs[LOCAL]['order17']['targets'].items()}==dict(A=119,B=111,C=108,D=127)
    assert all(F(r['max_excess_over_9n'])<=0 for r in inputs[LOCAL]['finite_cost']['rows'])
    assert{r['dimension']for r in inputs[LOCAL]['finite_cost']['rows']}>=set(range(92,98))
    assert inputs[MIXED]['mixed_replacement_complete']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in(BASE,BASE+'_review',LOCAL,MIXED)for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths +=[ROOT/f'results/{s}.md'for s in('fresh_cusp_link_shape_curvature_20261017','fresh_cusp_factor_payment_framework_20261017','fresh_distance14_unrestricted_dimension156_20261017','fresh_selective_distance9_dimension108_20261017')]
    out.update(status='PASS_FINITE_CROSS_IDEAL7_DIMENSION91',sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        dimension_bound=91,inherited_dimension_bound=97,excluded_dimensions=list(range(92,98)),
        unrestricted_finite_volume=True,arithmetic_hypothesis=False,cusp_restriction=None,
        finite_vertex_budget='9n',ideal_normalized_budget='sum_i min(15/2,r_i/2)',convex_order_cap=15,
        unchanged_compact_facet_ambient_bound=79,
        scope='Complete dimension transfer from the bound local face-payment/cost theorem, inherited actual cusp factor multiplicity, all-rank convex-order/shape proof and completed97 tail. External review, formalization and a new construction are not asserted.')
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--scalar-only',action='store_true');args=ap.parse_args()
    assert not(args.verify and args.scalar_only)
    d=run(args.scalar_only)
    if args.scalar_only:print(json.dumps(d,indent=2))
    else:
        p=ROOT/f'results/{TAG}.json'
        if args.verify:assert d==json.loads(p.read_text())
        else:p.write_text(json.dumps(d,indent=2)+'\n')
        print(d['status'],len(d['reduced_profiles']),'profiles',flush=True)
