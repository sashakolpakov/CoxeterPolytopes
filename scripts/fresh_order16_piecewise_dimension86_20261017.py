#!/usr/bin/env python3
"""Exact dimension86 transfer with a piecewise cusp-factor envelope."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb,prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order16_piecewise_dimension86_20261017'
BASE='fresh_ideal7_cross_dimension91_20261017'
LOCAL='fresh_order16_finite_cross_payments_20261017'
EMBED='fresh_order16_preferred_embeddings_20261017'
CYCLE='fresh_cycle16_half_coupon_transfer_20261017'

def digest(p):
    h=sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def budget(r):
    return F(r,2)if r<=14 else F(7)if r==15 else F(127,17)if r==16 else F(217,29)

def profile(ranks):
    d=sum(ranks);h=[1]
    for r in ranks:
        out=[0]*(len(h)+r)
        for i,a in enumerate(h):
            for j in range(r+1):out[i+j]+=a
        h=out
    V=prod(r+1 for r in ranks);assert sum(h)==V and h==h[::-1]
    G=sum(a*min(j,d-j)for j,a in enumerate(h));curvature=F(5*G-d*V,2*V)
    cost=sum(map(budget,ranks),F());moment=0;s=len(ranks)
    for k in range(s+1):
        for chosen in combinations(ranks,k):
            shift=sum(r+1 for r in chosen);m=(d-1)//2-shift
            if m>=0:moment+=(-1)**k*((d-2*shift)*comb(m+s,s)-2*s*comb(m+s,s+1))
    assert curvature==F(3*d,4)-F(5*moment,2*V)
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=V,h_polynomial=h,G=G,
        normalized_curvature=str(curvature),normalized_cost=str(cost),margin=str(curvature-cost))

def scalars():
    rows=[];minima=[]
    for d in range(86,91):
        shapes=set()
        for a in range(d//14+1):
            for b in range((d-14*a)//16+1):
                for c in range((d-14*a-16*b)//17+1):
                    for middle in(0,15):
                        for low in range(14):
                            high=d-14*a-16*b-17*c-middle-low
                            if high==0 or high>=18:
                                shapes.add(tuple(sorted([14]*a+[16]*b+[17]*c+([middle]if middle else[])+([low]if low else[])+([high]if high else[]))))
        group=[profile(rs)for rs in sorted(shapes)];assert all(F(r['margin'])>0 for r in group)
        rows+=group;n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>9*n
        smallest=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=smallest['margin'],minimizing_ranks=smallest['ranks'],finite_slack=str(beta-9*n)))
    assert len(rows)==3355 and min(F(r['margin'])for r in rows)==F(146417,911250)
    boundary=profile((1,)+(14,)*6);assert F(boundary['margin'])==F(-112409,1518750)
    # This cost increase is in the helpful direction: curvature-minus-cost
    # cannot increase when two ranks15 are replaced by14 and16.
    assert budget(14)+budget(16)-2*budget(15)==F(8,17)>0
    return dict(reduced_profiles=rows,dimension_minima=minima,boundary_at_n86=boundary,
        concentration=dict(repeated_cap_ranks=[14,16,17],at_most_one_low_part=list(range(1,14)),
            at_most_one_middle_part=15,at_most_one_high_part_minimum=18,
            double15_budget_increase='8/17',cost_direction='preserved or increased; curvature decreases'))

def run(scalar_only=False):
    out=scalars()
    if scalar_only:return dict(status='PASS_CONDITIONAL_PIECEWISE_86_SCALARS',profile_count=len(out['reduced_profiles']),dimension_minima=out['dimension_minima'],boundary=out['boundary_at_n86']['margin'])
    inputs={}
    for stem in(BASE,BASE+'_review',LOCAL,EMBED,CYCLE):
        d=json.loads((ROOT/f'results/{stem}.json').read_text());assert d['status'].startswith('PASS')
        for rel,h in d['sources'].items():assert digest(ROOT/rel)==h,(stem,rel)
        inputs[stem]=d
    assert inputs[BASE]['dimension_bound']==inputs[BASE+'_review']['dimension_bound']==91
    assert inputs[LOCAL]['order16_costs']==dict(A=112,B=107,C=103,D=111)
    assert inputs[LOCAL]['ideal_rank_envelope']=={'1_to_14':'r/2','15':'7','16':'127/17','17_and_larger':'217/29'}
    assert all(F(r['maximum_excess'])<=0 for r in inputs[LOCAL]['finite_cost']['rows'])
    assert {r['dimension']for r in inputs[LOCAL]['finite_cost']['rows']}>=set(range(87,92))
    assert inputs[CYCLE]['equality_core']['ideal_factor_cost']==112
    assert len(inputs[EMBED]['patterns'])==260
    stems=(BASE,BASE+'_review',LOCAL,EMBED,CYCLE)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths+=[ROOT/f'results/{s}.md'for s in('fresh_cusp_link_shape_curvature_20261017','fresh_cusp_factor_payment_framework_20261017','fresh_distance14_unrestricted_dimension156_20261017')]
    out.update(status='PASS_ORDER16_PIECEWISE_IDEAL_DIMENSION86',sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        dimension_bound=86,inherited_dimension_bound=91,excluded_dimensions=list(range(87,92)),
        unrestricted_finite_volume=True,arithmetic_hypothesis=False,cusp_restriction=None,
        finite_vertex_budget='9n',ideal_rank_envelope=inputs[LOCAL]['ideal_rank_envelope'],
        unchanged_compact_facet_ambient_bound=79,
        scope='Complete piecewise-envelope dimension transfer from the bound local order16 payment/cost theorem, actual cusp-factor multiplicity, all-rank convex-order/shape theorem and completed91 tail. Independent review and external/formal verification have separately stated scopes.')
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--scalar-only',action='store_true');args=ap.parse_args()
    assert not(args.verify and args.scalar_only);d=run(args.scalar_only)
    if args.scalar_only:print(json.dumps(d,indent=2))
    else:
        p=ROOT/f'results/{TAG}.json'
        if args.verify:assert d==json.loads(p.read_text())
        else:p.write_text(json.dumps(d,indent=2)+'\n')
        print(d['status'],len(d['reduced_profiles']),'profiles',flush=True)
