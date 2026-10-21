#!/usr/bin/env python3
"""Exact dimension80 transfer with a piecewise cusp-factor envelope."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from math import comb,prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ideal6_piecewise_dimension80_20261017'
BASE='fresh_order16_piecewise_dimension86_20261017'
LOCAL='fresh_ideal6_combined_payments_20261017'
EMBED='fresh_ideal6_small_factor_frontier_20261017'
D14='fresh_affine_d14_long_pair_transfer_20261017'
FINITE='fresh_broad_finite_half_coupons_20261017'

def digest(p):
    h=sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def budget(r):
    return F(r,2)if r<=13 else F(20,3)if r==14 else F(111,16)if r==15 else F(127,17)if r==16 else F(95,13)

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
    for d in range(80,86):
        shapes=set()
        for a in range(d//13+1):
            for b in range((d-13*a)//16+1):
                for c in range((d-13*a-16*b)//17+1):
                    for m14,m15,low in product((0,14),(0,15),range(13)):
                        high=d-13*a-16*b-17*c-m14-m15-low
                        if high==0 or high>=18:
                            shapes.add(tuple(sorted([13]*a+[16]*b+[17]*c+([m14]if m14 else[])+([m15]if m15 else[])+([low]if low else[])+([high]if high else[]))))
        group=[profile(rs)for rs in sorted(shapes)];assert all(F(r['margin'])>0 for r in group)
        rows+=group;n=d+1;beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>9*n
        smallest=min(group,key=lambda r:F(r['margin']))
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=smallest['margin'],minimizing_ranks=smallest['ranks'],finite_slack=str(beta-9*n)))
    assert [r['profiles']for r in minima]==[742,776,813,852,892,933]
    assert len(rows)==5008 and min(F(r['margin'])for r in rows)==F(940675,11294304)
    boundary=profile((1,)+(13,)*6);assert F(boundary['margin'])==F(-1131729,7529536)
    assert budget(13)+budget(15)-2*budget(14)==F(5,48)>0
    assert budget(14)+budget(16)-2*budget(15)==F(107,408)>0
    return dict(reduced_profiles=rows,dimension_minima=minima,boundary_at_n80=boundary,
        concentration=dict(repeated_cap_ranks=[13,16,17],at_most_one_low_part=list(range(1,13)),
            at_most_one_each_middle_part=[14,15],at_most_one_high_part_minimum=18,
            double14_budget_increase='5/48',double15_budget_increase='107/408',
            cost_direction='preserved or increased; curvature decreases',termination='Every move strictly increases the integer sum of squared ranks at fixed total rank.'))

def run(scalar_only=False):
    out=scalars()
    if scalar_only:return dict(status='PASS_CONDITIONAL_PIECEWISE_80_SCALARS',profile_count=len(out['reduced_profiles']),dimension_minima=out['dimension_minima'],boundary=out['boundary_at_n80']['margin'])
    inputs={}
    stems=(BASE,BASE+'_review',LOCAL,EMBED,D14,FINITE)
    for stem in stems:
        d=json.loads((ROOT/f'results/{stem}.json').read_text());assert d['status'].startswith('PASS')
        for rel,h in d['sources'].items():assert digest(ROOT/rel)==h,(stem,rel)
        inputs[stem]=d
    assert inputs[BASE]['dimension_bound']==inputs[BASE+'_review']['dimension_bound']==86
    assert inputs[LOCAL]['ideal_rank_envelope']=={'1_to_13':'r/2','14':'20/3','15':'111/16','16':'127/17','17_and_larger':'95/13'}
    assert inputs[LOCAL]['finite_budget_dimension_maximum']==inputs[FINITE]['dimension_maximum']==108
    assert len(inputs[EMBED]['patterns'])==437
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths+=[ROOT/f'results/{s}.md'for s in('fresh_cusp_link_shape_curvature_20261017','fresh_cusp_factor_payment_framework_20261017','fresh_distance14_unrestricted_dimension156_20261017')]
    out.update(status='PASS_IDEAL6_PIECEWISE_DIMENSION80',sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        dimension_bound=80,inherited_dimension_bound=86,excluded_dimensions=list(range(81,87)),
        unrestricted_finite_volume=True,arithmetic_hypothesis=False,cusp_restriction=None,
        finite_vertex_budget='9n',ideal_rank_envelope=inputs[LOCAL]['ideal_rank_envelope'],
        unchanged_compact_facet_ambient_bound=79,
        scope='Complete dimension transfer from the combined finite9/ideal6 local payments, actual cusp-factor multiplicity, all-rank concentration/shape theorem and completed86 tail. Independent arithmetic review and external/formal review have separately stated scopes.')
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
