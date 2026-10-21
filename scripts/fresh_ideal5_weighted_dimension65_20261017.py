#!/usr/bin/env python3
"""Dimension66 exclusion from reviewed local transfers and exact budgets."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ideal5_weighted_dimension65_20261017'
BASE='fresh_ideal5_weighted_dimension66_20261017'
BC='fresh_order12_bc_unit_coupon_transfer_20261017'
B14='fresh_order14_affine_b_fork_projection_transfer_20261017'
LARGE='fresh_large_finite_d_leaf_half_transfer_20261017'
FINITE='fresh_a2_d11_full_coupon_forest_budget_20261017'
SCALAR='fresh_order14_b84_scalar_frontier_20261017'
SPECIAL={11:F(16,3),12:F(75,13),13:F(6),14:F(32,5),15:F(107,16),16:F(120,17)}


def cost(r):
    return F(r,2) if r<=10 else SPECIAL[r] if r<=16 else F(163,23)


def terminal_shapes(total):
    def visit(left,minimum,prefix,low,middle,high):
        if not left:
            yield prefix;return
        for r in range(minimum,left+1):
            if (r<10 and low) or (r==13 and middle) or (r>=18 and high):
                continue
            yield from visit(left-r,r,prefix+(r,),low or r<10,
                             middle or r==13,high or r>=18)
    return visit(total,1,(),False,False,False)


def profile(ranks):
    h=[1]
    for r in ranks:
        new=[0]*(len(h)+r)
        for j,v in enumerate(h):
            for k in range(r+1):new[j+k]+=v
        h=new
    d=sum(ranks);V=prod(r+1 for r in ranks)
    assert len(h)==d+1 and h==h[::-1] and sum(h)==V
    G=sum(v*min(j,d-j) for j,v in enumerate(h))
    curvature=F(5*G-d*V,2*V);budget=sum((cost(r) for r in ranks),F())
    return dict(ranks=list(ranks),link_dimension=d,link_vertices=V,
                h_polynomial=h,G=G,normalized_curvature=str(curvature),
                normalized_cost=str(budget),margin=str(curvature-budget))


def digest(p):return sha256(p.read_bytes()).hexdigest()


def run():
    stems=(BASE,BASE+'_review',BC,B14,LARGE,FINITE,SCALAR)
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert digest(ROOT/p)==h,(s,p)
    assert data[BASE]['dimension_bound']==data[BASE+'_review']['dimension_bound']==66
    assert data[BASE]['unrestricted_finite_volume'] and not data[BASE]['arithmetic_hypothesis']
    assert data[BASE]['cusp_restriction'] is None
    assert data[BC]['new_order12_costs']==dict(A=60,B=64,C=63,D=62)
    assert data[BC]['new_coupon_increment_per_A2_D11_component_pair']=='2'
    assert data[BC]['existing_A1_D11_full_upgrade_required'] and not data[BC]['large_D_D10_coupon_used']
    assert data[B14]['new_order14_costs']==dict(A=70,B=84,C=84,D=82)
    assert data[B14]['finite_weights_unchanged'] and not data[B14]['added_finite_coupons']
    assert data[LARGE]['routing']==dict(whole_component_at_most17=2531,
        protected_simple_leaf_neighbor_degree2=18,two_protected_leaf9_corners=40,
        new_small_component_D10_half=16,adaptive_simple_leaf_neighbor_degree2=9)
    assert all(data[s]['triangle_minimum']==2 and data[s]['quadrilateral_minimum']==1
               for s in (BC,LARGE))
    maxima=data[FINITE]['certificate']['dimension_maxima']
    assert data[FINITE]['doubled_representative_cross_formula']==(
        'a*(2*d10+3*tD11to16+2*a17)+2*a1*(d11+d12)+2*a4*d13+4*a2*d11')
    assert data[FINITE]['maximum_rank']>=66
    assert maxima[66]['doubled_cost']==990 and maxima[66]['witness_components']==['D66']
    assert F(maxima[66]['finite_slack'])==0
    assert cost(12)+cost(14)-2*cost(13)==F(11,65)
    assert cost(11)+cost(13)-2*cost(12)==-F(8,39)
    rows=[profile(rs) for rs in terminal_shapes(65)]
    assert len(rows)==len({tuple(r['ranks']) for r in rows})==1659
    assert all(F(r['margin'])>0 for r in rows)
    minimum=min(rows,key=lambda r:F(r['margin']))
    assert (minimum['margin'],minimum['ranks'])==('13367/495924',[10,11,12,16,16])
    # Independently written former conditional scalar interface now has its
    # actual local BC theorem and finite-cost accounting supplied above.
    old_rows={tuple(r['ranks']):r for r in data[SCALAR]['common_domain_rows']}
    assert {tuple(r['ranks']) for r in rows}=={tuple(p) for p in data[SCALAR]['new_terminal_ranks']}
    for r in rows:
        old=old_rows[tuple(r['ranks'])]
        assert F(r['margin'])==F(old['new_margin'])+F(r['ranks'].count(11),12)
    boundary_rows=[profile(rs) for rs in terminal_shapes(64)]
    bm=min(boundary_rows,key=lambda r:F(r['margin']))
    assert len(boundary_rows)==1528
    assert (bm['margin'],bm['ranks'])==('-78444/485537',[12,12,12,12,16])
    assert sum(F(r['margin'])<0 for r in boundary_rows)==16
    assert maxima[65]['doubled_cost']==974 and maxima[65]['witness_components']==['D65']
    small=dict(data[BASE]['small_order_costs'])
    small['12']=data[BC]['new_order12_costs'];small['14']=data[B14]['new_order14_costs']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in stems
                           for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/f'results/{s}.md' for s in (
        'fresh_cusp_factor_payment_framework_20261017',
        'fresh_cusp_link_shape_curvature_20261017',
        'fresh_distance14_unrestricted_dimension156_20261017')]
    return dict(status='PASS_UNRESTRICTED_WEIGHTED_IDEAL5_DIMENSION65',dimension_bound=65,
        inherited_dimension_bound=66,excluded_dimensions=[66],unrestricted_finite_volume=True,
        arithmetic_hypothesis=False,cusp_restriction=None,
        finite_budget_source=FINITE,finite_budget_dimension_maximum=data[FINITE]['maximum_rank'],
        ideal_rank_envelope={str(r):str(v) for r,v in SPECIAL.items()},
        small_order_costs=small,reduced_profiles=rows,
        dimension_minima=[dict(ambient_dimension=66,profiles=len(rows),
            minimum_margin=minimum['margin'],minimizing_ranks=minimum['ranks'],
            finite_maximum='495',finite_curvature_coefficient='495',finite_slack='0',
            finite_maximum_witness=['D66'],strict_global_slack=325)],
        concentration=dict(repeated_cap_ranks=[10,11,12,14,15,16,17],
            at_most_one_low_part=list(range(1,10)),at_most_one_each_middle_part=[13],
            at_most_one_high_part_minimum=18,double13_budget_increase='11/65',
            invalid_double12_budget_increase='-8/39',
            cost_direction='Preserved or increased; normalized curvature decreases.',
            termination='Every nontrivial move increases the bounded integer sum of squared ranks.'),
        boundary_at_n65=dict(profile=bm,profiles=len(boundary_rows),negative_profiles=16,
            finite_maximum='487',finite_curvature_coefficient='480',finite_slack='-7',
            finite_maximum_witness=['D65'],new_bound64_asserted=False),
        unchanged_compact_arithmetic_bound=21,unchanged_arbitrary_compact_bound=29,
        effective_compact_facet_ambient_bound=65,older_auxiliary_compact_facet_bound=79,
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
        scope='Large finite D short-leaf halves plus small/D10 coupons, A2/D11 full upgrade, order12 B64/C63 and order14 B84. Exact finite maximum495 at66, complete1659 positive ideal comparisons under the corrected concentration, positive global constant325, and completed66 higher tail. No example, external refereeing or formalization asserted.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['dimension_minima'],flush=True)
