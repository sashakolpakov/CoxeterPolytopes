#!/usr/bin/env python3
"""Finite8 and the single D17 pair transfer: exact comparison for n<=75."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_ideal5_piecewise_dimension77_20261017 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ideal5_piecewise_dimension75_20261017'
BASE=B.TAG
LOCAL='fresh_finite8_ideal5_combined_payments_20261017'
PAIR='fresh_d17_singleton_pair_coupon_20261017'
FINITE='fresh_finite8_a4_d13_neighbor_coupon_budget_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
ENVELOPE={'1_to_12':'r/2','13':'87/14','14':'20/3','15':'111/16',
          '16':'126/17','17_and_larger':'163/23'}


def digest(p):return sha256(p.read_bytes()).hexdigest()


def profile(ranks):
    r=B.profile(ranks)
    decrease=F(ranks.count(16),17)
    r['normalized_cost']=str(F(r['normalized_cost'])-decrease)
    r['margin']=str(F(r['normalized_curvature'])-F(r['normalized_cost']))
    return r


def terminal_shapes(d):
    shapes=set()
    for a in range(d//12+1):
        for b in range((d-12*a)//14+1):
            for c in range((d-12*a-14*b)//16+1):
                for e in range((d-12*a-14*b-16*c)//17+1):
                    for p,q,low in product((0,13),(0,15),range(12)):
                        high=d-12*a-14*b-16*c-17*e-p-q-low
                        if high==0 or high>=18:
                            shapes.add(tuple(sorted([12]*a+[14]*b+[16]*c+[17]*e+
                                ([p]if p else[])+([q]if q else[])+([low]if low else[])+([high]if high else[]))))
    return sorted(shapes)


def scalars():
    rows=[];minima=[]
    for n in (76,77):
        group=[profile(rs)for rs in terminal_shapes(n-1)]
        assert all(F(r['margin'])>0 for r in group)
        smallest=min(group,key=lambda r:F(r['margin']))
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert beta>8*n
        minima.append(dict(ambient_dimension=n,profiles=len(group),minimum_margin=smallest['margin'],
            minimizing_ranks=smallest['ranks'],finite_slack=str(beta-8*n),strict_global_slack=5*(n-1)))
        rows+=group
    assert [r['profiles']for r in minima]==[1053,1111]
    assert [r['minimum_margin']for r in minima]==['11615/448188','41558/485537']
    assert [r['finite_slack']for r in minima]==['57','68']
    boundary=profile((12,12,12,12,12,14))
    assert F(boundary['margin'])==F(-16873,171366)
    double13=F(12,2)+F(20,3)-2*F(87,14)
    double15=F(20,3)+F(126,17)-2*F(111,16)
    assert double13==F(5,21)>0 and double15==F(83,408)>0
    return dict(reduced_profiles=rows,dimension_minima=minima,boundary_at_n75=boundary,
        concentration=dict(repeated_cap_ranks=[12,14,16,17],at_most_one_low_part=list(range(1,12)),
            at_most_one_each_middle_part=[13,15],at_most_one_high_part_minimum=18,
            double13_budget_increase=str(double13),double15_budget_increase=str(double15),
            cost_direction='preserved or increased; curvature decreases',
            termination='Every move strictly increases the bounded integer sum of squared parts at fixed total rank.'))


def run():
    out=scalars()
    stems=(BASE,BASE+'_review',LOCAL,PAIR,FINITE,SMALL)
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in stems}
    for stem,d in inputs.items():
        assert d['status'].startswith('PASS'),stem
        for name,h in d['sources'].items():assert digest(ROOT/name)==h,(stem,name)
    assert inputs[BASE]['dimension_bound']==inputs[BASE+'_review']['dimension_bound']==77
    assert inputs[LOCAL]['finite_vertex_budget']=='8n'
    assert inputs[LOCAL]['finite_budget_dimension_maximum']==80
    pair=inputs[PAIR];finite=inputs[FINITE]
    assert pair['diagram_order']==17 and pair['affine_type']=='D16'
    assert pair['old_cost']==127 and pair['new_cost']==126
    assert pair['removed_pair']==[1,11] and pair['affected_overlap_cases']==[2518]
    old=next(r for r in inputs[SMALL]['records']if r['order']==17 and r['kind']=='D')
    oldset=set(map(tuple,old['retained_paid_pairs']));newset=set(map(tuple,pair['new_paid_pairs']))
    assert len(oldset)==127 and len(newset)==126 and oldset-newset=={(1,11)}
    assert newset<=oldset and set(map(tuple,pair['triangle_pairs_retained']))<=newset
    for p in inputs[LOCAL]['short8_reduction']:
        if p['order']==17 and p['kind']=='D':assert not p['unpaid_short8_pairs']
    # Reconstruct the removed pair's distance rather than infer it from its name.
    target=old['target'];dist={1:0};queue=[1]
    for i in queue:
        for j,v in enumerate(target[i]):
            if v and j not in dist:dist[j]=dist[i]+1;queue.append(j)
    assert dist[11]==10
    core=pair['literal_core'];assert core['signature']==[17,1,1]
    assert core['omitted_radical_coefficients']=={'3':-2,'13':2}
    assert sum(F(c['half_coupon'])for c in core['corners'])==1
    paid_corners=[c for c in core['corners']if F(c['half_coupon'])]
    assert len(paid_corners)==2
    for corner in paid_corners:
        assert corner['signature']==[17,0,0] and F(corner['half_coupon'])==F(1,2)
        assert {(c['type'],c['role'])for c in corner['components']}=={
            ('A4','endpoint'),('D13','long_arm_neighbor_of_fork')}
    assert finite['dimension_maximum']==80 and finite['finite_vertex_budget']=='8n'
    assert finite['component_count']['A4_endpoints']==2 and finite['component_count']['D13_new_marks']==1
    assert finite['component_count']['half_weight']=='1/2'
    assert all(F(r['excess'])<=0 for r in finite['dimension_endpoints'])
    costs={m:dict(v)for m,v in inputs[LOCAL]['small_order_costs'].items()}
    assert costs['17']['D']==127;costs['17']['D']=126
    assert F(max(costs['17'].values()),17)==F(126,17)
    expected=dict(inputs[LOCAL]['ideal_rank_envelope']);expected['16']='126/17'
    assert expected==ENVELOPE
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}'for s in stems
        for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/f'results/{s}.md'for s in ('fresh_cusp_link_shape_curvature_20261017',
        'fresh_cusp_factor_payment_framework_20261017','fresh_distance14_unrestricted_dimension156_20261017')]
    out.update(status='PASS_IDEAL5_PIECEWISE_DIMENSION75',
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        dimension_bound=75,inherited_dimension_bound=77,excluded_dimensions=[76,77],
        unrestricted_finite_volume=True,arithmetic_hypothesis=False,cusp_restriction=None,
        finite_vertex_budget='8n',finite_budget_dimension_maximum=80,
        ideal_rank_envelope=ENVELOPE,small_order_costs=costs,
        local_change=dict(diagram_order=17,affine_type='D16',removed_pair=[1,11],
            previous_cost=127,new_cost=126,affected_overlap_case=2518,
            finite_coupon=pair['new_finite_coupon'],triangle_and_short8_pairs_retained=True),
        unchanged_compact_arithmetic_bound=21,unchanged_arbitrary_compact_bound=29,
        effective_compact_facet_ambient_bound=75,older_auxiliary_compact_facet_bound=79,
        scope='Complete finite8/ideal5 local theorem plus the one D17 singleton transfer, budgeted whole A4/D13-neighbor half coupons, all-rank concentration and shape inequality, and the completed77 higher tail. No arithmetic or cusp restriction; no new construction or external/formal validation claim.')
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--scalar-only',action='store_true');args=ap.parse_args()
    assert not(args.verify and args.scalar_only)
    if args.scalar_only:
        d=scalars();print('PASS_CONDITIONAL_DIMENSION75_SCALARS',d['dimension_minima'],flush=True)
    else:
        d=run();p=ROOT/f'results/{TAG}.json'
        if args.verify:assert d==json.loads(p.read_text())
        else:p.write_text(json.dumps(d,indent=2)+'\n')
        print(d['status'],flush=True)
