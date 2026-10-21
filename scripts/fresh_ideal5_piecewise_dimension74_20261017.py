#!/usr/bin/env python3
"""Three order15 B pair deletions and the exact unrestricted dimension74 bound."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_ideal5_piecewise_dimension75_20261017 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ideal5_piecewise_dimension74_20261017'
LOCAL='fresh_order15_affine_b_pair_reduction_20261017'
SMALL='fresh_distance8_ideal5_small_factor_20261017'
PAYMENTS='fresh_finite8_ideal5_combined_payments_20261017'
ENVELOPE=dict(B.ENVELOPE,**{'14':'97/15'})


def digest(p):return sha256(p.read_bytes()).hexdigest()


def profile(ranks):
    r=B.profile(ranks)
    r['normalized_cost']=str(F(r['normalized_cost'])-F(ranks.count(14),5))
    r['margin']=str(F(r['normalized_curvature'])-F(r['normalized_cost']))
    return r


def scalars():
    rows=[profile(rs)for rs in B.terminal_shapes(74)]
    assert len(rows)==999 and all(F(r['margin'])>0 for r in rows)
    low=min(rows,key=lambda r:F(r['margin']))
    assert low['ranks']==[10,12,12,12,12,16]and F(low['margin'])==F(63,74698)
    beta=F(74*70,8);finite=beta-8*75
    assert finite==F(95,2)>0
    d13=F(6)+F(97,15)-2*F(87,14)
    d15=F(97,15)+F(126,17)-2*F(111,16)
    assert d13==F(4,105)>0 and d15==F(7,2040)>0
    boundary=profile((1,)+(12,)*6)
    assert F(boundary['margin'])==F(-168143,742586)
    return dict(reduced_profiles=rows,dimension_minima=[dict(ambient_dimension=75,profiles=999,
        minimum_margin=low['margin'],minimizing_ranks=low['ranks'],finite_slack=str(finite),strict_global_slack=370)],
        boundary_at_n74=boundary,concentration=dict(repeated_cap_ranks=[12,14,16,17],
            at_most_one_low_part=list(range(1,12)),at_most_one_each_middle_part=[13,15],
            at_most_one_high_part_minimum=18,double13_budget_increase=str(d13),double15_budget_increase=str(d15),
            cost_direction='preserved or increased; normalized curvature decreases',
            termination='Every move strictly increases the bounded integer sum of squared parts at fixed total rank.'))


def run():
    out=scalars();stems=(B.TAG,B.TAG+'_review',LOCAL,SMALL,PAYMENTS)
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in stems}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for name,h in d['sources'].items():assert digest(ROOT/name)==h,(s,name)
    old=docs[B.TAG];local=docs[LOCAL]
    assert old['dimension_bound']==docs[B.TAG+'_review']['dimension_bound']==75
    assert old['finite_vertex_budget']=='8n'and old['unrestricted_finite_volume']and not old['arithmetic_hypothesis']
    oldrow=next(r for r in docs[SMALL]['records']if r['order']==15 and r['kind']=='B')
    paid=set(map(tuple,local['new_paid_pairs']));before=set(map(tuple,oldrow['retained_paid_pairs']))
    drop={(0,11),(11,13),(11,14)}
    assert local['diagram_order']==15 and local['affine_type']=='B14'
    assert local['old_cost']==100 and local['new_cost']==len(paid)==97
    assert paid<before and before-paid==drop==set(map(tuple,local['removed_pairs']))
    target=oldrow['target'];distances=[]
    for a,b in sorted(drop):
        ds={a:0};todo=[a]
        for i in todo:
            for j,v in enumerate(target[i]):
                if i!=j and v and j not in ds:ds[j]=ds[i]+1;todo.append(j)
        assert len(ds)==15 and ds[b]>9;distances.append([a,b,ds[b]])
    affected=[r['case']for r in oldrow['case_embeddings']if set(map(tuple,r['marked_images']))&drop]
    assert affected==local['affected_cases']==[1284,1623,1641,1646,1666,1933,2326]
    assert [r['index']for r in local['cases']]==affected
    assert all(r['finite_payment']==1 and r['outcome']['maximum_order']<=14 for r in local['cases'])
    assert local['new_finite_coupons']==[]
    required=next(r['required_pairs']for r in docs[PAYMENTS]['triangle_pair_retention']if r['order']==15 and r['kind']=='B')
    assert set(map(tuple,required))<=paid
    assert next(r for r in docs[PAYMENTS]['short8_reduction']if r['order']==15 and r['kind']=='B')['unpaid_short8_pairs']==[]
    costs={m:dict(v)for m,v in old['small_order_costs'].items()};costs['15']['B']=97
    assert max(costs['15'].values())==97 and F(97,15)==F(ENVELOPE['14'])
    assert {k:v for k,v in old['ideal_rank_envelope'].items()if k!='14'}=={k:v for k,v in ENVELOPE.items()if k!='14'}
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths +=[ROOT/f'results/{s}.md'for s in('fresh_cusp_link_shape_curvature_20261017',
        'fresh_cusp_factor_payment_framework_20261017','fresh_distance14_unrestricted_dimension156_20261017')]
    out.update(status='PASS_UNRESTRICTED_IDEAL5_DIMENSION74',dimension_bound=74,inherited_dimension_bound=75,
        excluded_dimensions=[75],unrestricted_finite_volume=True,arithmetic_hypothesis=False,cusp_restriction=None,
        finite_vertex_budget='8n',finite_budget_dimension_maximum=80,ideal_rank_envelope=ENVELOPE,small_order_costs=costs,
        local_change=dict(affine_type='B14',diagram_order=15,old_cost=100,new_cost=97,removed_pairs=sorted(map(list,drop)),
            removed_pair_distances=distances,affected_overlap_cases=affected,additional_finite_cost='0'),
        unchanged_compact_arithmetic_bound=21,unchanged_arbitrary_compact_bound=29,effective_compact_facet_ambient_bound=74,
        older_auxiliary_compact_facet_bound=79,sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        scope='Complete local seven-case replacement plus the inherited finite8/ideal5 theorem and all999 reduced cusp profiles exclude dimension75. The completed dimension75 theorem supplies all higher dimensions. No new finite-volume example.')
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['reduced_profiles']),out['dimension_minima'],flush=True)
