#!/usr/bin/env python3
"""Complete finite rank17 refinement excludes dimensions67 and68."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_ideal5_weighted_dimension68_20261017 as B

ROOT=B.ROOT;TAG='fresh_ideal5_weighted_dimension66_20261017'
LOCAL='fresh_finite_d17_long_pair_support_20261017'
FINITE='fresh_rank17_reduced_forest_frontier_20261017'


def digest(p):return sha256(p.read_bytes()).hexdigest()


def local_interface(local):
    assert local['new_payment_theorem'] and local['ideal_weights_unchanged']
    assert local['cross_component_coupons_unchanged']
    assert local['selected_domain_routing']==dict(whole_component_at_most16=2391,
        protected_leaf9=74,adaptive_leaf9_or_order16=9,eligible_bound17=140)
    out=[]
    for r in local['records']:
        paid=set(map(tuple,r['retained_pairs']));removed=set(map(tuple,r['removed_pairs']))
        assert paid.isdisjoint(removed) and paid|removed==set(combinations(range(17),2))
        assert len(paid)==r['new_cost']=={'A':114,'B':105,'D':127}[r['kind']]
        assert set(map(tuple,r['basic_pairs']))<=paid
        assert set(map(tuple,r['selected_high_pairs']))<=paid
        assert len(r['selected_high_contexts'])==140
        for c in r['selected_high_contexts']:
            assert all(tuple(m['marked_images'])in paid for m in c['maps'])
        assert all(tuple(sorted(a[i]for i in p))in paid for a in r['automorphisms']for p in paid)
        out.append(dict(kind=r['kind'],rank=17,old_cost=r['old_cost'],new_cost=r['new_cost'],
            basic_pairs=r['basic_pairs'],retained_pairs=r['retained_pairs'],removed_pairs=r['removed_pairs'],
            selected_contexts=sum(bool(c['maps'])for c in r['selected_high_contexts']),
            selected_maps=sum(len(c['maps'])for c in r['selected_high_contexts'])))
    totals=local['simultaneous_literal_totals']
    assert len(totals)==174 and all(F(r['retained_total'])>=1 for r in totals)
    debited=[r for r in totals if any(c['debit']for c in r['corners'])]
    assert len(debited)==1 and debited[0]['case']==662 and debited[0]['retained_total']=='2'
    return out


def run():
    stems=(B.TAG,B.TAG+'_review',LOCAL,FINITE,FINITE+'_review')
    data={}
    for stem in stems:
        d=json.loads((ROOT/f'results/{stem}.json').read_text())
        assert d['status'].startswith('PASS'),stem
        for name,h in d['sources'].items():assert digest(ROOT/name)==h,(stem,name)
        data[stem]=d
    base=data[B.TAG]
    assert base['dimension_bound']==data[B.TAG+'_review']['dimension_bound']==68
    assert base['unrestricted_finite_volume'] and not base['arithmetic_hypothesis'] and base['cusp_restriction']is None
    interface=local_interface(data[LOCAL])
    finite=data[FINITE]['variants'][1];review=data[FINITE+'_review']['variants'][1]
    assert(finite['assumed_A17_base'],finite['assumed_D17_base'])==(114,127)
    assert finite['assumed_B17_base_upper_bound']==127>=105
    assert finite['dimension_maxima']==review['dimension_maxima']
    groups=[];rows=[]
    for n,number,margin,shape in((67,1761,'851/165888',[11]*6),(68,1906,'411499/3234816',[11]*5+[12])):
        group=[B.profile(rs)for rs in B.terminal_shapes(n-1)]
        assert len(group)==number and all(F(r['margin'])>0 for r in group)
        m=min(group,key=lambda r:F(r['margin']))
        assert(m['margin'],m['ranks'])==(margin,shape)
        f=finite['dimension_maxima'][n];cost=F(f['doubled_cost'],2)
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert cost==8*n-32 and beta-cost>0
        assert f['witness_components']==[f'D{n}']
        rows+=group
        groups.append(dict(ambient_dimension=n,profiles=number,minimum_margin=margin,minimizing_ranks=shape,
            finite_maximum=str(cost),finite_curvature_coefficient=str(beta),finite_slack=str(beta-cost),
            finite_maximum_witness=f['witness_components'],strict_global_slack=5*(n-1)))
    boundary=B.profile((13,)*5)
    assert boundary['link_dimension']==65 and boundary['margin']=='-1565/10976'
    assert finite['dimension_maxima'][66]['doubled_cost']==992
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_UNRESTRICTED_WEIGHTED_IDEAL5_DIMENSION66',dimension_bound=66,
        inherited_dimension_bound=68,excluded_dimensions=[67,68],unrestricted_finite_volume=True,
        arithmetic_hypothesis=False,cusp_restriction=None,
        finite_rank17_interfaces=interface,finite_budget_dimension_maximum=75,
        finite_budget_at_excluded_dimensions={str(r['ambient_dimension']):r['finite_maximum']for r in groups},
        finite_budget_source_variant=1,ideal_weights_unchanged=True,cross_component_coupons_unchanged=True,
        ideal_rank_envelope=base['ideal_rank_envelope'],small_order_costs=base['small_order_costs'],
        concentration=base['concentration'],reduced_profiles=rows,dimension_minima=groups,
        boundary_at_n66=dict(profile=boundary,finite_maximum='496',finite_curvature_coefficient='495',
            finite_slack='-1',finite_maximum_witness=['D66'],new_bound65_asserted=False),
        unchanged_compact_arithmetic_bound=21,unchanged_arbitrary_compact_bound=29,
        effective_compact_facet_ambient_bound=66,older_auxiliary_compact_facet_bound=79,
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        scope='Complete simultaneous finite A17/B17/D17 deletions with all coupons/ideal weights unchanged. Independent exact finite forest certificate, all1761+1906 reduced profiles at67/68, positive finite/ideal/global margins, and completed68 higher tail. Internal proof; no example, outside refereeing or formalization claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['dimension_minima'],flush=True)
