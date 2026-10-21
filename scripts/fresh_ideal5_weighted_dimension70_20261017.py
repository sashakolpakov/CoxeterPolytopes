#!/usr/bin/env python3
"""Complete weighted order16 composition and exact unrestricted ceiling70."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_ideal5_weighted_dimension71_20261017 as B

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_ideal5_weighted_dimension70_20261017'
LOCAL = 'fresh_d16_partial_weight_transfer_20261017'
SMALL = 'fresh_distance8_ideal5_small_factor_20261017'
FINITE = 'fresh_finite8_d12_b10_coupon_budget_20261017'
ENVELOPE = dict(B.ENVELOPE, **{'15': '107/16'})


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def budget(rank):
    return F(rank, 2) if rank <= 11 else F(ENVELOPE[str(rank)]) if rank <= 16 else F(163, 23)


def profile(ranks):
    row = B.profile(ranks)
    cost = sum(map(budget, ranks), F())
    assert F(row['normalized_cost'])-cost == F(ranks.count(15), 4)
    row['normalized_cost'] = str(cost)
    row['margin'] = str(F(row['normalized_curvature'])-cost)
    return row


def scalar_comparison():
    # The previous complete terminal domain remains valid. It need not be minimal.
    rows = [profile(ranks) for ranks in B.terminal_shapes(70)]
    assert len(rows) == 1432 and all(F(r['margin']) > 0 for r in rows)
    lo = min(rows, key=lambda r: F(r['margin']))
    assert lo['margin'] == '23641/338688' and lo['ranks'] == [11, 11, 11, 11, 13, 13]
    beta = F(70*66, 8)
    assert beta-8*71 == F(19, 2)
    assert budget(11)+budget(13)-2*budget(12) == F(16, 91) > 0
    boundary = profile((11,)*5+(14,))
    assert boundary['link_dimension'] == 69 and F(boundary['margin']) == -F(1637, 51840)
    return dict(reduced_profiles=rows, dimension_minima=[dict(ambient_dimension=71,
        profiles=len(rows), minimum_margin=lo['margin'], minimizing_ranks=lo['ranks'],
        finite_slack=str(beta-8*71), strict_global_slack=350)], boundary_at_n70=boundary,
        concentration=dict(repeated_cap_ranks=[11, 13, 14, 15, 16, 17],
            at_most_one_low_part=list(range(1, 11)), at_most_one_each_middle_part=[12],
            at_most_one_high_part_minimum=18, double12_budget_increase='16/91',
            cost_direction='preserved or increased; normalized curvature decreases',
            termination='Each move strictly increases the bounded integer sum of squared parts at fixed total rank.',
            scope='The unchanged complete terminal domain is retained, without claiming it is minimal.'))


def local_interface(local, small):
    old = next(r for r in small['records'] if r['order'] == 16 and r['kind'] == 'D')
    old_paid = set(map(tuple, old['retained_paid_pairs']))
    half = set(map(tuple, local['half_weight_pairs']))
    weights = {tuple(r['pair']): F(r['weight']) for r in local['weighted_pairs']}
    assert local['affine_type'] == 'D15' and local['diagram_order'] == 16
    assert local['inherited_cost'] == len(old_paid) == 111 and len(half) == 10
    assert local['target'] == old['target']
    assert set(weights) == set(combinations(range(16), 2))
    assert weights == {p: F(1, 2) if p in half else F(int(p in old_paid)) for p in weights}
    assert sum(weights.values()) == F(local['new_cost']) == 106
    assert set(map(tuple, local['inherited_zero_pairs'])) == set(weights)-old_paid
    assert not local['added_finite_coupons']
    for key in ('triangle_pairs_retained_at_weight_one', 'short8_pairs_retained_at_weight_one'):
        assert all(weights[tuple(p)] == 1 for p in local[key])
    support = [r['case'] for r in old['case_embeddings'] if set(map(tuple, r['marked_images'])) & half]
    assert support == local['affected_cases'] == [r['index'] for r in local['cases']]
    assert len(support) == 54
    completions = [c for r in local['cases'] for c in r['completions']]
    assert len(completions) == 66
    assert sum(len(c['embeddings']) for c in completions) == local['affected_map_count'] == 328
    assert local['completion_outcomes'] == dict(excluded_lorentz_index=22,
        excluded_multiple_core_contact_roots=20, excluded_projection_above_norm=14,
        sealed_completion=10)
    for c in completions:
        assert all(weights[tuple(e['marked_image'])] == F(1, 2) for e in c['embeddings'])
        if c['outcome'] == 'sealed_completion':
            assert not c['new_roots']
            assert F(c['ideal_total']) == F(1, 2) and F(c['finite_payment_total']) >= F(1, 2)
            assert F(c['total_payment']) == F(c['ideal_total'])+F(c['finite_payment_total']) >= 1
        elif c['outcome'] == 'excluded_lorentz_index':
            assert c['signature'][1] > 1
        elif c['outcome'] == 'excluded_multiple_core_contact_roots':
            assert len(c['required_core_contacts']) > 1
        else:
            assert c['outcome'] == 'excluded_projection_above_norm' and F(c['projection_square']) > 2
    return dict(type='D15', order=16, old_cost=111, new_cost='106', half_weight_pairs=local['half_weight_pairs'],
        affected_cases=support, affected_maps=328, completion_outcomes=local['completion_outcomes'],
        new_finite_cost='0', total_payment_lower_bound='1',
        scope='All affected selected-factor embeddings are covered by the complete local theorem. Short8 and elementary triangle pairs retain unit weight; the inherited nine zero pairs and other local routes are unchanged.')


def run():
    stems = (B.TAG, B.TAG+'_review', LOCAL, SMALL, FINITE)
    inputs = {}
    for stem in stems:
        data = json.loads((ROOT/f'results/{stem}.json').read_text())
        assert data['status'].startswith('PASS'), stem
        for name, checksum in data['sources'].items():
            assert digest(ROOT/name) == checksum, (stem, name)
        inputs[stem] = data
    base = inputs[B.TAG]
    assert base['dimension_bound'] == inputs[B.TAG+'_review']['dimension_bound'] == 71
    assert base['unrestricted_finite_volume'] and not base['arithmetic_hypothesis']
    assert base['cusp_restriction'] is None
    finite = inputs[FINITE]
    assert finite['finite_vertex_budget'] == '8n' and finite['dimension_maximum'] == 75
    assert len(finite['dimension_endpoints']) == 13034
    interface = local_interface(inputs[LOCAL], inputs[SMALL])
    costs = {order: dict(row) for order, row in base['small_order_costs'].items()}
    assert costs['16'] == {'A': 80, 'B': 107, 'C': 103, 'D': 111}
    costs['16']['D'] = 106
    for order, row in costs.items():
        assert F(max(row.values()), int(order)) == budget(int(order)-1)
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}' for stem in stems
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    out = scalar_comparison()
    out.update(status='PASS_UNRESTRICTED_WEIGHTED_IDEAL5_DIMENSION70', dimension_bound=70,
        inherited_dimension_bound=71, excluded_dimensions=[71], unrestricted_finite_volume=True,
        arithmetic_hypothesis=False, cusp_restriction=None, finite_vertex_budget='8n',
        finite_budget_dimension_maximum=75, ideal_rank_envelope=ENVELOPE,
        small_order_costs=costs, weighted_order16=interface,
        unchanged_compact_arithmetic_bound=21, unchanged_arbitrary_compact_bound=29,
        effective_compact_facet_ambient_bound=70, older_auxiliary_compact_facet_bound=79,
        sources={str(path.relative_to(ROOT)): digest(path) for path in paths},
        scope='Complete weighted order16 local composition, unchanged finite8 budget, all1432 reduced cusp profiles at71 and the completed71 higher tail. No new finite-volume example, external refereeing or proof-assistant formalization claimed.')
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], len(result['reduced_profiles']), result['dimension_minima'], flush=True)
