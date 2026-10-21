#!/usr/bin/env python3
"""Order13 transfers, weighted D17 payments and the exact dimension71 bound."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_ideal5_piecewise_dimension74_20261017 as B
import fresh_order13_affine_b_singleton_transfer_20261017 as L

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_ideal5_weighted_dimension71_20261017'
DIAG = 'fresh_order13_payment_support_20261017'
D12 = 'fresh_affine_d12_leaf_pair_transfer_20261017'
CYCLE = 'fresh_cycle13_ideal5_transfer_20261017'
WEIGHTED = 'fresh_d17_partial_weight_transfer_20261017'
FINITE = 'fresh_finite8_d12_b10_coupon_budget_20261017'
PAYMENTS = 'fresh_finite8_ideal5_combined_payments_20261017'
ENVELOPE = {'1_to_11': 'r/2', '12': '75/13', '13': '87/14', '14': '97/15',
            '15': '111/16', '16': '120/17', '17_and_larger': '163/23'}


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def budget(r):
    return F(r, 2) if r <= 11 else F(ENVELOPE[str(r)]) if r <= 16 else F(163, 23)


def profile(ranks):
    row = B.profile(ranks)
    cost = sum(map(budget, ranks), F())
    assert F(row['normalized_cost'])-cost == F(3*ranks.count(12), 13)+F(6*ranks.count(16), 17)
    row['normalized_cost'] = str(cost)
    row['margin'] = str(F(row['normalized_curvature'])-cost)
    return row


def terminal_shapes(d):
    shapes = set()
    caps = (11, 13, 14, 15, 16, 17)
    def visit(i, rest, parts):
        if i == len(caps):
            for middle in (0, 12):
                for low in range(11):
                    high = rest-middle-low
                    if high == 0 or high >= 18:
                        shapes.add(tuple(sorted(parts+([middle] if middle else [])+
                            ([low] if low else [])+([high] if high else []))))
            return
        r = caps[i]
        for count in range(rest//r+1):
            visit(i+1, rest-count*r, parts+[r]*count)
    visit(0, d, [])
    assert all(sum(rs) == d for rs in shapes)
    return sorted(shapes)


def scalars():
    rows, minima = [], []
    for n in (72, 73, 74):
        group = [profile(rs) for rs in terminal_shapes(n-1)]
        assert all(F(r['margin']) > 0 for r in group)
        lo = min(group, key=lambda r: F(r['margin']))
        beta = F(n*(n-6), 8) if n % 2 == 0 else F((n-1)*(n-5), 8)
        assert beta-8*n > 0
        minima.append(dict(ambient_dimension=n, profiles=len(group), minimum_margin=lo['margin'],
            minimizing_ranks=lo['ranks'], finite_slack=str(beta-8*n), strict_global_slack=5*(n-1)))
        rows += group
    assert [r['profiles'] for r in minima] == [1530, 1635, 1748]
    assert [r['minimum_margin'] for r in minima] == ['6875/49152', '41401/145152', '135455/304128']
    assert [r['finite_slack'] for r in minima] == ['18', '28', '37']
    increment = budget(11)+budget(13)-2*budget(12)
    assert increment == F(16, 91) > 0
    # The old repeated-rank15 elimination is deliberately not used.
    assert budget(14)+budget(16)-2*budget(15) < 0
    boundary = profile((11,)*5+(15,))
    assert F(boundary['margin']) == -F(5597, 248832)
    return dict(reduced_profiles=rows, dimension_minima=minima, boundary_at_n71=boundary,
        concentration=dict(repeated_cap_ranks=[11, 13, 14, 15, 16, 17],
            at_most_one_low_part=list(range(1, 11)), at_most_one_each_middle_part=[12],
            at_most_one_high_part_minimum=18, double12_budget_increase=str(increment),
            cost_direction='preserved or increased; normalized curvature decreases',
            termination='Every move strictly increases the bounded integer sum of squared parts at fixed total rank.'))


def order13(inputs):
    diag, local, cycle, d12 = (inputs[x] for x in (DIAG, L.TAG, CYCLE, D12))
    assert diag['eligible_patterns'] == 801
    assert local['old_cost'] == 78 and local['new_cost'] == 77
    bdrop = {(0, 10), (10, 11), (10, 12)}
    long = local['conditional_long_orbit']
    assert set(map(tuple, long['pairs'])) == bdrop-{(0, 10)} and not long['removed']
    assert [r['index'] for r in long['cases']] == [166, 170, 1744, 2137]
    promoted = []
    for case in long['cases']:
        if case['finite_unit_total']:
            assert case['index'] in (166, 1744)
            promoted.append(dict(case=case['index'], payment='whole finite B13 unit'))
        else:
            assert case['index'] in (170, 2137)
            corner = next(r for r in case['corners'] if r['old_half'] == '1/2')
            selected = [corner['components'][i] for i in corner['marked_component_indices']]
            assert {r['type'] for r in selected} == {'A1', 'D12'}
            d = next(r for r in selected if r['type'] == 'D12')
            mark = next(m for m in corner['marks'] if m in d['nodes'])
            assert d['roles'][str(mark)] == 'short_fork_leaf'
            assert corner['signature'] == [13, 0, 0]
            promoted.append(dict(case=case['index'], payment='full whole A1/D12 short-leaf coupon'))
    assert all(F(r['normalized_determinant_lower_bound']) > 0
               for r in local['immediate_D14_strict_gap']['rows'])
    records = []
    for rec in diag['records']:
        kind, target = rec['kind'], rec['target']
        before = set(combinations(range(13), 2))
        required = set(map(tuple, rec['required_union']))
        paid = set(map(tuple, cycle['cycle_pair_rule']['retained'])) if kind == 'A' else (
            before-bdrop if kind == 'B' else required if kind == 'C' else set(map(tuple, d12['new_paid_pairs'])))
        if kind == 'A':
            assert paid == required and before-paid == set(map(tuple, cycle['cycle_pair_rule']['transferred']))
        if kind == 'B':
            assert paid == set(map(tuple, local['new_paid_pairs']))-set(map(tuple, long['pairs']))
        if kind == 'D':
            assert before-paid == set(map(tuple, d12['removed_pairs']))
        _, ds, features = L.P.distance_data(target)
        triangles = {p for p in before if ds[p[0]][p[1]] <= 5 or
            (ds[p[0]][p[1]] <= 8 and any(ds[i][f] <= 5 for i in p for f in features))}
        assert triangles <= paid
        missing_short = {p for p in before-paid if ds[p[0]][p[1]] <= 8}
        assert not missing_short if kind != 'A' else missing_short == before-paid
        affected = [r['case'] for r in rec['case_embeddings'] if set(map(tuple, r['marked_images']))-paid]
        if kind == 'B':
            assert affected == [160, 166, 170, 1206, 1744, 2137]
        elif kind == 'D':
            assert len(affected) == 127
        else:
            assert not affected
        immediate_missing = set(map(tuple, rec['immediate_pairs']))-paid
        assert not immediate_missing if kind in 'AC' else immediate_missing <= before-paid
        records.append(dict(kind=kind, order=13, cost=len(paid),
            retained_pairs=[list(p) for p in sorted(paid)], removed_pairs=[list(p) for p in sorted(before-paid)],
            elementary_triangle_pairs=[list(p) for p in sorted(triangles)],
            short8_pairs_needing_universal_cycle=[list(p) for p in sorted(missing_short)],
            selected_affected_cases=affected, immediate_affected_pairs=[list(p) for p in sorted(immediate_missing)]))
    assert [r['cost'] for r in records] == [65, 75, 73, 74]
    return dict(records=records, promoted_B_distance11_cases=promoted,
        scope='Complete order13 selected-pattern domain, universal cycle13 short-path replacement, B inner4 strict gaps and universal D12 leaf transfer; every other inherited immediate factor has order<=12 or an unchanged finite/short-path payment.')


def run():
    stems = (B.TAG, B.TAG+'_review', DIAG, L.TAG, D12, CYCLE, WEIGHTED, FINITE, PAYMENTS)
    inputs = {}
    for stem in stems:
        d = json.loads((ROOT / f'results/{stem}.json').read_text())
        assert d['status'].startswith('PASS'), stem
        for name, h in d['sources'].items():
            assert digest(ROOT / name) == h, (stem, name)
        inputs[stem] = d
    assert inputs[B.TAG]['dimension_bound'] == inputs[B.TAG+'_review']['dimension_bound'] == 74
    assert inputs[B.TAG]['unrestricted_finite_volume'] and not inputs[B.TAG]['arithmetic_hypothesis']
    finite = inputs[FINITE]
    assert finite['finite_vertex_budget'] == '8n' and finite['dimension_maximum'] == 75
    assert len(finite['dimension_endpoints']) == 13034 and all(F(r[-1]) <= 0 for r in finite['dimension_endpoints'])
    weighted = inputs[WEIGHTED]
    assert weighted['inherited_cost'] == 126 and F(weighted['new_cost']) == 120
    assert not weighted['added_finite_coupons'] and len(weighted['cases']) == 12
    assert all(F(r['payment_total']) >= 1 for r in weighted['cases'])
    assert sum(F(r['weight']) for r in weighted['weighted_pairs']) == 120
    weights = {tuple(r['pair']): F(r['weight']) for r in weighted['weighted_pairs']}
    for key in ('triangle_pairs_retained_at_weight_one', 'short8_pairs_retained_at_weight_one'):
        assert all(weights[tuple(p)] == 1 for p in weighted[key])
    interface = order13(inputs)
    costs = {m: dict(v) for m, v in inputs[B.TAG]['small_order_costs'].items()}
    costs['13'] = {r['kind']: r['cost'] for r in interface['records']}
    costs['17']['D'] = 120
    for order, row in costs.items():
        assert F(max(row.values()), int(order)) == budget(int(order)-1)
    out = scalars()
    paths = [Path(__file__)] + [ROOT / f'{folder}/{stem}.{ext}' for stem in stems
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths += [ROOT / f'results/{stem}.md' for stem in ('fresh_cusp_link_shape_curvature_20261017',
        'fresh_cusp_factor_payment_framework_20261017', 'fresh_distance14_unrestricted_dimension156_20261017')]
    out.update(status='PASS_UNRESTRICTED_WEIGHTED_IDEAL5_DIMENSION71', dimension_bound=71,
        inherited_dimension_bound=74, excluded_dimensions=[72, 73, 74],
        unrestricted_finite_volume=True, arithmetic_hypothesis=False, cusp_restriction=None,
        finite_vertex_budget='8n', finite_budget_dimension_maximum=75,
        ideal_rank_envelope=ENVELOPE, small_order_costs=costs, order13_interfaces=interface,
        weighted_order17=dict(type='D16', cost='120', half_weight_pairs=weighted['half_weight_pairs'],
            affected_cases=weighted['affected_cases'], new_finite_cost='0'),
        unchanged_compact_arithmetic_bound=21, unchanged_arbitrary_compact_bound=29,
        effective_compact_facet_ambient_bound=71, older_auxiliary_compact_facet_bound=79,
        sources={str(p.relative_to(ROOT)): digest(p) for p in paths},
        scope='Complete order13 and weighted order17 local composition, the enlarged finite8 budget through75, all4913 exact reduced cusp profiles at72..74, and the completed74 higher tail. No new finite-volume example, external review or proof-assistant formalization claimed.')
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    d = json.loads(json.dumps(run()))
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert d == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], len(d['reduced_profiles']), d['dimension_minima'], flush=True)
