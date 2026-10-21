#!/usr/bin/env python3
"""Compose the universal deep corridor, current short7 interfaces and rank65 budget.

The accompanying written proof supplies the geometric case split. This binder
checks its finite interfaces, exact scalar conclusion and frozen dependencies.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_deep_finite8_dimension64_20261018'
UNIVERSAL = 'fresh_universal_deep_finite8_corridor_20261018'
BUDGET = 'fresh_depth17_rank65_conditional_budget_20261018'
CYCLE = 'fresh_deep_finite8_cycle_interface_20261017'
TRIANGLE = 'fresh_finite7_ideal5_triangle_payments_20261017'
CURRENT = 'fresh_finite_d16_d17_long_pair_half_transfer_20261017'
IDEAL = 'fresh_ambient65_ideal_margin_closure_20261017'
TAIL = 'fresh_ideal5_weighted_dimension65_20261017'
BASE = 'fresh_ideal5_small_factor_payments_20261017'
SUP12 = 'fresh_order12_payment_support_20261017'
SUP13 = 'fresh_order13_payment_support_20261017'
BC12 = 'fresh_order12_bc_unit_coupon_transfer_20261017'
D12 = 'fresh_affine_d11_leaf_pair_transfer_20261017'
B13 = 'fresh_order13_affine_b_fork_projection_transfer_20261017'
D13 = 'fresh_affine_d12_leaf_pair_transfer_20261017'
B14 = 'fresh_order14_affine_b_fork_projection_transfer_20261017'
B15 = 'fresh_order15_affine_b_residual_terminal_transfer_20261017'
C15 = 'fresh_order15_affine_c_terminal_projection_transfer_20261017'
D15 = 'fresh_order15_affine_d_half_coupon_transfer_20261017'
D16 = 'fresh_d16_partial_weight_transfer_20261017'
D17 = 'fresh_order17_affine_d_three_quarter_transfer_20261017'
STEMS = (UNIVERSAL, BUDGET, CYCLE, TRIANGLE, CURRENT, IDEAL, TAIL,
         IDEAL+'_review', TAIL+'_review', BASE, SUP12, SUP13, BC12, D12,
         B13, D13, B14, B15, C15, D15, D16, D17,
         'fresh_finite8_ideal5_combined_payments_20261017',
         'fresh_finite8_ideal5_local_replacement_20261017',
         'fresh_cycle12_ideal5_transfer_20261017',
         'fresh_cycle13_ideal5_transfer_20261017',
         'fresh_cycle14_ideal5_transfer_20261017',
         'fresh_finite8_cycle_replacement_20261017',
         'fresh_ideal5_weighted_dimension71_20261017')


def classical(kind, n):
    g = [[0]*n for _ in range(n)]
    if kind == 'A':
        edges = [(i, (i+1) % n, 3) for i in range(n)]
    elif kind == 'B':
        edges = [(i, i+1, 4 if i == n-4 else 3) for i in range(n-3)]
        edges += [(0, n-2, 3), (0, n-1, 3)]
    elif kind == 'C':
        edges = [(i, i+1, 4 if i in (0, n-2) else 3) for i in range(n-1)]
    else:
        assert kind == 'D'
        edges = [(i, i+1, 3) for i in range(n-5)]
        edges += [(0, n-4, 3), (0, n-3, 3), (n-5, n-2, 3), (n-5, n-1, 3)]
    for i, j, label in edges:
        g[i][j] = g[j][i] = label
    distances = []
    for start in range(n):
        row = [None]*n
        row[start] = 0
        queue = deque([start])
        while queue:
            i = queue.popleft()
            for j, label in enumerate(g[i]):
                if label and row[j] is None:
                    row[j] = row[i]+1
                    queue.append(j)
        assert all(x is not None for x in row)
        distances.append(row)
    return g, distances


def unit_set(pairs):
    return {tuple(p): F(1) for p in pairs}


def short7_interfaces(data):
    entries = {}
    for row in data[BASE]['records']:
        entries[row['order'], row['kind']] = (row['target'], unit_set(row['paid_pairs']), BASE)
    for order, stem in ((12, SUP12), (13, SUP13)):
        for row in data[stem]['records']:
            g, dist = classical(row['kind'], order)
            assert g == row['target']
            # The completed cycle12/13 transfers leave exactly short5 units.
            pairs = [p for p in combinations(range(order), 2) if dist[p[0]][p[1]] <= 5]
            if row['kind'] != 'A':
                pairs = list(combinations(range(order), 2))
            entries[order, row['kind']] = (g, unit_set(pairs), stem)

    replacements = [
        (12, 'B', BC12, data[BC12]['B']['retained_unit_pairs']),
        (12, 'C', BC12, data[BC12]['C']['retained_unit_pairs']),
        (12, 'D', D12, data[D12]['retained_pairs']),
        (13, 'B', B13, data[B13]['retained_unit_pairs']),
        # The completed dimension71 theorem explicitly adopts this C12
        # basic-plus-selected set; SUP13 alone was only a diagnostic.
        (13, 'C', SUP13, next(r for r in data[SUP13]['records'] if r['kind'] == 'C')['required_union']),
        (13, 'D', D13, data[D13]['new_paid_pairs']),
        (14, 'B', B14, data[B14]['retained_pairs']),
        (15, 'B', B15, data[B15]['retained_unit_pairs']),
        (15, 'C', C15, data[C15]['retained_unit_pairs']),
    ]
    for order, kind, stem, pairs in replacements:
        g = entries[order, kind][0]
        if 'target' in data[stem]:
            assert data[stem]['target'] == g
        entries[order, kind] = (g, unit_set(pairs), stem)
    for order, stem in ((15, D15), (16, D16), (17, D17)):
        g, _, _ = entries[order, 'D']
        assert g == data[stem]['target']
        weights = {tuple(r['pair']): F(r['weight']) for r in data[stem]['weighted_pairs']}
        entries[order, 'D'] = (g, weights, stem)

    expected = {12: [60, 64, 63, 62], 13: [65, 73, 73, 74],
                14: [70, 84, 84, 82], 15: [75, 94, 91, 94],
                16: [80, 107, 103, 106], 17: [85, 111, 108, 119]}
    rows = []
    for (order, kind), (g, weights, stem) in sorted(entries.items()):
        graph, distances = classical(kind, order)
        assert g == graph
        assert sum(weights.values(), F()) == expected[order]['ABCD'.index(kind)]
        short = [p for p in combinations(range(order), 2) if distances[p[0]][p[1]] <= 7]
        missing = [p for p in short if weights.get(p, F()) != 1]
        if kind == 'A':
            assert missing == [p for p in short if distances[p[0]][p[1]] in (6, 7)]
        else:
            assert not missing, (order, kind, missing)
        rows.append(dict(order=order, kind=kind, target=graph,
                         final_weight_source=stem, total_weight=str(sum(weights.values(), F())),
                         weighted_pairs=[dict(pair=p, weight=str(w)) for p, w in sorted(weights.items())],
                         short7_pairs=short, unpaid_short7_pairs=missing,
                         transfer='retained cycle6/7 interface' if missing else 'literal units'))
    assert len(rows) == 24
    return rows


def run():
    data = {stem: json.loads((ROOT/f'results/{stem}.json').read_text()) for stem in STEMS}
    for stem, packet in data.items():
        assert packet['status'].startswith('PASS'), stem
        for name, digest in packet['sources'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest() == digest, (stem, name)
    corridor = data[UNIVERSAL]
    assert corridor['parameters']['feature_depth'] == 17
    assert F(corridor['inequalities']['strict_determinant_lower_bound']) == F(13674, 2645) > 0
    assert corridor['minimal_core_input']['order_bound'] == 10
    budget = data[BUDGET]
    assert budget['exact_modified_maximum'] == '479' and budget['curvature'] == '480'
    assert budget['geometric_payment_coverage_proved_here'] is False
    cycle = next(r for r in data[CYCLE]['depth_interfaces'] if r['depth'] == 17)
    assert cycle['tail']['Q_lower'] == 18 and cycle['tail']['minimum'] == ['29/144', '0']
    triangle = data[TRIANGLE]
    assert triangle['finite_distance_at_most'] == 7 and triangle['triangle_minimum'] == 2
    assert triangle['whole_ideal_factor_order_at_most'] == 9
    ideal = data[IDEAL]
    assert ideal['complete_terminal_profiles'] == 1650
    assert ideal['minimum_absolute_margin']['absolute_margin'] == '445/23'
    assert ideal['strengthened_global_additive_constant'] == '5280'
    assert data[TAIL]['dimension_bound'] == 65 and data[TAIL]['unrestricted_finite_volume'] is True
    interfaces = short7_interfaces(data)
    finite_slack = F(480)-F(479)
    ideal_slack = F(445, 23)
    constant = F(5280)
    assert finite_slack > 0 and ideal_slack > 0 and constant > 0
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
              for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    return dict(status='PASS_DEEP_FINITE8_UNRESTRICTED_DIMENSION64_COMPOSITION',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        dimension_bound=64, excluded_dimension=65, inherited_higher_tail_bound=65,
        unrestricted_finite_volume=True, arithmetic_hypothesis=False,
        finite_rule='Current finite rule with deep17 whole A/B/D distance8 units replaced by halves; all other weights and coupons unchanged.',
        ideal_weights_unchanged=True, triangle_minimum=2, quadrilateral_minimum=1,
        short7_small_factor_interfaces=interfaces,
        short7_tree_targets=18, short7_cycle_targets=6,
        quadrilateral_case_split=['No changed actual finite corner: all weights agree with the completed old rule.',
                                 'At least one changed actual finite corner: universal corridor supplies a full-star short7 path; finite and retained ideal units or the compatible cycle/tree transfer supply total>=1.'],
        finite_maximum='479', finite_curvature='480', finite_surplus='1',
        all_cusp_absolute_surplus='445/23', strict_global_constant='5280',
        impossible_inequality='0 >= v_f + (445/23)*v_i + 5280',
        unchanged_compact_arithmetic_bound=21, unchanged_arbitrary_compact_bound=29,
        scope='Complete written local-to-global composition with exact current short7 pair inclusion and scalar checks. Inherited geometric classification, ideal scalar census and finite budget retain their separate verification scopes. No new finite-volume example, external refereeing or formal verification.')


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
    print(result['status'], 'bound', result['dimension_bound'],
          'factor_interfaces', len(result['short7_small_factor_interfaces']))
