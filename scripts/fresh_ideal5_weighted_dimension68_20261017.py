#!/usr/bin/env python3
"""Compose the complete order12 transfers and exact finite budget: ceiling68.

Local geometric classifications are explicit frozen inputs. The new scalar
domain is enumerated from scratch, with direct integer polynomial products.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_ideal5_weighted_dimension68_20261017'
BASE = 'fresh_ideal5_weighted_dimension69_20261017'
SUPPORT = 'fresh_order12_payment_support_20261017'
A = 'fresh_cycle12_ideal5_transfer_20261017'
B = 'fresh_order12_affine_b_half_transfer_20261017'
C = 'fresh_order12_affine_c_endpoint_reduction_20261017'
D = 'fresh_affine_d11_leaf_pair_transfer_20261017'
A17 = 'fresh_finite_a17_endpoint_payment_support_20261017'
FOREST = 'fresh_a17_endpoint_d11_upgrade_forest_budget_20261017'
FINITE = 'fresh_a11_cycle_coupon_domination_20261017'
TRI = 'fresh_ideal5_finite8_triangle_payments_20261017'
COMBINED = 'fresh_finite8_ideal5_combined_payments_20261017'
SHAPE = 'fresh_cusp_link_shape_curvature_20261017'
ENVELOPE = {'11': '65/12', '12': '75/13', '13': '43/7',
            '14': '32/5', '15': '107/16', '16': '120/17'}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def budget(rank):
    return F(rank, 2) if rank <= 10 else F(ENVELOPE[str(rank)]) if rank <= 16 else F(163, 23)


def terminal_shapes(total):
    """Sorted partitions, with bounds enforced while adjoining each part."""
    def visit(left, minimum, prefix, low_used, twelve_used, high_used):
        if not left:
            yield prefix
            return
        for r in range(minimum, left+1):
            low, twelve, high = r <= 9, r == 12, r >= 18
            if (low and low_used) or (twelve and twelve_used) or (high and high_used):
                continue
            yield from visit(left-r, r, prefix+(r,), low_used or low,
                             twelve_used or twelve, high_used or high)
    yield from visit(total, 1, (), False, False, False)


def profile(ranks):
    h = [1]
    for r in ranks:
        new = [0]*(len(h)+r)
        for j, value in enumerate(h):
            for k in range(r+1):
                new[j+k] += value
        h = new
    d = sum(ranks); v = prod(r+1 for r in ranks)
    assert len(h) == d+1 and h == h[::-1] and sum(h) == v
    g = sum(value*min(j, d-j) for j, value in enumerate(h))
    curvature = F(5*g-d*v, 2*v)
    cost = sum(map(budget, ranks), F())
    return dict(ranks=list(ranks), link_dimension=d, link_vertices=v,
                h_polynomial=h, G=g, normalized_curvature=str(curvature),
                normalized_cost=str(cost), margin=str(curvature-cost))


def distances(g):
    n = len(g)
    dist = [[0 if i == j else 1 if g[i][j] else n+1 for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                dist[i][j] = min(dist[i][j], dist[i][k]+dist[k][j])
    assert max(map(max, dist)) < n
    return dist


def local_interfaces(data):
    support = data[SUPPORT]
    assert len(support['high_patterns']) == 758
    assert data[COMBINED]['high_overlap_contexts'] == 2614
    assert data[COMBINED]['low_overlap_contexts'] == 66
    triangle = data[TRI]['rule']
    assert triangle == dict(finite_distance_at_most=8, ideal_distance_at_most=5,
        ideal_feature_distance_at_most=5, ideal_feature_pair_distance_at_most=8,
        ideal_whole_factor_order_at_most=9)
    all_pairs = set(combinations(range(12), 2))
    weights = {
        'A': {p: F(p in set(map(tuple, data[A]['pair_rule']['retained']))) for p in all_pairs},
        'B': {tuple(r['pair']): F(r['weight']) for r in data[B]['weighted_pairs']},
        'C': {p: F(p in set(map(tuple, data[C]['new_paid_pairs']))) for p in all_pairs},
        'D': {p: F(p in set(map(tuple, data[D]['retained_pairs']))) for p in all_pairs},
    }
    rows = []
    for rec in support['records']:
        kind = rec['kind']; g = rec['target']; dist = distances(g)
        features = [i for i in range(12) if sum(bool(x) for x in g[i]) != 2
                    or any(x > 3 for x in g[i])]
        elementary = {p for p in all_pairs if dist[p[0]][p[1]] <= 5 or
            (dist[p[0]][p[1]] <= 8 and features and
             min(dist[a][f] for a in p for f in features) <= 5)}
        assert all(weights[kind][p] == 1 for p in elementary)
        assert set(weights[kind]) == all_pairs
        assert all(0 <= w <= 1 for w in weights[kind].values())
        cost = sum(weights[kind].values(), F())
        expected = {'A':60, 'B':65, 'C':65, 'D':62}[kind]
        assert cost == expected
        changed = sorted(p for p in all_pairs if weights[kind][p] != 1)
        if kind in ('B', 'C'):
            assert all(weights[kind][p] == 1 for p in map(tuple, rec['short8_pairs']))
            assert all(weights[kind][p] == 1 for p in map(tuple, rec['basic_pairs']))
        if kind == 'A':
            assert len(changed) == 6 and all(dist[i][j] == 6 for i,j in changed)
            assert changed == sorted(map(tuple, data[A]['pair_rule']['transferred']))
        if kind == 'D':
            assert len(changed) == 4 and all(dist[i][j] == 9 for i,j in changed)
            assert changed == sorted(map(tuple, data[D]['removed_pairs']))
        rows.append(dict(kind=kind, order=12, cost=int(cost), changed_pairs=changed,
            changed_distances=[dist[i][j] for i,j in changed],
            elementary_triangle_pairs=sorted(elementary),
            weighted_pairs=[dict(pair=p,weight=str(weights[kind][p])) for p in sorted(all_pairs)]))
    assert data[B]['affected_maps'] == 28 and len(data[B]['affected_cases']) == 14
    assert data[B]['low_actual_variants'] == 69 and not data[B]['low_affected_cases']
    assert not data[B]['added_finite_coupons'] and not data[C]['added_finite_coupons']
    assert data[C]['complete_high_patterns'] == 758 and len(data[C]['high_case_embeddings']) == 15
    assert data[C]['low_interface']['affine_order12_families'] == {'B':2, 'D':10}
    assert data[A17]['candidate_removed_pair'] == [0,16]
    assert data[A17]['old_same_component_cost'] == 136
    assert data[A17]['candidate_same_component_cost'] == 135
    assert data[A17]['legacy_replacements']['remaining_total_at_affected_case'] == '2'
    assert data[FINITE]['added_coupons'] == data[A]['finite_coupons']
    assert data[FOREST]['variants'][1]['D11_upgrade']
    assert data[FOREST+'_review']['variants'][1]['rank69_cost'] == '544'
    maxima = data[FINITE]['dimension_maxima']
    assert maxima == data[FOREST]['variants'][1]['dimension_maxima']
    assert maxima[69]['doubled_cost'] == 1088
    return rows


def run():
    stems = (BASE, BASE+'_review', SUPPORT, A, B, C, D, A17, FOREST,
             FOREST+'_review', FINITE, TRI, COMBINED, SHAPE)
    data = {}
    for stem in stems:
        doc = json.loads((ROOT/f'results/{stem}.json').read_text())
        assert doc['status'].startswith('PASS'), stem
        for name, checksum in doc['sources'].items():
            assert digest(ROOT/name) == checksum, (stem,name)
        data[stem] = doc
    base = data[BASE]
    assert base['dimension_bound'] == data[BASE+'_review']['dimension_bound'] == 69
    assert base['unrestricted_finite_volume'] and not base['arithmetic_hypothesis']
    assert base['cusp_restriction'] is None
    interfaces = local_interfaces(data)
    costs = {k:dict(v) for k,v in base['small_order_costs'].items()}
    costs['12'] = {r['kind']:r['cost'] for r in interfaces}
    assert all(F(max(row.values()), int(order)) == budget(int(order)-1)
               for order,row in costs.items())
    shapes = list(terminal_shapes(68))
    assert len(shapes) == len(set(shapes)) == 2060
    rows = [profile(rs) for rs in shapes]
    assert all(F(r['margin']) > 0 for r in rows)
    minimum = min(rows, key=lambda r:F(r['margin']))
    assert minimum['margin'] == '3641/13608' and minimum['ranks'] == [11]*5+[13]
    gain = budget(11)+budget(13)-2*budget(12)
    assert gain == F(23,1092) > 0
    # Low concentration only through rank10: rank11 is no longer linear.
    assert all(budget(r) == F(r,2) for r in range(1,11))
    n = 69; beta = F((n-1)*(n-5),8)
    finite = F(data[FINITE]['dimension_maxima'][n]['doubled_cost'],2)
    assert finite == beta == 544
    boundary = profile((11,)*5+(12,))
    assert boundary['margin'] == '411499/3234816'
    lower_finite = F(data[FINITE]['dimension_maxima'][68]['doubled_cost'],2)
    assert lower_finite-F(68*62,8) == 17
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}' for stem in stems
             for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_UNRESTRICTED_WEIGHTED_IDEAL5_DIMENSION68', dimension_bound=68,
        inherited_dimension_bound=69, excluded_dimensions=[69], unrestricted_finite_volume=True,
        arithmetic_hypothesis=False, cusp_restriction=None, finite_budget_dimension_maximum=75,
        finite_vertex_budget_at_excluded_dimension='544', ideal_rank_envelope=ENVELOPE,
        small_order_costs=costs, order12_interfaces=interfaces,
        finite_rule_changes=dict(A17_endpoint_deleted=True, A1_D11_short_leaf_full_upgrade=True,
                                 new_A11_cycle_coupons=data[A]['finite_coupons']),
        finite_maximum_witness=data[FINITE]['dimension_maxima'][69]['witness_components'],
        reduced_profiles=rows,
        dimension_minima=[dict(ambient_dimension=69, profiles=len(rows),
            minimum_margin=minimum['margin'], minimizing_ranks=minimum['ranks'],
            finite_maximum=str(finite), finite_curvature_coefficient=str(beta),
            finite_slack='0', strict_global_slack=340)],
        concentration=dict(repeated_cap_ranks=[10,11,13,14,15,16,17],
            at_most_one_low_part=list(range(1,10)), at_most_one_each_middle_part=[12],
            at_most_one_high_part_minimum=18, double12_budget_increase=str(gain),
            cost_direction='Preserved or increased; normalized curvature decreases.',
            termination='Every nontrivial move increases the bounded integer sum of squared ranks.'),
        boundary_at_n68=dict(profile=boundary, finite_maximum=str(lower_finite),
            finite_curvature_coefficient='527', finite_slack='-17',
            meaning='Finite obstruction remains; no ceiling67 is asserted.'),
        unchanged_compact_arithmetic_bound=21, unchanged_arbitrary_compact_bound=29,
        effective_compact_facet_ambient_bound=68, older_auxiliary_compact_facet_bound=79,
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
        scope='Complete A/B/C/D order12 local composition and A17 deletion with the full revised finite coupon budget; all2060 reduced ideal profiles at69; positive global constant340; completed69 supplies n>=70. Internal computer-assisted proof, no new example, external refereeing or formalization claimed.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['dimension_minima'],flush=True)
