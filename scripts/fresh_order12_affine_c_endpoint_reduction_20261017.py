#!/usr/bin/env python3
"""Remove only the endpoint pair of affine C11, composing the frozen69 rule.

The inherited all-label star classifications are not rerun. All induced maps
of their758 eligible high patterns into the C path are reconstructed here.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_order12_affine_c_endpoint_reduction_20261017'
SUPPORT = 'fresh_order12_payment_support_20261017'
BASE = 'fresh_ideal5_weighted_dimension69_20261017'
COMBINED = 'fresh_finite8_ideal5_combined_payments_20261017'
TRI = 'fresh_ideal5_finite8_triangle_payments_20261017'
D13 = 'fresh_affine_d13_long_pair_transfer_20261017'
D14 = 'fresh_affine_d14_long_pair_transfer_20261017'
D12 = 'fresh_affine_d12_leaf_pair_transfer_20261017'
B13 = 'fresh_order13_affine_b_singleton_transfer_20261017'
B15 = 'fresh_order15_affine_b_pair_reduction_20261017'
SINGLE = 'fresh_d17_singleton_pair_coupon_20261017'
HALF17 = 'fresh_d17_partial_weight_transfer_20261017'
HALF16 = 'fresh_d16_partial_weight_transfer_20261017'
LATEST = 'fresh_b14_b15_singleton_completion_transfer_20261017'
CYCLE12 = 'fresh_cycle12_ideal5_transfer_20261017'
D11 = 'fresh_affine_d11_leaf_pair_transfer_20261017'
A17 = 'fresh_finite_a17_endpoint_payment_support_20261017'
CYCLES = ['fresh_finite8_cycle_replacement_20261017',
          'fresh_finite8_ideal5_local_replacement_20261017',
          'fresh_cycle14_ideal5_transfer_20261017',
          'fresh_cycle13_ideal5_transfer_20261017']


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(stem):
    doc = json.loads((ROOT / f'results/{stem}.json').read_text())
    assert doc['status'].startswith('PASS'), stem
    for name, value in doc['sources'].items():
        assert digest(ROOT / name) == value, (stem, name)
    return doc


def target():
    g = [[0] * 12 for _ in range(12)]
    for i in range(11):
        g[i][i+1] = g[i+1][i] = 4 if i in (0, 10) else 3
    return g


def path_maps(pattern, g):
    """Complete induced embeddings into a path, by contiguous segments."""
    m = len(pattern)
    assert 1 <= m <= len(g)
    adj = [[j for j in range(m) if pattern[i][j]] for i in range(m)]
    assert all(pattern[i][i] == 0 for i in range(m))
    assert all(pattern[i][j] == pattern[j][i] for i in range(m) for j in range(m))
    seen, todo = set(), [0]
    while todo:
        i = todo.pop()
        if i not in seen:
            seen.add(i)
            todo += adj[i]
    assert len(seen) == m, 'marked component must be connected'
    if any(len(a) > 2 for a in adj) or sum(map(len, adj)) != 2*(m-1):
        return []
    if m == 1:
        return [[i] for i in range(len(g))]
    ends = [i for i, a in enumerate(adj) if len(a) == 1]
    assert len(ends) == 2
    order = [min(ends)]
    while len(order) < m:
        order += [j for j in adj[order[-1]] if j not in order]
    maps = []
    for start in range(len(g)-m+1):
        for vertices in (list(range(start, start+m)), list(range(start, start+m))[::-1]):
            image = [None] * m
            for i, j in zip(order, vertices):
                image[i] = j
            if all(pattern[i][j] == g[image[i]][image[j]] for i in range(m) for j in range(m)):
                maps.append(image)
    return sorted(maps)


def finite_corner(c):
    return c['signature'][1:] == [0, 0]


def later_interfaces(inputs):
    """Recount literal terminal outputs; analytic short-path branches are inherited."""
    out = []
    for stem, key, unit in ((D13, 'extreme_equality', 'unit'),
                            (D14, 'equality', 'finite_small17_unit')):
        rows = inputs[stem][key]['rows']
        assert len(rows) == 72
        for row in rows:
            total = F(0)
            for c in row['corners']:
                value = F(c[unit]) + F(c['half_coupon'])
                assert not value or finite_corner(c)
                total += value
            assert total >= 1
        out.append(dict(source=stem, literal_equality_rows=len(rows),
                        positive_terminal_payments='finite only',
                        other_ideal_outputs='direct/shared-neighbor short paths or explicit whole B/D immediate factors'))
    d = inputs[D12]
    for key, expected in (('nonextreme_zero_cases', 48), ('arbitrary_one_sided_cases', 32)):
        assert len(d[key]) == expected
        for r in d[key]:
            assert finite_corner(r['chosen']) and r['chosen']['unit'] == 1
    assert len(d['two_A1_spectators']['rows']) == 16
    for r in d['two_A1_spectators']['rows']:
        assert all(finite_corner(c) for c in r['corners'])
        assert sum((F(c['unit']) + F(c['old_half']) for c in r['corners']), F()) >= 1
    out.append(dict(source=D12, checked_finite_configurations=96, other_ideal_output_maximum_distance=4))
    b = inputs[B13]
    cases = b['selected_support']['cases'] + b['conditional_long_orbit']['cases']
    for r in cases:
        for c in r['corners']:
            value = F(c['finite_unit']) + F(c['old_half']) + F(c['new_half'])
            assert not value or finite_corner(c)
    assert len(cases) == 6
    out.append(dict(source=B13, checked_literal_cases=6,
                    other_ideal_outputs='short paths; immediate D14 inner4 excluded by strict Schur gap'))
    b = inputs[B15]
    assert len(b['cases']) == 7 and not b['new_finite_coupons']
    for r in b['cases']:
        assert r['finite_payment'] == 1 and r['outcome']['maximum_order'] <= 14
        assert r['outcome']['method'] in ('no_direct_contact_to_marked_component', 'zero_residual_endpoint_extension')
    out.append(dict(source=B15, checked_finite_cases=7, ideal_terminal_outputs=[]))
    single = inputs[SINGLE]['literal_core']['corners']
    assert sum((F(c['finite_unit']) + F(c['half_coupon']) for c in single), F()) >= 1
    assert all(not (F(c['finite_unit']) + F(c['half_coupon'])) or finite_corner(c) for c in single)
    out.append(dict(source=SINGLE, checked_literal_corners=4, other_ideal_output_maximum_distance=2))
    for stem, order in ((HALF17, 17), (HALF16, 16)):
        d = inputs[stem]
        rows = d['cases'] if order == 17 else [c for r in d['cases'] for c in r['completions'] if c['outcome'] == 'sealed_completion']
        for r in rows:
            total = F(0)
            for c in r['corners']:
                iw = F(c['ideal_weight']); fw = F(c['finite_unit']) + F(c['finite_half'])
                assert not fw or finite_corner(c)
                if iw:
                    assert c['kind'] == f'affineD{order-1}_order{order}' and iw == F(1, 2)
                total += iw + fw
            assert total >= 1
        out.append(dict(source=stem, checked_sealed_cases=len(rows),
                        positive_ideal_outputs=f'only the same affine D factor of order{order}'))
    latest = inputs[LATEST]
    assert not latest['added_finite_coupons']
    assert [r['order'] for r in latest['records']] == [14, 15]
    for r in latest['records']:
        assert r['map_outcomes'] == {'projection_exclusion': r['affected_map_count']}
    out.append(dict(source=LATEST, excluded_maps=sum(r['affected_map_count'] for r in latest['records']),
                    ideal_terminal_outputs=[]))
    for stem in CYCLES:
        out.append(dict(source=stem,
                        analytic_interface='Cycle alternatives: finite units/coupons or ideal path at most2. Tree alternatives: ideal marked distance at most9 with a persistent nearby feature.',
                        removed_pair_distance=11))
    # These two newer independent local transfers are not needed for the69
    # composition, but their outputs are also safe for the present deletion.
    cyc = inputs[CYCLE12]
    assert cyc['pair_rule']['cost'] == 60 and len(cyc['literal_cases']) == 15
    for r in cyc['literal_cases']:
        for c in r['corners']:
            if c['coupons']:
                assert c['inertia']['negative'] == c['inertia']['zero'] == 0
    out.append(dict(source=CYCLE12, interface_only=True, checked_literal_cases=15,
                    other_ideal_output_maximum_distance=4,
                    finite_coupon_budget='Separate obligation; no coupon introduced by the C deletion.'))
    d11 = inputs[D11]
    assert len(d11['arbitrary_one_sided_corners']) == len(d11['two_A1_A3_spectators']['rows']) == 32
    assert all(r['selected_corner']['unit'] == 1 for r in d11['arbitrary_one_sided_corners'])
    for r in d11['two_A1_A3_spectators']['rows']:
        assert F(r['unit_total']) + F(r['half_total']) >= 1
        assert all(not c['type'].startswith('affine') for v in r['corners'] for c in v['components'])
    out.append(dict(source=D11, interface_only=True, checked_finite_configurations=64,
                    other_ideal_output_maximum_distance=5,
                    finite_coupon_budget='Separate obligation; no coupon introduced by the C deletion.'))
    a17 = inputs[A17]
    assert a17['legacy_replacements']['affected_case'] == 662
    assert F(a17['legacy_replacements']['remaining_total_at_affected_case']) == 2
    changed = next(r for r in a17['legacy_replacements']['cases'] if r['case'] == 662)
    positive = [c for c in changed['corners'] if F(c['retained_weight'])]
    assert len(positive) == 2 and all(c['shared_finite_component_type'] == 'D17' for c in positive)
    assert all(r['A17_endpoint_units'] == 0 for r in a17['successor_interfaces']['typed_unit_interfaces'])
    out.append(dict(source=A17, interface_only=True,
                    only_debited_legacy_case=662, remaining_finite_D17_units=2,
                    new_ideal_output_pairs=[]))
    return out


def run():
    stems = [SUPPORT, BASE, COMBINED, TRI, D13, D14, D12, B13, B15, SINGLE, HALF17, HALF16, LATEST,
             CYCLE12, D11, A17] + CYCLES
    inputs = {s: read(s) for s in stems}
    support = inputs[SUPPORT]; base = inputs[BASE]; combined = inputs[COMBINED]
    assert base['dimension_bound'] == 69 and base['finite_budget_dimension_maximum'] == 75
    assert combined['high_overlap_contexts'] == 2614 and combined['low_overlap_contexts'] == 66
    g = target(); row = next(r for r in support['records'] if r['kind'] == 'C')
    assert row['target'] == g
    automorphisms = path_maps(g, g)
    assert automorphisms == [list(range(12)), list(range(12))[::-1]]
    removed = (0, 11); paid = set(combinations(range(12), 2)) - {removed}
    assert len(paid) == 65
    assert {tuple(sorted(a[i] for i in removed)) for a in automorphisms} == {removed}
    # Scaling the norm2 C Gram by its positive kernel gives a path Laplacian.
    laplacian = [[0]*12 for _ in range(12)]
    for i in range(11):
        laplacian[i][i] += 1; laplacian[i+1][i+1] += 1
        laplacian[i][i+1] = laplacian[i+1][i] = -1
    assert all(sum(r) == 0 for r in laplacian)
    basic = {(a, b) for a, b in combinations(range(12), 2)
             if b-a <= 5 or (b-a <= 9 and min(abs(a-f) for f in (0,1,10,11)) <= 5)
             or (b-a <= 9 and min(abs(b-f) for f in (0,1,10,11)) <= 5)}
    short8 = {(a, b) for a, b in combinations(range(12), 2) if b-a <= 8}
    triangle_rule = inputs[TRI]['rule']
    assert triangle_rule == dict(finite_distance_at_most=8, ideal_distance_at_most=5,
        ideal_feature_distance_at_most=5, ideal_feature_pair_distance_at_most=8,
        ideal_whole_factor_order_at_most=9)
    assert len(basic) == 63 and basic == set(map(tuple, row['basic_pairs']))
    assert short8 <= basic <= paid
    patterns = support['high_patterns']; assert len(patterns) == 758
    saved_maps = {r['case']: r for r in row['high_case_embeddings']}
    high = []; images = set()
    for p in patterns:
        maps = path_maps(p['labels'], g)
        pairs = {tuple(sorted(a[i] for i in p['marks'])) for a in maps}
        assert pairs <= paid
        if maps:
            saved = saved_maps[p['index']]
            assert maps == sorted(saved['embeddings']) and pairs == set(map(tuple, saved['marked_images']))
            high.append(dict(case=p['index'],embeddings=maps,marked_images=sorted(pairs)))
            images |= pairs
        else:
            assert p['index'] not in saved_maps
    assert len(high) == 15 and sum(len(r['embeddings']) for r in high) == 30
    assert images == set(map(tuple, row['selected_high_pairs']))
    assert basic | images == paid
    low = support['low_corner_domain']['cases']; assert len(low) == 66
    variants = [v for r in low for v in r['variants']]; assert len(variants) == 69
    affine12 = [v for v in variants if v['signature'][2] and v['order'] == 12]
    assert Counter(v['affine_kind'] for v in affine12) == Counter(B=2,D=10)
    assert all('C' not in v['order12_target_maps'] for v in variants)
    immediate = support['immediate_D13_alternatives']
    assert len(immediate) == 6 and {r['target_kind'] for r in immediate} == {'B','D'}
    assert {r['order'] for r in inputs[D14]['profile_formulas']['small_affine_rows']} == {13}
    assert not row['selected_low_pairs'] and not row['immediate_pairs']
    interfaces = later_interfaces(inputs)
    paths = [Path(__file__)] + [ROOT/f'{folder}/{s}.{ext}' for s in stems
              for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_ORDER12_AFFINE_C_ENDPOINT_PAYMENT_REDUCTION',
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
        affine_type='C11',diagram_order=12,target=g,positive_kernel_scaling_squares=['1/2']+['1']*10+['1/2'],
        scaled_Gram_path_Laplacian=laplacian,automorphisms=automorphisms,
        old_cost=66,new_cost=65,removed_pairs=[removed],removed_pair_distance=11,
        new_paid_pairs=sorted(paid),basic_pairs=sorted(basic),short8_pairs=sorted(short8),
        elementary_triangle_pairs=sorted(short8),triangle_exceptional_maximum_order=9,
        complete_high_patterns=758,high_case_embeddings=high,selected_high_pairs=sorted(images),
        low_interface=dict(chosen_cases=66,actual_variants=69,affine_order12_families=dict(Counter(v['affine_kind'] for v in affine12)),C_images=[]),
        immediate_interfaces=dict(D13=immediate,D14_whole_factor_orders=[13]),
        later_transfer_interfaces=interfaces,added_finite_coupons=[],
        inherited_complete_rule=BASE,triangle_minimum=2,quadrilateral_minimum=1,
        finite_budget_unchanged=True,finite_budget_dimension_maximum=75,new_dimension_bound=None,
        proof_scope='Complete local composition of the frozen69 rule with only C order12 endpoint(0,11) removed. Inherited all-label core classifications, support closure and analytic transfer proofs remain inputs; no new global theorem or future unreviewed transfer composition is asserted.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['old_cost'],out['new_cost'],len(out['high_case_embeddings']),flush=True)
