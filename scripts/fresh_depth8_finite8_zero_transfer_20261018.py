#!/usr/bin/env python3
"""Exact certificates for the depth-eight finite-eight zero transfer.

The report supplies the all-label minimal-core and actual-corner arguments.
No forest optimization, geometric census or dimension conclusion occurs here.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_depth8_finite8_zero_transfer_20261018'
STEMS = (
    'fresh_depth11_finite8_zero_transfer_20261018',
    'fresh_universal_deep_finite8_corridor_20261018',
    'fresh_deep_finite8_cycle_interface_20261017',
    'fresh_distance8_ideal_corner_replacement_20261017',
    'fresh_finite8_cycle_replacement_20261017',
    'fresh_finite8_ideal5_local_replacement_20261017',
    'fresh_deep_finite8_dimension64_20261018',
    'fresh_finite7_ideal5_triangle_payments_20261017',
    'fresh_cycle11_half_pair_transfer_20261018',
    'fresh_order11_bcd_half_pair_transfer_20261018',
    'fresh_order11_order17_dimension63_20261018',
    'fresh_sparse_coupon_order11_13_dimension62_20261018',
    'fresh_order17_e6_unit_coupon_zero_transfer_20261018',
    'fresh_order14_ambient62_ideal_margin_closure_20261018',
    'fresh_rank15_rank16_dimension61_20261018',
)


def inverse(n, i, j):
    return F((min(i, j)+1)*(n-max(i, j)), n+1)


def path_check(n):
    det = [F(1), F(2)]
    for k in range(2, n+1):
        det.append(2*det[-1]-det[-2])
    assert det == list(map(F, range(1, n+2)))
    for i, j in product(range(n), repeat=2):
        value = 2*inverse(n, i, j)
        if i:
            value -= inverse(n, i-1, j)
        if i+1 < n:
            value -= inverse(n, i+1, j)
        assert value == int(i == j) and inverse(n, i, j) > 0
    return dict(order=n, inverse_product_entries=n*n,
                leading_determinants=list(map(str, det[1:])))


def finite_matrix(kind, n):
    # The B terminal is scaled by sqrt(2), an invertible congruence:
    # its diagonal is4 and incident entry-2. This is not a norm-two Gram.
    g = [[F(2 if i == j else 0) for j in range(n)] for i in range(n)]
    edges = [(i, i+1) for i in range(n-1)] if kind != 'D' else (
        [(i, i+1) for i in range(n-3)] + [(n-3, n-2), (n-3, n-1)])
    for i, j in edges:
        g[i][j] = g[j][i] = -1
    if kind == 'B':
        g[-1][-1] = 4
        g[-2][-1] = g[-1][-2] = -2
    h = [row[:] for row in g]
    pivots = []
    for i in range(n):
        p = h[i][i]
        assert p > 0
        pivots.append(p)
        for j in range(i+1, n):
            for k in range(i+1, n):
                h[j][k] -= h[j][i]*h[i][k]/p
    determinant = F(1)
    for p in pivots:
        determinant *= p
    assert determinant == (n+1 if kind == 'A' else 4)
    return dict(kind=kind, order=n, rational_congruence_gram=[list(map(str, r)) for r in g],
                positive_pivots=list(map(str, pivots)), determinant=str(determinant))


def retained_contact():
    gate = [dict(outward_distance=j, inverse_diagonal=str(inverse(8, j-1, j-1)))
            for j in range(3, 7)]
    assert [F(r['inverse_diagonal']) for r in gate] == [F(2), F(20, 9), F(20, 9), F(2)]
    assert all(F(r['inverse_diagonal']) >= 2 for r in gate)
    j7, j8 = inverse(8, 6, 6), inverse(8, 7, 7)
    both = j7+j8+2*inverse(8, 6, 7)
    assert j7 == F(14, 9) and j8 == F(8, 9)
    assert 2*j7 > 2 and F(5, 2)*j8 > 2 and both == 4 > 2
    # For labels m>=5, theta_m^2 >= (3+sqrt(5))/2 >5/2.
    assert F(5) > F(2)**2
    cases = [
        dict(only_contact=8, B_order=9, extra_white=False,
             contact_label=[3, 4], finite_outcomes=['A10', 'B10']),
        dict(only_contact=7, B_order=9, ninth_B_node='one inward step',
             contact_label=[3], finite_outcomes=['A10']),
        dict(only_contact=7, B_order=9, ninth_B_node='outward step8',
             contact_label=[3], finite_outcomes=['D10']),
        dict(only_contact=7, B_order=8, extra_white=False,
             contact_label=[3], finite_outcomes=['A9']),
        dict(only_contact=7, B_order=8, extra_white=True,
             extra_white_contacts_only='c', extra_bond_label=[3, 4],
             finite_outcomes=['A10', 'B10']),
    ]
    return dict(feature_depth=8, retained_outward_whites=6,
        minimal_indefinite_order_bound=10, B_order_at_most=9,
        edge_gates=gate, strict_PD_excludes_equality=True,
        contact7_inverse=str(j7), contact8_inverse=str(j8),
        contact7_nonsimple_projection_lower=str(2*j7),
        contact8_label5_projection_strict_lower=str(F(5, 2)*j8),
        simultaneous7_and8_projection_lower=str(both),
        disconnected_K_minus_c_extra_white_gate=str(F(5, 2)*j8),
        complete_missed_retained_contact_cases=cases,
        excluded_finite_matrices=[finite_matrix('A', 9), finite_matrix('A', 10),
                                  finite_matrix('B', 10), finite_matrix('D', 10)])


def corridor():
    n, a, b = 21, 6, 14
    left, right = list(range(4, 8)), list(range(13, 17))
    alpha = F(17, 352)
    blocks, cross = [], []
    for side, ids in [('left', left), ('right', right)]:
        for i, j in product(ids, repeat=2):
            u, v = (i+1, j+1) if side == 'left' else (n-i, n-j)
            c = inverse(n, i, j)-alpha*u*v
            assert c >= 0
            if i == j:
                assert c >= 2
            blocks.append(dict(side=side, i=i, j=j, factors=[u, v],
                inverse_entry=str(inverse(n, i, j)), residual_entry=str(c)))
    for i, j in product(left, right):
        u, v = i+1, n-j
        assert inverse(n, i, j) == F(u*v, n+1)
        cross.append(dict(i=i, j=j, factors=[u, v], entry=str(inverse(n, i, j))))
    coefficient = alpha**2-F(1, n+1)**2
    assert alpha > F(1, n+1) and coefficient == F(33, 123904)
    minimum_u = min(i+1 for i in left)
    minimum_v = min(n-j for j in right)
    determinant_minimum = coefficient*minimum_u**2*minimum_v**2
    assert determinant_minimum == F(20625, 123904) > 0
    return dict(path_order=n, marked_positions=[a, b], left_support=left, right_support=right,
        inverse_checks=[path_check(8), path_check(n)], block_factor=str(alpha),
        residual_blocks=blocks, cross_factorizations=cross,
        inequalities=['q-2 >= alpha*U^2', 'r-2 >= alpha*V^2', 's=U*V/22'],
        determinant_coefficient=str(coefficient), minimum_U=minimum_u, minimum_V=minimum_v,
        strict_determinant_lower_bound=str(determinant_minimum))


# (u_max*l1, u_max*l2, c0*u_max^2), all rational even for endpoint4.
PROFILES = {
    'isolated': (F(0), F(0), F(2)),
    'end3_0': (F(1), F(0), F(1)),
    'end3_1': (F(0), F(1), F(1)),
    'end4_0': (F(1), F(0), F(0)),
    'end4_1': (F(0), F(1), F(0)),
    'next3_0': (F(2), F(0), F(0)),
    'next3_1': (F(0), F(2), F(0)),
    'both_ends3': (F(1, 2), F(1, 2), F(0)),
}


def surd_multiply(x, y):
    a, b = map(F, x)
    c, d = map(F, y)
    return a*c+2*b*d, a*d+b*c


def profile_check(data):
    saved = data['fresh_distance8_ideal_corner_replacement_20261017']['cycle_profiles']
    assert [p['name'] for p in saved] == list(PROFILES)
    rows = []
    for p in saved:
        moments = [surd_multiply(p['umax'], m) for m in p['moment']]
        norm = surd_multiply(p['c0'], surd_multiply(p['umax'], p['umax']))
        assert all(v[1] == 0 for v in moments+[norm])
        actual = tuple(v[0] for v in moments+[norm])
        assert actual == PROFILES[p['name']]
        assert surd_multiply(p['umax'], p['minimum_height']) == (F(1), F(0))
        rows.append(dict(name=p['name'], position=p['position'],
                         normalized_moments=list(map(str, actual))))
    return rows


def classify(name, q):
    if name in ('isolated', 'both_ends3'):
        return 'unrestricted', None
    side = int(name[-1])
    distance = (8 if side == 0 else q)-(2 if name.startswith('next') else 1)+1
    if name == 'end3_0':
        assert distance == 8
        return 'changed_deep8_requires_neighbor', distance
    if distance <= 8:
        return 'retained_nondeep_finite_unit', distance
    if distance == 9:
        if name.startswith(('end4', 'next3')):
            return 'retained_sealed_leaf9_unit', distance
        return 'old_leaf9_requires_neighbor', distance
    return 'unrestricted', distance


def cycle_region(q):
    allowed = ('unrestricted', 'changed_deep8_requires_neighbor', 'old_leaf9_requires_neighbor')
    classes = [dict(name=p, classification=classify(p, q)[0],
                    distance=classify(p, q)[1]) for p in PROFILES]
    keep = [r for r in classes if r['classification'] in allowed]
    expected = {9: ['isolated', 'end3_0', 'end3_1', 'both_ends3'],
                10: ['isolated', 'end3_0', 'end3_1', 'end4_1', 'both_ends3'],
                11: ['isolated', 'end3_0', 'end3_1', 'end4_1', 'next3_1', 'both_ends3']}
    assert [r['name'] for r in keep] == expected[q]
    a, b = F(1, 8), F(1, q)
    rows, pairs = [], []
    for cr, dr in product(keep, repeat=2):
        c, d = cr['name'], dr['name']
        side = ('c' if cr['classification'].endswith('requires_neighbor') else
                'd' if dr['classification'].endswith('requires_neighbor') else None)
        constants = [PROFILES[c][2], PROFILES[d][2]]
        if side is not None:
            constants[0 if side == 'c' else 1] -= F(1, 2)
        assert min(constants) >= 0
        values = []
        for u, v in product((0, 1), repeat=2):
            l = u*PROFILES[c][0]+v*PROFILES[d][0]
            m = u*PROFILES[c][1]+v*PROFILES[d][1]
            norm = u*PROFILES[c][2]+v*PROFILES[d][2]
            raw = 1-2*(a*l+b*m)-(a+b)*norm-a*b*(l-m)**2
            correction = (a+b)*(F(u, 2) if side == 'c' else F(v, 2) if side == 'd' else F())
            lower = raw+correction
            assert lower > 0
            values.append(lower)
            rows.append(dict(c=c, d=d, upper_reciprocals=[u, v], fixed_cauchy_side=side,
                raw=str(raw), correction=str(correction), determinant_lower_bound=str(lower)))
        pairs.append(dict(c=c, d=d, fixed_cauchy_side=side,
            corrected_normalized_constants=list(map(str, constants)), minimum=str(min(values))))
    minimum = min(F(r['minimum']) for r in pairs)
    assert minimum == {9: F(1, 18), 10: F(1, 10), 11: F(1, 11)}[q]
    return dict(P=8, Q_lower=q, all_Q_at_least=q == 11, classifications=classes,
                retained_profiles=[r['name'] for r in keep], ordered_pairs=pairs,
                rectangle_values=rows, minimum=str(minimum))


def run():
    data = {s: json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem, packet in data.items():
        assert packet['status'].startswith('PASS'), stem
        for p, digest in packet['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == digest, (stem, p)
    current = data['fresh_rank15_rank16_dimension61_20261018']
    assert current['dimension_bound'] == 61 and not current['arithmeticity_required']
    old = data['fresh_depth11_finite8_zero_transfer_20261018']
    assert old['triangle_minimum'] == 2 and old['quadrilateral_minimum'] == 1
    protected = [dict(arcs=list(pq), fixed_leaf_distance=min(pq)-1) for pq in ((6, 8), (7, 8), (8, 8))]
    assert all(r['fixed_leaf_distance'] < 8 for r in protected)
    regions = [cycle_region(q) for q in (9, 10, 11)]
    assert sum(len(r['ordered_pairs']) for r in regions) == 77
    assert sum(len(r['rectangle_values']) for r in regions) == 308
    # The current changes to short ideal roles require the actual-corner
    # priorities proved in the report. The arithmetic below is not that proof.
    inherited = data['fresh_order11_order17_dimension63_20261018']['retained_small_tree_short8_interfaces']
    assert len(inherited) == 18
    paths = [Path(__file__)] + [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
             for folder, ext in [('scripts', 'py'), ('results', 'json'), ('results', 'md')]]
    return dict(status='PASS_DEPTH8_FINITE8_ZERO_LOCAL_TRANSFER',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        retained_contact=retained_contact(), corridor=corridor(),
        inherited_normalized_profiles=profile_check(data), protected_cycle_cases=protected,
        cycle_regions=regions, cycle_ordered_pairs=77, cycle_rectangle_values=308,
        current_baseline='fresh_rank15_rank16_dimension61_20261018',
        feature_depth=8, finite_distance=8, new_weight='0', previous_feature_depth=11,
        triangle_minimum=2, quadrilateral_minimum=1, ideal_weights_unchanged=True,
        all_other_finite_weights_and_coupons_unchanged=True,
        local_theorem='In the completed dimension61 simultaneous rule delete precisely all remaining whole finite A/B/D distance8 payments whose marks both have feature depth>=8. Every triangle still receives>=2 and every quadrilateral>=1.',
        scope='Written universal local theorem and exact matrix/profile inequalities. No new finite budget, ideal reduction, global dimension bound or finite-volume example. Published classification and inherited actual-face geometry remain explicit inputs.')


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
    print(result['status'], 'corridor', result['corridor']['strict_determinant_lower_bound'],
          'cycle_minima', [r['minimum'] for r in result['cycle_regions']], flush=True)
