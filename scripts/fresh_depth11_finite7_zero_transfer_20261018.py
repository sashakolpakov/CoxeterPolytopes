#!/usr/bin/env python3
"""Exact corridor and triangle interfaces for deep finite-distance-seven zeros.

The report supplies the universal geometry and pointwise rule composition.
No finite-forest optimization or new dimension theorem is computed here.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import fresh_finite7_ideal5_triangle_payments_20261017 as T

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_depth11_finite7_zero_transfer_20261018'
BASE = 'fresh_order11_order16_dimension60_20261018'
STEMS = (
    BASE, 'fresh_depth8_finite8_zero_transfer_20261018',
    'fresh_universal_deep_finite8_corridor_20261018', T.TAG,
    'fresh_finite8_ideal5_local_replacement_20261017',
    'fresh_cycle11_half_pair_transfer_20261018',
    'fresh_cycle12_ideal5_transfer_20261017',
    'fresh_cycle13_ideal5_transfer_20261017',
    'fresh_cycle14_ideal5_transfer_20261017',
    'fresh_order11_order17_dimension63_20261018',
    'fresh_sparse_coupon_order11_13_dimension62_20261018',
    'fresh_order14_ambient62_ideal_margin_closure_20261018',
    'fresh_order16_d_distance9_half_transfer_20261018',
    'fresh_order11_bc_zero_transfer_20261018',
    'fresh_cusp_factor_payment_framework_20261017',
)


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def inverse(n, i, j):
    return F((min(i, j) + 1) * (n - max(i, j)), n + 1)


def path_check(n):
    determinants = [F(1), F(2)]
    for k in range(2, n + 1):
        determinants.append(2 * determinants[-1] - determinants[-2])
    assert determinants == list(map(F, range(1, n + 2)))
    for i, j in product(range(n), repeat=2):
        v = 2 * inverse(n, i, j)
        if i:
            v -= inverse(n, i - 1, j)
        if i + 1 < n:
            v -= inverse(n, i + 1, j)
        assert v == int(i == j) and inverse(n, i, j) > 0
    return dict(order=n, inverse_entries=n*n,
                leading_determinants=list(map(str, determinants[1:])))


def corridor():
    n, a, b = 26, 9, 16
    assert b - a == 7 and a == n - 1 - b == 9
    arm = [dict(outward_distance=j, inverse_diagonal=str(inverse(11, j-1, j-1)))
           for j in range(3, 10)]
    assert min(F(r['inverse_diagonal']) for r in arm) == F(9, 4) > 2
    left, right = list(range(7, 11)), list(range(15, 19))
    alpha = F(122, 3267)
    assert alpha == min((inverse(n, i, i)-2)/F((i+1)**2) for i in left)
    entries, cross = [], []
    for side, ids in (('left', left), ('right', right)):
        for i, j in product(ids, repeat=2):
            u, v = (i+1, j+1) if side == 'left' else (n-i, n-j)
            residual = inverse(n, i, j)-alpha*u*v
            assert residual >= 0
            if i == j:
                assert residual >= 2
            entries.append(dict(side=side, i=i, j=j, factors=[u, v],
                inverse_entry=str(inverse(n, i, j)), residual_entry=str(residual)))
    for i, j in product(left, right):
        u, v = i+1, n-j
        assert inverse(n, i, j) == F(u*v, n+1)
        cross.append(dict(i=i, j=j, factors=[u, v], value=str(inverse(n, i, j))))
    coefficient = alpha**2-F(1, n+1)**2
    assert alpha > F(1, n+1) and coefficient == F(1, 43923)
    minimum = coefficient*8**4
    assert minimum == F(4096, 43923) > 0
    return dict(feature_depth=11, retained_outward_whites=9,
        minimal_indefinite_order_bound=10, retained_core_contact_radius=8,
        path_order=n, marked_positions=[a, b], left_support=left, right_support=right,
        strict_edge_arm_gates=arm, inverse_checks=[path_check(11), path_check(n)],
        block_factor=str(alpha), residual_blocks=entries, cross_factorizations=cross,
        determinant_coefficient=str(coefficient), minimum_U=8, minimum_V=8,
        strict_determinant_lower_bound=str(minimum),
        conclusion='A new deep finite7 corner forces a white-interior path of length<=6 between an adjacent black pair.')


def triangles(saved):
    fresh = T.run()
    assert fresh == saved
    protected = []
    for row in fresh['rows']:
        for c in row['corners']:
            if not c['elementary_payment'] or c['length'] < 7:
                continue
            assert c['length'] == 7
            w = c['feature_witness']
            assert w is not None and w['endpoint'] in c['pair']
            assert len(w['path'])-1 <= 5 < 11
            assert w['path'][0] == w['endpoint'] and w['path'][-1] == w['feature']
            assert (len(w['feature_neighbors']) >= 3 or
                    any(m >= 4 for _, m in w['incident_labels']))
            protected.append(dict(diagram=row['diagram'], marks=row['marks'],
                pair=c['pair'], feature_witness=w))
    assert len(protected) == 14
    exceptional = fresh['exceptional_resolution']['rows']
    assert [len(r['alternate_path'])-1 for r in exceptional] == [6, 6, 6, 7]
    e = exceptional[-1]
    assert e['marks'] == [1, 2, 9] and e['alternate_pair'] == [1, 9]
    assert e['fixed_simple_leaf'] == 1 and e['fixed_degree_two_neighbor'] == 0
    assert e['fixed_simple_leaf'] in e['alternate_pair']
    return dict(complete_high_core_marked_triples=len(fresh['rows']),
        retained_distance7_elementary_payments=protected,
        exceptional_alternate_lengths=[6, 6, 6, 7],
        exceptional_distance7_marked_leaf=e,
        low_core_reason='All-label order<=6 Steiner argument already gives two paths of length<=5.',
        triangle_minimum=2)


def run():
    docs = {s: json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s, d in docs.items():
        assert d['status'].startswith('PASS'), (s, d['status'])
        for p, expected in d['sources'].items():
            assert digest(ROOT/p) == expected, (s, p)
    assert docs[BASE]['dimension_bound'] == 60
    assert docs['fresh_depth8_finite8_zero_transfer_20261018']['feature_depth'] == 8
    small = docs['fresh_order11_order17_dimension63_20261018']['retained_small_tree_short8_interfaces']
    assert {(r['order'], r['kind']) for r in small} == set(product(range(12, 18), 'BCD'))
    assert all(r['short8_unit_pairs'] > 0 for r in small)
    tri = triangles(docs[T.TAG])
    paths = [Path(__file__), Path(T.__file__), Path(T.T.__file__), Path(T.T.T.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    paths = list(dict.fromkeys(paths))
    return dict(status='PASS_DEPTH11_FINITE7_ZERO_LOCAL_TRANSFER',
        sources={str(p.relative_to(ROOT)): digest(p) for p in paths},
        current_baseline=BASE, current_completed_dimension_bound=60,
        finite_distance=7, feature_depth=11, new_weight='0',
        corridor=corridor(), triangle_interface=tri,
        triangle_minimum=2, quadrilateral_minimum=1,
        ideal_weights_unchanged=True, all_other_finite_weights_and_coupons_unchanged=True,
        current_short6_ideal_interface='Cycle11 distance5 first; all remaining short5 units and small-tree short6 units retained. Unpaid distance6 cycles12/13/14 or>=15 and feature-free classical trees use their inherited universal proofs. Their finite7/8 outputs have a fixed nearby leaf or branch, while sealed components are too small for the new deletion.',
        whole_face_composition='If any newly zero finite7 actual corner occurs, the universal corridor gives an adjacent short6 actual corner. The current short6 ideal replacements retain a complete unit under the new rule. If none occurs, every actual corner weight agrees pointwise with the completed dimension60 rule.',
        new_dimension_bound=None, new_finite_volume_example=False,
        scope='Universal written local theorem, exact A11/A26 inverse and Schur inequalities, fresh1025 marked-triple reconstruction and protected finite7 outputs. Shared inherited triangle routines. No new core classification, finite-forest optimization, ideal convolution, ambient dimension bound or finite-volume example.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        assert not path.exists()
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['corridor']['strict_determinant_lower_bound'],
          result['triangle_interface']['complete_high_core_marked_triples'], flush=True)
