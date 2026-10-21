#!/usr/bin/env python3
"""A local cycle interface for prospective deep finite distance-eight halves.

No complete changed payment system or dimension bound is asserted.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path

import fresh_finite8_cycle_replacement_20261017 as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_deep_finite8_cycle_interface_20261017'
A, B = R.A, R.B


def endpoint_geometry(depth):
    # Endpoint3 has marked distance L. Its opposite old endpoint is a
    # whole-component leaf, at distance other_arc-1 from the old mark.
    rows = []
    assert depth in (17, 20)
    for P, Q in [(6, 8), (7, 8)] + [(8, q) for q in range(8, depth+1)]:
        other = P if Q == 8 else Q
        if P == 8:
            other = Q
        assert other - 1 < depth
        rows.append(dict(P=P, Q=Q, attaching_arc=8,
                         fixed_leaf_distance_from_old_mark=other-1,
                         outcome='not_a_deep_pair'))
    assert len(rows) == depth-5
    return rows


def cycle_tail(depth):
    P, Q = 8, depth+1
    profiles = B.cycle_profiles()
    keep = []
    classification = []
    for p in profiles:
        old, distance = R.classification(p, P, Q)
        kind = old
        if p['name'] == 'end3_0':
            assert p['position'] == [1, 0] and old == 'finite8' and distance == 8
            assert p['c0'] == A.Q(1)
            kind = 'deep_pair_requires_spectator_neighbor'
        classification.append(dict(name=p['name'], classification=kind,
                                   marked_distance_at_Q_lower=distance,
                                   profile=B.wire(p)))
        if kind in ('unrestricted', 'deep_pair_requires_spectator_neighbor'):
            keep.append((p, kind))
    assert [p['name'] for p, _ in keep] == [
        'isolated', 'end3_0', 'end3_1', 'end4_1', 'next3_1', 'both_ends3']
    rows, pairs, all_values = [], [], []
    for (c, ck), (d, dk) in product(keep, repeat=2):
        special = 'deep_pair_requires_spectator_neighbor'
        side = 'c' if ck == special else 'd' if dk == special else None
        cc, dd = dict(c), dict(d)
        if side == 'c':
            cc['c0'] -= A.Q(F(1, 2))
        if side == 'd':
            dd['c0'] -= A.Q(F(1, 2))
        assert cc['c0'].sign() >= 0 and dd['c0'].sign() >= 0
        assert all(x.sign() >= 0 for x in cc['moment'] + dd['moment'])
        values = []
        for ub, vb in product((0, 1), repeat=2):
            u = c['umax'] if ub else A.Q()
            v = d['umax'] if vb else A.Q()
            raw = A.value(A.Q(F(1, P)), A.Q(F(1, Q)), c, d, u, v)
            correction = A.Q(F(1, P) + F(1, Q)) * (
                u*u/2 if side == 'c' else v*v/2 if side == 'd' else A.Q())
            value = A.value(A.Q(F(1, P)), A.Q(F(1, Q)), cc, dd, u, v)
            assert value == raw + correction and value.sign() > 0
            values.append(value)
            rows.append(dict(c=c['name'], d=d['name'], fixed_cauchy_side=side,
                             upper_height_reciprocals=[ub, vb],
                             raw=raw.wire(), correction=correction.wire(),
                             determinant_lower_bound=value.wire()))
        minimum = R.minimum(values)
        pairs.append(dict(c=c['name'], d=d['name'], fixed_cauchy_side=side,
                          corrected_norm_constants=[cc['c0'].wire(), dd['c0'].wire()],
                          minimum=minimum.wire()))
        all_values += values
    minimum = R.minimum(all_values)
    assert len(pairs) == 36 and len(rows) == 144
    assert minimum == A.Q({17: F(29, 144), 20: F(73, 336)}[depth])
    assert [(p['c'], p['d']) for p in pairs if p['minimum'] == minimum.wire()] == [
        ('end3_0', 'end3_0')]
    return dict(P=P, Q_lower=Q, all_Q_at_least=True,
                profile_classification=classification, pairs=pairs, rectangles=rows,
                minimum=minimum.wire(),
                monotone_tail=f'The six-profile list and height bounds are constant for Q>={Q}; nonnegative corrected norm constants and moments make the lower bound nonincreasing in 1/Q and concave in reciprocal heights.')


def run():
    paths = [Path(__file__)]
    for stem in ('fresh_finite8_cycle_replacement_20261017',
                 'fresh_finite8_ideal5_local_replacement_20261017',
                 'fresh_cycle14_ideal5_transfer_20261017',
                 'fresh_distance8_ideal_corner_replacement_20261017',
                 'fresh_simple_cusp_affine_a_quadrilateral_20261017'):
        paths += [ROOT / f'{folder}/{stem}.{ext}'
                  for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    for stem in ('fresh_cycle12_ideal5_transfer_20261017',
                 'fresh_cycle13_ideal5_transfer_20261017',
                 'fresh_finite_d16_d17_long_pair_half_transfer_20261017',
                 'fresh_large_finite_d_leaf_half_transfer_20261017',
                 'fresh_finite7_ideal5_triangle_payments_20261017',
                 'fresh_full_a43_corridor_high_domain_20261017'):
        paths.append(ROOT / f'results/{stem}.md')
    return dict(status='PASS_DEEP_FINITE8_LOCAL_CYCLE_INTERFACE',
                sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                         for p in paths},
                depth_interfaces=[dict(depth=depth,
                                       protected_endpoint_cases=endpoint_geometry(depth),
                                       tail=cycle_tail(depth)) for depth in (17, 20)],
                counts=dict(depth_thresholds=2, protected_endpoint_cases=27,
                            residual_profiles_per_depth=6,
                            ordered_pairs=72, rectangle_values=288),
                theorem='For either depth17 or20, the inherited ideal-five cycle/tree transfers still apply if whole finite A/B/D distance8 pairs with both marks at least that feature depth are reduced from units to halves, with all other current payments unchanged. The sole newly residual family is P=8,Q>=depth+1, with strict corrected lower bound29/144 at depth17 and73/336 at depth20.',
                scope='Local interface only. No complete modified quadrilateral coverage, finite-forest budget, global dimension bound, new polytope or independent arithmetic implementation.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], data['counts'], 'minima',
          [(x['depth'], x['tail']['minimum']) for x in data['depth_interfaces']], flush=True)
