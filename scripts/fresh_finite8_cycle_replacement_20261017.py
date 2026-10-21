#!/usr/bin/env python3
"""Finite8/leaf9 repair of the ideal6 cycle and feature-free tree transfers.

No high-core overlap census or ambient dimension conclusion is asserted.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path

import fresh_distance8_ideal_corner_replacement_20261017 as B
import fresh_cycle16_half_coupon_transfer_20261017 as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_finite8_cycle_replacement_20261017'
A = B.A


def minimum(values):
    out = values[0]
    for v in values[1:]:
        if (v-out).sign() < 0:
            out = v
    return out


def classification(p, P, Q):
    if p['name'] in ('isolated', 'both_ends3'):
        return 'unrestricted', None
    side = next(i for i, x in enumerate(p['position']) if x)
    dist = (P, Q)[side]-p['position'][side]+1
    if dist <= 8:
        return 'finite8', dist
    if dist == 9:
        if p['name'].startswith(('end4', 'next3')):
            return 'sealed_finite_leaf9', dist
        assert p['name'].startswith('end3')
        return 'requires_spectator_contact', dist
    return 'unrestricted', dist


def single_profiles():
    rows = []
    # Q=11 also represents the entire monotone tail Q>=11 below.
    for P in (7, 8, 9):
        for Q in range(max(P, 15-P), 12):
            m = P+Q
            for p in B.cycle_profiles():
                if p['name'] in ('isolated', 'both_ends3'):
                    continue
                side = next(i for i, x in enumerate(p['position']) if x)
                pos = p['position'][side]
                contact = pos if side == 0 else m-pos
                w2 = 2 if p['name'].startswith('end4') else 1
                ids = list(range(1, m))+[m]
                edges = [(i, i+1, 1) for i in range(1, m-1)]+[(m, contact, w2)]
                adj = {i: set() for i in ids}
                for i, j, _ in edges:
                    adj[i].add(j)
                    adj[j].add(i)
                leaf = m-1 if side == 0 else 1
                near = next(iter(adj[leaf]))
                assert len(adj[leaf]) == 1 and len(adj[near]) == 2
                kind, distance = classification(p, P, Q)
                assert C.distances(adj, ids, m)[P] == distance
                # Schur at c uses the A_(m-1) endpoint/next-node inverse.
                old_inverse = F(m-1, m) if pos == 1 else F(2*(m-2), m)
                c_inverse = 1/(2-w2*old_inverse)
                expected = F(m, 2) if w2 == 2 else F(m, 4) if pos == 2 else F(m, m+1)
                assert c_inverse == expected
                if kind == 'sealed_finite_leaf9':
                    assert c_inverse > 2
                rows.append(dict(P=P, Q=Q, profile=p['name'], classification=kind,
                    marked_distance=distance, positive_pivots=C.tree_pd(ids, edges),
                    fixed_leaf=leaf, fixed_degree_two_neighbor=near,
                    new_root_inverse_diagonal=str(c_inverse)))
    assert len(rows) == 66
    return rows


def rectangles():
    rows = []
    summaries = []
    profiles = B.cycle_profiles()
    counts = Counter()
    for P in (7, 8, 9):
        for Q in range(max(P, 15-P), 12):
            keep = [(p, classification(p, P, Q)[0]) for p in profiles]
            keep = [(p, k) for p, k in keep if k in ('unrestricted', 'requires_spectator_contact')]
            vals = []
            for (c, ck), (d, dk) in product(keep, repeat=2):
                # One fixed coordinate Cauchy bound only. Never add two such
                # lower bounds, and never switch sides at a rectangle vertex.
                side = 'c' if ck == 'requires_spectator_contact' else 'd' if dk == 'requires_spectator_contact' else None
                cc, dd = dict(c), dict(d)
                if side == 'c':
                    cc['c0'] -= A.Q(F(1, 2))
                if side == 'd':
                    dd['c0'] -= A.Q(F(1, 2))
                assert cc['c0'].sign() >= 0 and dd['c0'].sign() >= 0
                pair_vals = []
                for ub, vb in product((0, 1), repeat=2):
                    u, v = c['umax'] if ub else A.Q(), d['umax'] if vb else A.Q()
                    raw = A.value(A.Q(F(1, P)), A.Q(F(1, Q)), c, d, u, v)
                    correction = A.Q(F(1, P)+F(1, Q))*(u*u/2 if side == 'c' else v*v/2 if side == 'd' else A.Q())
                    value = A.value(A.Q(F(1, P)), A.Q(F(1, Q)), cc, dd, u, v)
                    assert value == raw+correction
                    pair_vals.append(value)
                    rows.append(dict(P=P, Q_lower=Q, all_Q_at_least=Q == 11,
                        c=c['name'], d=d['name'], fixed_cauchy_side=side,
                        upper_height_reciprocals=[ub, vb], raw=raw.wire(),
                        correction_lower_bound=correction.wire(), determinant_lower_bound=value.wire()))
                lo = minimum(pair_vals)
                special = c['name'] == d['name'] == 'isolated' and (P, Q) in ((7, 8), (7, 9), (8, 8))
                if special:
                    assert lo == A.Q(1-4*(F(1, P)+F(1, Q)))
                    counts['isolated_special_pairs'] += 1
                else:
                    assert lo.sign() > 0
                    vals += pair_vals
                    counts['strict_ordered_pairs'] += 1
            summaries.append(dict(P=P, Q_lower=Q, all_Q_at_least=Q == 11,
                retained_profiles=[dict(name=p['name'],kind=k) for p,k in keep],
                strict_pair_minimum=minimum(vals).wire()))
    assert len(rows) == 740
    assert counts == dict(isolated_special_pairs=3, strict_ordered_pairs=182)
    return dict(rows=rows,regions=summaries,counts=dict(counts))


def isolated_cases():
    rows = []
    for P, Q in ((7, 8), (7, 9), (8, 8)):
        alpha = F(1, P)+F(1, Q)
        gate = 1-3*alpha
        required = 4-1/alpha
        assert gate > 0 and 0 <= required < F(1, 2)
        rows.append(dict(P=P,Q=Q,one_height_at_least_sqrt2_gap=str(gate),
            required_unscaled_spectator_difference=str(required),
            outcome='sealed_cycle16_two_half_coupons' if required == 0 else 'excluded_by_disjoint_support_gap'))
    assert [r['required_unscaled_spectator_difference'] for r in rows] == ['4/15','1/16','0']
    return rows


def tree_checks():
    rows = []
    for P in (7, 8, 9):
        for c, d in product(B.tree_profiles(), repeat=2):
            if P <= 8 and any(p['name'] == 'inner_next3' for p in (c,d)):
                rows.append(dict(P=P,c=c['name'],d=d['name'],kind='immediate_payment',
                    length=P-1,finite_payment='distance<=8',
                    ideal_payment='distance<=6' if P == 7 else 'distance7 with a mark within1 of a persistent fork'))
                continue
            costs = [2*p['moment']/p['minimum_height']+p['c0']/p['minimum_height']**2 for p in (c,d)]
            value = A.Q(1)-sum(costs,A.Q())/P
            lower = {7:F(1,7),8:F(1,4),9:F(1,9)}[P]
            assert (value-A.Q(lower)).sign() >= 0
            rows.append(dict(P=P,c=c['name'],d=d['name'],kind='positive_gap',gap=value.wire()))
    assert Counter(r['kind'] for r in rows) == dict(immediate_payment=30,positive_gap=162)
    return rows


def run():
    data = dict(single_profile_checks=single_profiles(),cycle=rectangles(),
        isolated_special_cases=isolated_cases(),tree_checks=tree_checks())
    paths = [Path(__file__)]
    for stem in ('fresh_distance8_ideal_corner_replacement_20261017',
                 'fresh_ideal_distance6_mixed_replacement_20261017',
                 'fresh_cycle16_half_coupon_transfer_20261017'):
        paths += [ROOT/f'{folder}/{stem}.{ext}' for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    data.update(status='PASS_FINITE8_LEAF9_IDEAL6_LOCAL_REPLACEMENT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        counts=dict(cycle_regions=11,cycle_rectangle_values=740,single_profile_checks=66,
            strict_ordered_cycle_pairs=182,isolated_special_pairs=3,tree_cases=192),
        theorem='Actual ideal cycle corners of order>=15 and marked distance7,8,9, and classical tree corners at those distances outside feature radius6, transfer to finite distance<=8, finite distance9 at a whole-component leaf, or the ideal6/feature6 rule; the cycle16 antipodal equality transfers to two sealed A1/A15-middle half-coupons.',
        hypotheses=dict(actual_quadrilateral=True,arbitrary_supporting_whites=True,
            simple_cusps_required=False,arithmetic_required=False),
        scope='Local transfer only. The full distance8 overlap census, all preferred small-factor obligations, total finite budgets and any ambient dimension improvement remain separate.')
    return data


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['counts'],flush=True)
