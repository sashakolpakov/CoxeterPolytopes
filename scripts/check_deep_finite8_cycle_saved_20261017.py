#!/usr/bin/env python3
"""Independent rational check of the two saved deep finite8 cycle tables.

No project imports. Checks the complete scalar domains and arithmetic, not
the inherited geometric classification or the complete payment system.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / 'results/fresh_deep_finite8_cycle_interface_20261017.json'

# At a maximum reciprocal height, list (u*l1, u*l2, c0*u^2).
# For end4_1, u=1/sqrt(2), l2=sqrt(2), c0=0, so all three
# entries are rational. Both-end3 has u=1/2 and moments (1,1).
PROFILES = {
    'isolated': (F(0), F(0), F(2)),
    'end3_0': (F(1), F(0), F(1)),
    'end3_1': (F(0), F(1), F(1)),
    'end4_1': (F(0), F(1), F(0)),
    'next3_1': (F(0), F(2), F(0)),
    'both_ends3': (F(1, 2), F(1, 2), F(0)),
}


def rational(pair):
    assert isinstance(pair, list) and len(pair) == 2
    assert F(pair[1]) == 0, ('expected rational scalar', pair)
    return F(pair[0])


def main():
    data = json.loads(PATH.read_text())
    assert data['status'] == 'PASS_DEEP_FINITE8_LOCAL_CYCLE_INTERFACE'
    assert [x['depth'] for x in data['depth_interfaces']] == [17, 20]
    totals = {'protected': 0, 'pairs': 0, 'values': 0}
    for interface in data['depth_interfaces']:
        depth = interface['depth']
        expected_endpoints = [(6, 8), (7, 8)] + [(8, q) for q in range(8, depth + 1)]
        endpoint_rows = interface['protected_endpoint_cases']
        assert [(r['P'], r['Q']) for r in endpoint_rows] == expected_endpoints
        for row in endpoint_rows:
            distance = (row['P'] if row['P'] < 8 else row['Q']) - 1
            assert row['attaching_arc'] == 8
            assert row['fixed_leaf_distance_from_old_mark'] == distance < depth
            assert row['outcome'] == 'not_a_deep_pair'
        totals['protected'] += len(endpoint_rows)

        tail = interface['tail']
        assert tail['P'] == 8 and tail['Q_lower'] == depth + 1
        assert tail['all_Q_at_least'] is True
        a, b = F(1, 8), F(1, depth + 1)
        expected_pairs = list(product(PROFILES, repeat=2))
        pairs = tail['pairs']
        assert [(r['c'], r['d']) for r in pairs] == expected_pairs
        rows = tail['rectangles']
        expected_keys = [(c, d, u, v) for c, d in expected_pairs for u, v in product((0, 1), repeat=2)]
        assert [(r['c'], r['d'], *r['upper_height_reciprocals']) for r in rows] == expected_keys
        values_by_pair = {pair: [] for pair in expected_pairs}
        for row in rows:
            c, d = row['c'], row['d']
            ub, vb = row['upper_height_reciprocals']
            side = 'c' if c == 'end3_0' else 'd' if d == 'end3_0' else None
            assert row['fixed_cauchy_side'] == side
            l = ub * PROFILES[c][0] + vb * PROFILES[d][0]
            m = ub * PROFILES[c][1] + vb * PROFILES[d][1]
            norm = ub * PROFILES[c][2] + vb * PROFILES[d][2]
            raw = 1 - 2 * (a * l + b * m) - (a + b) * norm - a * b * (l - m)**2
            correction = (a + b) * (F(ub, 2) if side == 'c' else F(vb, 2) if side == 'd' else 0)
            value = raw + correction
            assert rational(row['raw']) == raw
            assert rational(row['correction']) == correction
            assert rational(row['determinant_lower_bound']) == value > 0
            values_by_pair[c, d].append(value)
        for row in pairs:
            key = row['c'], row['d']
            assert rational(row['minimum']) == min(values_by_pair[key])
            side = 'c' if key[0] == 'end3_0' else 'd' if key[1] == 'end3_0' else None
            assert row['fixed_cauchy_side'] == side
            # These are unnormalized c0 constants, before multiplying by u^2.
            constants = [F(2) if name == 'isolated' else F(1) if name.startswith('end3') else F(0) for name in key]
            if side is not None:
                constants[0 if side == 'c' else 1] -= F(1, 2)
            assert [rational(v) for v in row['corrected_norm_constants']] == constants
            assert min(constants) >= 0
        minimum = min(min(v) for v in values_by_pair.values())
        assert minimum == F(5, 16) - F(2, depth + 1)
        assert minimum == {17: F(29, 144), 20: F(73, 336)}[depth]
        assert rational(tail['minimum']) == minimum
        assert [pair for pair, vals in values_by_pair.items() if min(vals) == minimum] == [('end3_0', 'end3_0')]
        totals['pairs'] += len(pairs)
        totals['values'] += len(rows)
    assert totals == {'protected': 27, 'pairs': 72, 'values': 288}
    assert data['counts'] == dict(depth_thresholds=2, protected_endpoint_cases=27,
                                 residual_profiles_per_depth=6, ordered_pairs=72,
                                 rectangle_values=288)
    print('PASS_INDEPENDENT_DEEP_FINITE8_RATIONAL_SCALARS', totals,
          'certificate_sha256', sha256(PATH.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
