#!/usr/bin/env python3
"""Universal depth11 short7 corridor and compatible zero-weight face rule.

Exact rational identities only; geometric classification and composition are
written in the companion report. No new dimension conclusion is asserted.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_depth11_finite8_zero_transfer_20261018'
STEMS = ('fresh_universal_deep_finite8_corridor_20261018',
         'fresh_deep_finite8_cycle_interface_20261017',
         'fresh_deep_finite8_dimension64_20261018',
         'fresh_finite7_ideal5_triangle_payments_20261017',
         'fresh_finite8_ideal5_local_replacement_20261017',
         'fresh_finite8_cycle_replacement_20261017')
DEPTH, OUTWARD, N, LEFT, RIGHT = 11, 9, 27, 9, 17


def inverse(n, i, j):
    return F((min(i, j)+1)*(n-max(i, j)), n+1)


def corridor():
    assert 10-2 <= OUTWARD and RIGHT-LEFT == 8
    checks = []
    for n in (DEPTH, N):
        determinants = [F(1), F(2)]
        for k in range(2, n+1):
            determinants.append(2*determinants[-1]-determinants[-2])
        assert determinants == [F(i+1) for i in range(n+1)]
        for i in range(n):
            for j in range(n):
                v = 2*inverse(n, i, j)
                if i:
                    v -= inverse(n, i-1, j)
                if i+1 < n:
                    v -= inverse(n, i+1, j)
                assert v == int(i == j) and inverse(n, i, j) > 0
        checks.append(dict(order=n, inverse_product_entries=n*n,
                           leading_principal_determinants=list(map(str, determinants[1:]))))
    gates = [dict(distance=j, inverse_diagonal=str(inverse(DEPTH, j-1, j-1)))
             for j in range(3, OUTWARD+1)]
    assert min(F(r['inverse_diagonal']) for r in gates) == F(9, 4) > 2
    left, right = list(range(7, 11)), list(range(16, 20))
    assert left == list(range(LEFT-2, LEFT+2)) and right == list(range(RIGHT-1, RIGHT+3))
    alpha = F(17, 308)
    blocks = []
    for side, positions in [('left', left), ('right', right)]:
        for i, j in product(positions, repeat=2):
            u, v = (i+1, j+1) if side == 'left' else (N-i, N-j)
            entry = inverse(N, i, j)
            assert entry >= alpha*u*v > 0
            blocks.append(dict(side=side, i=i, j=j, entry=str(entry),
                               lower_bound=str(alpha*u*v)))
    cross = []
    for i, j in product(left, right):
        v = F((i+1)*(N-j), N+1)
        assert inverse(N, i, j) == v
        cross.append(dict(i=i, j=j, entry=str(v), factors=[i+1, N-j]))
    q0 = min(inverse(N, i, i) for i in left+right)
    rho2 = 1/(F((N+1)**2)*alpha**2)
    coefficients = [((1-rho2)*q0*q0-4*q0+4), (1-rho2)*q0-2, 1-rho2]
    assert q0 == F(40, 7) and rho2 == F(121, 289)
    assert coefficients == [F(36, 289), F(382, 289), F(168, 289)]
    assert min(coefficients) > 0
    return dict(feature_depth=DEPTH, retained_outward_whites=OUTWARD,
        path_order=N, marked_positions=[LEFT, RIGHT],
        left_contact_support=left, right_contact_support=right,
        minimal_core_order_bound=10, retained_component_radius_bound=8,
        inverse_checks=checks, edge_gates=gates, block_bounds=blocks,
        cross_factorizations=cross, block_factor=str(alpha), q_minimum=str(q0),
        correlation_squared=str(rho2), shifted_determinant_coefficients=list(map(str, coefficients)),
        strict_determinant_minimum=str(coefficients[0]))


# Normalized endpoint moments (u*l1,u*l2,c0*u^2) at maximal reciprocal
# height. For endpoint4 the square root cancels exactly. These are the
# six complete tail profiles from the frozen cycle-interface theorem.
PROFILES = {
    'isolated': (F(0), F(0), F(2)),
    'end3_0': (F(1), F(0), F(1)),
    'end3_1': (F(0), F(1), F(1)),
    'end4_1': (F(0), F(1), F(0)),
    'next3_1': (F(0), F(2), F(0)),
    'both_ends3': (F(1, 2), F(1, 2), F(0)),
}


def cycles():
    protected = [dict(arcs=[p, q], fixed_leaf_distance=(p if p<8 else q)-1)
                 for p, q in [(6, 8), (7, 8)] + [(8, q) for q in range(8, DEPTH+1)]]
    assert len(protected) == 6 and all(r['fixed_leaf_distance'] < DEPTH for r in protected)
    a, b = F(1, 8), F(1, DEPTH+1)
    rows, pairs = [], []
    for c, d in product(PROFILES, repeat=2):
        side = 'c' if c == 'end3_0' else 'd' if d == 'end3_0' else None
        constants = [F(2) if p == 'isolated' else F(1) if p.startswith('end3') else F(0)
                     for p in (c, d)]
        if side is not None:
            constants[0 if side == 'c' else 1] -= F(1, 2)
        assert min(constants) >= 0
        values = []
        for u, v in product((0, 1), repeat=2):
            l = u*PROFILES[c][0] + v*PROFILES[d][0]
            m = u*PROFILES[c][1] + v*PROFILES[d][1]
            norm = u*PROFILES[c][2] + v*PROFILES[d][2]
            raw = 1-2*(a*l+b*m)-(a+b)*norm-a*b*(l-m)**2
            correction = (a+b)*(F(u, 2) if side == 'c' else F(v, 2) if side == 'd' else F())
            value = raw+correction
            assert value > 0
            values.append(value)
            rows.append(dict(c=c, d=d, upper_reciprocals=[u, v], fixed_cauchy_side=side,
                raw=str(raw), correction=str(correction), determinant_lower_bound=str(value)))
        pairs.append(dict(c=c, d=d, corrected_norm_constants=list(map(str, constants)), minimum=str(min(values))))
    minimum = min(F(r['determinant_lower_bound']) for r in rows)
    assert len(pairs) == 36 and len(rows) == 144
    assert minimum == F(5, 16)-F(2, 12) == F(7, 48)
    assert [(r['c'], r['d']) for r in pairs if F(r['minimum']) == minimum] == [('end3_0', 'end3_0')]
    return dict(protected_endpoints=protected, P=8, Q_at_least=12,
                normalized_profiles={k:list(map(str, v)) for k,v in PROFILES.items()},
                pairs=pairs, rectangle_values=rows, minimum=str(minimum),
                argument='Same six-profile list for all Q>=12; fixed-side Cauchy correction, nonnegative corrected constants, reciprocal-height concavity and monotonicity in1/Q.')


def run():
    packets = {s: json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem, data in packets.items():
        assert data['status'].startswith('PASS'), stem
        for path, expected in data['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == expected, (stem, path)
    base = packets['fresh_deep_finite8_dimension64_20261018']
    assert base['short7_tree_targets'] == 18 and base['short7_cycle_targets'] == 6
    assert base['ideal_weights_unchanged']
    c, cyc = corridor(), cycles()
    paths = [Path(__file__), ROOT/'scripts/check_deep_finite8_cycle_saved_20261017.py']
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_DEPTH11_FINITE8_ZERO_LOCAL_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        corridor=c, cycle_interface=cyc, changed_weight='0', previous_weights=['1','1/2'],
        triangle_minimum=2, quadrilateral_minimum=1, ideal_weights_unchanged=True,
        local_theorem='Delete precisely the finite distance8 pairs in whole A/B/D components when both marks have feature depth>=11. Retain every other current finite and ideal weight and coupon. Triangles still receive>=2 and quadrilaterals>=1.',
        scope='Complete written universal local transfer, with exact corridor and cycle endpoint inequalities. No new finite budget, ideal reduction, global dimension bound or finite-volume example.')


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
    print(result['status'], 'corridor', result['corridor']['strict_determinant_minimum'],
          'cycle', result['cycle_interface']['minimum'])
