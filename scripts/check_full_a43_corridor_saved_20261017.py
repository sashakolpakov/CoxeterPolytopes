#!/usr/bin/env python3
"""Independent saved-template arithmetic, without importing project modules."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_full_a43_corridor_high_domain_20261017'
ZERO = (F(), F())


def decode(value):
    return tuple(map(F, value))


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[0], -x[1]


def multiply(x, y):
    coefficients = [F(), F(), F()]
    for i, a in enumerate(x):
        for j, b in enumerate(y):
            coefficients[i+j] += a*b
    return coefficients[0]+2*coefficients[2], coefficients[1]


def sign(x):
    """Use successively tighter rational enclosures, not a squared-norm test."""
    a, b = x
    if not b:
        return (a > 0)-(a < 0)
    bits = 8
    while True:
        denominator = 1 << bits
        floor = isqrt(2*denominator*denominator)
        low, high = F(floor, denominator), F(floor+1, denominator)
        lo, hi = sorted((a+b*low, a+b*high))
        if lo > 0:
            return 1
        if hi < 0:
            return -1
        bits *= 2


def inverse(matrix):
    """Generic rational Gauss-Jordan inversion, independent of the path formula."""
    n = len(matrix)
    rows = [[F(x) for x in row]+[F(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        divisor = rows[j][j]
        rows[j] = [x/divisor for x in rows[j]]
        for i in range(n):
            if i == j or not rows[i][j]:
                continue
            multiplier = rows[i][j]
            rows[i] = [x-multiplier*y for x, y in zip(rows[i], rows[j])]
    assert all(rows[i][j] == int(i == j) for i in range(n) for j in range(n))
    return [row[n:] for row in rows]


def projection(x, y, matrix):
    total = ZERO
    for i, a in x.items():
        for j, b in y.items():
            value = multiply(a, b)
            total = add(total, (matrix[i][j]*value[0], matrix[i][j]*value[1]))
    return total


def main():
    packet = json.loads((ROOT/f'results/{TAG}.json').read_text())
    for name, digest in packet['sources'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    assert packet['status'] == 'PASS_FULL_A43_CORRIDOR_ALTERNATIVE_ON_COMPLETE_6703_HIGH_DOMAIN'
    cartan = [[2 if i == j else -1 if abs(i-j) == 1 else 0 for j in range(43)] for i in range(43)]
    full_inverse = inverse(cartan)
    white = [i for i in range(43) if i not in (17, 25)]
    white_inverse_small = inverse([[cartan[i][j] for j in white] for i in white])
    white_inverse = [[F() for _ in range(43)] for _ in range(43)]
    for i, a in enumerate(white):
        for j, b in enumerate(white):
            white_inverse[a][b] = white_inverse_small[i][j]
    counts, outcome_counts = Counter(), {}
    for name, template in packet['templates'].items():
        sides = []
        for constraints, mark, slots in zip(template['fixed_slot_constraints'], (17, 25), ((16, 15, 18), (26, 27, 24))):
            fixed = [None if x is None else decode(x) for x in constraints]
            domains = []
            for k, value in enumerate(fixed[1:]):
                allowed = [ZERO, (F(1), F())] + ([] if k == 1 else [(F(), F(1))])
                assert value is None or value in allowed
                domains.append(allowed if value is None else [value])
            side = []
            for values in product(*domains):
                if values[0] != ZERO and values[1] != ZERO:
                    continue
                vector = {mark: fixed[0], **dict(zip(slots, values))}
                side.append({i: x for i, x in vector.items() if x != ZERO})
            sides.append(side)
        for side, saved in zip(sides, template['side_profiles']):
            assert [[[i, list(map(str, value))] for i, value in sorted(vector.items())] for vector in side] == saved
        expected_choices = list(product(range(len(sides[0])), range(len(sides[1]))))
        assert [tuple(r['side_choices']) for r in template['joint_profiles']] == expected_choices
        outcomes = Counter()
        for record in template['joint_profiles']:
            i, j = record['side_choices']
            x, y = sides[0][i], sides[1][j]
            slacks = [add((F(2), F()), neg(projection(v, v, white_inverse))) for v in (x, y)]
            assert [list(map(str, s)) for s in slacks] == record['strict_edge_slacks']
            if any(sign(s) <= 0 for s in slacks):
                reason = 'strict_actual_edge_failure'
            else:
                a, b = [add((F(2), F()), neg(projection(v, v, full_inverse))) for v in (x, y)]
                cross = neg(projection(x, y, full_inverse))
                determinant = add(multiply(a, b), neg(multiply(cross, cross)))
                assert [list(map(str, z)) for z in (a, b)] == record['schur_diagonals']
                assert list(map(str, cross)) == record['schur_cross']
                assert list(map(str, determinant)) == record['determinant']
                assert sign(a) < 0 and sign(determinant) > 0
                reason = 'two_negative_directions'
            assert record['reason'] == reason
            outcomes[reason] += 1
            counts['joint_profiles'] += 1
        assert dict(outcomes) == template['profile_outcomes']
        outcome_counts[name] = outcomes
        counts['templates'] += 1
    actual_uses, weighted, statuses = Counter(), Counter(), Counter()
    for case in packet['cases']:
        statuses[case['status']] += 1
        for orientation in case['orientations']:
            counts['orientations'] += 1
            if orientation['status'] == 'complete_contact_template':
                name = orientation['template']
                actual_uses[name] += 1
                weighted.update(outcome_counts[name])
            else:
                assert orientation['status'] == 'fixed_contact_exclusion' and orientation['witnesses']
                weighted['no_contact_profile'] += 1
    assert dict(actual_uses) == packet['template_multiplicities']
    assert dict(weighted) == packet['weighted_profile_outcomes']
    assert dict(statuses) == packet['corner_outcomes']
    assert counts == dict(templates=379, joint_profiles=4983, orientations=13644)
    assert len(packet['cases']) == 26812
    print(json.dumps(dict(status='PASS_INDEPENDENT_SAVED_A43_TEMPLATE_ARITHMETIC', counts=dict(counts),
        weighted_profile_outcomes=dict(weighted), source_bindings=len(packet['sources']),
        scope='No project imports. Generic rational inverses for A43 and white41; independent rational enclosure signs for sqrt2; complete saved379-slot-template domains and4983 edge/Schur checks; saved orientation multiplicities. No independent native-core reconstruction, old-block embedding completeness, all-label geometry or inherited overlap census replay.'), sort_keys=True))


if __name__ == '__main__':
    main()
