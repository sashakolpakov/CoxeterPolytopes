#!/usr/bin/env python3
"""Exact polynomial ledger for the analytic nine-direction frame proof.

This checks scalar identities and constants, not the geometric frame premises.
Standard library only; no searches or optimization. Write is exclusive.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'results/mixed_eleven_nine_direction_global_20261020.json'
PINS = [
    'results/right_angled_mixed_eleven_nine_direction_global_proof_20261020.md',
    'scripts/check_mixed_eleven_nine_direction_global_20261020.py',
]


def trim(p):
    p = [F(c) for c in p]
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p


def add(a, b):
    p = [F(0)] * max(len(a), len(b))
    for i, c in enumerate(a):
        p[i] += c
    for i, c in enumerate(b):
        p[i] += c
    return trim(p)


def scale(a, c):
    return trim([F(c) * x for x in a])


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    p = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            p[i + j] += x * y
    return trim(p)


def power(a, n):
    p = [F(1)]
    for _ in range(n):
        p = mul(p, a)
    return p


def derivative(a):
    return trim([i * a[i] for i in range(1, len(a))] or [0])


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def check():
    x, one, mu = [0, 1], [1], [1, 1]
    # Global frame and clipping constants, using integer square comparisons.
    assert 139**2 > 80**2 * 3       # r>(11/20).
    assert 100 * 121 > 324 * 37    # 10/sqrt37>18/11.
    assert F(15, 11) - F(18, 11) == -F(3, 11)
    clip = F(11**6, 128 * 18**3)
    assert clip < 3
    K = F(20, 11)**2 + 4 * F(20, 11)
    assert K == F(1280, 121)
    assert K + 3 == F(1643, 121)
    assert 11 > 9                  # 1/(2sqrt11)<1/6.
    # Exact reciprocal remainder after d/mu=x.
    taylor = add(mul([1, -2, 3, -4], power(mu, 2)),
                 mul([0, 0, 0, 0, 1], [5, 4]))
    assert taylor == one
    # The fourth-power tail inequality has a manifestly nonnegative factor.
    tail = [256, -256, 0, 0, 27]
    assert tail == mul(power([-4, 3], 2), [16, 8, 3])
    # Unordered-three-subset inclusion probabilities give these moments.
    p1, p2, p3 = F(1, 3), F(1, 12), F(1, 84)
    assert p1 - p2 == F(1, 4)
    assert p1 - 3 * p2 + 2 * p3 == F(3, 28)
    assert 3 * p2 - 6 * p3 == F(5, 28)
    assert p1 - 7 * p2 + 12 * p3 == -F(3, 28)
    assert F(5, 28) - F(3, 28) / 9 == F(1, 6)
    # The zero-sum cubic proof: 7M/8 with M^2=8V/9.
    assert F(7, 8)**2 * F(8, 9) == F(49, 72)
    assert F(3, 7)**2 * F(49, 72) == F(1, 8)
    # V upper envelope and maximum 2/11 by completing its square.
    v = [0, 2, -F(11, 2)]
    assert sub([F(2, 11)], v) == scale(power([-F(2, 11), 1], 2), F(11, 2))
    q = [4, -11]
    product = mul(power(q, 2), power(mu, 5))
    assert derivative(product) == mul(mul(q, power(mu, 4)), [-2, -77])
    # Clear the (1+x)^5 denominator of 1-x/6-envelope(v).
    lhs = mul([1, -F(1, 6)], power(mu, 5))
    lhs = sub(lhs, power(mu, 3))
    lhs = sub(lhs, scale(mul(v, mu), F(3, 4)))
    lhs = sub(lhs, scale(v, F(1, 6)))
    lhs = sub(lhs, scale(mul(power(v, 2), power(mu, 5)), F(1643, 726)))
    R = sub([F(233, 24), F(275, 24), F(10, 3), F(1, 6), -F(1, 6)],
            scale(product, F(1643, 2904)))
    assert lhs == mul(power(x, 2), R)
    margin = F(233, 24) - F(1643, 2904) * 16
    assert margin == F(635, 968) and margin > 0
    assert F(4, 11) < 1
    assert 84 * F(63, 100) == F(1323, 25)
    assert F(756, 25) + F(27, 100) * 84 == F(1323, 25)
    return dict(
        schema=1, status='PASS_EXACT_GLOBAL_NINE_DIRECTION_IDENTITIES',
        pins={p: digest(ROOT / p) for p in PINS},
        ratio_lower='11/20', clipping_coefficient_upper='3',
        reciprocal_remainder_upper=str(K),
        polynomial_margin=str(margin),
        cleared_difference_coefficients=[str(c) for c in lhs],
        R_coefficients=[str(c) for c in R],
        contract='mean positive part <=63/100-delta/18',
        equality='Nine selected entries1/3, other two0',
        scope='Scalar identities plus written analytic frame proof; no geometry or dimension theorem',
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    assert not (args.write and args.verify)
    result = check()
    if args.write:
        with CERT.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    else:
        assert json.loads(CERT.read_text()) == result, 'Certificate mismatch'
    print(json.dumps({k: v for k, v in result.items() if k != 'pins'}))


if __name__ == '__main__':
    main()
