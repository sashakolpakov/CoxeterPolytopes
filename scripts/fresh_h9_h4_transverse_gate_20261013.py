#!/usr/bin/env python3
"""Standalone rational-pair arithmetic for an all-label H4 transverse gate."""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h4_transverse_gate_20261013'


def q(a=0, b=0):
    return F(a), F(b)


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def mul(x, y):
    a, b = x
    c, d = y
    return a*c + 5*b*d, a*d + b*c


def neg(x):
    return -x[0], -x[1]


def inv(x):
    a, b = x
    den = a*a - 5*b*b
    assert den
    return a/den, -b/den


def sign(x):
    a, b = x
    if b == 0:
        return (a > 0) - (a < 0)
    if a >= 0 and b > 0:
        return 1
    if a <= 0 and b < 0:
        return -1
    value = a*a - 5*b*b if a > 0 else 5*b*b - a*a
    return (value > 0) - (value < 0)


def total(xs):
    value = q()
    for x in xs:
        value = add(value, x)
    return value


def enc(x):
    return list(map(str, x))


def run(verify=False):
    G = [[q(2*int(i == j)) for j in range(4)] for i in range(4)]
    for i, x in enumerate([q('-1/2', '-1/2'), q(-1), q(-1)]):
        G[i][i+1] = G[i+1][i] = x
    B = [
        [q(14, 6), q('33/2', '15/2'), q(11, 5), q('11/2', '5/2')],
        [q('33/2', '15/2'), q(21, 9), q(14, 6), q(7, 3)],
        [q(11, 5), q(14, 6), q(10, 4), q(5, 2)],
        [q('11/2', '5/2'), q(7, 3), q(5, 2), q(3, 1)],
    ]
    assert all(total(mul(G[i][k], B[k][j]) for k in range(4)) == q(int(i == j))
               for i in range(4) for j in range(4))
    assert all(sign(x) > 0 for row in B for x in row)
    T = [row[:] for row in G]
    pivots = []
    for k in range(4):
        pivot = T[k][k]
        assert sign(pivot) > 0
        pivots.append(pivot)
        for i in range(k+1, 4):
            for j in range(k+1, 4):
                T[i][j] = add(T[i][j], neg(mul(mul(T[i][k], inv(pivot)), T[k][j])))
    minimum = q(3, 1)
    assert all(sign(add(B[i][i], neg(minimum))) >= 0 for i in range(4))
    assert B[3][3] == minimum and sign(add(minimum, q(-2))) > 0
    out = {
        'status': 'PASS_EXACT_ALL_LABEL_H4_TRANSVERSE_SCHUR_GATE',
        'arithmetic': 'Rational pairs a+b*sqrt(5), positive real embedding',
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'norm_two_H4_Gram': [[enc(x) for x in row] for row in G],
        'inverse': [[enc(x) for x in row] for row in B],
        'positive_ldl_pivots': list(map(enc, pivots)),
        'minimum_inverse_diagonal': enc(minimum),
        'projected_norm_lower_bound_for_nonright_coxeter_cap': enc(minimum),
        'perpendicular_component_norm_upper_bound': enc(add(q(2), neg(minimum))),
        'scope': 'For a norm-two real Coxeter normal acute to all four norm-two H4 simple roots, '
                 'either all four pairings vanish or at least one has magnitude at least one. '
                 'The positive inverse makes its H4 projection square at least 3+sqrt(5)>2. '
                 'The five-wall Gram cannot be positive semidefinite, and its perpendicular '
                 'component is timelike. Hence any wall meeting the whole H4 face must be '
                 'orthogonal to H4. Includes every finite angle label, parallel and dotted '
                 'relations. Does not exclude caps away from the face or any compact H9 completion.',
    }
    path = ROOT/'results'/f'{TAG}.json'
    if verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    run(ap.parse_args().verify)
