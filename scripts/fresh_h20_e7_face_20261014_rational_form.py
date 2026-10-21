#!/usr/bin/env python3
"""Rational similarity of the E7 face form; exact Hilbert symbols and PARI replay.
An existence certificate for a rational change of basis, not an explicit map
or a matching corner/chamber in the target reflection lattice.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt, prod
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_face_20261014_rational_form'
SOURCE = 'results/fresh_h20_e7_face_20261014.json'


def unit(x, p):
    a, b, v = x.numerator, x.denominator, 0
    while a % p == 0: a //= p; v += 1
    while b % p == 0: b //= p; v -= 1
    m = 8 if p == 2 else p
    return v, a*pow(b, -1, m) % m


def hilbert(a, b, p):
    r, u = unit(a, p); s, v = unit(b, p)
    if p == 2:
        exponent = ((u-1)*(v-1)//4 + r*(v*v-1)//8 + s*(u*u-1)//8) % 2
        return -1 if exponent else 1
    value = (-1)**((r*s*((p-1)//2)) % 2)
    if s % 2: value *= -1 if pow(u, (p-1)//2, p) == p-1 else 1
    if r % 2: value *= -1 if pow(v, (p-1)//2, p) == p-1 else 1
    return value


def run():
    d = json.loads((ROOT/SOURCE).read_text())
    for p, h in d['sources'].items(): assert sha256((ROOT/p).read_bytes()).hexdigest() == h
    A = [[F(x) for x in row] for row in d['intrinsic_normalization']['form']]
    M = [r[:] for r in A]; diag = []
    for i in range(14):
        pivot = M[i][i]; assert pivot; diag.append(pivot)
        for j in range(i+1, 14):
            for k in range(i+1, 14): M[j][k] -= M[j][i]*M[i][k]/pivot
    assert diag == [F(k+1, k) for k in range(1, 12)]+[F(2), F(3, 2), F(-2)]
    diagonals = [diag, [F(1)]*13+[F(-2)], [F(3)]*13+[F(-6)]]
    primes = [2, 3, 5, 7, 11]
    for D in diagonals:
        for x in D:
            v = abs(x.numerator*x.denominator)
            for p in primes:
                while v % p == 0: v //= p
            assert v == 1
    values = [[prod(hilbert(a, b, p) for a, b in combinations(D, 2)) for p in primes] for D in diagonals]
    assert values == [[-1, -1, 1, 1, 1], [1]*5, [-1, -1, 1, 1, 1]]
    assert all(sum(x > 0 for x in D) == 13 and sum(x < 0 for x in D) == 1 for D in diagonals)
    assert all(prod(row) == 1 for row in values)  # Real Hasse symbols are+1.
    ratio = prod(diagonals[0])/prod(diagonals[2])
    assert ratio > 0 and isqrt(ratio.numerator)**2 == ratio.numerator and isqrt(ratio.denominator)**2 == ratio.denominator
    gp = shutil.which('gp'); assert gp, 'PARI/GP is required for independent replay'
    statements = []
    for D in diagonals:
        for p in primes:
            expression = '*'.join(f'hilbert({a},{b},{p})' for a, b in combinations(D, 2))
            statements.append(f'print({expression});')
    result = subprocess.run([gp, '-q', '-f'], input='\n'.join(statements)+'\nquit;\n', text=True, capture_output=True, check=True)
    assert not result.stderr, result.stderr
    pari = list(map(int, result.stdout.split())); assert pari == sum(values, [])
    paths = [SOURCE, f'scripts/{TAG}.py']
    return dict(status='PASS_E7_FACE_RATIONALLY_SIMILAR_TO_F2_H13_WITH_SCALE3', sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        diagonalizations=[list(map(str, D)) for D in diagonals], finite_places=primes, hasse_invariants=values, signature=[13, 1],
        determinant_ratio_to_scaled_target=str(ratio), determinant_ratio_square_root=str(F(isqrt(ratio.numerator), isqrt(ratio.denominator))),
        independent_PARI_hilbert_products=pari, omitted_places='All other odd-prime coefficients are units, hence Hilbert symbols+1. At infinity there is only one negative coefficient, hence Hasse+1.',
        formula_reference='https://math.mit.edu/classes/18.782/2013fa/LectureNotes10.pdf',
        scope='Rational quadratic-form classification gives diag(A11,A2,-2) isometric overQ to3*diag(1^13,-2), but not to the unscaled latter form. No explicit rational isometry, integral lattice identification, reflection-host transfer or matching A11+A2 corner is supplied. The separate nineteen-wall source-corner obstruction is unchanged.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); a = ap.parse_args()
    d = run(); p = ROOT/f'results/{TAG}.json'
    if a.verify: assert d == json.loads(p.read_text())
    else: p.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], d['hasse_invariants'])
