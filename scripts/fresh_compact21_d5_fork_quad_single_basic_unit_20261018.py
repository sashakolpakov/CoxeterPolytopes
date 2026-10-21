#!/usr/bin/env python3
"""Tiny exact certificate for all one-basic-unit D5 fork quad patterns.

The completeness and all-distance exclusion arguments are in the companion
written proof. This program expands determinant identities, checks every
finite endpoint pattern in that proof, and checks its algebraic obstructions.
It neither searches opposite distances nor solves a global charging problem.
"""
import argparse
import fcntl
from fractions import Fraction as Q
from hashlib import file_digest
from itertools import permutations, product
import json
from math import isqrt
from pathlib import Path
import resource
import signal
import sys
import fresh_compact29_small_faces_20261015 as P

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_d5_fork_quad_single_basic_unit_20261018'


class Poly:
    """Sparse rational polynomials in x,y; no interpolation or float arithmetic."""
    def __init__(self, value=0):
        self.t = ({k: Q(v) for k, v in value.items() if v}
                  if isinstance(value, dict) else ({(0, 0): Q(value)} if value else {}))

    def __add__(self, other):
        other = other if isinstance(other, Poly) else Poly(other)
        terms = self.t.copy()
        for k, v in other.t.items():
            terms[k] = terms.get(k, Q())+v
        return Poly(terms)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.t.items()})

    def __sub__(self, other):
        return self+-Poly(other) if not isinstance(other, Poly) else self+-other

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        other = other if isinstance(other, Poly) else Poly(other)
        terms = {}
        for (a, b), u in self.t.items():
            for (c, d), v in other.t.items():
                k = (a+c, b+d)
                terms[k] = terms.get(k, Q())+u*v
        return Poly(terms)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        assert isinstance(exponent, int) and exponent >= 0
        value = Poly(1)
        for _ in range(exponent):
            value *= self
        return value

    def __eq__(self, other):
        return self.t == (other.t if isinstance(other, Poly) else Poly(other).t)

    def encode(self):
        return [{'x_degree': a, 'y_degree': b, 'coefficient': str(c)}
                for (a, b), c in sorted(self.t.items())]


def det_poly(matrix):
    n = len(matrix)
    total = Poly()
    for perm in permutations(range(n)):
        inversions = sum(perm[j] > perm[i] for i in range(n) for j in range(i))
        term = Poly((-1)**inversions)
        for i, j in enumerate(perm):
            term *= matrix[i][j]
        total += term
    return total


def sigma(x):
    a, b, c, d = x.v
    assert not b and not d
    return P.K((a, 0, -c, 0))


def integral(x):
    a, b, c, d = x.v
    assert not b and not d
    return (2*c).denominator == 1 and (a-c).denominator == 1


def rational_square(q):
    q = Q(q)
    return q >= 0 and isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator


def certify_nonsquare(x):
    a, b, c, d = x.v
    assert not b and not d
    norm = x*sigma(x)
    assert not norm.v[1] and not norm.v[2] and not norm.v[3]
    if not c:
        assert not rational_square(a) and not rational_square(a/5)
        reason = 'rational value is neither a rational square nor five times one'
    else:
        assert not rational_square(norm.v[0])
        reason = 'nonsquare rational field norm'
    return {'value': P.enc(x), 'norm': P.enc(norm), 'obstruction': reason}


def basic_payment(matrix, i, j):
    frontier, seen = {i}, {i}
    for _ in range(2):
        frontier = {v for u in frontier for v in range(len(matrix))
                    if v not in seen and matrix[u][v]}
        if j in frontier:
            return 1
        seen |= frontier
    return 0


def run():
    one = P.K(1)
    sqrt2, sqrt5 = P.K((0, 1, 0, 0)), P.K((0, 0, 1, 0))
    phi = (1+sqrt5)/2
    labels = {3: one, 4: sqrt2}
    # These seven endpoint types are complete by the rank-five spherical
    # classification containing A4. The proof does not sample large labels.
    endpoints = [(None, 0), (0, 3), (0, 4), (1, 3), (2, 3), (3, 3), (3, 4)]
    zero = []
    for contact, label in endpoints:
        g = [[P.K(2*(i == j)) for j in range(5)] for i in range(5)]
        for i in range(3):
            g[i][i+1] = g[i+1][i] = -one
        if contact is not None:
            g[contact][4] = g[4][contact] = -labels[label]
        assert P.psd_status(g) == 1
        if basic_payment(g, 3, 4) == 0:
            zero.append((contact, label))
    assert zero == [(None, 0), (0, 3), (0, 4), (1, 3)]
    permitted = []
    for a, b in product(zero, repeat=2):
        g = [[P.K(2*(i == j)) for j in range(5)] for i in range(5)]
        for i, j in ((0, 1), (1, 2)):
            g[i][j] = g[j][i] = -one
        for new, (contact, label) in ((3, a), (4, b)):
            if contact is not None:
                g[contact][new] = g[new][contact] = -labels[label]
        if P.psd_status(g) == 1 and not basic_payment(g, 3, 4):
            permitted.append((a, b))
    expected = [((None, 0), (None, 0))]
    for a in zero[1:]:
        expected += [((None, 0), a), (a, (None, 0))]
    expected += [((0, 3), (1, 3)), ((1, 3), (0, 3))]
    assert set(permitted) == set(expected) and len(permitted) == 9

    x, y = Poly({(1, 0): 1}), Poly({(0, 1): 1})
    # Here x=X is a squared opposite product; c has already been eliminated.
    r03 = [[Q(5, 4)-x*Q(1, 2), Q(-3, 4), Q(-1, 4)],
           [Q(-3, 4), Q(5, 4), y-Q(1, 4)],
           [Q(-1, 4), y-Q(1, 4), Q(5, 4)]]
    r13 = [[Q(5, 4)-x*Q(1, 2), Q(-3, 4), Q(-1, 2)],
           [Q(-3, 4), Q(5, 4), y-Q(1, 2)],
           [Q(-1, 2), y-Q(1, 2), 1]]
    f03 = (5-2*x)*(3+y-2*y*y)+3*y-7
    f13 = (5-2*x)*(1+y-y*y)+3*y-5
    assert 8*det_poly(r03) == f03
    assert 4*det_poly(r13) == f13
    # The rescaled label-four root has norm4 and contact -2. Its residual
    # matrix is exactly r13 with z=G(d,old opposite) in place of y.
    r01 = [[Q(5, 4), Q(-3, 4), x-Q(1, 4), Q(-1, 2)],
           [Q(-3, 4), Q(5, 4), Q(-1, 4), y-Q(1, 2)],
           [x-Q(1, 4), Q(-1, 4), Q(5, 4), Q(-1, 2)],
           [Q(-1, 2), y-Q(1, 2), Q(-1, 2), 1]]
    f01 = (2*x-3)*((x+1)*y*y-(x+2)*y)-2*(x-1)**2
    assert 2*det_poly(r01) == f01

    finite_x = [phi*phi, P.K(3), (5+sqrt5)/2]
    obstructions = {}
    for name, disc in (
            ('N-P03', lambda a: 25*a*a-50*a+9),
            ('N-P13_or_N-P04', lambda a: 5*a*a-14*a+9)):
        obstructions[name] = []
        for xx in finite_x:
            a = 5-2*xx
            record = certify_nonsquare(disc(a))
            record['finite_X'] = P.enc(xx)
            obstructions[name].append(record)
    # Finite y alternatives in the two norm-two attached-root branches.
    yphi_x03 = (25+12*phi)/11
    assert (5-2*yphi_x03)*(3-phi-2*phi*phi)-3*phi-7 == 0
    assert not integral(yphi_x03)
    yminus1_x13 = P.K(Q(13, 2))
    yphi_x13 = 2+5*phi/4
    assert (5-2*yminus1_x13)*(-1)-8 == 0
    assert (5-2*yphi_x13)*(1-phi-phi*phi)-3*phi-5 == 0
    assert not integral(yminus1_x13) and not integral(yphi_x13)
    # In the rescaled label-four branch finite z can only be -2.
    # Labels3,5 require sqrt2 in the field; label6 requires sqrt6;
    # label10 requires sqrt(5+sqrt5). Each is excluded exactly.
    label4_field_obstructions = [certify_nonsquare(v) for v in (P.K(2), P.K(6), 5+sqrt5)]
    finite_z_x = P.K(Q(18, 5))
    assert (5-2*finite_z_x)*(-5)-11 == 0 and not integral(finite_z_x)

    # The P03-P13 branch: x< -1, so a finite x is only -phi.
    # A finite y is -1 or -phi. All three quadratic discriminants fail.
    obstructions['P03-P13'] = [certify_nonsquare(v) for v in
                              (73+89*phi, P.K(104), 153+224*phi)]
    # Direct coefficient checks after substituting the finite side.
    xx = -phi
    a, b, c = (2*xx-3)*(xx+1), -(2*xx-3)*(xx+2), -2*(xx-1)*(xx-1)
    assert b*b-4*a*c == 73+89*phi
    assert 3*3-4*(4*phi)*(-9*phi-5) == 153+224*phi
    assert 2*det_poly(r01) == f01

    sources = [f'scripts/{TAG}.py', 'scripts/fresh_compact29_small_faces_20261015.py',
               'results/fresh_compact21_d5_fork_quad_local_witness_20261018.json',
               'results/fresh_compact21_single_field_20261016.json']
    pins = {}
    for path in sources:
        with (ROOT/path).open('rb') as stream:
            pins[path] = file_digest(stream, 'sha256').hexdigest()
    return {'status': 'PASS_COMPLETE_ONE_BASIC_UNIT_D5_FORK_QUAD_PATTERN',
            'sources': pins, 'field': 'Q(sqrt5)',
            'finite_A4_endpoint_types': endpoints,
            'zero_basic_endpoint_types': zero,
            'ordered_zero_extra_quad_patterns': permitted,
            'unordered_patterns': 5, 'excluded_patterns': 4,
            'corank_polynomials': {'N-P03': f03.encode(), 'N-P13_or_N-P04': f13.encode(),
                                   'P03-P13': f01.encode()},
            'finite_opposite_nonsquare_obstructions': obstructions,
            'finite_y_nonintegral_X': P.enc([yphi_x03, yminus1_x13, yphi_x13]),
            'rescaled_label4_field_obstructions': label4_field_obstructions,
            'rescaled_label4_finite_z_nonintegral_X': P.enc(finite_z_x),
            'remaining_pattern': 'N-N, classified by the pinned earlier certificate and written proof',
            'conditional_conclusion': 'Every one-basic-unit D5 fork quad excludes all designated E8 triangles',
            'global_charge_saving_proved': False, 'new_dimension_bound': None}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.alarm(59)
    rss = lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    assert rss() < 1024**3
    lock_path = ROOT/'formal/.proof-batches/exclusive.lock'
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = run()
        assert rss() < 1024**3
        target = ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text()) == json.loads(json.dumps(result))
        else:
            assert not target.exists()
            target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'endpoint_types': 7,
                      'ordered_quad_patterns': 9, 'excluded_unordered_patterns': 4,
                      'peak_rss_bytes': rss(), 'new_dimension_bound': None}))
