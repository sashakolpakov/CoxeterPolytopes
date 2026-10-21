#!/usr/bin/env python3
"""Fixed exact Q(sqrt13) two-edge algebra; use outer --parent-lock.

No check runs on import. Written norm-pool/classification/germ arguments
remain inputs. No optimizer, cone census, height or Pell search.
"""
import argparse
import hashlib
import json
import signal
from dataclasses import dataclass
from fractions import Fraction as F
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/compact_h9_sqrt13_two_edges_20261020.json'
PINS = {
    'results/compact_h9_sqrt13_b9_two_edge_route_20261020.md': 'b922647902d32beccd98fbe5a2781ddcbb9cafaeb10a12902eb3ebbce809547f',
    'results/review_compact_h9_sqrt13_two_edges_resource_20261020.md': '20b3381d2d20ff7efb0d4df7b0642d1b8dc0eecb1924005bc0502878467f96de',
    'scripts/run_bounded_research_check_20261018.py': '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py': '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py': 'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py': '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}


@dataclass(frozen=True, eq=False)
class Q:
    """Exact a+b*omega, omega^2=omega+3; integrality means a,b in Z."""
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    @staticmethod
    def cast(value):
        return value if isinstance(value, Q) else Q(value)

    def __eq__(self, other):
        try:
            other = Q.cast(other)
        except (TypeError, ValueError):
            return NotImplemented
        return self.a == other.a and self.b == other.b

    def __add__(self, other):
        other = Q.cast(other)
        return Q(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.a, -self.b)

    def __sub__(self, other):
        return self+-Q.cast(other)

    def __rsub__(self, other):
        return Q.cast(other)+-self

    def __mul__(self, other):
        other = Q.cast(other)
        return Q(self.a*other.a+3*self.b*other.b,
                 self.a*other.b+self.b*other.a+self.b*other.b)

    __rmul__ = __mul__

    def conjugate(self):
        return Q(self.a+self.b, -self.b)

    def norm(self):
        return self.a*self.a+self.a*self.b-3*self.b*self.b

    def __truediv__(self, other):
        other = Q.cast(other)
        assert other.norm() != 0
        c = self*other.conjugate()
        return Q(c.a/other.norm(), c.b/other.norm())

    def __rtruediv__(self, other):
        return Q.cast(other)/self

    def sign(self):
        a, b = 2*self.a+self.b, self.b
        s = lambda t: (t > 0)-(t < 0)
        if not a:
            return s(b)
        if not b or s(a) == s(b):
            return s(a)
        return s(a)*s(a*a-13*b*b)

    def integral(self):
        return self.a.denominator == self.b.denominator == 1


def enc(x):
    x = Q.cast(x)
    return [str(x.a), str(x.b)]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def add(x, y):
    return [a+b for a, b in zip(x, y)]


def scale(a, x):
    return [a*b for b in x]


def mv(matrix, x):
    return [sum(a*b for a, b in zip(row, x)) for row in matrix]


def transpose(matrix):
    return list(map(list, zip(*matrix)))


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def positive_definite(matrix):
    a = [row[:] for row in matrix]
    for k in range(len(a)):
        pivot = a[k][k]
        assert pivot.sign() > 0
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                a[i][j] = a[i][j]-a[i][k]*a[k][j]/pivot


def check():
    paths = list(PINS)+[str(Path(__file__).resolve().relative_to(ROOT))]
    sources = {p: digest(ROOT/p) for p in paths}
    assert all(sources[p] == value for p, value in PINS.items())
    omega, one = Q(0, 1), Q(1)
    eta = omega+1
    sq = lambda x: x*x
    assert sq(omega) == omega+3 and sq(eta) == 3*eta+1
    assert eta.norm() == -1 and 1/eta == eta-3
    assert omega.integral() and not (omega/2).integral()
    assert not Q(F(1, 2), F(1, 2)).integral()
    assert omega.conjugate() == 1-omega and eta.conjugate() == -1/eta
    assert (omega-2).sign() > 0 and (omega-3).sign() < 0
    assert (Q(11)*Q(11)-9*Q(13)).sign() > 0
    basis = [[Q(int(i == j)) for i in range(10)] for j in range(10)]
    gram = [[Q() for _ in range(10)] for _ in range(10)]
    for i in range(10):
        gram[i][i] = -eta if i == 0 else one
    dot = lambda x, y: sum(a*b for a, b in zip(x, mv(gram, y)))
    alpha = [add(scale(-1, basis[i]), basis[i+1]) for i in range(1, 9)]
    alpha.append(scale(-1, basis[9]))
    p = basis[0]
    v = add(p, scale(omega, basis[1]))
    q = add(scale(omega, p), scale(eta, basis[1]))
    u = add(scale(eta, p), scale(eta+1, add(basis[1], basis[2])))
    d = add(u, basis[3])
    w = add(scale(2, u), v)
    H = add(scale(3*eta, p), add(scale(4*eta+1, basis[1]), scale(3*eta+1, basis[2])))
    roots = alpha+[v, d]
    assert [dot(x, x) for x in roots] == [Q(2)]*8+[one, Q(2), Q(2)]
    assert dot(p, p) == -eta and dot(q, q) == -2*eta and dot(H, H) == -eta

    # Three first-edge profile identities.
    first = []
    for name, norm, extra, R in [('short', 1, 0, 1), ('long', 2, 0, 2),
                                 ('B9', 2, 1, 1)]:
        const = scale(extra, basis[2])
        assert dot(const, const)-norm == -R
        assert dot(p, const) == dot(basis[1], const) == 0
        first.append({'name': name, 'norm': norm, 'R': R})
    assert sq(omega)-eta == 2 and dot(v, alpha[0]) == -omega
    tau_v = eta/omega
    assert (eta-sq(tau_v)).sign() > 0
    # Physical x<=1 and |x'|<=sqrt(2eta)<sqrt13-1 force integer x.
    assert (Q(16)-2*Q(0, 2)-2*eta).sign() > 0  # (sqrt13-1)^2-2eta
    assert dot(q, v) == 0 and dot(q, p).sign() < 0
    assert all(dot(q, r).sign() <= 0 for r in alpha+[v])
    assert [i for i, r in enumerate(alpha+[v]) if dot(q, r) == Q()] == list(range(1, 10))
    assert sq(omega) == eta+2

    # Six second-edge norm identities.
    profiles = [('isolated_short', 1, 0, 0, 4), ('isolated_long', 2, 0, 0, 8),
                ('B2_v', 1, 1, 0, 2), ('A2_v', 2, 1, 0, 6),
                ('B8', 2, 0, 1, 4), ('B9_bridge', 2, 1, 1, 2)]
    second = []
    for name, norm, ell, t, R in profiles:
        const = add(scale(F(-ell, 2), v), scale(t, basis[3]))
        ux, uz = scale(F(1, 2), q), scale(F(1, 2), basis[2])
        co = [dot(const, const)-norm, 2*dot(const, ux), 2*dot(const, uz),
              dot(ux, ux), dot(uz, uz), 2*dot(ux, uz)]
        assert co == [Q(F(-R, 4)), Q(), Q(), -eta/2, Q(F(1, 4)), Q()]
        second.append({'name': name, 'norm': norm, 'ell': ell, 't': t,
                       'R': R, 'norm_coefficients': list(map(enc, co))})
    assert dot(u, u) == one and dot(u, v) == -one
    tau_u = sq(eta)/(eta+1)
    assert (2*eta-sq(tau_u)).sign() > 0
    assert 4*sq(eta+1)-2*eta*sq(eta) == 2
    # Check binary parity on all four residue classes in O/(2)=F4.
    for a in range(2):
        for b in range(2):
            X = Q(a, b)
            for ell in range(2):
                time = (omega*X-ell)/2
                z1 = (eta*X-omega*ell)/2
                assert (time.integral() and z1.integral()) == (a == ell and b == ell)

    # Search z=Y=X/2 for even rows, z=X for odd rows.
    # Each has 0<z<6, |z'|<7/2, hence b=0,1,2 and -4<=a<=5.
    assert (36-3*sq(eta)).sign() > 0 and (Q(F(49, 4))-3*eta).sign() > 0
    assert (Q(13)-Q(F(49, 4))).sign() > 0
    assert (9*Q(13)-Q(F(361, 4))).sign() > 0
    assert (5-2*omega).sign() > 0
    cases = [('R4_even', 4, 0, 2, [(1, 0)]),
             ('R8_even', 8, 0, 2, [(1, 0), (0, 1), (1, 1)]),
             ('R6_odd', 6, 1, 1, [(-1, 1), (1, 1), (3, 1)]),
             ('R2_odd', 2, 1, 1, [(1, 1)])]
    domains = {}
    for name, R, ell, divisor, expected in cases:
        found = []
        for b in range(3):
            for a in range(-4, 6):
                z = Q(a, b)
                X = divisor*z
                if z.sign() <= 0 or (R*sq(eta)/2-sq(X)).sign() < 0:
                    continue
                if (R*eta/2-sq(X.conjugate())).sign() < 0:
                    continue
                if ell and (a % 2 != 1 or b % 2 != 1):
                    continue
                found.append((a, b))
        assert found == expected
        domains[name] = found
    # Fixed square tables and nonsquare norms.
    radicands = [(2*eta+1, 3), (2*eta+2, 12),
                 (eta+1, 3), (5*eta+3, 29), (13*eta+5, 51)]
    for rad, expected in radicands:
        assert rad.norm() == expected and isqrt(expected)**2 != expected
    assert 2*eta*sq(omega)+2 == 2*sq(eta+1)
    assert 2*eta*sq(eta)+2 == 4*sq(eta+1)
    # sqrt2 is not in Q(sqrt13): squaring a+b*sqrt13 has ab=0,
    # reducing to the rational nonsquares2 or26.
    assert all(isqrt(n)**2 != n for n in (2, 26))
    for X, rad in [(eta-2, eta+1), (eta, 5*eta+3), (eta+2, 13*eta+5)]:
        assert (2*eta*sq(X)+6)/4 == rad
    assert (2*eta*sq(eta)+2)/4 == sq(eta+1)
    # Survivors.
    def reconstruct(X, Z, ell, t):
        return add(add(scale(X/2, q), scale(F(-ell, 2), v)),
                   add(scale(Z/2, basis[2]), scale(t, basis[3])))
    assert reconstruct(eta, 2*(eta+1), 1, 0) == u
    assert reconstruct(eta, 2*(eta+1), 1, 1) == d
    assert reconstruct(2*eta, 4*(eta+1), 0, 0) == w
    A = scale(-1, basis[3])
    assert A == [sum(alpha[j][i] for j in range(2, 9)) for i in range(10)]
    assert u == add(d, A) and w == add(add(scale(2, d), v), scale(2, A))

    # All55 actual pairs and integral reflections of the fixed eleven walls.
    pair_rows = []
    special = {(0, 9): -omega, (1, 10): -eta, (2, 10): Q(-1), (9, 10): Q(-1)}
    for i, r in enumerate(roots):
        assert all(x.integral() for x in r) and any(x.norm() in (1, -1) for x in r)
        qr = dot(r, r)
        assert all((2*x*y/qr).integral() for x in r for y in mv(gram, r))
        for j in range(i+1, 11):
            product = dot(r, roots[j])
            expected = Q(-1) if i < 8 and j == i+1 else Q()
            expected = special.get((i, j), expected)
            assert product == expected
            ratio4 = 4*sq(product)/(qr*dot(roots[j], roots[j]))
            label = 'ultraparallel' if (i, j) in ((0, 9), (1, 10)) else (
                2 if product == Q() else 4 if ratio4 == Q(2) else 3)
            assert label == 'ultraparallel' and (ratio4-4).sign() > 0 or (
                label != 'ultraparallel' and ratio4 == Q({2: 0, 3: 1, 4: 2}[label]))
            pair_rows.append([i, j, enc(product), label])
    old = [[dot(x, y) for y in alpha] for x in alpha]
    new_roots = [v, d]+alpha[2:]
    new = [[dot(x, y) for y in new_roots] for x in new_roots]
    assert old == new
    positive_definite(old)
    middle_roots = alpha[1:]+[v]
    positive_definite([[dot(x, y) for y in middle_roots] for x in middle_roots])
    assert all(dot(p, r) == Q() for r in alpha)
    assert [i for i, r in enumerate(roots) if dot(H, r) == Q()] == list(range(2, 11))
    assert [dot(H, roots[i]) for i in (0, 1)] == [-eta, -3*eta-1]
    assert dot(p, H).sign() < 0 and all(dot(H, r).sign() <= 0 for r in roots)

    # Original-lattice columns, exact inverse, involution and time orientation.
    col1 = add(scale(-(eta+1), p), add(scale(-2*eta, basis[1]), scale(-(eta+1), basis[2])))
    columns = [H, col1, scale(-1, u)]+basis[3:]
    matrix = transpose(columns)
    assert [[dot(x, y) for y in columns] for x in columns] == gram
    assert all(x.integral() for col in columns for x in col)
    I = [[Q(int(i == j)) for j in range(10)] for i in range(10)]
    gram_inverse = [row[:] for row in gram]
    gram_inverse[0][0] = -1/eta
    inverse = mm(mm(gram_inverse, transpose(matrix)), gram)
    assert mm(matrix, inverse) == I and inverse == matrix and mm(matrix, matrix) == I
    assert all(x.integral() for row in inverse for x in row)
    assert mv(matrix, p) == H and mv(matrix, alpha[0]) == v
    assert mv(matrix, alpha[1]) == d and mv(matrix, v) == alpha[0]
    assert mv(matrix, d) == alpha[1]
    assert all(mv(matrix, x) == x for x in alpha[2:])
    trace = sum(matrix[i][i] for i in range(10))
    assert trace == Q(6)
    center = [one]+[Q(F(i, 100)) for i in range(9, 0, -1)]
    assert dot(center, center).sign() < 0 and dot(center, p).sign() < 0
    assert all(dot(center, r).sign() < 0 for r in roots)
    assert sources == {path: digest(ROOT/path) for path in paths}, 'source changed during check'
    return {'status': 'PASS_FIXED_TWO_EDGES_AND_CHAMBER_INVOLUTION_ALGEBRA',
            'sources': sources, 'owns_lock': False, 'field_encoding': 'a+b*omega; rational strings',
            'scope': 'written primitive alphabet, finite classification and full actual-germ arguments are inputs; no volume or reflectivity conclusion',
            'first_profiles': first, 'second_profiles': second, 'second_domains': domains,
            'nonsquare_norms': [n for _, n in radicands],
            'roots': [list(map(enc, r)) for r in roots], 'pairs': pair_rows,
            'old_new_corner_gram': [[enc(x) for x in row] for row in old],
            'first_endpoint': list(map(enc, q)), 'second_endpoint': list(map(enc, H)),
            'strict_common_center': list(map(enc, center)),
            'involution_columns': [list(map(enc, r)) for r in columns],
            'involution_trace': enc(trace), 'infinite_order': False, 'new_example': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.write and CERT.exists():
        raise SystemExit('refusing to overwrite certificate')
    signal.alarm(45)
    result = json.loads(json.dumps(check(), sort_keys=True))
    payload = json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n'
    assert len(payload.encode()) <= 15360, 'fixed certificate output limit'
    assert result['sources'] == {p: digest(ROOT/p) for p in result['sources']}
    if args.verify:
        assert CERT.stat().st_size <= 15360
        assert result == json.loads(CERT.read_text()), 'saved certificate differs'
    else:
        with CERT.open('x') as stream:
            stream.write(payload)
    signal.alarm(0)
    print(json.dumps({'status': result['status'], 'mode': 'replay' if args.verify else 'generation',
                      'pair_count': 55, 'second_profiles': 6, 'involution_trace': result['involution_trace']}))


if __name__ == '__main__':
    main()
