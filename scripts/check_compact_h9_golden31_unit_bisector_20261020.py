#!/usr/bin/env python3
"""Fixed golden31 identities. Root-owned bounded wrapper only; no search."""
import argparse
import hashlib
import json
import os
import resource
import signal
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'compact_h9_golden31_unit_bisector_20261020'
CERT = ROOT / ('results/' + STEM + '.json')
MAX_CERT = 128 * 1024
PINS = {
    'results/compact_h9_golden31_unit_bisector_preparation_20261020.md':
        '7259c05e994fd48d8bb5f70c172e98f9c0120bbfa13fc58f3a6b08174e5f9ea1',
    'results/compact_h9_golden31_fixed_checker_preparation_20261020.md':
        '550a64e17a4008434060c46766c3cdeade08be575a9dbef86591e5a545711294',
    'results/review_compact_h9_golden31_unit_bisector_root_20261020.md':
        '4cf76e054d1b8461a46d86fce9421a0b7139b741f6f45938761a6a15a71a80f9',
    'results/bugaenko_facet_binary_extension_and_collar_audit_20261020.md':
        '774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227',
    'results/compact_h9_golden_discriminant59_edge_route_20261020.md':
        '725d8775cebec51d8c18b4b0a0d3be0255e8b1472621adb4d3cefa656adb429e',
    'data/known_examples/bugaenko_compact_h8.json':
        'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
    'data/known_examples/bugaenko_vinberg_transcript.json':
        'c4178fc4678bacb48c430b2cec18e16577078824fdd6dea39375e63ed6ae015a',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
}


class Q:
    """Exactly a+b*phi; phi^2=phi+1. No floating arithmetic."""
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        self.a, self.b = F(a), F(b)

    @staticmethod
    def of(value):
        return value if isinstance(value, Q) else Q(value)

    def __add__(self, other):
        other = Q.of(other)
        return Q(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-Q.of(other))

    def __rsub__(self, other):
        return Q.of(other) - self

    def __mul__(self, other):
        other = Q.of(other)
        return Q(self.a*other.a + self.b*other.b,
                 self.a*other.b + self.b*other.a + self.b*other.b)

    __rmul__ = __mul__

    def conjugate(self):
        return Q(self.a + self.b, -self.b)

    def norm(self):
        return self.a*self.a + self.a*self.b - self.b*self.b

    def __truediv__(self, other):
        other = Q.of(other)
        assert other.norm() != 0
        c = other.conjugate()
        return self * Q(c.a/other.norm(), c.b/other.norm())

    def __rtruediv__(self, other):
        return Q.of(other) / self

    def __eq__(self, other):
        other = Q.of(other)
        return self.a == other.a and self.b == other.b

    def sign(self):
        # Twice the value is A+B*sqrt(5). Opposite signs are decided
        # by squared rational magnitudes, never a decimal approximation.
        a, b = 2*self.a + self.b, self.b
        if not a:
            return (b > 0) - (b < 0)
        if not b or a*b > 0:
            return (a > 0) - (a < 0)
        difference = a*a - 5*b*b
        assert difference != 0
        v = a if difference > 0 else b
        return (v > 0) - (v < 0)

    def integral(self):
        return self.a.denominator == self.b.denominator == 1

    def encoded(self):
        return [str(self.a), str(self.b)]


Z, O, PHI = Q(), Q(1), Q(0, 1)


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            result.update(block)
    return result.hexdigest()


def bounded_certificate_read():
    with CERT.open('rb') as stream:
        saved = stream.read(MAX_CERT + 1)
    assert len(saved) <= MAX_CERT, 'saved certificate exceeds cap'
    return saved


def transpose(a):
    return [list(row) for row in zip(*a)]


def identity(n):
    return [[Q(int(i == j)) for j in range(n)] for i in range(n)]


def mm(a, b):
    assert len(a[0]) == len(b)
    return [[sum((x*y for x, y in zip(row, col)), Z)
             for col in zip(*b)] for row in a]


def mv(a, v):
    return [sum((x*y for x, y in zip(row, v)), Z) for row in a]


def inverse(a):
    # Only the positive definite eight-node Cartan matrix is inverted.
    n = len(a)
    rows = [a[i][:] + identity(n)[i] for i in range(n)]
    for j in range(n):
        pivot = rows[j][j]
        assert pivot != Z
        rows[j] = [x/pivot for x in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [x-factor*y for x, y in zip(rows[i], rows[j])]
    assert [row[:n] for row in rows] == identity(n)
    return [row[n:] for row in rows]


def positive_pivots(a):
    work = [row[:] for row in a]
    pivots = []
    for j in range(len(work)):
        pivot = work[j][j]
        assert pivot.sign() > 0
        pivots.append(pivot)
        for i in range(j+1, len(work)):
            for k in range(j+1, len(work)):
                work[i][k] = work[i][k] - work[i][j]*work[j][k]/pivot
    return pivots


def scalar(a, v):
    return [a*x for x in v]


def add(*vectors):
    return [sum((v[i] for v in vectors), Z) for i in range(len(vectors[0]))]


def outer(a, b):
    return [[x*y for y in b] for x in a]


def matrix_sum(*matrices):
    return [[sum((a[i][j] for a in matrices), Z)
             for j in range(len(matrices[0][0]))]
            for i in range(len(matrices[0]))]


def matrix_scale(a, matrix):
    return [scalar(a, row) for row in matrix]


def trace(matrix):
    return sum((matrix[i][i] for i in range(len(matrix))), Z)


def enc(value):
    if isinstance(value, Q):
        return value.encoded()
    if isinstance(value, F):
        return str(value)
    if isinstance(value, list):
        return [enc(v) for v in value]
    if isinstance(value, dict):
        return {k: enc(v) for k, v in value.items()}
    return value


def check(owner):
    before = {path: digest(ROOT/path) for path in PINS}
    assert before == PINS, 'pinned input changed'
    own = Path(__file__).resolve().relative_to(ROOT).as_posix()
    before[own] = digest(ROOT/own)
    assert os.getppid() == owner
    phi = PHI
    assert phi*phi == phi+1 and phi.norm() == -1
    assert phi.sign() == 1 and phi.conjugate().sign() == -1
    assert (2*phi-3-Q(F(1, 5))).sign() > 0
    assert Q(1, -2).sign() < 0 and Q(-1, 2).sign() > 0
    assert Q(2, -1).sign() > 0 and Q(-2, 1).sign() < 0
    assert (1/phi).integral() and 1/phi == phi-1

    c8 = matrix_scale(Q(2), identity(8))
    edges = [(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)]
    for i, j in edges:
        c8[i][j] = c8[j][i] = Q(-1)
    ci = inverse(c8)
    assert mm(c8, ci) == mm(ci, c8) == identity(8)
    assert all(v.integral() for row in ci for v in row)
    pivots8 = positive_pivots(c8)
    det8 = O
    for value in pivots8:
        det8 = det8*value
    assert det8 == O
    weights = {i: [ci[j][i] for j in range(8)] for i in (0, 1, 7)}
    assert weights[0] == [Q(x) for x in (4, 5, 7, 10, 8, 6, 4, 2)]
    assert weights[1] == [Q(x) for x in (5, 8, 10, 15, 12, 9, 6, 3)]
    assert weights[7] == [Q(x) for x in (2, 3, 4, 6, 5, 4, 3, 2)]

    b, delta = 5+2*phi, 2+7*phi
    assert delta == phi*b and b.norm() == 31 and delta.norm() == -31
    assert all(31 % d for d in (2, 3, 5))
    assert b.sign() > 0 and b.conjugate().sign() > 0
    assert (-delta).sign() < 0 and (-delta.conjugate()).sign() > 0
    G = [[Z for _ in range(10)] for _ in range(10)]
    G[0][0], G[0][1], G[1][0], G[1][1] = -2*phi, phi+1, phi+1, Q(2)
    for i in range(8):
        for j in range(8):
            G[i+2][j+2] = c8[i][j]
    assert G[0][0]*G[1][1]-G[0][1]*G[1][0] == -delta
    assert (Q(2)-G[0][1]*G[1][0]/G[0][0]).sign() > 0
    assert (-2*phi.conjugate()).sign() > 0
    assert (2-(phi.conjugate()+1)*(phi.conjugate()+1)
            /(-2*phi.conjugate())).sign() > 0
    Gi = [[Z for _ in range(10)] for _ in range(10)]
    Gi[0][0], Gi[1][1] = -2/delta, 2*phi/delta
    Gi[0][1] = Gi[1][0] = (phi+1)/delta
    for i in range(8):
        for j in range(8):
            Gi[i+2][j+2] = ci[i][j]
    I = identity(10)
    assert mm(Gi, G) == mm(G, Gi) == I

    def unit(i):
        return [Q(int(i == j)) for j in range(10)]

    def dot(u, v):
        return sum((x*y for x, y in zip(u, mv(G, v))), Z)

    def lift(v):
        assert len(v) == 9
        return [v[0], Z] + v[1:]

    fixture = json.loads((ROOT/'data/known_examples/bugaenko_compact_h8.json').read_text())
    decode = lambda entry: Q(F(entry['a']), F(entry['b']))
    oldG = [[decode(x) for x in row] for row in fixture['lattice_gram']]
    assert oldG == [[G[i][j] for j in [0]+list(range(2, 10))]
                   for i in [0]+list(range(2, 10))]
    roots = [lift([decode(x) for x in row]) for row in fixture['simple_roots']]
    assert len(roots) == 11
    gram = [[dot(u, v) for v in roots] for u in roots]
    assert gram == [[decode(x) for x in row] for row in fixture['wall_gram']]
    assert all(dot(v, v) == 2 for v in roots)
    for index, native in ((1, 0), (2, 2), (3, 3), (4, 4),
                          (5, 5), (6, 6), (7, 7), (9, 1)):
        assert roots[index] == scalar(Q(-1), unit(native+2))
    omega = {i: [Z, Z]+weights[i] for i in weights}
    assert roots[8] == add(unit(0), scalar(phi, omega[7]))
    assert roots[0] == add(scalar(phi, unit(0)), scalar(phi, omega[0]))

    r = add(scalar(phi, unit(0)), scalar(Q(2), unit(1)))
    c, g = scalar(Q(-1), unit(1)), roots[8]
    t = add(c, g)
    p0, p1 = unit(0), add(scalar(Q(2), unit(0)), scalar(-(phi+1), unit(1)))
    assert mv(G, r) == [Z, b]+[Z]*8
    assert dot(r, r) == 2*b and dot(c, c) == dot(g, g) == 2
    assert dot(r, c) == -b and dot(c, g) == -(phi+1)
    assert dot(t, t) == -2/phi and dot(p0, t) == -1-3*phi
    assert dot(t, c) == dot(t, g) == -1/phi
    assert dot(p0, p0) == -2*phi and dot(p1, p1) == -2*delta
    assert dot(p0, r) == dot(p1, c) == Z
    assert dot(p0, c).sign() < 0 and dot(p1, r).sign() < 0
    assert r[0]*c[1]-r[1]*c[0] == -phi
    for vector in roots+[r, c]:
        assert all(x.integral() for x in scalar(2/dot(vector, vector), mv(G, vector)))

    # Fixed complete finite-angle arithmetic seals for the inherited
    # prime-isolation proof, including non-native arrangement angles.
    angles = [O, Q(2), phi+1, Q(3), 2+phi, 2-phi, 3-phi]
    angle_seals = []
    for x in angles:
        assert x.sign() > 0 and (4-x).sign() > 0
        assert x.norm().denominator == (x-4).norm().denominator == 1
        assert int(x.norm()) % 31 and int((x-4).norm()) % 31
        angle_seals.append([x, x.norm(), (x-4).norm()])
    first_bound = 4/(31*(4-b.conjugate()))
    assert (4-b.conjugate()-Q(F(1, 5))).sign() > 0
    assert (1-first_bound).sign() > 0
    # The two binary norm forms in the proof and their conjugate Schur
    # coefficient use the same scalar b*(4-b)/4.
    assert b-b*b/4 == b*(4-b)/4
    assert (b.conjugate()*(4-b.conjugate())).sign() > 0

    p = p1
    assert dot(p, c) == Z and dot(p, g) == -delta
    assert all(dot(unit(i+2), g) == (phi if i == 7 else Z) for i in range(8))
    assert dot(omega[7], omega[7]) == 2
    e7 = [scalar(Q(-1), unit(i+2)) for i in range(7)]
    e7gram = [[dot(u, v) for v in e7] for u in e7]
    pivots7, det7 = positive_pivots(e7gram), O
    for value in pivots7:
        det7 = det7*value
    assert det7 == 2
    assert all(dot(v, t) == dot(v, c) == dot(v, g) == Z for v in e7)
    assert all(dot(v, omega[7]) == Z for v in e7)
    k = b-4*phi
    assert k == 5-2*phi and k.norm() == 11
    assert k.sign() > 0 and k.conjugate().sign() > 0
    short_norm, long_norm = F(16, 341), F(16, 11)
    assert short_norm == F(16, 1)/(b.norm()*k.norm()) < 1
    assert 1 < long_norm < 2
    assert ((phi+1)-4/k).sign() > 0 and (1-4/k.conjugate()).sign() > 0
    assert k*(phi+1) == phi+3
    assert (4/k).norm() == long_norm
    # Ordinary E7 identification uses the written golden even-E8 norm2
    # coefficient argument, not a root enumeration.

    sx = add(scalar(9+16*phi, unit(0)), scalar(6+10*phi, omega[1]))
    sy = add(unit(1), scalar(-phi/2, add(sx, scalar(Q(-1), unit(0)))))
    cols = [sx, sy, unit(9), scalar(Q(-1), roots[10]), unit(8),
            unit(7), unit(6), unit(5), unit(4), unit(2)]
    S = transpose(cols)
    P8 = [[Q(int(i == j and i >= 2)) for j in range(10)] for i in range(10)]
    ow = outer(omega[7], mv(G, omega[7]))
    ot = outer(t, mv(G, t))
    W = matrix_sum(I, matrix_scale(Q(-2), P8), ow)
    Jt = matrix_sum(matrix_scale(Q(-1), I), matrix_scale(-phi, ot))
    H = mm(Jt, W)
    assert H == matrix_sum(matrix_scale(Q(-1), I), matrix_scale(Q(2), P8),
                           matrix_scale(Q(-1), ow), matrix_scale(-phi, ot))
    SH, HS = mm(S, H), mm(H, S)
    for matrix in (S, W, Jt, H, SH):
        assert all(x.integral() for row in matrix for x in row)
        assert mm(mm(transpose(matrix), G), matrix) == G
        adjoint = mm(mm(Gi, transpose(matrix)), G)
        assert mm(adjoint, matrix) == mm(matrix, adjoint) == I
        assert all(x.integral() for row in adjoint for x in row)
    for matrix in (S, W, Jt, H):
        assert mm(matrix, matrix) == I
    assert mm(SH, HS) == mm(HS, SH) == I
    assert mv(S, r) == r
    for i in range(11):
        target = 8-i if i <= 8 else 19-i
        assert mv(S, roots[i]) == roots[target]
    assert mv(W, t) == mv(Jt, t) == mv(H, t) == t
    assert all(mv(W, v) == scalar(Q(-1), v) and mv(H, v) == v for v in e7)
    assert dot(p0, mv(S, p0)).sign() < 0
    assert dot(p0, t).sign() < 0 and dot(p0, mv(SH, p0)).sign() < 0
    st = add(scalar(8+13*phi, unit(0)), scalar(Q(-1), unit(1)),
             scalar(5+8*phi, omega[1]), scalar(phi, omega[0]))
    assert mv(S, t) == st
    assert dot(omega[7], omega[1]) == 3 and dot(omega[7], omega[0]) == 2
    assert trace(S) == 2 and trace(mm(S, P8)) == -8-16*phi
    assert dot(omega[7], mv(S, omega[7])) == -10-18*phi
    assert dot(t, st) == -20-36*phi
    assert trace(H) == 6
    assert trace(SH) == (-trace(S)+2*trace(mm(S, P8))
                         -dot(omega[7], mv(S, omega[7]))-phi*dot(t, st))
    assert trace(SH) == 28+42*phi and (trace(SH)-10).sign() > 0

    assert os.getppid() == owner
    assert {path: digest(ROOT/path) for path in before} == before, 'source changed during check'
    return enc({
        'status': 'PASS_FIXED_GOLDEN31_UNIT_BISECTOR',
        'encoding': '[a,b] means a+b*phi; all entries are exact rational strings',
        'scope': 'literal golden31 only; silver31 and general binary family are distinct',
        'written_premises': ['golden PID and unit classification; primitive norm-ideal argument',
                             'complete finite-angle list and prime isolation geometry',
                             'actual full Bugaenko facet; actual edge and E7 germ proof',
                             'local chamber-germ theorem; coset and arithmetic covolume bridge'],
        'source_sha256': before,
        'G': G, 'E8_edges': [list(edge) for edge in edges], 'E8_inverse': ci,
        'E8_pivots': pivots8, 'E7_gram': e7gram, 'E7_pivots': pivots7,
        'fixture_wall_gram': gram, 'b': b, 'delta': delta,
        'norm_pool': [Q(2), 2*b], 'finite_angle_prime_seals': angle_seals,
        'r': r, 'c': c, 'g': g, 't': t, 'p0': p0, 'p1': p1,
        'first_edge_norm_bound': first_bound,
        'germ_short_norm_bound': short_norm, 'germ_long_norm_bound': long_norm,
        'germ_long_physical_square_bound': 4/k,
        'germ_long_conjugate_square_bound': 4/k.conjugate(),
        'S': S, 'W': W, 'H': H, 'SH': SH,
        'trace_S': trace(S), 'trace_SP8': trace(mm(S, P8)),
        'omega7_Somega7': dot(omega[7], mv(S, omega[7])),
        't_St': dot(t, st), 'trace_SH': trace(SH),
        'completed_fixed_checks': True,
        'no_height_pell_or_cone_search': True,
    })


def main():
    if not __debug__:
        raise RuntimeError('assertions must be enabled')
    assert __debug__, 'assertions must be enabled'
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', action='store_true')
    actions.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    assert sys.dont_write_bytecode, 'invoke with -B'
    owner = int(os.environ.get('COXETER_CHECK_OWNER', '0'))
    assert owner > 0 and os.getppid() == owner, 'bounded parent owner required'
    signal.alarm(45)
    resource.setrlimit(resource.RLIMIT_CPU, (45, 45))
    saved = None
    if args.write:
        assert not CERT.exists(), 'refusing certificate overwrite'
    else:
        assert CERT.is_file()
        saved = bounded_certificate_read()
    result = check(owner)
    payload = (json.dumps(result, indent=2, sort_keys=True)+'\n').encode()
    assert len(payload) <= MAX_CERT
    assert os.getppid() == owner
    if args.write:
        with CERT.open('xb') as stream:
            stream.write(payload)
    else:
        assert saved == payload, 'fresh replay differs from certificate snapshot'
        assert bounded_certificate_read() == saved, 'saved certificate changed during replay'
    signal.alarm(0)
    print(json.dumps({'status': result['status'], 'mode': 'write' if args.write else 'verify',
                      'certificate_bytes': len(payload),
                      'certificate_sha256': hashlib.sha256(payload).hexdigest()}))


if __name__ == '__main__':
    main()
