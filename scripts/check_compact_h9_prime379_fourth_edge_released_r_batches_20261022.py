#!/usr/bin/env python3
"""PREPARED, UNREVIEWED, UNEXECUTED: norm379 fourth edge, releasing r.

Standalone stdlib source; helpers copied from frozen56d8294b, never imported.
Five profiles x three phases, then separately authenticated combination.
Nine long selectors are empty by the pinned reviewed proof; six short
selectors stream at most324 scalar pairs. WHOLE original HALF seeds,
all15 original known walls, exact positive self recognition, every tie
and a unique full nine-wall cone are mandatory. No outcome is presumed.
Root alone owns eventual reviewed, guarded execution.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SELF = Path(__file__).resolve()
SPEC = "results/compact_h9_prime379_fourth_edge_released_r_checker_spec_20261022.md"
PINS = {
    'results/compact_h9_prime379_fourth_edge_released_r_checker_spec_20261022.md':
        'c107eb4b71c51bd3f6ded30fb5d38088a44e36bcea5487ffb19680441f875dc2',
    'results/compact_h9_prime379_fourth_edge_integral_long_stabilizer_and_half_folded_domain_candidate_20261022.md':
        'df2bc3a297bd2fa1322914fe4fac532f5aa7abbe36c79a39f689df2fd7cfa336',
    'results/review_prime379_fourth_edge_integral_long_stabilizer_half_domain_root_acceptance_20261022.md':
        'ccf44144471e7db2e330e1d532be3c5b17ab760cf76ff521c8befa836225364b',
    'results/review_prime379_fourth_edge_integral_long_stabilizer_half_domain_full_scope_compact_reviewer_20261022.md':
        '8ff25dad59574e5a4f60f6621e384c0166d5fdec91f3e112ddb1f380dfe12636',
    'results/compact_h9_prime379_third_edge_and_e7_two_a1_germ_candidate_20261022.md':
        'd52b2c15ef72ddff43c12a30fd5424c5ce5aa0ed9677a005f05ec3597fd7968d',
    'results/review_prime379_third_edge_and_e7_two_a1_germ_root_acceptance_20261022.md':
        '870bc701ced864f1e0bf6bbf1e97e91380bbc7436b693c0d0196926eee85fc42',
    'results/prime379_third_released_m_runtime_manifest_20261022.json':
        'ec8eb9bb4558703db8a2d9253bd7c6c425f71c2d0820426e6e565494593b2824',
    'results/compact_h9_prime379_first_binary_edge_and_quarter_folded_transverse_domain_candidate_20261022.md':
        '7e7ffd616d0e025aeaff84abae1c9fcfb7d930c8a3879a7bb43734f308cacd3a',
    'results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md':
        'b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce',
    'data/known_examples/bugaenko_compact_h8.json':
        'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
    'results/prime379_third_released_m_combine_20261022.json':
        '6d7627bbf2c2b7224aaa539eff2678da30ed061d5890e01a8a0f0d1cdee4917e',
    'scripts/check_compact_h9_prime379_third_edge_released_m_batches_20261022.py':
        '56d8294bc78fa037de7413ac076799bf572da09659bac1ae5ca86301abe2bd7d',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}
PROFILES = ("norm2_orthogonal", "norm2_e8", "norm2b_orthogonal",
            "norm2b_simple", "norm2b_golden")
PHASES = tuple(range(3))
TIE_LIMIT = 128
FILE_LIMIT = 2*1024**2
COUNT_KEYS = ("scalar_rows", "inside_bounds", "folded_window", "field_roots",
              "half_seed", "original_integral", "admitted")


class F:
    """Exact a+b*phi; phi^2=phi+1, with rational coefficients."""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a, self.b = Fraction(a), Fraction(b)

    @staticmethod
    def of(value):
        return value if isinstance(value, F) else F(value)

    def __add__(self, other):
        q = F.of(other)
        return F(self.a + q.a, self.b + q.b)

    __radd__ = __add__

    def __neg__(self):
        return F(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-F.of(other))

    def __rsub__(self, other):
        return F.of(other) + (-self)

    def __mul__(self, other):
        q = F.of(other)
        return F(self.a*q.a + self.b*q.b,
                 self.a*q.b + self.b*q.a + self.b*q.b)

    __rmul__ = __mul__

    def conjugate(self):
        return F(self.a + self.b, -self.b)

    def norm(self):
        return self.a*self.a + self.a*self.b - self.b*self.b

    def __truediv__(self, other):
        q = F.of(other)
        n = q.norm()
        if not n:
            raise ZeroDivisionError("zero field denominator")
        z = self*q.conjugate()
        return F(z.a/n, z.b/n)

    def __eq__(self, other):
        if not isinstance(other, (F, int, Fraction)):
            return NotImplemented
        q = F.of(other)
        return self.a == q.a and self.b == q.b

    def __hash__(self):
        return hash((self.a, self.b)) if self.b else hash(self.a)

    def sign(self):
        a, b = self.a + self.b/2, self.b/2  # a+b*sqrt(5)
        if not a:
            return (b > 0) - (b < 0)
        if not b or (a > 0) == (b > 0):
            return (a > 0) - (a < 0)
        n = a*a - 5*b*b
        s = (n > 0) - (n < 0)
        return s if a > 0 else -s

    def integral(self):
        return self.a.denominator == self.b.denominator == 1

    def half(self):
        return (2*self).integral()


ZERO, ONE, PHI = F(), F(1), F(0, 1)
D, B = ONE + 20*PHI, F(20, 21)
H = tuple(F(x) for x in (2, 3, 4, 6, 5, 4, 3, 2))
J = tuple(F(x) for x in (4, 6, 8, 12, 10, 8, 6, 3))


def fenc(q):
    return [q.a.numerator, q.a.denominator, q.b.numerator, q.b.denominator]


def fdec(row):
    if not isinstance(row, list) or len(row) != 4 or any(type(x) is not int for x in row):
        raise ValueError("invalid field encoding")
    if row[1] <= 0 or row[3] <= 0:
        raise ValueError("invalid rational denominator")
    return F(Fraction(row[0], row[1]), Fraction(row[2], row[3]))


def venc(v):
    if not all(q.integral() for q in v):
        raise ValueError("nonintegral original vector")
    return [[int(q.a), int(q.b)] for q in v]


def vdec(rows):
    if not isinstance(rows, list) or len(rows) != 10:
        raise ValueError("invalid original vector dimension")
    if any(not isinstance(r, list) or len(r) != 2 or
           any(type(x) is not int for x in r) for r in rows):
        raise ValueError("invalid original vector encoding")
    return tuple(F(*r) for r in rows)


def vkey(v):
    return tuple((q.a, q.b) for q in v)


def rational_sqrt(q):
    if q < 0:
        return None
    n, d = isqrt(q.numerator), isqrt(q.denominator)
    return Fraction(n, d) if n*n == q.numerator and d*d == q.denominator else None


def field_roots(z):
    """All exact Q(phi) square roots, including pure sqrt(5) and zero."""
    aa, bb = z.a + z.b/2, z.b/2
    roots = set()
    if not bb:
        x = rational_sqrt(aa)
        y = rational_sqrt(aa/5)
        if x is not None:
            roots.update((F(x), F(-x)))
        if y is not None:
            roots.update((F(-y, 2*y), F(y, -2*y)))
    else:
        n = rational_sqrt(aa*aa - 5*bb*bb)
        if n is None:
            return ()
        for nn in (n, -n):
            x = rational_sqrt((aa + nn)/2)
            if x is None or not x:
                continue
            y = bb/(2*x)
            q = F(x-y, 2*y)
            if q*q == z:
                roots.update((q, -q))
    if any(q*q != z for q in roots):
        raise AssertionError("field-square reconstruction error")
    return tuple(sorted(roots, key=lambda q: (q.a, q.b)))


def floor_field(q):
    bounds = (q.a + q.b*Fraction(8, 5), q.a + q.b*Fraction(81, 50))
    low, high = min(bounds).__floor__(), max(bounds).__floor__() + 1
    while high - low > 1:
        mid = (low + high)//2
        if (q-mid).sign() >= 0:
            low = mid
        else:
            high = mid
    return low


def ceil_field(q):
    return -floor_field(-q)


def dot(v, w, gram):
    return sum((v[i]*gram[i][j]*w[j]
                for i in range(len(v)) for j in range(len(v))
                if v[i] != ZERO and gram[i][j] != ZERO and w[j] != ZERO), ZERO)


def mv(matrix, vector):
    return tuple(sum((a*b for a, b in zip(row, vector)), ZERO) for row in matrix)


def mm(left, right):
    return tuple(tuple(sum((left[i][k]*right[k][j] for k in range(len(right))), ZERO)
                       for j in range(len(right[0]))) for i in range(len(left)))


def identity(n):
    return tuple(tuple(F(i == j) for j in range(n)) for i in range(n))


def reflection(v, gram):
    norm = dot(v, v, gram)
    cov = tuple(sum((gram[i][j]*v[j] for j in range(len(v))), ZERO)
                for i in range(len(v)))
    return tuple(tuple(F(i == j)-2*v[i]*cov[j]/norm for j in range(len(v)))
                 for i in range(len(v)))


def add(v, w, factor=ONE):
    return tuple(a + factor*b for a, b in zip(v, w))


def scale(v, q):
    return tuple(q*a for a in v)


def original(v3, epsilon):
    h = v3[2] - Fraction(3, 2)*epsilon
    return tuple(v3[:2]) + tuple(h*x + epsilon*y for x, y in zip(H, J))


def full_basis_images(matrix, complement=ONE):
    """Full original images; active e7 projection is H/2, not H."""
    native = [tuple(F(k == j+2) for k in range(10)) for j in range(8)]
    hh = (ZERO, ZERO, ONE)
    correction = original(scale(add(mv(matrix, hh), hh, -complement),
                                F(Fraction(1, 2))), 0)
    return [original(mv(matrix, (ONE, ZERO, ZERO)), 0),
            original(mv(matrix, (ZERO, ONE, ZERO)), 0)] + \
           [scale(v, complement) for v in native[:7]] + \
           [add(scale(native[7], complement), correction)]


def read_json_bounded(path):
    with path.open("rb") as stream:
        raw = stream.read(FILE_LIMIT+1)
    if len(raw) > FILE_LIMIT:
        raise ValueError("bounded JSON input exceeded2MiB")
    return json.loads(raw)


def geometry():
    """Direct accepted long-wall skew boost on the WHOLE original module."""
    g3 = ((-2*PHI, ONE, ZERO), (ONE, F(10), ZERO), (ZERO, ZERO, F(2)))
    gram = [[ZERO for _ in range(10)] for _ in range(10)]
    gram[0][:2], gram[1][:2] = [-2*PHI, ONE], [ONE, F(10)]
    for i in range(8):
        gram[i+2][i+2] = F(2)
    for i, j in ((0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)):
        gram[i+2][j+2] = gram[j+2][i+2] = F(-1)
    gram = tuple(tuple(row) for row in gram)
    native = [tuple(F(k == j+2) for k in range(10)) for j in range(8)]
    hh, xx = (ZERO, ZERO, ONE), (ONE, ZERO, ZERO)
    assert B == PHI*D and B.norm() == 379
    assert D.sign() > 0 and D.conjugate().sign() < 0
    assert B.sign() > 0 and B.conjugate().sign() > 0
    assert dot((ZERO, ZERO)+H, (ZERO, ZERO)+H, gram) == 2
    assert dot((ZERO, ZERO)+J, (ZERO, ZERO)+J, gram) == 6
    assert dot((ZERO, ZERO)+H, (ZERO, ZERO)+J, gram) == 3
    for i, v in enumerate(native):
        assert dot((ZERO, ZERO)+H, v, gram) == F(i == 7)
        assert dot((ZERO, ZERO)+J, v, gram) == F(i == 6)

    calibration = read_json_bounded(ROOT/"data/known_examples/bugaenko_compact_h8.json")
    def number(row):
        return F(Fraction(row["a"]), Fraction(row["b"]))
    calgram = tuple(tuple(number(x) for x in row) for row in calibration["lattice_gram"])
    assert len(calgram) == 9 and all(len(row) == 9 for row in calgram)
    indexes = (0,)+tuple(range(2, 10))
    assert all(calgram[i][j] == gram[indexes[i]][indexes[j]]
               for i in range(9) for j in range(9))
    base = []
    for row in calibration["simple_roots"]:
        assert len(row) == 9
        base.append((number(row[0]), ZERO)+tuple(number(x) for x in row[1:]))
    assert len(base) == 11 and len(set(base)) == 11
    assert all(integral_reflector(wall, ONE, gram) for wall in base)
    assert all(scale(v, -1) in base for v in native)

    old_q = (F(18, 25), F(-12, -17), ZERO)
    m3 = (F(2, 5), F(-2, -3), ZERO)
    r3 = add(add(scale(old_q, F(2572, Fraction(8323, 2))), m3, -PHI/2),
             hh, F(14856, 24036))
    w3 = (F(28490145, 46098027), F(-19256316, -31157373),
          F(2815372, 4555369))
    r0 = (ONE, 2*PHI)+(ZERO,)*8
    m, r, w = original(m3, 0), original(r3, 0), original(w3, 0)
    assert dot(old_q, old_q, g3) == -2*D
    assert dot(m3, m3, g3) == dot(r3, r3, g3) == 2
    assert dot(w3, w3, g3) == 2*B and dot(r3, w3, g3) == ZERO
    assert dot(m3, r3, g3) == -PHI
    assert all(integral_reflector(v, d, gram)
               for v, d in ((r0, B), (m, ONE), (r, ONE), (w, B)))
    assert all((c/B).integral() for c in mv(gram, w))
    known = base+[r0, m, r, w]
    assert len(set(known)) == 15

    old = read_json_bounded(ROOT/"results/prime379_third_released_m_combine_20261022.json")
    assert old["kind"] == "prime379_released_m_complete_domain"
    assert old["status"] == "PASS_COMPLETE_DOMAIN"
    assert old["profile_count"] == 5 and old["selector_count"] == 55
    assert old["germ"]["certified"] is True and old["germ"]["cone_candidates"] == 1
    assert old["known_wall_count"] == 14
    assert fdec(old["collision"]) == F(Fraction(1154937, 745919), Fraction(1948905, 745919))
    assert vdec(old["germ"]["new_simple_wall"]) == w
    assert len(old["minimizers"]) == 1 and vdec(old["minimizers"][0]["root"]) == w
    old_basis = [scale(v, -1) for v in native[:7]]+[r, w]
    old_gram = tuple(tuple(dot(x, y, gram) for y in old_basis) for x in old_basis)
    assert tuple(tuple(fdec(x) for x in row) for row in old["germ"]["gram"]) == old_gram
    old_point = tuple(fdec(x) for x in old["germ"]["endpoint_original_ray"])
    assert len(old_point) == 10 and dot(old_point, old_point, gram).sign() < 0
    assert dot(old_point, (ONE, ZERO)+(ZERO,)*8, gram).sign() < 0
    assert all(dot(old_point, wall, gram).sign() <= 0 for wall in known)
    assert {vkey(v) for v in known if dot(old_point, v, gram) == ZERO} == \
        {vkey(v) for v in old_basis}

    # Active B=E-P_w; the extension is ZERO on E7. Old e7 projects to H/2.
    basis3 = identity(3)
    pwcols = [scale(w3, dot(v, w3, g3)/(2*B)) for v in basis3]
    projector = tuple(tuple(F(i == j)-pwcols[j][i] for j in range(3))
                      for i in range(3))
    a, y, h = w3
    skew = ((h/D, 10*h/D, (-10*y-a)/D),
            (2*PHI*h/D, -h/D, (y-2*PHI*a)/D),
            (-y/2, a/2, ZERO))
    assert h/D == F(84392, 136549)
    assert (-10*y-a)/D == F(4918155, 7957743)
    assert (y-2*PHI*a)/D == F(-3340830, -5405577)
    assert mm(projector, projector) == projector
    assert mv(projector, w3) == (ZERO, ZERO, ZERO)
    assert mv(projector, r3) == r3
    assert mv(skew, w3) == (ZERO, ZERO, ZERO)
    assert mm(skew, skew) == tuple(tuple(PHI*x for x in row) for row in projector)
    assert mm(skew, projector) == skew == mm(projector, skew)
    gs = mm(g3, skew)
    assert all(gs[i][j] == -gs[j][i] for i in range(3) for j in range(3))
    for operator in (projector, skew):
        assert all((2*x).integral() for image in full_basis_images(operator, ZERO)
                   for x in image)
    n3 = mv(skew, r3)
    if dot(n3, xx, g3).sign() > 0:
        skew = tuple(tuple(-x for x in row) for row in skew)
        n3 = mv(skew, r3)
    assert dot(n3, xx, g3).sign() < 0 and dot(n3, n3, g3) == -2*PHI
    assert dot(n3, r3, g3) == ZERO == dot(n3, w3, g3)
    assert mv(skew, n3) == scale(r3, PHI)
    assert positive_proportional(original(n3, 0), old_point)

    matrix = tuple(tuple(F(i == j)+(PHI-1)*projector[i][j]-skew[i][j]
                         for j in range(3)) for i in range(3))
    inverse = tuple(tuple(F(i == j)+(PHI-1)*projector[i][j]+skew[i][j]
                          for j in range(3)) for i in range(3))
    assert mm(matrix, inverse) == identity(3) == mm(inverse, matrix)
    assert mv(matrix, w3) == w3 == mv(inverse, w3)
    assert mv(matrix, n3) == add(scale(n3, PHI), r3, -PHI)
    assert mv(matrix, r3) == add(scale(r3, PHI), n3, F(-1))
    assert dot(n3, mv(matrix, n3), g3) == -2*PHI*PHI
    powers, inverse_powers = [identity(3)], [identity(3)]
    for _ in range(2):
        powers.append(mm(matrix, powers[-1]))
        inverse_powers.append(mm(inverse, inverse_powers[-1]))
    for family in (powers, inverse_powers):
        for j in PHASES:
            images = full_basis_images(family[j])
            assert all(x.half() for image in images for x in image)
            assert all(dot(images[i], images[k], gram) == gram[i][k]
                       for i in range(10) for k in range(10))
            assert all(images[i+2] == native[i] for i in range(7))
            if j == 2:
                assert all(x.integral() for image in images for x in image)
            assert mm(family[j], (inverse_powers if family is powers else powers)[j]) == identity(3)
    for sign, operator in ((-1, powers[2]), (1, inverse_powers[2])):
        expected = tuple(tuple(F(i == j)+2*PHI*projector[i][j]+sign*2*PHI*skew[i][j]
                               for j in range(3)) for i in range(3))
        assert operator == expected

    line_q, line_h = original(n3, 0), scale(r, F(-1))
    retained = [scale(v, -1) for v in native[:7]]+[w]
    assert all(dot(line_q, v, gram) == ZERO == dot(line_h, v, gram) for v in retained)
    assert dot(line_h, r, gram) == -2
    assert all(dot(line_q, v, gram).sign() <= 0 for v in known)
    assert {vkey(v) for v in known if dot(line_q, v, gram) == ZERO} == \
        {vkey(v) for v in old_basis}
    upper = 2*PHI*PHI/(1+2*PHI)
    upper3 = scale(mv(powers[2], r3), F(-1))
    upper_root = original(upper3, 0)
    assert upper3 == add(scale(n3, 2*PHI), r3, -(1+2*PHI))
    assert integral_reflector(upper_root, ONE, gram)
    assert upper.sign() > 0 and (PHI-upper*upper).sign() > 0
    assert dot(upper_root, line_q, gram) == -4*PHI*PHI
    assert dot(upper_root, line_h, gram) == 2*(1+2*PHI)
    assert dot(upper_root, add(line_q, line_h, upper), gram) == ZERO
    assert positive_proportional(original(mv(powers[2], n3), 0),
                                 add(line_q, line_h, upper))

    # Rational seals for the reviewed empty-long proof, NOT an observed census.
    assert (F(Fraction(13, 8))-PHI).sign() > 0
    assert (B-F(Fraction(268, 5))).sign() > 0
    assert (B.conjugate()-F(Fraction(55, 8))).sign() > 0
    assert Fraction(845, 34304) < Fraction(16, 625)
    assert Fraction(8, 55) < Fraction(1521, 10000)
    assert Fraction(11, 5)**2 < 5 and Fraction(16, 25) < 1
    assert Fraction(43, 5)**2 < 16*5 and Fraction(63, 10)**2 < 9*5
    return {"g3": g3, "gram": gram, "q": n3, "r3": r3, "w3": w3,
            "A": PHI, "powers": powers, "inverse_powers": inverse_powers,
            "line_q": line_q, "line_h": line_h, "start": ZERO,
            "start_ray": line_q, "starting_wall": r, "upper": upper,
            "upper_root": upper_root, "known": known, "retained": retained}


def record(path):
    digest, size = hashlib.sha256(), 0
    with path.open("rb") as stream:
        while chunk := stream.read(1024*1024):
            digest.update(chunk)
            size += len(chunk)
    return {"path": str(path.relative_to(ROOT)), "bytes": size, "sha256": digest.hexdigest()}


def bindings():
    rows = {}
    for name, expected in PINS.items():
        row = record(ROOT/name)
        if row["sha256"] != expected:
            raise RuntimeError("changed pinned input: "+name)
        rows[name] = row
    rows[str(SELF.relative_to(ROOT))] = record(SELF)
    return rows


def output_path(raw):
    path = (ROOT/raw).resolve()
    if path.parent != ROOT/"results" or path.suffix != ".json" or path.exists():
        raise ValueError("output must be a NEW direct results/*.json destination")
    return path


def write_new(path, data):
    text = json.dumps(data, indent=2, sort_keys=True)+"\n"
    if len(text.encode("utf-8")) > FILE_LIMIT:
        raise RuntimeError("output storage bound exceeded; no truncation")
    with path.open("x", encoding="utf-8") as stream:
        stream.write(text)


def profile(name):
    return {"norm2_orthogonal": (ONE, ZERO, 0),
            "norm2_e8": (ONE, ZERO, 1),
            "norm2b_orthogonal": (B, ZERO, 0),
            "norm2b_simple": (B, B, 0),
            "norm2b_golden": (B, B*PHI, 0)}[name]


def scalar_rows(name):
    if name not in ("norm2_orthogonal", "norm2_e8"):
        raise ValueError("long selectors have a proved empty domain, not a scalar loop")
    limit, radius = (3 if name == "norm2_orthogonal" else 2), F(4)
    for e in range(-limit, limit+1):
        center = e*(PHI-1)
        for c in range(ceil_field(center-radius), floor_field(center+radius)+1):
            yield F(c, e)/4


def inside_scalar_bounds(name, scalar):
    if name not in ("norm2_orthogonal", "norm2_e8"):
        return False
    bound = F(Fraction(23, 20) if name == "norm2_orthogonal" else Fraction(23, 40))
    return (bound-scalar).sign() > 0 and (bound+scalar).sign() > 0 and \
        (ONE-scalar.conjugate()).sign() >= 0 and \
        (ONE+scalar.conjugate()).sign() >= 0


def integral_reflector(v, d, gram):
    if len(v) != 10 or not all(x.integral() for x in v):
        return False
    if d not in (ONE, B) or dot(v, v, gram) != 2*d:
        return False
    # Evenness and prime norm379 make these normalized vectors primitive:
    # a squared nonunit divisor cannot divide1 or b. Check FULL reflection.
    cov = mv(gram, v)
    return all((x*y/d).integral() for x in v if x != ZERO
               for y in cov if y != ZERO)


def positive_proportional(vector, wall):
    """Only full TEN-coordinate equality and positive ratio recognize self."""
    if len(vector) != 10 or len(wall) != 10:
        raise ValueError("self recognition requires all ten original coordinates")
    pivot = next((i for i, x in enumerate(wall) if x != ZERO), None)
    if pivot is None:
        raise ValueError("zero known wall")
    ratio = vector[pivot]/wall[pivot]
    return ratio.sign() > 0 and all(x == ratio*y for x, y in zip(vector, wall))


def known_self_exemptions(root, geo):
    return [i for i, wall in enumerate(geo["known"]) if positive_proportional(root, wall)]


def validate_original(v, d, geo):
    if not integral_reflector(v, d, geo["gram"]):
        return None
    numerator = -dot(v, geo["line_q"], geo["gram"])
    derivative = dot(v, geo["line_h"], geo["gram"])
    if derivative.sign() <= 0:
        return None
    tau = numerator/derivative
    if (tau-geo["start"]).sign() <= 0 or (geo["upper"]-tau).sign() < 0:
        return None
    if dot(v, geo["start_ray"], geo["gram"]).sign() >= 0:
        return None
    for wall in geo["known"]:
        if positive_proportional(v, wall):
            continue
        if dot(v, wall, geo["gram"]).sign() > 0:
            return None
    return tau


def validate_part_root(root, phase, name, geo):
    d, kappa, epsilon = profile(name)
    if d == B:
        return None  # Each long selector is empty by the pinned complete proof.
    gram = geo["gram"]
    native6 = tuple(F(k == 8) for k in range(10))
    native7 = tuple(F(k == 9) for k in range(10))
    h = dot(root, native7, gram)
    if dot(root, native6, gram) != epsilon:
        return None
    expected = tuple(root[:2])+tuple(h*x+epsilon*y for x, y in zip(H, J))
    if expected != root:
        return None
    if dot(root, original(geo["w3"], 0), gram) != -kappa:
        return None
    active = tuple(root[:2])+(h+Fraction(3, 2)*epsilon,)
    seed = mv(geo["inverse_powers"][phase], active)
    if not all(x.half() for x in original(seed, epsilon)):
        return None
    z = dot(seed, geo["r3"], geo["g3"])/2
    if not (4*z).integral() or not inside_scalar_bounds(name, z):
        return None
    if abs((4*z).b) > (3 if name == "norm2_orthogonal" else 2):
        return None
    kk = 2*d-kappa*kappa/(2*B)-F(Fraction(3, 2))*epsilon*epsilon
    if (kk*(geo["A"]+1)/4-z*z).sign() < 0:
        return None
    if dot(seed, geo["w3"], geo["g3"]) != -kappa:
        return None
    time = -dot(seed, geo["q"], geo["g3"])/(2*PHI)
    rebuilt = add(add(scale(geo["q"], time), geo["w3"], -kappa/(2*B)),
                  geo["r3"], z)
    if rebuilt != seed:
        return None
    if dot(seed, seed, geo["g3"])+Fraction(3, 2)*epsilon*epsilon != 2*d:
        return None
    if time*time != (2*z*z+kappa*kappa/(2*B)+Fraction(3, 2)*epsilon*epsilon-2*d)/(2*PHI):
        return None
    return validate_original(root, d, geo)


def keep_minimum(minimum, ties, tau, root, phase, name):
    if minimum is None or (tau-minimum).sign() < 0:
        minimum, ties = tau, {}
    if tau == minimum:
        key = vkey(root)
        if key not in ties:
            ties[key] = (root, phase, name)
        if len(ties) > TIE_LIMIT:
            raise RuntimeError("complete minimum tie storage bound exceeded; no truncation")
    return minimum, ties


def encode_ties(ties, geo):
    return [{"root": venc(ties[key][0]), "phase": ties[key][1], "profile": ties[key][2],
             "known_self_exemptions": known_self_exemptions(ties[key][0], geo)}
            for key in sorted(ties)]


def unit(name, geo):
    d, kappa, epsilon = profile(name)
    kk = 2*d-kappa*kappa/(2*B)-F(Fraction(3, 2))*epsilon*epsilon
    assert kk.sign() > 0
    minimum, ties, selectors = None, {}, []
    for phase in PHASES:
        counts = {key: 0 for key in COUNT_KEYS}
        if d == B:
            # The reviewed quarter-reflection argument gives 4z in bO.
            # Its physical/conjugate seals force z=0, then a^2<0.
            # No observed absence, integral folded reflection or omitted loop.
            assert (-kk/(2*PHI)).sign() < 0
            selectors.append({"phase": phase, "profile": name, "counts": counts,
                              "domain_status": "PROVED_EMPTY_LONG"})
            continue
        for z in scalar_rows(name):
            counts["scalar_rows"] += 1
            if not inside_scalar_bounds(name, z):
                continue
            counts["inside_bounds"] += 1
            if (kk*(geo["A"]+1)/4-z*z).sign() < 0:
                continue
            counts["folded_window"] += 1
            square = (2*z*z+kappa*kappa/(2*B)+Fraction(3, 2)*epsilon*epsilon-2*d)/(2*PHI)
            for time in field_roots(square):
                counts["field_roots"] += 1
                seed = add(add(scale(geo["q"], time), geo["w3"], -kappa/(2*B)),
                           geo["r3"], z)
                if not all(x.half() for x in original(seed, epsilon)):
                    continue
                counts["half_seed"] += 1
                assert dot(seed, seed, geo["g3"])+Fraction(3, 2)*epsilon*epsilon == 2*d
                assert dot(seed, geo["w3"], geo["g3"]) == -kappa
                root = original(mv(geo["powers"][phase], seed), epsilon)
                if not all(x.integral() for x in root):
                    continue
                counts["original_integral"] += 1
                tau = validate_original(root, d, geo)
                if tau is None:
                    continue
                counts["admitted"] += 1
                minimum, ties = keep_minimum(minimum, ties, tau, root, phase, name)
        cap = 63 if name == "norm2_orthogonal" else 45
        assert counts["scalar_rows"] <= cap
        assert counts["field_roots"] <= 2*counts["folded_window"]
        selectors.append({"phase": phase, "profile": name, "counts": counts,
                          "domain_status": "FINITE_SCALAR_BOX"})
    return {"kind": "prime379_released_r_profile_unit", "status": "PASS_DECLARED_DOMAIN",
            "profile": name, "selectors": selectors,
            "starting_collision": fenc(geo["start"]),
            "upper_collision": fenc(geo["upper"]), "known_wall_count": len(geo["known"]),
            "collision": None if minimum is None else fenc(minimum),
            "minimizers": encode_ties(ties, geo), "germ_certified": False,
            "scope": "Fourth-edge ONE profile/three phases only; long emptiness uses pinned proof; no germ claim"}


def positive_definite(matrix):
    a = [list(row) for row in matrix]
    for k in range(len(a)):
        if a[k][k].sign() <= 0:
            return False
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                a[i][j] = a[i][j]-a[i][k]*a[k][j]/a[k][k]
    return True


def solve(matrix, rhs):
    n = len(matrix)
    a = [list(row)+[b] for row, b in zip(matrix, rhs)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if a[i][k] != ZERO), None)
        if pivot is None:
            return None
        a[k], a[pivot] = a[pivot], a[k]
        p = a[k][k]
        a[k] = [q/p for q in a[k]]
        for i in range(n):
            if i != k and a[i][k] != ZERO:
                c = a[i][k]
                a[i] = [x-c*y for x, y in zip(a[i], a[k])]
    return tuple(row[-1] for row in a)


def germ(minimum, ties, geo):
    point = add(geo["line_q"], geo["line_h"], minimum)
    gram = geo["gram"]
    assert dot(point, point, gram).sign() < 0
    assert dot(point, geo["start_ray"], gram).sign() < 0
    assert (minimum-geo["start"]).sign() > 0
    assert (geo["upper"]-minimum).sign() >= 0
    if any(dot(point, wall, gram).sign() > 0 for wall in geo["known"]):
        raise RuntimeError("complete minimum violates a known actual wall")
    chosen = []
    for root, _, _ in ties.values():
        basis = geo["retained"]+[root]
        gg = tuple(tuple(dot(x, y, gram) for y in basis) for x in basis)
        if not positive_definite(gg) or any(dot(point, wall, gram) != ZERO for wall in basis):
            continue
        all_inside = True
        for other, _, _ in ties.values():
            coeff = solve(gg, tuple(dot(wall, other, gram) for wall in basis))
            if coeff is None or any(c.sign() < 0 for c in coeff):
                all_inside = False
                break
            rebuilt = tuple(sum((coeff[j]*basis[j][i] for j in range(9)), ZERO)
                            for i in range(10))
            if rebuilt != other:
                all_inside = False
                break
        if all_inside:
            chosen.append((root, gg))
    unique = len(chosen) == 1
    return {"certified": unique, "cone_candidates": len(chosen),
            "new_simple_wall": venc(chosen[0][0]) if unique else None,
            "gram": [[fenc(x) for x in row] for row in chosen[0][1]] if unique else None,
            "endpoint_original_ray": [fenc(x) for x in point],
            "scope": "Actual full finite germ ONLY with unique cone and all15 authenticated selectors"}


def valid_counts(counts, name):
    if not isinstance(counts, dict) or set(counts) != set(COUNT_KEYS):
        return False
    if any(type(x) is not int or x < 0 for x in counts.values()):
        return False
    if profile(name)[0] == B:
        return all(x == 0 for x in counts.values())
    cap = 63 if name == "norm2_orthogonal" else 45
    return 0 < counts["scalar_rows"] <= cap and \
        counts["inside_bounds"] <= counts["scalar_rows"] and \
        counts["folded_window"] <= counts["inside_bounds"] and \
        counts["field_roots"] <= 2*counts["folded_window"] and \
        counts["half_seed"] <= counts["field_roots"] and \
        counts["original_integral"] <= counts["half_seed"] and \
        counts["admitted"] <= counts["original_integral"]


def combine(parts, before, geo):
    required = {(phase, name) for name in PROFILES for phase in PHASES}
    seen, profiles, receipts = set(), set(), []
    totals = {key: 0 for key in COUNT_KEYS}
    minimum, ties = None, {}
    for raw, expected_sha in parts:
        if len(expected_sha) != 64 or any(c not in "0123456789abcdef" for c in expected_sha):
            raise ValueError("--part needs its exact lowercase SHA-256")
        path = (ROOT/raw).resolve()
        if path.parent != ROOT/"results" or not path.is_file() or path.stat().st_size > FILE_LIMIT:
            raise ValueError("part must be a bounded existing direct results file")
        rec = record(path)
        if rec["sha256"] != expected_sha or any(rec["path"] == r["path"] for r in receipts):
            raise ValueError("changed or repeated part path")
        data = read_json_bounded(path)
        name = data.get("profile")
        if data.get("kind") != "prime379_released_r_profile_unit" or \
                data.get("status") != "PASS_DECLARED_DOMAIN" or data.get("bindings") != before:
            raise ValueError("not a completed current-source profile unit")
        if data.get("starting_collision") != fenc(geo["start"]) or \
                data.get("upper_collision") != fenc(geo["upper"]) or \
                data.get("known_wall_count") != 15 or data.get("germ_certified") is not False:
            raise ValueError("incorrect fourth-edge interval, known walls or partial-profile scope")
        if name not in PROFILES or name in profiles:
            raise ValueError("missing/duplicate/invalid profile")
        profiles.add(name)
        selectors = data.get("selectors")
        if not isinstance(selectors, list) or len(selectors) != 3:
            raise ValueError("each profile must contain all three selectors")
        for selector in selectors:
            phase = selector.get("phase")
            key = (phase, selector.get("profile"))
            if type(phase) is not int or key not in required or key in seen or key[1] != name:
                raise ValueError("invalid or repeated typed selector")
            status = "PROVED_EMPTY_LONG" if profile(name)[0] == B else "FINITE_SCALAR_BOX"
            if selector.get("domain_status") != status or not valid_counts(selector.get("counts"), name):
                raise ValueError("invalid selector counters or domain bound")
            for counter in COUNT_KEYS:
                totals[counter] += selector["counts"][counter]
            seen.add(key)
        tau = None if data["collision"] is None else fdec(data["collision"])
        rows = data["minimizers"]
        if not isinstance(rows, list) or len(rows) > TIE_LIMIT or (tau is None) != (len(rows) == 0):
            raise ValueError("inconsistent part minimum")
        admitted = sum(selector["counts"]["admitted"] for selector in selectors)
        if (admitted == 0) != (tau is None) or len(rows) > admitted:
            raise ValueError("part minimum disagrees with admitted-root counters")
        part_roots = set()
        for row in rows:
            phase = row.get("phase")
            if type(phase) is not int or phase not in PHASES or row.get("profile") != name:
                raise ValueError("invalid minimizer witness selector")
            if next(s["counts"]["admitted"] for s in selectors if s["phase"] == phase) == 0:
                raise ValueError("minimizer witness selector has no admitted roots")
            root = vdec(row["root"])
            key = vkey(root)
            if key in part_roots or validate_part_root(root, phase, name, geo) != tau:
                raise ValueError("invalid/duplicate original reflector, phase or collision")
            if row.get("known_self_exemptions") != known_self_exemptions(root, geo):
                raise ValueError("incorrect full-original known-wall self exemptions")
            part_roots.add(key)
            minimum, ties = keep_minimum(minimum, ties, tau, root, phase, name)
        receipts.append(rec)
    if seen != required or profiles != set(PROFILES):
        raise ValueError("all five profiles and all15 typed selectors are mandatory")
    if totals["scalar_rows"] > 324 or totals["field_roots"] > 648:
        raise ValueError("complete scalar or field-root domain bound exceeded")
    if minimum is None:
        raise RuntimeError("empty complete domain contradicts the proved finite upper")
    cone = germ(minimum, ties, geo)
    return {"kind": "prime379_released_r_complete_domain",
            "status": "PASS_COMPLETE_DOMAIN" if cone["certified"] else "COMPLETE_DOMAIN_GERM_UNRESOLVED",
            "parts": sorted(receipts, key=lambda row: row["path"]),
            "profile_count": len(profiles), "selector_count": len(seen),
            "proved_empty_selectors": 9,
            "counts": totals,
            "collision": fenc(minimum), "minimizers": encode_ties(ties, geo), "germ": cone,
            "starting_collision": fenc(geo["start"]),
            "starting_original_wall": venc(geo["starting_wall"]),
            "upper_collision": fenc(geo["upper"]), "upper_original_mirror": venc(geo["upper_root"]),
            "known_wall_count": len(geo["known"]),
            "scope": "Fourth edge releasing r ONLY with pinned domain/all15 selectors; no volume/period/example claim"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", action="store_true")
    parser.add_argument("--profile", choices=PROFILES)
    parser.add_argument("--combine", action="store_true")
    parser.add_argument("--part", nargs=2, action="append", default=[], metavar=("PATH", "SHA256"))
    parser.add_argument("--output")
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error("optimized Python is unsupported: mandatory geometry assertions must remain enabled")
    before = bindings()
    if args.plan:
        if args.combine or args.profile or args.part or args.output:
            parser.error("--plan is standalone")
        print(json.dumps({"status": "PLAN_ONLY", "bindings": before,
              "units": [{"profile": name, "phases": list(PHASES)} for name in PROFILES],
              "selectors": 15, "proved_empty_long_selectors": 9,
              "maximum_scalar_rows": 324, "executed_units": 0}, sort_keys=True))
        return
    if not os.environ.get("COXETER_CHECK_OWNER"):
        parser.error("use the unchanged capped runner with --parent-lock")
    if not args.output:
        parser.error("a NEW --output is required")
    destination = output_path(args.output)
    if args.combine:
        if args.profile or len(args.part) != 5:
            parser.error("--combine needs exactly five --part PATH SHA256 pairs and no profile")
    elif args.profile is None or args.part:
        parser.error("a unit needs --profile and no --part")
    geo = geometry()
    data = combine(args.part, before, geo) if args.combine else unit(args.profile, geo)
    if bindings() != before:
        raise RuntimeError("source or pinned input changed during execution")
    if args.combine and any(record(ROOT/row["path"]) != row for row in data["parts"]):
        raise RuntimeError("part changed during combine")
    data["bindings"] = before
    data["execution_boundary"] = "root must separately bind actual capped runner ledgers; source hashes are not runtime evidence"
    write_new(destination, data)
    print(json.dumps({"status": data["status"], "output_record": record(destination),
          "profile": data.get("profile"), "selector_count": len(data.get("selectors", [])) or data.get("selector_count"),
          "minimum_ties": len(data["minimizers"]),
          "germ_certified": data.get("germ", {}).get("certified", False)}, sort_keys=True))


if __name__ == "__main__":
    main()


