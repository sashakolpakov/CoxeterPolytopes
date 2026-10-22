#!/usr/bin/env python3
"""PREPARED, UNREVIEWED CODE, UNEXECUTED: composite551 transverse edge.

Seven disjoint profiles, all four phases, then authenticated combination.
Standalone stdlib only: copied exact helpers are local, never imported.
Whole original QUARTER seeds, both proper norms, all ten coordinates,
explicit primitivity/reflection, all13 known walls, every tied hyperplane
and the unique full cone are mandatory. No outcome is presumed.
Root alone owns eventual independently reviewed, guarded execution.
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
PINS = {
    'results/compact_h9_composite551_transverse_batches_checker_spec_20261022.md':
        '1252228722f70b11a2877d4804f9277c67f667e260c237e934a3bd9c9f160f04',
    'results/compact_h9_composite551_full_base_first_binary_edge_and_transverse_contract_candidate_20261022.md':
        '33058946a4dd33df9559d5a2d3d388c87d09b3cd94be4ede1a40161243e73f22',
    'results/review_compact_h9_composite551_full_base_first_edge_and_transverse_domain_full_scope_replay_reviewer_20261022.md':
        '21e7adfb3cbaef05326f7edee0001ac83994d557a5271e1dc56c8fa34aeb911b',
    'results/review_composite551_full_base_first_binary_edge_and_transverse_domain_root_acceptance_20261022.md':
        '0d538a3fb62ee24772c4c8651acce103afdc0b0a462723bdbd4f4623d09a5a20',
    'results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md':
        'b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce',
    'results/bugaenko_facet_binary_extension_and_collar_audit_20261020.md':
        '774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227',
    'results/golden_odd_prime_e8_long_corner_original_module_recognition_candidate_20261022.md':
        '503a332347bae12cf24c328764c085baf0462f6b14c66cef7fcbddd0bc5d7c1e',
    'results/golden_e8_binary_family_full_base_facet_20261020.md':
        'ea5d96e331ac91cb9efc457db2ddaa24bac3341e9e6781038150765895675afd',
    'data/known_examples/bugaenko_compact_h8.json':
        'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
    'scripts/check_compact_h9_prime379_fourth_edge_released_r_batches_20261022.py':
        'b3a31476176a2c6267e4ac9dcb1ab16eb2e57dd5c2f61b4c673f37085e970219',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}
PROFILES = ("norm2_orthogonal", "norm2_simple", "norm2_golden", "norm2_e8",
            "norm2p_orthogonal", "norm2q_orthogonal", "norm2b_orthogonal")
PHASES = (0, 1, 2, 3)
TIE_LIMIT = 128
FILE_LIMIT = 4*1024**2
COUNT_KEYS = ("scalar_rows", "inside_bounds", "folded_window", "field_roots",
              "quarter_seed", "original_integral", "primitive_reflector", "admitted")
SCOPE = "composite551_first_transverse_edge_only"


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
D, B = ONE + 24*PHI, F(24, 25)
P, QDIV = F(4, 1), F(5, 4)
H = tuple(F(x) for x in (2, 3, 4, 6, 5, 4, 3, 2))
K6 = tuple(F(x) for x in (4, 6, 8, 12, 10, 8, 6, 3))


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
    return tuple(v3[:2]) + tuple(h*x + epsilon*y for x, y in zip(H, K6))


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


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key")
        result[key] = value
    return result


def invalid_number(value):
    raise ValueError("noninteger/nonfinite JSON number: " + value)


def read_json_bounded(path):
    with path.open("rb") as stream:
        raw = stream.read(FILE_LIMIT+1)
    if len(raw) > FILE_LIMIT:
        raise ValueError("bounded JSON input exceeded4MiB")
    return json.loads(raw, object_pairs_hook=unique_object,
                      parse_float=invalid_number, parse_constant=invalid_number)


def record(path):
    digest, size = hashlib.sha256(), 0
    with path.open("rb") as stream:
        while chunk := stream.read(65536):
            digest.update(chunk)
            size += len(chunk)
    return {"path": str(path.relative_to(ROOT)), "bytes": size,
            "sha256": digest.hexdigest()}


def valid_sha(value):
    return isinstance(value, str) and len(value) == 64 and \
        all(c in "0123456789abcdef" for c in value)


def bindings(expected_self):
    if not valid_sha(expected_self):
        raise ValueError("mandatory externally frozen self SHA256")
    rows = {}
    for name, expected in PINS.items():
        row = record(ROOT/name)
        if row["sha256"] != expected:
            raise RuntimeError("changed pinned input: " + name)
        rows[name] = row
    me = record(SELF)
    if me["sha256"] != expected_self:
        raise RuntimeError("worker self SHA256 mismatch")
    rows[str(SELF.relative_to(ROOT))] = me
    return rows


def output_path(raw):
    path = (ROOT/raw).resolve()
    if path.parent != ROOT/"results" or path.suffix != ".json" or path.exists():
        raise ValueError("output must be a NEW direct results/*.json file")
    return path


def write_new(path, data):
    text = json.dumps(data, sort_keys=True, indent=2, allow_nan=False)+"\n"
    if len(text.encode()) > FILE_LIMIT:
        raise RuntimeError("output storage bound exceeded; no truncation")
    with path.open("x", encoding="utf-8") as stream:
        stream.write(text)


def ogcd(a, b):
    """Euclidean gcd in Z[phi]; exact nearest-coordinate division."""
    if not a.integral() or not b.integral():
        raise ValueError("integer-ring gcd input")
    while b != ZERO:
        quotient = a/b
        rounded = F((quotient.a+Fraction(1, 2)).__floor__(),
                    (quotient.b+Fraction(1, 2)).__floor__())
        remainder = a-rounded*b
        if not remainder.integral() or abs(remainder.norm()) >= abs(b.norm()):
            raise ArithmeticError("Euclidean norm failed to decrease")
        a, b = b, remainder
    return a


def primitive(v):
    common = ZERO
    for x in v:
        if not x.integral():
            return False
        common = ogcd(common, x)
        if abs(common.norm()) == 1:
            return True
    return False


def integral_reflector(v, d, gram):
    if len(v) != 10 or d not in (ONE, P, QDIV, B) or \
            not all(x.integral() for x in v) or dot(v, v, gram) != 2*d:
        return False
    if not primitive(v):
        return False
    cov = mv(gram, v)
    return all((x*y/d).integral() for x in v if x != ZERO
               for y in cov if y != ZERO)


def positive_proportional(vector, wall):
    if len(vector) != 10 or len(wall) != 10:
        raise ValueError("self recognition needs all ten original coordinates")
    pivot = next((i for i, x in enumerate(wall) if x != ZERO), None)
    if pivot is None:
        raise ValueError("zero known wall")
    ratio = vector[pivot]/wall[pivot]
    return ratio.sign() > 0 and all(x == ratio*y for x, y in zip(vector, wall))


def hyperplane_key(v):
    pivot = next((x for x in v if x != ZERO), None)
    if pivot is None:
        raise ValueError("zero mirror")
    return vkey(tuple(x/pivot for x in v))


def geometry():
    """Reconstruct the accepted first germ and ORIGINAL quarter boost."""
    g3 = ((-2*PHI, ONE, ZERO), (ONE, F(12), ZERO), (ZERO, ZERO, F(2)))
    gram = [[ZERO for _ in range(10)] for _ in range(10)]
    gram[0][:2], gram[1][:2] = [-2*PHI, ONE], [ONE, F(12)]
    for i in range(8):
        gram[i+2][i+2] = F(2)
    for i, j in ((0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)):
        gram[i+2][j+2] = gram[j+2][i+2] = F(-1)
    gram = tuple(tuple(row) for row in gram)
    native = [tuple(F(k == j+2) for k in range(10)) for j in range(8)]
    x3, y3, hh = (ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (ZERO, ZERO, ONE)
    q3, r3 = (F(11, 10), F(-5, -7), ZERO), (F(1, 2), F(-1, -1), ZERO)
    r0 = (ONE, 2*PHI)+(ZERO,)*8
    q, r, h = original(q3, 0), original(r3, 0), original(hh, 0)
    xx = original(x3, 0)
    assert B == P*QDIV == PHI*D and B.norm() == 551
    assert P.norm() == 19 and QDIV.norm() == 29
    assert all(d.sign() > 0 and d.conjugate().sign() > 0 for d in (P, QDIV, B))
    assert D.sign() > 0 and D.conjugate().sign() < 0
    assert dot(q3, q3, g3) == -2*D and dot(r3, r3, g3) == 2
    assert dot(q3, r3, g3) == ZERO
    assert dot(q, xx, gram) == F(-25, -49)
    assert dot(r0, r0, gram) == 2*B
    assert mv(gram, r0) == (ZERO, D)+(ZERO,)*8
    assert mv(gram, r) == (F(-5, -7), F(-11, -10))+(ZERO,)*8
    assert all(integral_reflector(v, d, gram) for v, d in ((r0, B), (r, ONE)))
    for i, v in enumerate(native):
        assert dot(h, v, gram) == F(i == 7)
        assert dot((ZERO, ZERO)+K6, v, gram) == F(i == 6)
    assert dot(h, h, gram) == 2
    assert dot((ZERO, ZERO)+K6, (ZERO, ZERO)+K6, gram) == 6
    assert dot(h, (ZERO, ZERO)+K6, gram) == 3
    assert add(scale(q3, F(1, 1)), r3, -F(5, 7)) == scale(x3, 2)
    assert add(scale(q3, F(1, 2)), r3, -F(11, 10)) == scale(y3, 2)

    calibration = read_json_bounded(ROOT/"data/known_examples/bugaenko_compact_h8.json")
    def number(row):
        return F(Fraction(row["a"]), Fraction(row["b"]))
    cg = tuple(tuple(number(x) for x in row) for row in calibration["lattice_gram"])
    indexes = (0,)+tuple(range(2, 10))
    assert len(cg) == 9 and all(len(row) == 9 for row in cg)
    assert all(cg[i][j] == gram[indexes[i]][indexes[j]]
               for i in range(9) for j in range(9))
    base = []
    for row in calibration["simple_roots"]:
        assert len(row) == 9
        base.append((number(row[0]), ZERO)+tuple(number(x) for x in row[1:]))
    assert len(base) == len(set(base)) == 11
    assert all(integral_reflector(v, ONE, gram) for v in base)
    assert all(scale(v, -1) in base for v in native)
    known = base+[r0, r]
    assert len(set(known)) == 13
    assert all(dot(q, wall, gram).sign() <= 0 for wall in known)
    start_germ = [scale(v, -1) for v in native]+[r]
    assert {vkey(v) for v in known if dot(q, v, gram) == ZERO} == \
        {vkey(v) for v in start_germ}

    # Active operators are zero on E7; retain old e7's H/2 projection.
    ecols, jcols = [], []
    for v in identity(3):
        aq, ah = -dot(q3, v, g3)/(2*D), dot(hh, v, g3)/2
        ecols.append(add(scale(q3, aq), hh, ah))
        jcols.append(add(scale(hh, D*aq), q3, ah))
    projection = tuple(tuple(ecols[j][i] for j in range(3)) for i in range(3))
    skew = tuple(tuple(jcols[j][i] for j in range(3)) for i in range(3))
    assert mm(projection, projection) == projection
    assert mv(projection, r3) == (ZERO, ZERO, ZERO)
    assert mv(skew, q3) == scale(hh, D) and mv(skew, hh) == q3
    assert mv(skew, r3) == (ZERO, ZERO, ZERO)
    assert mm(skew, skew) == tuple(tuple(D*x for x in row) for row in projection)
    assert mm(skew, projection) == skew == mm(projection, skew)
    gs = mm(g3, skew)
    assert all(gs[i][j] == -gs[j][i] for i in range(3) for j in range(3))
    for operator in (projection, skew):
        assert all((2*x).integral() for image in full_basis_images(operator, ZERO)
                   for x in image)
    aa, bb = F(Fraction(47, 2), 37), F(Fraction(7, 2), 6)
    assert aa*aa-D*bb*bb == ONE and (F(84)-aa).sign() > 0
    matrix = tuple(tuple(F(i == j)+(aa-1)*projection[i][j]+bb*skew[i][j]
                         for j in range(3)) for i in range(3))
    inverse = tuple(tuple(F(i == j)+(aa-1)*projection[i][j]-bb*skew[i][j]
                          for j in range(3)) for i in range(3))
    assert mm(matrix, inverse) == identity(3) == mm(inverse, matrix)
    powers, inverse_powers = [identity(3)], [identity(3)]
    for _ in range(3):
        powers.append(mm(matrix, powers[-1]))
        inverse_powers.append(mm(inverse, inverse_powers[-1]))
    for family, other in ((powers, inverse_powers), (inverse_powers, powers)):
        for j in PHASES:
            images = full_basis_images(family[j])
            assert all((4*x).integral() for image in images for x in image)
            assert all(dot(images[i], images[k], gram) == gram[i][k]
                       for i in range(10) for k in range(10))
            assert all(images[i+2] == native[i] for i in range(7))
            assert mv(family[j], r3) == r3
            assert mm(family[j], other[j]) == identity(3)
            assert dot(q3, mv(family[j], q3), g3).sign() < 0
            if j == 3:
                assert all(x.integral() for image in images for x in image)
    t, s = 2*aa, 2*bb
    assert t*t-1 == 4*F(1921, 3108)
    a3, b3 = (t*t*t-3*t)/2, s*(t*t-1)/2
    assert ((a3-1)/2).integral() and (b3/2).integral()
    assert mv(powers[3], q3) == add(scale(q3, a3), hh, D*b3)
    upper = D*b3/a3
    upper_root = full_basis_images(powers[3])[9]
    assert integral_reflector(upper_root, ONE, gram)
    assert dot(q, upper_root, gram) == -D*b3
    assert dot(h, upper_root, gram) == a3
    assert upper.sign() > 0 and (D-upper*upper).sign() > 0
    assert dot(add(q, h, upper), upper_root, gram) == ZERO

    # Known positive self mirrors cannot be a new cutter on this interval.
    positive_known = {F(25, 49)/2: 0, B/2: 0}
    for wall in known:
        derivative, value = dot(h, wall, gram), dot(q, wall, gram)
        if derivative.sign() > 0:
            tau = -value/derivative
            assert tau in positive_known and tau.sign() > 0
            assert (tau*tau-D).sign() > 0
            positive_known[tau] += 1
    assert positive_known == {F(25, 49)/2: 2, B/2: 1}
    retained = [scale(v, -1) for v in native[:7]]+[r]
    assert all(dot(q, v, gram) == ZERO == dot(h, v, gram) for v in retained)
    assert dot(h, scale(native[7], -1), gram) == -1
    return {"gram": gram, "g3": g3, "q3": q3, "r3": r3,
            "line_q": q, "line_h": h, "A": aa, "powers": powers,
            "inverse_powers": inverse_powers, "upper": upper,
            "upper_root": upper_root, "known": known, "retained": retained}


def profile(name):
    return {"norm2_orthogonal": (ONE, ZERO, 0),
            "norm2_simple": (ONE, ONE, 0),
            "norm2_golden": (ONE, PHI, 0),
            "norm2_e8": (ONE, ZERO, 1),
            "norm2p_orthogonal": (P, ZERO, 0),
            "norm2q_orthogonal": (QDIV, ZERO, 0),
            "norm2b_orthogonal": (B, ZERO, 0)}[name]


def box(name):
    d = profile(name)[0]
    return (12, 14) if d == ONE else (6, 6) if d == P else \
        (6, 5) if d == QDIV else (2, 2)


def scalar_rows(name):
    d = profile(name)[0]
    cmax, fmax = box(name)
    for c in range(-cmax, cmax+1):
        for f in range(-fmax, fmax+1):
            yield d*F(c, f)/4


def scalar_in_box(name, z):
    d = profile(name)[0]
    coeff = 4*z/d
    cmax, fmax = box(name)
    return coeff.integral() and abs(coeff.a) <= cmax and abs(coeff.b) <= fmax


def inside_scalar_bounds(name, z, geo):
    d = profile(name)[0]
    return (d*(geo["A"]+1)/2-z*z).sign() >= 0 and \
        (d.conjugate()-z.conjugate()*z.conjugate()).sign() >= 0


def validate_original(root, d, geo, checked=False):
    if not checked and not integral_reflector(root, d, geo["gram"]):
        return None
    value = dot(root, geo["line_q"], geo["gram"])
    derivative = dot(root, geo["line_h"], geo["gram"])
    if value.sign() >= 0 or derivative.sign() <= 0:
        return None
    tau = -value/derivative
    if tau.sign() <= 0 or (geo["upper"]-tau).sign() < 0:
        return None
    # geometry() proved that none of these known mirrors cuts the interval.
    if any(positive_proportional(root, wall) for wall in geo["known"]):
        return None
    if any(dot(root, wall, geo["gram"]).sign() > 0 for wall in geo["known"]):
        return None
    return tau


def validate_part_root(root, phase, name, geo):
    d, kappa, epsilon = profile(name)
    gram = geo["gram"]
    native = [tuple(F(k == j+2) for k in range(10)) for j in range(8)]
    if any(dot(root, native[i], gram) != F(epsilon if i == 6 else 0)
           for i in range(7)):
        return None
    h = dot(root, native[7], gram)
    expected = tuple(root[:2])+tuple(h*x+epsilon*y for x, y in zip(H, K6))
    if expected != root or dot(root, original(geo["r3"], 0), gram) != -kappa:
        return None
    active = tuple(root[:2])+(h+Fraction(3, 2)*epsilon,)
    seed = mv(geo["inverse_powers"][phase], active)
    if not all((4*x).integral() for x in original(seed, epsilon)):
        return None
    z = seed[2]
    kk = 2*d-kappa*kappa/2-F(Fraction(3, 2))*epsilon*epsilon
    if not scalar_in_box(name, z) or not inside_scalar_bounds(name, z, geo) or \
            (kk*(geo["A"]+1)/4-z*z).sign() < 0:
        return None
    if dot(seed, geo["r3"], geo["g3"]) != -kappa:
        return None
    time = -dot(seed, geo["q3"], geo["g3"])/(2*D)
    rebuilt = add(add(scale(geo["q3"], time), geo["r3"], -kappa/2),
                  (ZERO, ZERO, ONE), z)
    square = (2*z*z+kappa*kappa/2+Fraction(3, 2)*epsilon*epsilon-2*d)/(2*D)
    if rebuilt != seed or time not in field_roots(square):
        return None
    if dot(seed, seed, geo["g3"])+Fraction(3, 2)*epsilon*epsilon != 2*d:
        return None
    if original(mv(geo["powers"][phase], seed), epsilon) != root:
        return None
    return validate_original(root, d, geo)


def keep_minimum(minimum, ties, tau, root, phase, name):
    if minimum is None or (tau-minimum).sign() < 0:
        minimum, ties = tau, {}
    if tau == minimum:
        key = hyperplane_key(root)
        if key in ties and not positive_proportional(root, ties[key][0]):
            raise ArithmeticError("oppositely oriented tied hyperplane")
        if key not in ties:
            ties[key] = (root, phase, name)
        if len(ties) > TIE_LIMIT:
            raise RuntimeError("all-tie storage limit exceeded; no truncation")
    return minimum, ties


def encode_ties(ties):
    return [{"root": venc(ties[key][0]), "phase": ties[key][1],
             "profile": ties[key][2]} for key in sorted(ties)]


def unit(name, geo):
    d, kappa, epsilon = profile(name)
    kk = 2*d-kappa*kappa/2-F(Fraction(3, 2))*epsilon*epsilon
    assert kk.sign() > 0
    minimum, ties, selectors = None, {}, []
    cmax, fmax = box(name)
    for phase in PHASES:
        counts = {key: 0 for key in COUNT_KEYS}
        for z in scalar_rows(name):
            counts["scalar_rows"] += 1
            if not inside_scalar_bounds(name, z, geo):
                continue
            counts["inside_bounds"] += 1
            if (kk*(geo["A"]+1)/4-z*z).sign() < 0:
                continue
            counts["folded_window"] += 1
            square = (2*z*z+kappa*kappa/2+Fraction(3, 2)*epsilon*epsilon-2*d)/(2*D)
            for time in field_roots(square):
                counts["field_roots"] += 1
                seed = add(add(scale(geo["q3"], time), geo["r3"], -kappa/2),
                           (ZERO, ZERO, ONE), z)
                if not all((4*x).integral() for x in original(seed, epsilon)):
                    continue
                counts["quarter_seed"] += 1
                assert dot(seed, seed, geo["g3"])+Fraction(3, 2)*epsilon*epsilon == 2*d
                assert dot(seed, geo["r3"], geo["g3"]) == -kappa
                root = original(mv(geo["powers"][phase], seed), epsilon)
                if not all(x.integral() for x in root):
                    continue
                counts["original_integral"] += 1
                if not integral_reflector(root, d, geo["gram"]):
                    continue
                counts["primitive_reflector"] += 1
                tau = validate_original(root, d, geo, checked=True)
                if tau is not None:
                    counts["admitted"] += 1
                    minimum, ties = keep_minimum(minimum, ties, tau, root, phase, name)
        assert counts["scalar_rows"] == (2*cmax+1)*(2*fmax+1)
        assert counts["field_roots"] <= 2*counts["folded_window"]
        selectors.append({"profile": name, "phase": phase, "counts": counts})
    return {"kind": "composite551_transverse_profile_unit", "schema": 1,
            "status": "PASS_DECLARED_DOMAIN", "scope": SCOPE,
            "profile": name, "selectors": selectors, "starting_collision": fenc(ZERO),
            "upper_collision": fenc(geo["upper"]), "known_wall_count": 13,
            "collision": None if minimum is None else fenc(minimum),
            "minimizers": encode_ties(ties), "germ_certified": False}


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
        a[k] = [x/p for x in a[k]]
        for i in range(n):
            if i != k and a[i][k] != ZERO:
                c = a[i][k]
                a[i] = [x-c*y for x, y in zip(a[i], a[k])]
    return tuple(row[-1] for row in a)


def germ(minimum, ties, geo):
    point = add(geo["line_q"], geo["line_h"], minimum)
    gram = geo["gram"]
    if dot(point, point, gram).sign() >= 0 or \
            dot(point, geo["line_q"], gram).sign() >= 0:
        raise RuntimeError("minimum endpoint must be future timelike")
    if any(dot(point, wall, gram).sign() > 0 for wall in geo["known"]):
        raise RuntimeError("minimum violates a known actual wall")
    chosen = []
    for root, _, _ in ties.values():
        basis = geo["retained"]+[root]
        gg = tuple(tuple(dot(x, y, gram) for y in basis) for x in basis)
        if not positive_definite(gg) or any(dot(point, x, gram) != ZERO for x in basis):
            continue
        inside = True
        for other, _, _ in ties.values():
            coeff = solve(gg, tuple(dot(x, other, gram) for x in basis))
            if coeff is None or any(c.sign() < 0 for c in coeff):
                inside = False
                break
            rebuilt = tuple(sum((coeff[j]*basis[j][i] for j in range(9)), ZERO)
                            for i in range(10))
            if rebuilt != other:
                inside = False
                break
        if inside:
            chosen.append((root, gg))
    unique = len(chosen) == 1
    return {"certified": unique, "cone_candidates": len(chosen),
            "new_simple_wall": venc(chosen[0][0]) if unique else None,
            "gram": [[fenc(x) for x in row] for row in chosen[0][1]] if unique else None,
            "endpoint_original_ray": [fenc(x) for x in point],
            "scope": "conditional actual full finite germ; separate geometric interpretation required"}


def valid_counts(counts, name):
    if not isinstance(counts, dict) or set(counts) != set(COUNT_KEYS) or \
            any(type(v) is not int or v < 0 for v in counts.values()):
        return False
    cmax, fmax = box(name)
    return counts["scalar_rows"] == (2*cmax+1)*(2*fmax+1) and \
        counts["inside_bounds"] <= counts["scalar_rows"] and \
        counts["folded_window"] <= counts["inside_bounds"] and \
        counts["field_roots"] <= 2*counts["folded_window"] and \
        counts["quarter_seed"] <= counts["field_roots"] and \
        counts["original_integral"] <= counts["quarter_seed"] and \
        counts["primitive_reflector"] <= counts["original_integral"] and \
        counts["admitted"] <= counts["primitive_reflector"]


def combine(parts, before, geo):
    required = {(phase, name) for name in PROFILES for phase in PHASES}
    seen, profiles, receipts = set(), set(), []
    totals = {key: 0 for key in COUNT_KEYS}
    minimum, ties = None, {}
    for raw, expected_sha in parts:
        if not valid_sha(expected_sha):
            raise ValueError("each part needs its exact SHA256")
        path = (ROOT/raw).resolve()
        if path.parent != ROOT/"results" or not path.is_file() or \
                not 0 < path.stat().st_size <= FILE_LIMIT:
            raise ValueError("part must be a bounded direct results file")
        rec = record(path)
        if rec["sha256"] != expected_sha or any(rec["path"] == r["path"] for r in receipts):
            raise ValueError("changed or duplicate part")
        data = read_json_bounded(path)
        name = data.get("profile")
        if data.get("kind") != "composite551_transverse_profile_unit" or \
                type(data.get("schema")) is not int or data["schema"] != 1 or \
                data.get("status") != "PASS_DECLARED_DOMAIN" or \
                data.get("scope") != SCOPE or data.get("bindings") != before:
            raise ValueError("not a completed current-source profile unit")
        if name not in PROFILES or name in profiles:
            raise ValueError("missing/duplicate/invalid profile")
        profiles.add(name)
        if data.get("starting_collision") != fenc(ZERO) or \
                data.get("upper_collision") != fenc(geo["upper"]) or \
                type(data.get("known_wall_count")) is not int or \
                data["known_wall_count"] != 13 or data.get("germ_certified") is not False:
            raise ValueError("incorrect interval, known walls or partial-unit scope")
        selectors = data.get("selectors")
        if not isinstance(selectors, list) or len(selectors) != 4:
            raise ValueError("each profile must contain four phase selectors")
        for selector in selectors:
            phase = selector.get("phase")
            key = (phase, selector.get("profile"))
            if type(phase) is not int or key not in required or key in seen or \
                    key[1] != name or not valid_counts(selector.get("counts"), name):
                raise ValueError("invalid/repeated selector or counters")
            seen.add(key)
            for counter in COUNT_KEYS:
                totals[counter] += selector["counts"][counter]
        tau = None if data.get("collision") is None else fdec(data["collision"])
        rows = data.get("minimizers")
        if not isinstance(rows, list) or len(rows) > TIE_LIMIT or \
                (tau is None) != (len(rows) == 0):
            raise ValueError("inconsistent profile minimum")
        admitted = sum(s["counts"]["admitted"] for s in selectors)
        if (admitted == 0) != (tau is None) or len(rows) > admitted:
            raise ValueError("minimum disagrees with admitted-root counters")
        part_roots = set()
        for row in rows:
            phase = row.get("phase")
            if type(phase) is not int or phase not in PHASES or row.get("profile") != name:
                raise ValueError("invalid minimizer witness")
            if next(s["counts"]["admitted"] for s in selectors if s["phase"] == phase) == 0:
                raise ValueError("minimizer phase has no admission")
            root = vdec(row["root"])
            key = hyperplane_key(root)
            if key in part_roots or validate_part_root(root, phase, name, geo) != tau:
                raise ValueError("invalid/duplicate original root, phase or collision")
            part_roots.add(key)
            minimum, ties = keep_minimum(minimum, ties, tau, root, phase, name)
        if record(path) != rec:
            raise RuntimeError("part changed while reading")
        receipts.append(rec)
    if seen != required or profiles != set(PROFILES):
        raise ValueError("all seven profiles and all28 selectors are mandatory")
    if totals["scalar_rows"] != 12948 or totals["field_roots"] > 25896:
        raise ValueError("complete scalar/root count bound")
    if minimum is None:
        raise RuntimeError("empty domain contradicts the proved finite upper")
    cone = germ(minimum, ties, geo)
    return {"kind": "composite551_transverse_complete_domain", "schema": 1,
            "status": "PASS_COMPLETE_DOMAIN" if cone["certified"] else
                      "COMPLETE_DOMAIN_GERM_UNRESOLVED", "scope": SCOPE,
            "parts": sorted(receipts, key=lambda r: r["path"]),
            "profile_count": 7, "selector_count": len(seen), "counts": totals,
            "collision": fenc(minimum), "minimizers": encode_ties(ties), "germ": cone,
            "starting_collision": fenc(ZERO), "upper_collision": fenc(geo["upper"]),
            "upper_original_mirror": venc(geo["upper_root"]), "known_wall_count": 13}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-source-sha256", required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--profile", choices=PROFILES)
    group.add_argument("--combine", action="store_true")
    parser.add_argument("--part", nargs=2, action="append", default=[], metavar=("PATH", "SHA256"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error("optimized Python disables mandatory geometry assertions")
    owner = os.environ.get("COXETER_CHECK_OWNER", "")
    if not owner.isdigit() or os.getppid() != int(owner):
        parser.error("use the unchanged capped runner with --parent-lock")
    if args.combine and len(args.part) != 7 or not args.combine and args.part:
        parser.error("combine needs seven --part PATH SHA pairs; a profile accepts none")
    before = bindings(args.expected_source_sha256)
    destination = output_path(args.output)
    geo = geometry()
    data = combine(args.part, before, geo) if args.combine else unit(args.profile, geo)
    if bindings(args.expected_source_sha256) != before:
        raise RuntimeError("source or pinned input changed during execution")
    if args.combine and any(record(ROOT/r["path"]) != r for r in data["parts"]):
        raise RuntimeError("part changed during combination")
    data["bindings"] = before
    data["execution_boundary"] = "root must bind actual capped ledgers/stdout; source preparation is not execution"
    write_new(destination, data)
    print(json.dumps({"status": data["status"], "output_record": record(destination),
          "profile": data.get("profile"), "selector_count": len(data.get("selectors", [])) or
          data.get("selector_count"), "minimum_ties": len(data["minimizers"]),
          "germ_certified": data.get("germ", {}).get("certified", False)}, sort_keys=True))


if __name__ == "__main__":
    main()


