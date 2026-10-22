#!/usr/bin/env python3
"""PREPARED, UNREVIEWED, UNEXECUTED: norm379 third actual edge, releasing m.

Standalone stdlib source; helpers copied from frozen e631, never imported.
Five profiles x eleven phases, then separately authenticated combination.
Transport w=F V^j seed uses the WHOLE original quarter seed and ten coordinates.
ALL FOURTEEN original known walls, exact positive self recognition, strict
starting interval and all ties/full nine-wall cone are mandatory. No outcome,
period, volume or example is presumed. Root owns eventual guarded execution.
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
SPEC = "results/compact_h9_prime379_third_edge_released_m_checker_spec_20261022.md"
PINS = {
    'results/compact_h9_prime379_third_edge_released_m_checker_spec_20261022.md':
        '6b43db1fcd1d38e440288c371a48877c4e41c38c8e018df32e5df666353f42f0',
    'results/compact_h9_prime379_first_transverse_wall_and_e7_golden_germ_candidate_20261022.md':
        'b0fc4656d5a41f88bc7941ceb1d3efff8be95775427d5e382e6bdec123264179',
    'results/review_prime379_first_transverse_wall_full_germ_root_acceptance_20261022.md':
        '5728de2f30a10d35b12450c5851bfce01064d997c6f6569f70c1bc0b065420ec',
    'results/review_prime379_first_transverse_runtime_and_full_germ_compact_reviewer_20261022.md':
        '9cc01b56bef53f5ed17560771cb0e3328adc4a733f817a3ad7f3e0f757482d92',
    'results/prime379_transverse_runtime_manifest_20261022.json':
        'ae53f5fc76d9d62c28c6bc4a905ce094acb4c6165c483cea9a3c6f11357e8c20',
    'results/compact_h9_prime379_first_binary_edge_and_quarter_folded_transverse_domain_candidate_20261022.md':
        '7e7ffd616d0e025aeaff84abae1c9fcfb7d930c8a3879a7bb43734f308cacd3a',
    'results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md':
        'b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce',
    'data/known_examples/bugaenko_compact_h8.json':
        'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
    'results/compact_h9_prime379_third_edge_released_m_complete_folded_domain_candidate_20261022.md':
        '0ff2e60377f26528378c9803024765fc443c7503ae071cc797936d49db3f56ac',
    'results/review_prime379_third_edge_complete_domain_root_acceptance_20261022.md':
        'a1d2763bf66faa04d3f91a65633ee05a5708a0ca60936053f9196a3096ed67cb',
    'results/review_prime379_third_edge_released_m_complete_domain_full_scope_compact_reviewer_20261022.md':
        '5e578817b4b08f3b99720f8587c78528d3b3ebcfefd5ff68ffc5c6d422a4e2c3',
    'results/prime379_transverse_combine_20261022.json':
        '861e3b5385339946ab1691e70ed3d9bfbf37932a1b814b6deb9ddbb6d958fce3',
    'scripts/check_compact_h9_prime379_quarter_folded_transverse_batches_20261022.py':
        'e6310c0a345d682b7769a4d6d7feb9dedd52b3229582c67858300900a426b7cd',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}
PROFILES = ("norm2_orthogonal", "norm2_simple", "norm2_golden",
            "norm2_e8", "norm2b_orthogonal")
PHASES = tuple(range(11))
TIE_LIMIT = 128
FILE_LIMIT = 2*1024**2
COUNT_KEYS = ("scalar_rows", "inside_bounds", "folded_window", "field_roots",
              "quarter_seed", "original_integral", "admitted")


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

    def quarter(self):
        return (4*self).integral()


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
    """Reconstruct original L, accepted start, F and the NEW outgoing edge."""
    g3 = ((-2*PHI, ONE, ZERO), (ONE, F(10), ZERO), (ZERO, ZERO, F(2)))
    q = (F(18, 25), F(-12, -17), ZERO)
    m = (F(2, 5), F(-2, -3), ZERO)
    hh = (ZERO, ZERO, ONE)
    assert B == PHI*D and B.norm() == 379
    assert D.sign() > 0 and D.conjugate().sign() < 0
    assert B.sign() > 0 and B.conjugate().sign() > 0
    assert dot(q, q, g3) == -2*D
    assert dot(m, m, g3) == 2 and dot(q, m, g3) == ZERO
    assert dot(q, (ONE, ZERO, ZERO), g3).sign() < 0
    assert dot(hh, hh, g3) == 2

    # E projects onto(q,H), Jq=DH, JH=q; both are ZERO on E7.
    basis3 = identity(3)
    ecols = [add(scale(q, -dot(v, q, g3)/(2*D)), hh, v[2]) for v in basis3]
    jcols = [add(scale(q, v[2]), hh, -dot(v, q, g3)/2) for v in basis3]
    projector = tuple(tuple(ecols[j][i] for j in range(3)) for i in range(3))
    skew = tuple(tuple(jcols[j][i] for j in range(3)) for i in range(3))
    assert mm(projector, projector) == projector
    assert mv(projector, q) == q and mv(projector, hh) == hh
    assert mv(projector, m) == (ZERO, ZERO, ZERO)
    assert mv(skew, q) == scale(hh, D) and mv(skew, hh) == q
    assert mv(skew, m) == (ZERO, ZERO, ZERO)
    assert mm(skew, skew) == tuple(tuple(D*x for x in row) for row in projector)
    assert mm(skew, projector) == skew == mm(projector, skew)
    gs = mm(g3, skew)
    assert all(gs[i][j] == -gs[j][i] for i in range(3) for j in range(3))
    for operator in (projector, skew):
        assert all((2*x).integral() for image in full_basis_images(operator, ZERO)
                   for x in image)

    aa, bb = F(133, Fraction(433, 2)), F(23, Fraction(75, 2))
    assert aa*aa-D*bb*bb == ONE
    assert aa.sign() > 0 and bb.sign() > 0 and (F(484)-aa).sign() > 0
    matrix = tuple(tuple(F(i == j)+(aa-1)*projector[i][j]+bb*skew[i][j]
                         for j in range(3)) for i in range(3))
    inverse = tuple(tuple(F(i == j)+(aa-1)*projector[i][j]-bb*skew[i][j]
                          for j in range(3)) for i in range(3))
    assert mm(matrix, inverse) == identity(3) == mm(inverse, matrix)
    assert mv(matrix, m) == m and mv(inverse, m) == m
    columns = [tuple(matrix[i][j] for i in range(3)) for j in range(3)]
    assert all(dot(columns[i], columns[j], g3) == g3[i][j]
               for i in range(3) for j in range(3))
    assert dot(q, mv(matrix, q), g3) == -2*D*aa
    powers, inverse_powers = [identity(3)], [identity(3)]
    for _ in range(10):
        powers.append(mm(matrix, powers[-1]))
        inverse_powers.append(mm(inverse, inverse_powers[-1]))
    for family in (powers, inverse_powers):
        for j in PHASES:
            images = full_basis_images(family[j])
            assert all(x.quarter() for image in images for x in image)
            if j == 10:
                assert all(x.integral() for image in images for x in image)
    assert mm(powers[10], inverse_powers[10]) == identity(3)
    q10 = mv(powers[10], q)
    a10, b10 = -dot(q10, q, g3)/(2*D), q10[2]/D
    assert q10 == add(scale(q, a10), hh, D*b10)
    assert a10.sign() > 0 and b10.sign() > 0
    assert a10*a10-D*b10*b10 == ONE
    upper = D*b10/a10
    assert upper.sign() > 0 and (D-upper*upper).sign() > 0

    gram = [[ZERO for _ in range(10)] for _ in range(10)]
    gram[0][:2], gram[1][:2] = [-2*PHI, ONE], [ONE, F(10)]
    for i in range(8):
        gram[i+2][i+2] = F(2)
    for i, j in ((0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)):
        gram[i+2][j+2] = gram[j+2][i+2] = F(-1)
    gram = tuple(tuple(row) for row in gram)
    native = [tuple(F(k == j+2) for k in range(10)) for j in range(8)]
    assert dot((ZERO, ZERO)+H, (ZERO, ZERO)+H, gram) == 2
    assert dot((ZERO, ZERO)+J, (ZERO, ZERO)+J, gram) == 6
    assert dot((ZERO, ZERO)+H, (ZERO, ZERO)+J, gram) == 3

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
    assert all(all(x.integral() for x in wall) and dot(wall, wall, gram) == 2
               for wall in base)
    assert all(scale(v, -1) in base for v in native)
    r0 = (ONE, 2*PHI)+(ZERO,)*8
    original_m, original_q = original(m, 0), original(q, 0)
    assert dot(r0, r0, gram) == 2*B
    assert all(x.integral() for x in original_m)
    assert mv(gram, original_m)[:2] == (F(-12, -17), F(-18, -25))
    known = base+[r0, original_m]
    assert len(set(known)) == 13
    assert all(dot(original_q, wall, gram).sign() <= 0 for wall in known)
    upper_root = tuple(full_basis_images(powers[10])[9])
    assert integral_reflector(upper_root, ONE, gram)
    assert dot(upper_root, original_q, gram) == -D*b10
    assert dot(upper_root, (ZERO, ZERO)+H, gram) == a10
    assert dot(upper_root, original(add(q, hh, upper), 0), gram) == ZERO
    # NEW edge transport; no old filter or minimum is used for this interval.
    a0, z0 = F(2572, Fraction(8323, 2)), F(14856, 24036)
    tau0 = F(Fraction(9575, 8376), Fraction(3997, 1396))
    r3 = add(add(scale(q, a0), m, -PHI/2), hh, z0)
    r = original(r3, 0)
    assert dot(r3, r3, g3) == 2 and dot(m, r3, g3) == -PHI
    assert all(x.integral() for x in r) and integral_reflector(r, ONE, gram)
    assert tau0 == D*a0/z0 and (upper-tau0).sign() > 0
    c0 = z0*z0-D*a0*a0
    assert c0 == (3-PHI)/4 and (c0-F(Fraction(1, 4))).sign() > 0
    start3 = add(q, hh, tau0)
    start = original(start3, 0)
    assert dot(start3, start3, g3).sign() < 0
    assert dot(start3, q, g3).sign() < 0
    assert dot(start, r, gram) == ZERO == dot(start, original_m, gram)

    old = read_json_bounded(ROOT/"results/prime379_transverse_combine_20261022.json")
    assert old["kind"] == "prime379_complete_quarter_domain"
    assert old["status"] == "PASS_COMPLETE_DOMAIN"
    assert old["profile_count"] == 5 and old["selector_count"] == 55
    assert old["germ"]["certified"] is True and old["germ"]["cone_candidates"] == 1
    assert fdec(old["collision"]) == tau0
    assert vdec(old["germ"]["new_simple_wall"]) == r
    assert len(old["minimizers"]) == 1 and vdec(old["minimizers"][0]["root"]) == r
    old_basis = [scale(v, -1) for v in native[:7]]+[original_m, r]
    old_gram = tuple(tuple(dot(x, y, gram) for y in old_basis) for x in old_basis)
    assert tuple(tuple(fdec(x) for x in row) for row in old["germ"]["gram"]) == old_gram
    assert [fdec(x) for x in old["germ"]["endpoint_original_ray"]] == list(start)

    unit = 2-PHI
    u3, u = add(m, r3), add(original_m, r)
    assert unit*PHI*PHI == ONE
    assert dot(u3, u3, g3) == 2*unit
    transport = mm(reflection(r3, g3), reflection(u3, g3))
    transport_inverse = mm(reflection(u3, g3), reflection(r3, g3))
    full_transport = mm(reflection(r, gram), reflection(u, gram))
    full_inverse = mm(reflection(u, gram), reflection(r, gram))
    assert mm(transport, transport_inverse) == identity(3) == mm(transport_inverse, transport)
    assert mm(full_transport, full_inverse) == identity(10) == mm(full_inverse, full_transport)
    for operator, full in ((transport, full_transport), (transport_inverse, full_inverse)):
        images = full_basis_images(operator)
        assert all(x.integral() for image in images for x in image)
        assert all(images[j][i] == full[i][j] for i in range(10) for j in range(10))
        assert all(dot(images[i], images[j], gram) == gram[i][j]
                   for i in range(10) for j in range(10))
        assert mv(operator, start3) == start3 and mv(full, start) == start
        assert all(mv(full, v) == v for v in native[:7])
    assert mv(transport, m) == r3
    assert mv(transport, r3) == add(scale(m, -1), r3, -PHI)
    assert mv(transport_inverse, r3) == m
    assert mv(transport_inverse, m) == add(scale(m, -PHI), r3, F(-1))
    fifth, full_fifth = identity(3), identity(10)
    for _ in range(5):
        fifth, full_fifth = mm(transport, fifth), mm(full_transport, full_fifth)
    assert fifth == identity(3) and full_fifth == identity(10)

    line_q = mv(full_transport, original_q)
    line_h = mv(full_transport, (ZERO, ZERO)+H)
    assert add(line_q, line_h, tau0) == start
    assert dot(line_q, start, gram).sign() < 0
    assert dot(line_h, original_m, gram) == -2*z0
    retained = [scale(v, -1) for v in native[:7]]+[r]
    assert all(dot(line_q, w, gram) == ZERO == dot(line_h, w, gram) for w in retained)
    known.append(r)
    assert len(set(known)) == 14
    assert all(dot(start, w, gram).sign() <= 0 for w in known)
    assert {vkey(w) for w in known if dot(start, w, gram) == ZERO} == \
        {vkey(w) for w in old_basis}
    upper_root = mv(full_transport, upper_root)
    assert integral_reflector(upper_root, ONE, gram)
    assert dot(upper_root, line_h, gram) == a10
    assert dot(upper_root, line_q, gram) == -D*b10
    assert dot(upper_root, start, gram).sign() < 0
    assert dot(upper_root, add(line_q, line_h, upper), gram) == ZERO
    return {"g3": g3, "gram": gram, "q": q, "m": m, "A": aa,
            "powers": powers, "inverse_powers": inverse_powers,
            "transport": transport, "transport_inverse": transport_inverse,
            "full_transport": full_transport, "full_inverse": full_inverse,
            "line_q": line_q, "line_h": line_h, "start": tau0, "start_ray": start,
            "starting_wall": r, "upper": upper, "upper_root": upper_root,
            "known": known, "retained": retained}


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
            "norm2_simple": (ONE, ONE, 0),
            "norm2_golden": (ONE, PHI, 0),
            "norm2_e8": (ONE, ZERO, 1),
            "norm2b_orthogonal": (B, ZERO, 0)}[name]


def scalar_rows(name):
    large = name == "norm2b_orthogonal"
    limit, radius = (4, Fraction(8, 5)) if large else (30, Fraction(4))
    for e in range(-limit, limit+1):
        center = e*(PHI-1)
        for c in range(ceil_field(center-radius), floor_field(center+radius)+1):
            scalar = F(c, e)/4
            yield B*scalar if large else scalar, scalar


def inside_scalar_bounds(name, scalar):
    if name == "norm2b_orthogonal":
        return (F(Fraction(9, 4))-scalar).sign() > 0 and \
            (F(Fraction(9, 4))+scalar).sign() > 0 and \
            (F(Fraction(2, 5))-scalar.conjugate()).sign() > 0 and \
            (F(Fraction(2, 5))+scalar.conjugate()).sign() > 0
    return (F(16)-scalar).sign() > 0 and (F(16)+scalar).sign() > 0 and \
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
    gram = geo["gram"]
    native6 = tuple(F(k == 8) for k in range(10))
    native7 = tuple(F(k == 9) for k in range(10))
    # Profile and inverse phase belong to v=F^-1*w, never to w itself.
    inverse_root = mv(geo["full_inverse"], root)
    h = dot(inverse_root, native7, gram)
    if dot(inverse_root, native6, gram) != epsilon:
        return None
    expected = tuple(inverse_root[:2])+tuple(h*x+epsilon*y for x, y in zip(H, J))
    if expected != inverse_root:
        return None
    if dot(inverse_root, original(geo["m"], 0), gram) != -kappa:
        return None
    active = tuple(inverse_root[:2])+(h+Fraction(3, 2)*epsilon,)
    seed = mv(geo["inverse_powers"][phase], active)
    if not all(x.quarter() for x in original(seed, epsilon)):
        return None
    z = seed[2]
    scalar = z/B if d == B else z
    if not (4*scalar).integral() or not inside_scalar_bounds(name, scalar):
        return None
    kk = 2*d-kappa*kappa/2-F(Fraction(3, 2))*epsilon*epsilon
    if (kk*(geo["A"]+1)/4-z*z).sign() < 0:
        return None
    if dot(seed, geo["m"], geo["g3"]) != -kappa:
        return None
    time = -dot(seed, geo["q"], geo["g3"])/(2*D)
    rebuilt = add(add(scale(geo["q"], time), geo["m"], -kappa/2),
                  (ZERO, ZERO, ONE), z)
    if rebuilt != seed:
        return None
    if dot(seed, seed, geo["g3"])+Fraction(3, 2)*epsilon*epsilon != 2*d:
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
    kk = 2*d-kappa*kappa/2-F(Fraction(3, 2))*epsilon*epsilon
    assert kk.sign() > 0
    minimum, ties, selectors = None, {}, []
    for phase in PHASES:
        counts = {key: 0 for key in COUNT_KEYS}
        for z, scalar in scalar_rows(name):
            counts["scalar_rows"] += 1
            if not inside_scalar_bounds(name, scalar):
                continue
            counts["inside_bounds"] += 1
            if (kk*(geo["A"]+1)/4-z*z).sign() < 0:
                continue
            counts["folded_window"] += 1
            square = (2*z*z+kappa*kappa/2+Fraction(3, 2)*epsilon*epsilon-2*d)/(2*D)
            for time in field_roots(square):
                counts["field_roots"] += 1
                seed = add(add(scale(geo["q"], time), geo["m"], -kappa/2),
                           (ZERO, ZERO, ONE), z)
                # Only WHOLE original quarter membership; no restriction on time.
                if not all(x.quarter() for x in original(seed, epsilon)):
                    continue
                counts["quarter_seed"] += 1
                assert dot(seed, seed, geo["g3"])+Fraction(3, 2)*epsilon*epsilon == 2*d
                assert dot(seed, geo["m"], geo["g3"]) == -kappa
                old3 = mv(geo["powers"][phase], seed)
                candidate3 = mv(geo["transport"], old3)
                root = original(candidate3, epsilon)
                if not all(x.integral() for x in root):
                    continue
                counts["original_integral"] += 1
                tau = validate_original(root, d, geo)
                if tau is None:
                    continue
                counts["admitted"] += 1
                minimum, ties = keep_minimum(minimum, ties, tau, root, phase, name)
        cap = 36 if d == B else 549
        assert counts["scalar_rows"] <= cap
        assert counts["field_roots"] <= 2*counts["folded_window"]
        selectors.append({"phase": phase, "profile": name, "counts": counts})
    return {"kind": "prime379_released_m_profile_unit", "status": "PASS_DECLARED_DOMAIN",
            "profile": name, "selectors": selectors,
            "starting_collision": fenc(geo["start"]),
            "upper_collision": fenc(geo["upper"]), "known_wall_count": len(geo["known"]),
            "collision": None if minimum is None else fenc(minimum),
            "minimizers": encode_ties(ties, geo), "germ_certified": False,
            "scope": "Third-edge ONE profile over11 phases only; no complete edge/germ or volume claim"}


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
            "scope": "Actual full finite germ ONLY with unique cone and all55 authenticated selectors"}


def valid_counts(counts, large):
    if not isinstance(counts, dict) or set(counts) != set(COUNT_KEYS):
        return False
    if any(type(x) is not int or x < 0 for x in counts.values()):
        return False
    return 0 < counts["scalar_rows"] <= (36 if large else 549) and \
        counts["inside_bounds"] <= counts["scalar_rows"] and \
        counts["folded_window"] <= counts["inside_bounds"] and \
        counts["field_roots"] <= 2*counts["folded_window"] and \
        counts["quarter_seed"] <= counts["field_roots"] and \
        counts["original_integral"] <= counts["quarter_seed"] and \
        counts["admitted"] <= counts["original_integral"]


def combine(parts, before, geo):
    required = {(phase, name) for name in PROFILES for phase in PHASES}
    seen, profiles, receipts = set(), set(), []
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
        if data.get("kind") != "prime379_released_m_profile_unit" or \
                data.get("status") != "PASS_DECLARED_DOMAIN" or data.get("bindings") != before:
            raise ValueError("not a completed current-source profile unit")
        if data.get("starting_collision") != fenc(geo["start"]) or \
                data.get("upper_collision") != fenc(geo["upper"]) or \
                data.get("known_wall_count") != 14 or data.get("germ_certified") is not False:
            raise ValueError("incorrect third-edge interval, known walls or partial-profile scope")
        if name not in PROFILES or name in profiles:
            raise ValueError("missing/duplicate/invalid profile")
        profiles.add(name)
        selectors = data.get("selectors")
        if not isinstance(selectors, list) or len(selectors) != 11:
            raise ValueError("each profile must contain all11 selectors")
        for selector in selectors:
            phase = selector.get("phase")
            key = (phase, selector.get("profile"))
            if type(phase) is not int or key not in required or key in seen or key[1] != name:
                raise ValueError("invalid or repeated typed selector")
            if not valid_counts(selector.get("counts"), name == "norm2b_orthogonal"):
                raise ValueError("invalid selector counters or domain bound")
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
        raise ValueError("all five profiles and all55 typed selectors are mandatory")
    if minimum is None:
        raise RuntimeError("empty complete domain contradicts the proved finite upper")
    cone = germ(minimum, ties, geo)
    return {"kind": "prime379_released_m_complete_domain",
            "status": "PASS_COMPLETE_DOMAIN" if cone["certified"] else "COMPLETE_DOMAIN_GERM_UNRESOLVED",
            "parts": sorted(receipts, key=lambda row: row["path"]),
            "profile_count": len(profiles), "selector_count": len(seen),
            "collision": fenc(minimum), "minimizers": encode_ties(ties, geo), "germ": cone,
            "starting_collision": fenc(geo["start"]),
            "starting_original_wall": venc(geo["starting_wall"]),
            "upper_collision": fenc(geo["upper"]), "upper_original_mirror": venc(geo["upper_root"]),
            "known_wall_count": len(geo["known"]),
            "scope": "Third edge releasing m ONLY with pinned domain/all55 selectors; no volume/period/example claim"}


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
              "selectors": 55, "maximum_scalar_rows": 24552, "executed_units": 0}, sort_keys=True))
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


