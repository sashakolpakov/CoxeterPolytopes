#!/usr/bin/env python3
"""PREPARED, UNREVIEWED, UNEXECUTED: norm239 E7+2A1 released-S units.

One --phase/--profile per bounded job, then a separate --combine of ALL15
current units. This NEW 594-row domain uses g=seed.S and both field-square
time signs. The WHOLE original seed is HALF-integral, not its coefficients.
No old worker or scientific package is imported or executed.
Only ORIGINAL-integral reflectors earn candidate credit; no winner is fixed.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELF = Path(__file__).resolve()
SPEC = "results/compact_h9_prime239_e7_two_a1_released_s_checker_spec_20261022.md"
PINS = {
    SPEC: "4815a93b6e22c7b4e97f6ac0ca6b0c142410a80816564f2f960c95b8272b3a13",
    "results/compact_h9_prime239_e7_two_a1_released_s_small_half_domain_candidate_20261022.md":
        "97f61d7c3f4260107c37458c9ae51dd38cd85f90ac33423580bbe5d6295525c5",
    "results/compact_h9_prime239_fourth_edge_and_e7_two_a1_germ_candidate_20261022.md":
        "258a6e99703e4d65ed7d71020fc642edd7c36a644e59acf63bde94ab5e9d49ef",
    "results/prime239_released_r_runtime_manifest_20261022.json":
        "e3e6a54e34ddab2d4213bd355c1f21c6534b275e6496b54d8b0326d9ad3bd85f",
    "results/compact_h9_prime239_full_base_and_first_binary_edge_candidate_20261021.md":
        "879a3feba630072a4b127e68a3729c4da3fca4f3c44630d9aae623ddd899de02",
    "results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md":
        "b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce",
    "scripts/check_compact_h9_prime239_e7_a2_released_r_batches_20261022.py":
        "c2dc88dc5f355e631317c8c1f829c2ec2413e1f742a9ff1aebaab0b8f1e7d380",
}
PROFILES = ("norm2_orthogonal", "norm2_e8", "norm2b_orthogonal",
            "norm2b_simple", "norm2b_golden")
TIE_LIMIT = 128


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
D, B = ONE + 16*PHI, F(16, 17)
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


def geometry():
    """Derive the NEW skew boost, full old-basis powers and finite upper."""
    g3 = ((-2*PHI, ONE, ZERO), (ONE, F(8), ZERO), (ZERO, ZERO, F(2)))
    m, r0 = (F(3, 7), F(-3, -5), ZERO), (ONE, 2*PHI, ZERO)
    rr = (F(13, 22), F(-10, -17), F(1, 2))
    wall_s = (F(55, 90), F(-43, -69), F(6, 10))
    wall_t = (F(12529, 20274), F(-9655, -15623), F(1446, 2343))
    previous = (F(2206, 2368)/131, F(-1370, -2030)/131,
                F(-695, 830)/131)
    qq = add(previous, add(rr, wall_s, F(Fraction(1, 2))),
             -F(10500, -6310)/1441)
    assert qq == (F(21786, 873)/1441, F(-9285, -5310)/1441,
                  F(-5475, 5040)/1441)
    assert dot(wall_s, wall_s, g3) == 2
    assert dot(wall_t, wall_t, g3) == 2*B
    assert dot(wall_s, wall_t, g3) == ZERO
    assert dot(qq, wall_s, g3) == dot(qq, wall_t, g3) == ZERO
    qnorm = dot(qq, qq, g3)
    assert qnorm.sign() < 0 and dot(qq, previous, g3).sign() < 0
    assert B == PHI*D and B.norm() == 239
    assert B.sign() > 0 and B.conjugate().sign() > 0

    aa, yy, hh = wall_t
    skew = ((hh/D, 8*hh/D, (-8*yy-aa)/D),
            (2*PHI*hh/D, -hh/D, (yy-2*PHI*aa)/D),
            (-yy/2, aa/2, ZERO))
    cross = ((ZERO, -hh, yy), (hh, ZERO, -aa), (-yy, aa, ZERO))
    assert mm(g3, skew) == cross
    assert mv(skew, wall_t) == (ZERO, ZERO, ZERO)
    cov_t = mv(g3, wall_t)
    project_t = tuple(tuple(wall_t[i]*cov_t[j]/(2*B) for j in range(3))
                      for i in range(3))
    plane = tuple(tuple(F(i == j)-project_t[i][j] for j in range(3))
                  for i in range(3))
    assert mm(skew, skew) == tuple(tuple(PHI*x for x in row) for row in plane)
    assert mm(project_t, project_t) == project_t
    assert mm(skew, project_t) == tuple((ZERO, ZERO, ZERO) for _ in range(3))
    assert hh/D == F(54, 87)
    assert (-8*yy-aa)/D == F(2407, 3894)
    assert (yy-2*PHI*aa)/D == F(-1867, -3021)
    # FULL original basis, including e7: these operators are ZERO on E7.
    for operator in (identity(3), project_t, skew):
        assert all((2*c).integral()
                   for image in full_basis_images(operator, ZERO) for c in image)

    alpha = dot(mv(skew, qq), wall_s, g3)/2
    assert alpha != ZERO
    if alpha.sign() > 0:
        skew = tuple(tuple(-x for x in row) for row in skew)
    alpha = -dot(mv(skew, qq), wall_s, g3)/2
    assert alpha.sign() > 0
    assert mv(skew, qq) == scale(wall_s, -alpha)
    assert alpha*alpha == -PHI*qnorm/2
    assert mv(skew, wall_s) == scale(qq, -PHI/alpha)
    matrix = tuple(tuple(F(i == j)+(PHI-1)*plane[i][j]+skew[i][j]
                         for j in range(3)) for i in range(3))
    inverse = tuple(tuple(F(i == j)+(PHI-1)*plane[i][j]-skew[i][j]
                          for j in range(3)) for i in range(3))
    assert mm(matrix, inverse) == identity(3)
    assert mv(matrix, wall_t) == wall_t
    assert dot(qq, mv(matrix, qq), g3) == PHI*qnorm
    columns = [tuple(matrix[i][j] for i in range(3)) for j in range(3)]
    assert all(dot(columns[i], columns[j], g3) == g3[i][j]
               for i in range(3) for j in range(3))
    powers, inverse_powers = [identity(3)], [identity(3)]
    for _ in range(2):
        powers.append(mm(matrix, powers[-1]))
        inverse_powers.append(mm(inverse, inverse_powers[-1]))
    for family, sign in ((powers, 1), (inverse_powers, -1)):
        for j in range(3):
            images = full_basis_images(family[j])
            assert all(c.half() for image in images for c in image)
            if j == 2:
                assert all(c.integral() for image in images for c in image)
                expected = tuple(tuple(F(i == k)+2*PHI*plane[i][k]
                                       +sign*2*PHI*skew[i][k]
                                       for k in range(3)) for i in range(3))
                assert family[j] == expected

    pp = mv(powers[2], qq)
    assert pp == add(scale(qq, 1+2*PHI), wall_s, -2*PHI*alpha)
    assert dot(pp, pp, g3) == qnorm and dot(pp, qq, g3).sign() < 0
    upper = 2*PHI*alpha/(1+2*PHI)
    ideal2 = -qnorm/2
    assert upper.sign() > 0 and (ideal2-upper*upper).sign() > 0
    upper_root = scale(mv(powers[2], wall_s), -1)
    assert all(c.integral() for c in original(upper_root, 0))
    assert dot(upper_root, upper_root, g3) == 2
    assert dot(upper_root, pp, g3) == ZERO
    assert dot(upper_root, qq, g3).sign() < 0
    assert dot(upper_root, wall_s, g3).sign() < 0
    assert dot(upper_root, qq, g3)/dot(upper_root, wall_s, g3) == upper

    gram = [[ZERO for _ in range(10)] for _ in range(10)]
    gram[0][:2], gram[1][:2] = [-2*PHI, ONE], [ONE, F(8)]
    for i in range(8):
        gram[i+2][i+2] = F(2)
    for i, j in ((0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)):
        gram[i+2][j+2] = gram[j+2][i+2] = F(-1)
    gram = tuple(tuple(row) for row in gram)
    original_t = original(wall_t, 0)
    cov_full = mv(gram, original_t)
    assert all((c/B).integral() for c in cov_full)
    assert cov_full[:2] == (F(-50203, -81229), F(-64711, -104710))
    old_native = [tuple(F(-1 if k == j+2 else 0) for k in range(10))
                  for j in range(8)]
    known = old_native + [original(r0, 0), original(m, 0),
                          original(rr, 0), original(wall_s, 0), original_t]
    return {"g3": g3, "gram": gram, "Q": qq, "W": wall_s,
            "S": wall_s, "T": wall_t, "negative_square": -qnorm,
            "powers": powers, "inverse_powers": inverse_powers,
            "upper": upper, "half_window2": ideal2*(PHI-1)/(PHI+1),
            "known": known, "retained": old_native[:7] + [original_t]}


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
            raise RuntimeError("changed pinned input: " + name)
        rows[name] = row
    rows[str(SELF.relative_to(ROOT))] = record(SELF)
    return rows


def output_path(raw):
    path = (ROOT/raw).resolve()
    if path.parent != ROOT/"results" or path.suffix != ".json" or path.exists():
        raise ValueError("output must be a NEW direct results/*.json destination")
    return path


def write_new(path, data):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write("\n")


def profile(name):
    return {"norm2_orthogonal": (0, ZERO, ONE),
            "norm2_e8": (1, ZERO, ONE),
            "norm2b_orthogonal": (0, ZERO, B),
            "norm2b_simple": (0, B, B),
            "norm2b_golden": (0, B*PHI, B)}[name]


def scalar_rows(name):
    large = name.startswith("norm2b_")
    limit, radius = (1, Fraction(7, 4)) if large else (4, Fraction(4))
    for c in range(-limit, limit+1):
        center = c*(PHI-1)
        for aa in range(ceil_field(center-radius), floor_field(center+radius)+1):
            z = F(aa, c)/2
            yield (B*z if large else z), z


def inside_scalar_bounds(name, z):
    if name.startswith("norm2b_"):
        return (F(Fraction(1, 2))-z).sign() > 0 and \
            (F(Fraction(1, 2))+z).sign() > 0 and \
            (F(Fraction(7, 8))-z.conjugate()).sign() > 0 and \
            (F(Fraction(7, 8))+z.conjugate()).sign() > 0
    return (F(3)-z).sign() > 0 and (F(3)+z).sign() > 0 and \
        (F(2)-z.conjugate()).sign() >= 0 and (F(2)+z.conjugate()).sign() >= 0


def integral_reflector(v, d, gram):
    if not all(q.integral() for q in v) or dot(v, v, gram) != 2*d:
        return False
    if d == ONE:
        return True  # Original even integral module, exact square2.
    cov = tuple(sum((gram[i][j]*v[j] for j in range(10)), ZERO) for i in range(10))
    return all((x*y/d).integral() for x in v if x != ZERO for y in cov if y != ZERO)


def validate_original(v, d, geo):
    if not integral_reflector(v, d, geo["gram"]):
        return None
    qo, wo = original(geo["Q"], 0), original(geo["W"], 0)
    numerator, derivative = -dot(v, qo, geo["gram"]), -dot(v, wo, geo["gram"])
    if numerator.sign() <= 0 or derivative.sign() <= 0:
        return None
    z = numerator/derivative
    if (geo["upper"]-z).sign() < 0:
        return None
    if any(dot(v, old, geo["gram"]).sign() > 0 for old in geo["known"]):
        return None
    return z


def validate_part_root(root, phase, name, geo):
    """Recover and revalidate the NEW typed scalar and whole folded seed."""
    epsilon, kappa, d = profile(name)
    native6 = tuple(F(k == 8) for k in range(10))
    native7 = tuple(F(k == 9) for k in range(10))
    h = dot(root, native7, geo["gram"])
    if dot(root, native6, geo["gram"]) != epsilon:
        return None
    expected = tuple(root[:2]) + tuple(h*x+epsilon*y for x, y in zip(H, J))
    if expected != root:
        return None
    if dot(root, original(geo["T"], 0), geo["gram"]) != -kappa:
        return None
    final3 = tuple(root[:2]) + (h+Fraction(3, 2)*epsilon,)
    seed = mv(geo["inverse_powers"][phase], final3)
    if not all(q.half() for q in original(seed, epsilon)):
        return None
    g = dot(seed, geo["S"], geo["g3"])
    scalar = g/B if name.startswith("norm2b_") else g
    if not scalar.half() or not inside_scalar_bounds(name, scalar):
        return None
    if dot(seed, geo["T"], geo["g3"]) != -kappa:
        return None
    if dot(seed, seed, geo["g3"]) + Fraction(3, 2)*epsilon*epsilon != 2*d:
        return None
    qdot, wdot = dot(seed, geo["Q"], geo["g3"]), dot(seed, geo["W"], geo["g3"])
    if wdot != g:
        return None
    if wdot == ZERO or (geo["half_window2"]-(qdot/wdot)*(qdot/wdot)).sign() < 0:
        return None
    return validate_original(root, d, geo)


def keep_minimum(minimum, ties, z, root):
    if minimum is None or (z-minimum).sign() < 0:
        minimum, ties = z, {}
    if z == minimum:
        ties[vkey(root)] = root
        if len(ties) > TIE_LIMIT:
            raise RuntimeError("finite local positive-root tie bound exceeded; no truncation")
    return minimum, ties


def unit(phase, name, geo):
    epsilon, kappa, d = profile(name)
    counts = {k: 0 for k in ("scalar_rows", "inside_bounds", "field_roots",
              "half_seed", "folded_window", "original_integral", "admitted")}
    minimum, ties = None, {}
    for g, scalar in scalar_rows(name):
        counts["scalar_rows"] += 1
        if not inside_scalar_bounds(name, scalar):
            continue
        counts["inside_bounds"] += 1
        square = (g*g/2+kappa*kappa/(2*B)
                  +Fraction(3, 2)*epsilon*epsilon-2*d)/geo["negative_square"]
        for time in field_roots(square):
            counts["field_roots"] += 1
            # The rational coefficients need NOT be half-integral.
            # Only the WHOLE ORIGINAL seed is tested below.
            seed = add(scale(geo["Q"], time), geo["S"], g/2)
            seed = add(seed, geo["T"], -kappa/(2*B))
            if not all(q.half() for q in original(seed, epsilon)):
                continue
            counts["half_seed"] += 1
            assert dot(seed, seed, geo["g3"]) + Fraction(3, 2)*epsilon*epsilon == 2*d
            assert dot(seed, geo["S"], geo["g3"]) == g
            assert dot(seed, geo["T"], geo["g3"]) == -kappa
            denominator = g
            assert dot(seed, geo["W"], geo["g3"]) == denominator
            if denominator == ZERO:
                continue
            folded_z = dot(seed, geo["Q"], geo["g3"])/denominator
            if (geo["half_window2"]-folded_z*folded_z).sign() < 0:
                continue
            counts["folded_window"] += 1
            candidate3 = mv(geo["powers"][phase], seed)
            candidate = original(candidate3, epsilon)
            if not all(q.integral() for q in candidate):
                continue
            counts["original_integral"] += 1
            z = validate_original(candidate, d, geo)
            if z is None:
                continue
            counts["admitted"] += 1
            minimum, ties = keep_minimum(minimum, ties, z, candidate)
    cap = 12 if name.startswith("norm2b_") else 81
    assert counts["scalar_rows"] <= cap
    assert counts["field_roots"] <= 2*counts["inside_bounds"]
    return {"kind": "released_s_phase_profile_unit", "status": "PASS_DECLARED_DOMAIN",
            "phase": phase, "profile": name, "counts": counts,
            "collision": None if minimum is None else fenc(minimum),
            "minimizers": [venc(ties[k]) for k in sorted(ties)],
            "germ_certified": False,
            "scope": "One E7+2A1 released-S unit; no aggregate firstness or chamber conclusion"}


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
    gram = geo["gram"]
    p3 = add(geo["Q"], geo["W"], -minimum)
    assert dot(p3, p3, geo["g3"]).sign() < 0
    point = original(p3, 0)
    if any(dot(point, old, gram).sign() > 0 for old in geo["known"]):
        raise RuntimeError("minimum violates a known actual wall")
    chosen = []
    for root in ties.values():
        basis = geo["retained"] + [root]
        gg = tuple(tuple(dot(x, y, gram) for y in basis) for x in basis)
        if not positive_definite(gg):
            continue
        if any(dot(point, v, gram) != ZERO for v in basis):
            continue
        all_inside = True
        for other in ties.values():
            coeff = solve(gg, tuple(dot(x, other, gram) for x in basis))
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
    return {"certified": len(chosen) == 1,
            "cone_candidates": len(chosen),
            "new_simple_wall": venc(chosen[0][0]) if len(chosen) == 1 else None,
            "gram": [[fenc(q) for q in row] for row in chosen[0][1]] if len(chosen) == 1 else None,
            "endpoint_original_ray": [fenc(q) for q in point],
            "scope": "Conditional actual full germ, only if unique cone and all15 current units"}


def combine(paths, before, geo):
    required = {(j, name) for j in range(3) for name in PROFILES}
    seen, receipts, minimum, ties = set(), [], None, {}
    for raw in paths:
        path = (ROOT/raw).resolve()
        if path.parent != ROOT/"results" or not path.is_file() or path.stat().st_size > 2*1024**2:
            raise ValueError("part must be an existing small direct results file")
        rec = record(path)
        with path.open(encoding="utf-8") as stream:
            data = json.load(stream)
        key = (data.get("phase"), data.get("profile"))
        if data.get("kind") != "released_s_phase_profile_unit" or data.get("status") != "PASS_DECLARED_DOMAIN":
            raise ValueError("part is not a completed unit")
        if key not in required or key in seen or data.get("bindings") != before:
            raise ValueError("duplicate/missing/changed-source unit")
        seen.add(key)
        receipts.append(rec)
        z = None if data["collision"] is None else fdec(data["collision"])
        rows = data["minimizers"]
        if len(rows) > TIE_LIMIT or (z is None) != (len(rows) == 0):
            raise ValueError("inconsistent part minimum")
        for row in rows:
            root = vdec(row)
            if validate_part_root(root, key[0], key[1], geo) != z:
                raise ValueError("invalid recorded original reflector or collision")
            minimum, ties = keep_minimum(minimum, ties, z, root)
    if seen != required:
        raise ValueError("all15 phase/profile units are required")
    if minimum is None:
        raise RuntimeError("empty full domain contradicts the proved finite upper")
    cone = germ(minimum, ties, geo)
    return {"kind": "released_s_complete_half_domain", "status": "PASS_COMPLETE_DOMAIN",
            "parts": sorted(receipts, key=lambda q: q["path"]),
            "unit_count": len(seen), "collision": fenc(minimum),
            "minimizers": [venc(ties[k]) for k in sorted(ties)], "germ": cone,
            "scope": "First released-S collision conditional on pinned E7+2A1 germ/contract; no volume/period claim"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", action="store_true")
    parser.add_argument("--phase", type=int, choices=range(3))
    parser.add_argument("--profile", choices=PROFILES)
    parser.add_argument("--combine", action="store_true")
    parser.add_argument("--part", action="append", default=[])
    parser.add_argument("--output")
    args = parser.parse_args()
    before = bindings()
    if args.plan:
        if args.combine or args.phase is not None or args.profile or args.part or args.output:
            parser.error("--plan is standalone")
        print(json.dumps({"status": "PLAN_ONLY", "bindings": before,
              "units": [{"phase": j, "profile": p} for j in range(3) for p in PROFILES],
              "maximum_scalar_rows": 594, "executed_units": 0}, sort_keys=True))
        return
    if not os.environ.get("COXETER_CHECK_OWNER"):
        parser.error("use the bounded research runner with --parent-lock")
    if not args.output:
        parser.error("a NEW --output is required")
    destination = output_path(args.output)
    if args.combine:
        if args.phase is not None or args.profile or len(args.part) != 15:
            parser.error("--combine needs exactly15 --part paths and no unit selector")
    elif args.phase is None or args.profile is None or args.part:
        parser.error("a unit needs --phase and --profile, without --part")
    geo = geometry()
    data = combine(args.part, before, geo) if args.combine else unit(args.phase, args.profile, geo)
    if bindings() != before:
        raise RuntimeError("worker or pinned input changed during unit")
    if args.combine and any(record(ROOT/q["path"]) != q for q in data["parts"]):
        raise RuntimeError("part changed during combine")
    data["bindings"] = before
    data["execution_boundary"] = "caller must separately bind actual capped runner ledgers"
    write_new(destination, data)
    print(json.dumps({"status": data["status"], "output": str(destination.relative_to(ROOT)),
          "output_record": record(destination), "phase": data.get("phase"),
          "profile": data.get("profile"), "counts": data.get("counts"),
          "minimum_ties": len(data["minimizers"]),
          "germ_certified": data.get("germ", {}).get("certified", False)}, sort_keys=True))


if __name__ == "__main__":
    main()

