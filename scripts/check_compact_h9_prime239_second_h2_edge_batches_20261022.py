#!/usr/bin/env python3
"""PREPARED, UNEXECUTED: exact folded norm239 second-edge units.

Run one --phase/--profile unit through the existing bounded runner with
--parent-lock.  Only a separate --combine of all 20 current units can
certify the complete domain and attempt the actual simple-cone test.
No other worker or scientific package is imported.  No result is hardcoded.
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
SPEC = "results/compact_h9_prime239_folded_second_h2_edge_checker_spec_20261022.md"
PINS = {
    SPEC: "4e070627629fd0e56aa3dff34de6cdd1f424c25a3677271dd9a6f2f7b09a5a29",
    "results/compact_h9_prime239_first_e7_golden_transverse_edge_candidate_20261021.md":
        "2a207fb82414203e34625d323238abc4c297541b6e34f23c4d28b16b7c3efa15",
    "results/compact_h9_prime239_second_h2_edge_finite_upper_contract_candidate_20261021.md":
        "1d11ee290abb97a9da0aaab1d86496071de46d2cf74a02489457fb7f369f89b4",
    "results/compact_h9_prime239_full_base_and_first_binary_edge_candidate_20261021.md":
        "879a3feba630072a4b127e68a3729c4da3fca4f3c44630d9aae623ddd899de02",
    "results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md":
        "b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce",
}
PROFILES = ("orthogonal", "simple", "golden", "e8", "longnorm")
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

    def quarter(self):
        return (4*self).integral()


ZERO, ONE, PHI = F(), F(1), F(0, 1)
A, D, B = ONE + PHI, ONE + 16*PHI, F(16, 17)
E, FF = F(7, -4), F(3, -1)
TT = F(563, 914)
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


def geometry():
    g3 = ((-2*PHI, ONE, ZERO), (ONE, F(8), ZERO), (ZERO, ZERO, F(2)))
    s, u, a, vv = F(21, 33), F(17, 25), F(3, 7), F(3, 5)
    q, m, r0 = (s, -u, ZERO), (a, -vv, ZERO), (ONE, 2*PHI, ZERO)
    rr = (F(13, 22), F(-10, -17), F(1, 2))
    hh = (ZERO, ZERO, ONE)
    nn = add(rr, hh)
    ss = reflection(nn, g3)
    bm = mm(reflection(r0, g3), reflection(m, g3))
    matrix = mm(mm(ss, bm), ss)
    inverse = mm(mm(ss, mm(reflection(m, g3), reflection(r0, g3))), ss)
    powers, inverse_powers = [identity(3)], [identity(3)]
    for _ in range(3):
        powers.append(mm(matrix, powers[-1]))
        inverse_powers.append(mm(inverse, inverse_powers[-1]))
    qq = add(scale(q, 2*A), hh, D)
    ww = add(m, rr, PHI/2)
    assert dot(q, q, g3) == -2*D and dot(m, m, g3) == 2
    assert dot(rr, rr, g3) == 2 and dot(rr, m, g3) == -PHI
    assert dot(qq, qq, g3) == -2*D*E and dot(qq, ww, g3) == 0
    assert dot(ww, ww, g3) == FF/2 and mv(matrix, rr) == rr
    assert sum((matrix[i][i] for i in range(3)), ZERO) == TT + 1
    assert mm(matrix, inverse) == identity(3)
    # Original basis includes e7, not only the integral vector H=Omega7.
    for family in (powers, inverse_powers):
        for j in range(4):
            images = [mv(family[j], (ONE, ZERO, ZERO)),
                      mv(family[j], (ZERO, ONE, ZERO)),
                      scale(add(mv(family[j], hh), hh, F(-1)), F(Fraction(1, 2)))]
            assert all(c.quarter() for image in images for c in image)
            if j == 3:
                assert all(c.integral() for image in images for c in image)
    pp = mv(powers[3], qq)
    lam = -dot(pp, qq, g3)/(2*D*E)
    eta = dot(pp, ww, g3)/(FF/2)
    upper = -eta/lam
    assert lam.sign() > 0 and dot(pp, m, g3).sign() < 0
    assert upper.sign() > 0 and (4*D/A-upper*upper).sign() > 0
    gram = [[ZERO for _ in range(10)] for _ in range(10)]
    gram[0][:2], gram[1][:2] = [-2*PHI, ONE], [ONE, F(8)]
    for i in range(8):
        gram[i+2][i+2] = F(2)
    for i, j in ((0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)):
        gram[i+2][j+2] = gram[j+2][i+2] = F(-1)
    gram = tuple(tuple(row) for row in gram)
    old_native = [tuple(F(-1 if k == j+2 else 0) for k in range(10))
                  for j in range(8)]
    known = old_native + [original(r0, 0), original(m, 0), original(rr, 0)]
    return {"g3": g3, "gram": gram, "q": q, "Q": qq, "W": ww,
            "R": rr, "powers": powers, "inverse_powers": inverse_powers,
            "upper": upper, "half_window2": 4*D/A*(TT-2)/(TT+2),
            "known": known, "retained": old_native[:7] + [original(rr, 0)]}


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
    return {"orthogonal": (0, ZERO, ONE), "simple": (0, ONE, ONE),
            "golden": (0, PHI, ONE), "e8": (1, ZERO, ONE),
            "longnorm": (0, ZERO, B)}[name]


def scalar_rows(name):
    large = name == "longnorm"
    limit, radius = (97, Fraction(86, 25)) if large else (625, Fraction(8))
    for c in range(-limit, limit+1):
        center = c*(PHI-1)
        for aa in range(ceil_field(center-radius), floor_field(center+radius)+1):
            z = F(aa, c)/4
            yield (B*z if large else z), z


def inside_scalar_bounds(name, z):
    if name == "longnorm":
        return (F(53)-z).sign() > 0 and (F(53)+z).sign() > 0 and \
            (F(Fraction(43, 50))-z.conjugate()).sign() > 0 and \
            (F(Fraction(43, 50))+z.conjugate()).sign() > 0
    return (F(342)-z).sign() > 0 and (F(342)+z).sign() > 0 and \
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
    """Recheck the recorded typed folded root, not only its final norm."""
    epsilon, kappa, d = profile(name)
    native6 = tuple(F(k == 8) for k in range(10))
    native7 = tuple(F(k == 9) for k in range(10))
    h = dot(root, native7, geo["gram"])
    if dot(root, native6, geo["gram"]) != epsilon:
        return None
    expected = tuple(root[:2]) + tuple(h*x+epsilon*y for x, y in zip(H, J))
    if expected != root or dot(root, original(geo["R"], 0), geo["gram"]) != -kappa:
        return None
    final3 = tuple(root[:2]) + (h+Fraction(3, 2)*epsilon,)
    seed = mv(geo["inverse_powers"][phase], final3)
    if not all(q.quarter() for q in original(seed, epsilon)):
        return None
    seed_h = seed[2]-Fraction(3, 2)*epsilon
    scalar = seed_h/B if name == "longnorm" else seed_h
    if not scalar.quarter() or not inside_scalar_bounds(name, scalar):
        return None
    qdot, wdot = dot(seed, geo["Q"], geo["g3"]), dot(seed, geo["W"], geo["g3"])
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
              "quarter_seed", "folded_window", "original_integral", "admitted")}
    minimum, ties = None, {}
    for h, scalar in scalar_rows(name):
        counts["scalar_rows"] += 1
        if not inside_scalar_bounds(name, scalar):
            continue
        counts["inside_bounds"] += 1
        lh = 2*A*(2*h+3*epsilon) + 2*kappa/PHI
        m0 = 4*h*h + 12*epsilon*h + 12*epsilon*epsilon - 4*d
        square = D*(lh*lh-(D-1)*m0)
        for square_root in field_roots(square):
            counts["field_roots"] += 1
            t0 = (D*lh+square_root)/(D*(D-1))
            k0 = D*t0-lh
            if not t0.quarter() or not k0.quarter():
                continue
            seed = add(scale(geo["q"], t0/2),
                       (F(3, 7), -F(3, 5), ZERO), -k0/2)
            seed = add(seed, (ZERO, ZERO, ONE), h+Fraction(3, 2)*epsilon)
            seed_original = original(seed, epsilon)
            if not all(q.quarter() for q in seed_original):
                continue
            counts["quarter_seed"] += 1
            assert dot(seed, seed, geo["g3"]) + Fraction(3, 2)*epsilon*epsilon == 2*d
            assert dot(seed, geo["R"], geo["g3"]) == -kappa
            denominator = k0+PHI*kappa/2
            if denominator == ZERO:
                continue
            folded_z = D*(2*A*t0-2*h-3*epsilon)/denominator
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
    cap = 1365 if name == "longnorm" else 21267
    assert counts["scalar_rows"] <= cap
    assert counts["field_roots"] <= 2*counts["inside_bounds"]
    return {"kind": "phase_profile_unit", "status": "PASS_DECLARED_DOMAIN",
            "phase": phase, "profile": name, "counts": counts,
            "collision": None if minimum is None else fenc(minimum),
            "minimizers": [venc(ties[k]) for k in sorted(ties)],
            "germ_certified": False,
            "scope": "One folded unit only; no aggregate firstness or chamber conclusion"}


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
            "scope": "Conditional actual full germ, only if unique cone and all20 current units"}


def combine(paths, before, geo):
    required = {(j, name) for j in range(4) for name in PROFILES}
    seen, receipts, minimum, ties = set(), [], None, {}
    for raw in paths:
        path = (ROOT/raw).resolve()
        if path.parent != ROOT/"results" or not path.is_file() or path.stat().st_size > 2*1024**2:
            raise ValueError("part must be an existing small direct results file")
        rec = record(path)
        with path.open(encoding="utf-8") as stream:
            data = json.load(stream)
        key = (data.get("phase"), data.get("profile"))
        if data.get("kind") != "phase_profile_unit" or data.get("status") != "PASS_DECLARED_DOMAIN":
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
        raise ValueError("all20 phase/profile units are required")
    if minimum is None:
        raise RuntimeError("empty full domain contradicts the proved finite upper")
    cone = germ(minimum, ties, geo)
    return {"kind": "complete_folded_domain", "status": "PASS_COMPLETE_DOMAIN",
            "parts": sorted(receipts, key=lambda q: q["path"]),
            "unit_count": len(seen), "collision": fenc(minimum),
            "minimizers": [venc(ties[k]) for k in sorted(ties)], "germ": cone,
            "scope": "First collision conditional on pinned starting-germ/contract proofs; no volume/period claim"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", action="store_true")
    parser.add_argument("--phase", type=int, choices=range(4))
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
              "units": [{"phase": j, "profile": p} for j in range(4) for p in PROFILES],
              "maximum_scalar_rows": 345732, "executed_units": 0}, sort_keys=True))
        return
    if not os.environ.get("COXETER_CHECK_OWNER"):
        parser.error("use the bounded research runner with --parent-lock")
    if not args.output:
        parser.error("a NEW --output is required")
    destination = output_path(args.output)
    if args.combine:
        if args.phase is not None or args.profile or len(args.part) != 20:
            parser.error("--combine needs exactly20 --part paths and no unit selector")
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
