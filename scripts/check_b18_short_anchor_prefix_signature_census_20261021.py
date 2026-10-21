#!/usr/bin/env python3
"""Complete B18 short-anchor/prefix screen; cubic matches remain necessary only.

Use the bounded parent-lock runner. Exact rational Sturm, field norms and
quadratic K-square tests; no cubic field-isomorphism or chamber claim.
"""
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import permutations
import json
from math import isqrt
import os
from pathlib import Path

import check_a18_isolated_cap_domain_and_norm_screen_20261021 as exact

ROOT = Path(__file__).resolve().parents[1]
SPEC = "results/arithmetic_whole_b18_short_anchor_prefix_signature_census_spec_20261021.md"
PINS = {
    "results/arithmetic_whole_b18_short_cap_inventory_and_degree_three_reduction_20261021.md":
        "0d9cffcc58700e38985e31865380457bbbf885d0c329cb6c53ea69c089f675e1",
    "results/arithmetic_whole_b18_point48_weak_prefix_inventory_and_cap_census_domain_20261021.md":
        "04dfdc05a18b6aa2d7d7ba71a0e6b3824601dce14a797e113a038f209ff003cb",
    "scripts/check_a18_isolated_cap_domain_and_norm_screen_20261021.py":
        "29da7434634c80fb001784d9d5f0d69c3d6247a7b6de1e55e2a0eb6b7ff22932",
}


def poly(degree, t, n, s=None):
    return exact.trim((n, -t, 1) if degree == 2 else (-n, s, -t, 1))


def disc(p):
    if len(p) == 3:
        return int(p[1]*p[1]-4*p[0])
    m, s, t = -p[0], p[1], -p[2]
    return int(t*t*s*s-4*s**3-4*t**3*m-27*m*m+18*t*s*m)


def irreducible(p):
    if len(p) == 3:
        d = disc(p)
        return d > 0 and isqrt(d)**2 != d
    m = abs(int(p[0]))
    assert m
    return all(exact.value(p, Q(r)) and exact.value(p, Q(-r))
               for r in range(1, m+1) if m % r == 0)


def isolate(sequence, a, b):
    count = exact.roots(sequence, a, b)
    if not count:
        return []
    if count == 1:
        return [(a, b)]
    middle = (a+b)/2
    return isolate(sequence, a, middle)+isolate(sequence, middle, b)


def interval_value(p, a, b):
    lo = hi = Q(0)
    for c in reversed(p):
        products = (lo*a, lo*b, hi*a, hi*b)
        lo, hi = min(products)+c, max(products)+c
    return lo, hi


def sign_at(p, sequence, box, expression):
    q = exact.trim(exact.remainder(exact.trim(expression), p))
    if not q:
        return 0, box
    a, b = box
    while True:
        lo, hi = interval_value(q, a, b)
        if lo > 0:
            return 1, (a, b)
        if hi < 0:
            return -1, (a, b)
        middle = (a+b)/2
        if exact.roots(sequence, a, middle):
            b = middle
        else:
            a = middle


def domain_roots(p, physical, other, tests=()):
    sequence = exact.sturm(p)
    d = len(p)-1
    if exact.roots(sequence, *physical) != 1 or exact.roots(sequence, *other) != d-1:
        return None
    boxes = [physical]+isolate(sequence, *other)
    for which, expression, required in tests:
        targets = [0] if which == "physical" else list(range(1, d))
        for j in targets:
            sign, boxes[j] = sign_at(p, sequence, boxes[j], expression)
            if sign != required:
                return None
    return boxes


def add_poly(a, b):
    n = max(len(a), len(b))
    return exact.trim((a[j] if j < len(a) else 0)+(b[j] if j < len(b) else 0)
                      for j in range(n))


def scale_poly(a, c):
    return exact.trim(c*x for x in a)


def mul_poly(a, b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return exact.trim(out)


def spatial_time(projection, norm):
    total = ()
    for i, a in projection:
        for j, b in projection:
            total = add_poly(total, scale_poly(mul_poly(a, b), min(i, j)))
    return add_poly(total, scale_poly(norm, -1))


def permutation_det(a):
    total, n = Q(0), len(a)
    for perm in permutations(range(n)):
        value = Q(-1 if sum(perm[i] > perm[j] for i in range(n) for j in range(i+1, n)) % 2 else 1)
        for i in range(n):
            value *= a[i][perm[i]]
        total += value
    return total


def norm(p, expression):
    d = len(p)-1
    cols = [exact.reduce_mod((Q(0),)*j+tuple(expression), p) for j in range(d)]
    matrix = [[cols[j][i] for j in range(d)] for i in range(d)]
    answer = exact.determinant(matrix)
    assert answer == permutation_det(matrix)
    return answer


def closed_norm(p, kind, index=None):
    d, t = len(p)-1, -p[-2]
    n = p[0] if d == 2 else -p[0]
    s = p[1] if d == 3 else None
    if kind == "prefix_isolated":
        i = index
        return n*(i*i*n-2*i*t+4) if d == 2 else n*(i**3*n-2*i*i*s+4*i*t-8)
    if kind == "short_isolated":
        return n*(81*n-18*t+4)/4 if d == 2 else n*(729*n-162*s+36*t-8)/8
    if kind in ("short_far3", "short_near3"):
        e = 1 if kind == "short_far3" else -1
        if d == 2:
            return (81*n*n+18*e*t*n+40*n-18*t*t-4*e*t+4)/4
        return (729*n*n+162*e*n*s+360*t*n+116*e*n-162*s*s
                -36*e*t*s-80*s+36*t*t+8*e*t-8)/8
    return None


def make_model(kind, index, p, boxes, time, raw_norm, projection):
    p, time = tuple(p), exact.trim(time)
    assert exact.reduce_mod(spatial_time(projection, raw_norm), p) == exact.reduce_mod(time, p)
    n = norm(p, time)
    expected = closed_norm(p, kind, index)
    assert expected is None or n == expected
    assert n and (n > 0) == (len(p) == 4)
    # For rational n/d its squareclass is squarefree(n*d), not numerator alone.
    part, square = exact.squarefree(abs(n.numerator)*n.denominator)
    d = disc(p)
    assert d > 0
    field_part, field_factor = exact.squarefree(d)
    return {"kind": kind, "index": index, "degree": len(p)-1,
            "polynomial": [int(x) for x in p], "discriminant": d,
            "discriminant_squarefree": field_part, "discriminant_square_factor": field_factor,
            "root_boxes_physical_first": [[str(a), str(b)] for a, b in boxes],
            "time_polynomial": [str(x) for x in time], "time_norm": str(n),
            "time_norm_signed_squarefree": part if n > 0 else -part,
            "time_norm_square_factor": square,
            "raw_norm_polynomial": [str(x) for x in raw_norm],
            "raw_projection": [[i, [str(x) for x in coeff]] for i, coeff in projection]}


def anchor_domain(degrees, counts):
    out = []
    for d in degrees:
        boxes = [("short_isolated", t, n, s) for t in range(5, 38)
                 for s in (range(1, 17) if d == 3 else (None,))
                 for n in ((1,) if d == 3 else range(1, 9))]
        boxes += [("short_far3", t, n, s) for t in (range(2, 7) if d == 3 else range(3, 7))
                  for s in (range(-7, 5) if d == 3 else (None,))
                  for n in ((-2, -1, 1, 2) if d == 3 else (-3, -2, -1, 1, 2))]
        boxes += [("short_near3", t, n, s) for t in range(5, 8)
                  for s in (range(-4, 9) if d == 3 else (None,))
                  for n in ((-2, -1, 1, 2) if d == 3 else (-2, -1, 1, 2, 3))]
        assert len(boxes) == (299 if d == 2 else 924)
        for kind, t, n, s in boxes:
            key = f"anchors:{d}:{kind}"
            counts[key+":visited"] += 1
            p = poly(d, t, n, s)
            if not irreducible(p):
                counts[key+":reducible_or_nonreal_quad"] += 1
                continue
            if kind == "short_isolated":
                roots = domain_roots(p, (Q(4), Q(338, 9)), (Q(0), Q(2, 9)))
                time, raw_norm, projection = (0, -1, Q(9, 2)), (0, 1), [(18, (0, Q(1, 2)))]
            elif kind == "short_far3":
                roots = domain_roots(p, (Q(2), Q(25, 4)), (Q(-3, 5), Q(3, 8)), (
                    ("physical", (-8, 0, 1), 1),
                    ("physical", (-3210, 18, 81), -1),
                    ("other", (-1, 1, Q(9, 2)), -1)))
                time, raw_norm, projection = (-1, 1, Q(9, 2)), (2,), [(18, (0, Q(1, 2))), (1, (1,))]
            else:
                roots = domain_roots(p, (Q(4), Q(13, 2)), (Q(-3, 8), Q(3, 5)), (
                    ("physical", (-4, -4, 1), 1),
                    ("physical", (-3210, -18, 81), -1),
                    ("other", (-1, -1, Q(9, 2)), -1)))
                time, raw_norm, projection = (-1, -1, Q(9, 2)), (2,), [(18, (-1, Q(1, 2))), (17, (1,))]
            if roots is None:
                counts[key+":outside_exact_strips"] += 1
                continue
            out.append(make_model(kind, 18, p, roots, time, raw_norm, projection))
            counts[key+":retained"] += 1
    if 2 in degrees:
        p = poly(2, 0, -2)
        out.append(make_model("short_finite_silver", 18, p,
            [(Q(1), Q(2)), (Q(-2), Q(-1))], (24, 17), (2,),
            [(17, (1,)), (18, (0, Q(1, 2)))]))
        counts["anchors:2:finite_typed:visited"] = 1
    return out


def prefix_domain(degrees, counts):
    out = []
    for d in degrees:
        for i in (range(5, 17) if d == 2 else range(5, 8)):
            for t in range(5, 17):
                for s in (range(1, {5: 12, 6: 10, 7: 9}[i]+1) if d == 3 else (None,)):
                    for n in (range(1, (32-1)//i+1) if d == 2 else range(1, 3 if i == 5 else 2)):
                        counts[f"prefix:{d}:isolated:visited"] += 1
                        p = poly(d, t, n, s)
                        if not irreducible(p):
                            counts[f"prefix:{d}:isolated:reducible_or_nonreal_quad"] += 1
                            continue
                        roots = domain_roots(p, (Q(4), Q(16)), (Q(0), Q(2, i)))
                        if roots is None:
                            counts[f"prefix:{d}:isolated:outside_exact_strips"] += 1
                            continue
                        out.append(make_model("prefix_isolated", i, p, roots,
                            (0, -2, i), (0, 2), [(i, (0, 1))]))
                        counts[f"prefix:{d}:isolated:retained"] += 1
        assert counts[f"prefix:{d}:isolated:visited"] == (420 if d == 2 else 516)
    if 2 in degrees:
        def typed(kind, i, p, time, norm, projection):
            sequence = exact.sturm(p)
            roots = [(Q(1), Q(2)), (Q(-2), Q(-1))] if p[1] == 0 else [(Q(1), Q(2)), (Q(-1), Q(0))]
            assert all(exact.roots(sequence, *b) == 1 for b in roots)
            out.append(make_model(kind, i, p, roots, time, norm, projection))
            counts["prefix:2:finite_typed:visited"] += 1
        silver, golden, bronze = poly(2, 0, -2), poly(2, 1, -1), poly(2, 0, -3)
        typed("prefix_own4", 5, silver, (14, 10), (2,), [(5, (0, 1)), (6, (1,))])
        for i in range(5, 11):
            typed("prefix_far3", i, silver, (3*i+1, 2*(i+1)), (2,), [(i, (1, 1)), (1, (1,))])
        for i in (5, 6):
            typed("prefix_golden_far_D", i, golden, (i, i+4), (2,), [(i, (0, 1)), (2, (1,))])
        typed("prefix_golden_near_D", 5, golden, (6, 11), (2,), [(5, (0, 1)), (3, (1,))])
        typed("prefix_golden_far3", 5, golden, (4, 7), (2,), [(5, (0, 1)), (1, (1,))])
        typed("prefix_short12", 16, bronze, (Q(55, 2), 16), (1,),
              [(16, (Q(1, 2), Q(1, 2))), (18, (Q(1, 2),))])
        assert counts["prefix:2:finite_typed:visited"] == 12
    return out


@dataclass(frozen=True)
class K:
    D: int
    a: Q = Q(0)
    b: Q = Q(0)

    def __post_init__(self):
        object.__setattr__(self, "a", Q(self.a)); object.__setattr__(self, "b", Q(self.b))

    def coerce(self, y):
        if not isinstance(y, K):
            return K(self.D, y)
        assert y.D == self.D
        return y

    def __add__(self, y):
        y = self.coerce(y); return K(self.D, self.a+y.a, self.b+y.b)

    __radd__ = __add__

    def __neg__(self):
        return K(self.D, -self.a, -self.b)

    def __mul__(self, y):
        y = self.coerce(y)
        return K(self.D, self.a*y.a+self.D*self.b*y.b, self.a*y.b+self.b*y.a)

    __rmul__ = __mul__

    def norm(self):
        return self.a*self.a-self.D*self.b*self.b

    def __truediv__(self, y):
        y = self.coerce(y); n = y.norm(); assert n
        z = self*K(self.D, y.a, -y.b)
        return K(self.D, z.a/n, z.b/n)

    def sign(self):
        s = lambda x: (x > 0)-(x < 0)
        if not self.a: return s(self.b)
        if not self.b or s(self.a) == s(self.b): return s(self.a)
        return s(self.a)*s(self.norm())

    def data(self):
        return [str(self.a), str(self.b)]


def rational_sqrt(q):
    q = Q(q)
    if q < 0: return None
    a, b = isqrt(q.numerator), isqrt(q.denominator)
    return Q(a, b) if a*a == q.numerator and b*b == q.denominator else None


def field_sqrt(x):
    assert x.sign() > 0
    n = rational_sqrt(x.norm())
    if n is None: return None
    for r in (n, -n):
        a = rational_sqrt((x.a+r)/2)
        if a is None: continue
        if a:
            z = K(x.D, a, x.b/(2*a))
        elif not x.b:
            b = rational_sqrt(x.a/x.D)
            if b is None: continue
            z = K(x.D, 0, b)
        else: continue
        if z*z == x:
            return z if z.sign() > 0 else -z
    return None


def quadratic_time(model):
    D = model["discriminant_squarefree"]
    theta = K(D, Q(-model["polynomial"][1], 2), Q(model["discriminant_square_factor"], 2))
    out = K(D)
    for coefficient in reversed(model["time_polynomial"]):
        out = out*theta+Q(coefficient)
    assert out.sign() > 0 and out.norm() == Q(model["time_norm"])
    return out


def signature(model):
    return (model["degree"], model["discriminant_squarefree"], model["time_norm_signed_squarefree"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", choices=("2", "3", "both"), default="both")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    assert os.environ.get("COXETER_CHECK_OWNER"), "Use bounded parent-lock runner"
    output = (ROOT/args.out).resolve()
    assert output.parent == ROOT/"results" and not output.exists()
    paths = [ROOT/p for p in PINS]+[ROOT/SPEC, Path(__file__).resolve()]
    pins = [{"path": str(p.relative_to(ROOT)), "sha256": exact.digest(p), "bytes": p.stat().st_size} for p in paths]
    assert all(p["sha256"] == PINS[p["path"]] for p in pins if p["path"] in PINS)
    degrees = (2, 3) if args.degree == "both" else (int(args.degree),)
    counts = Counter()
    anchors, prefixes = anchor_domain(degrees, counts), prefix_domain(degrees, counts)
    for models in (anchors, prefixes):
        for i, model in enumerate(models): model["id"] = i
    grouped = defaultdict(list)
    for model in prefixes: grouped[signature(model)].append(model)
    matches = []
    for anchor in anchors:
        candidates = grouped.get(signature(anchor), [])
        retained, rejected = [], []
        for prefix in candidates:
            row = {"prefix": prefix["id"]}
            if anchor["degree"] == 2:
                fa, fb = quadratic_time(anchor), quadratic_time(prefix)
                root = field_sqrt(fb/fa)
                if root is None:
                    rejected.append(prefix["id"])
                    continue
                assert root.sign() > 0 and root*root == fb/fa
                row["positive_time_ratio_K_sqrt"] = root.data()
            else:
                row["status"] = "NECESSARY_SIGNATURE_MATCH_CUBIC_FIELD_MAP_AND_K_SQUARE_UNRESOLVED"
            retained.append(row)
        matches.append({"anchor": anchor["id"], "signature_targets": [p["id"] for p in candidates],
                        "quadratic_nonsquare_targets": rejected, "retained_targets": retained})
    if 2 in degrees:
        # A genuine necessary pair must survive; an empty screen is not expected.
        a = next(m for m in anchors if m["kind"] == "short_finite_silver")
        b = next(m for m in prefixes if m["kind"] == "prefix_far3" and m["index"] == 5)
        fixture = next(m for m in matches[a["id"]]["retained_targets"] if m["prefix"] == b["id"])
        assert fixture["positive_time_ratio_K_sqrt"] == ["-2", "2"]
    result = {"status": "COMPLETE_NECESSARY_B18_ANCHOR_PREFIX_SCREEN_WITH_EXACT_QUADRATIC_TIME_TEST",
        "scope": "All reviewed degree2/3 anchor-prefix domains covered. Cubic matches are necessary discriminant/norm signatures ONLY. Quadratic matches verify same-field time squares ONLY. No mutual Gram/native-angle test, actual chamber or whole-B18 exclusion.",
        "input_pins": pins, "degrees": list(degrees), "counts": dict(sorted(counts.items())),
        "anchors": anchors, "prefixes": prefixes, "matches": matches,
        "anchor_counts": dict(sorted(Counter(str(m["degree"]) for m in anchors).items())),
        "prefix_counts": dict(sorted(Counter(str(m["degree"]) for m in prefixes).items())),
        "retained_pair_counts": {str(d): sum(len(m["retained_targets"]) for m in matches if anchors[m["anchor"]]["degree"] == d) for d in degrees},
        "anchors_with_retained_targets": [m["anchor"] for m in matches if m["retained_targets"]]}
    for path, pin in zip(paths, pins):
        assert exact.digest(path) == pin["sha256"] and path.stat().st_size == pin["bytes"]
    with output.open("x") as f: json.dump(result, f, separators=(",", ":"), sort_keys=True); f.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("anchors", "prefixes", "matches", "input_pins")}))
    print(json.dumps({"output": str(output.relative_to(ROOT)), "sha256": exact.digest(output), "bytes": output.stat().st_size}))


if __name__ == "__main__":
    main()
