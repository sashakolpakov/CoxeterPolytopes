#!/usr/bin/env python3
"""Complete necessary pair graph for the 44 reviewed silver A18 weak-cap models.

Exact Fraction arithmetic in Q(sqrt(2)); use the bounded parent-lock runner.
This is not a chamber construction, weak-height feasibility test, or exclusion.
"""

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
from itertools import combinations
import json
from math import isqrt
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = "results/arithmetic_whole_a18_silver_weak_cap_pair_graph_spec_20261021.md"
PINS = {
    "results/arithmetic_whole_a18_quadratic_weak_cap_inventory_and_remaining_interface_20261021.md":
        "4710025219a3d8a2d5a485de27d3a9f6962fa61289ef9e6120ee880ee4eec618",
    "results/arithmetic_whole_a18_actual_silver_field_reduction_20261021.md":
        "fb185ba3e154863520a62729348339ebde51443f5f94c77ff0abaa528c96485a",
}


@dataclass(frozen=True)
class K:
    a: Q = Q(0)
    b: Q = Q(0)

    def __post_init__(self):
        object.__setattr__(self, "a", Q(self.a))
        object.__setattr__(self, "b", Q(self.b))

    def __add__(self, y):
        y = field(y)
        return K(self.a + y.a, self.b + y.b)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, y):
        return self + -field(y)

    def __rsub__(self, y):
        return field(y) + -self

    def __mul__(self, y):
        y = field(y)
        return K(self.a*y.a + 2*self.b*y.b, self.a*y.b + self.b*y.a)

    __rmul__ = __mul__

    def norm(self):
        return self.a*self.a - 2*self.b*self.b

    def conj(self):
        return K(self.a, -self.b)

    def __truediv__(self, y):
        y = field(y)
        denominator = y.norm()
        assert denominator
        z = self*y.conj()
        return K(z.a/denominator, z.b/denominator)

    def sign(self):
        if not self.a:
            return rational_sign(self.b)
        if not self.b or rational_sign(self.a) == rational_sign(self.b):
            return rational_sign(self.a)
        return rational_sign(self.a)*rational_sign(self.norm())

    def integral(self):
        return self.a.denominator == self.b.denominator == 1

    def data(self):
        return [str(self.a), str(self.b)]


def field(x):
    return x if isinstance(x, K) else K(x)


def rational_sign(x):
    return (x > 0) - (x < 0)


def rational_sqrt(x):
    x = Q(x)
    if x < 0:
        return None
    p, q = isqrt(x.numerator), isqrt(x.denominator)
    return Q(p, q) if p*p == x.numerator and q*q == x.denominator else None


def positive_field_sqrt(x):
    """Return the physical positive K-square root, or None, exhaustively."""
    assert x.sign() > 0
    r = rational_sqrt(x.norm())
    if r is None:
        return None
    for signed_norm in (r, -r):
        u = rational_sqrt((x.a + signed_norm)/2)
        if u is None:
            continue
        if u:
            z = K(u, x.b/(2*u))
        elif not x.b:
            v = rational_sqrt(x.a/2)
            if v is None:
                continue
            z = K(0, v)
        else:
            continue
        if z*z == x:
            z = z if z.sign() > 0 else -z
            assert z.sign() > 0 and z*z == x
            return z
    return None


def inverse(i, j):
    assert 1 <= i <= 18 and 1 <= j <= 18
    return Q(min(i, j)*(19-max(i, j)), 19)


def spatial(u, v):
    return sum((a*b*inverse(i, j) for i, a in u for j, b in v), K())


def models():
    result = []

    def add(kind, i, projection, norm, helper=None, parameter=None):
        time = spatial(projection, projection)-norm
        assert time.sign() > 0 and time.conj().sign() < 0
        assert norm.sign() > 0 and norm.conj().sign() > 0
        result.append({"id": len(result), "kind": kind, "index": i,
                       "helper": helper, "parameter": parameter,
                       "projection": projection, "norm": norm, "time": time})

    for i in range(5, 15):
        for j in (1, 18):
            add("far3", i, [(i, K(1, 1)), (j, K(1))], K(2), j, K(1, 1))
    for i, j in ((5, 6), (6, 7), (14, 13), (13, 12)):
        add("near4", i, [(i, K(0, 1)), (j, K(1))], K(2), j, K(0, 1))
    for beta in (K(3, 2), K(6, 4)):
        for i in range(5, 15):
            add("isolated_dot", i, [(i, beta)], 2*beta, parameter=beta)
    assert len(result) == 44
    return result


def pair(a, b):
    row = {"models": [a["id"], b["id"]], "indices": [a["index"], b["index"]]}
    if a["index"] == b["index"]:
        return dict(row, reason="same_index")
    ratio = b["time"]/a["time"]
    root = positive_field_sqrt(ratio)
    if root is None:
        return dict(row, reason="time_ratio_not_square_in_K")
    product = a["time"]*root
    assert product.sign() > 0 and product*product == a["time"]*b["time"]
    gram = spatial(a["projection"], b["projection"])-product
    row.update(time_ratio_root=root.data(), positive_time_product=product.data(),
               raw_mutual_product=gram.data())
    if gram.sign() > 0:
        return dict(row, reason="positive_mutual_product")
    if not gram.integral():
        return dict(row, reason="nonintegral_mutual_cycle")
    beta = 4*gram*gram/(a["norm"]*b["norm"])
    row["pair_beta"] = beta.data()
    if not beta.integral():
        return dict(row, reason="nonintegral_two_cycle")
    if gram.sign() == 0:
        return dict(row, reason="compatible", pair_type="orthogonal")
    if beta == K(4):
        return dict(row, reason="parallel_forbidden_by_compactness")
    if (beta-4).sign() > 0:
        return dict(row, reason="compatible", pair_type="dotted")
    angles = {K(1): 3, K(2): 4, K(3): 6, K(2, 1): 8}
    if beta in angles:
        return dict(row, reason="compatible", pair_type="finite", label=angles[beta])
    return dict(row, reason="forbidden_finite_angle")


def components(count, rows):
    neighbors = [set() for _ in range(count)]
    for row in rows:
        if row["reason"] == "compatible":
            i, j = row["models"]
            neighbors[i].add(j)
            neighbors[j].add(i)
    unseen, result = set(range(count)), []
    while unseen:
        start = min(unseen)
        reached, stack = set(), [start]
        while stack:
            i = stack.pop()
            if i not in reached:
                reached.add(i)
                stack.extend(neighbors[i]-reached)
        unseen -= reached
        result.append(sorted(reached))
    return result


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    assert os.environ.get("COXETER_CHECK_OWNER"), "Use bounded parent-lock runner"
    out = (ROOT/args.out).resolve()
    assert out.parent == ROOT/"results" and not out.exists()
    paths = [ROOT/p for p in PINS] + [ROOT/SPEC, Path(__file__).resolve()]
    pins = [{"path": str(p.relative_to(ROOT)), "sha256": digest(p),
             "bytes": p.stat().st_size} for p in paths]
    assert all(p["sha256"] == PINS[p["path"]] for p in pins if p["path"] in PINS)
    caps = models()
    rows = [pair(a, b) for a, b in combinations(caps, 2)]
    assert len(rows) == 946
    # Independently specified genuine local pair: a false empty graph must fail.
    witness = next(r for r in rows if r["models"] == [1, 4])
    assert witness["reason"] == "compatible" and witness["pair_type"] == "dotted"
    assert witness["raw_mutual_product"] == ["-2", "-1"]
    serialized = []
    for cap in caps:
        serialized.append({k: (v.data() if isinstance(v, K) else
            [[i, c.data()] for i, c in v] if k == "projection" else v)
            for k, v in cap.items()})
    result = {"status": "EXACT_COMPLETE_NECESSARY_SILVER_WEAK_CAP_PAIR_GRAPH",
              "scope": "Necessary local arithmetic pair constraints only. No weak-height feasibility, full chamber, whole-A18 exclusion or new dimension ceiling.",
              "input_pins": pins, "models": serialized, "pairs": rows,
              "model_count": len(caps), "pair_count": len(rows),
              "reason_counts": dict(sorted(Counter(r["reason"] for r in rows).items())),
              "compatible_edges": [r["models"] for r in rows if r["reason"] == "compatible"],
              "connected_components": components(len(caps), rows),
              "genuine_silver_pair_preserved": True}
    for path, pin in zip(paths, pins):
        assert digest(path) == pin["sha256"] and path.stat().st_size == pin["bytes"]
    with out.open("x") as stream:
        json.dump(result, stream, separators=(",", ":"), sort_keys=True)
        stream.write("\n")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("models", "pairs", "input_pins")}))
    print(json.dumps({"output": str(out.relative_to(ROOT)), "sha256": digest(out),
                      "bytes": out.stat().st_size}))


if __name__ == "__main__":
    main()
