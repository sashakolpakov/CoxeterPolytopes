#!/usr/bin/env python3
"""Exact scalar certificate only. Root schedules all execution under its guard.

No optimizer, geometry census, scientific imports, or inner flock. The parent
wrapper owns the shared lock/resource gates. --write refuses existing output.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import signal

ROOT = Path(__file__).resolve().parents[1]
PINS = {
    "results/right_angled_mixed_positive_part_opposition_cap_research_20261019.md": "19ec3d814759cd2d9b9b1df8ef128f1df7d7afee68590aa3310479c954a65c9f",
    "results/right_angled_mixed_opposition_moments_research_20261018.md": "a95079d23556e82bcfbb369abccfb751fe96655ff03d95682b0a7e0c8f6bdd7e",
    "results/right_angled_mixed_opposition_power_two_research_20261018.md": "ab7ed2a9e9c04947538589aeee11accb238551d20fca017e58aa9bd892701a2f",
    "results/right_angled_mixed_polygon_three_face_coupling_research_20261018.md": "f2949232322ec28ba5c2d10a97a30a0a47a519c3dc6061e88dd95c4cf0b31a41",
    "results/right_angled_mixed_weighted_triangles_research_20261018.md": "2e9fac6d7c82493e8e7e1e0172235cc3b910f5f3b05d7ded0d7056bdae217a60",
}
DEFAULT = ROOT / "results/right_angled_mixed_positive_part_opposition_cap_research_20261019.json"


def digest(path):
    h = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def add(*polys):
    out = [F(0)] * max(map(len, polys))
    for p in polys:
        for i, value in enumerate(p):
            out[i] += value
    while len(out) > 1 and not out[-1]:
        out.pop()
    return out


def scale(p, value):
    return [value * x for x in p]


def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return add(out)


def deriv(p, order=1):
    for _ in range(order):
        p = [i * p[i] for i in range(1, len(p))] or [F(0)]
    return p


def at(p, x):
    out = F(0)
    for value in reversed(p):
        out = out * x + value
    return out


def bernstein(p, a, b):
    m, h = len(p) - 1, b - a
    shifted = [sum(p[k] * comb(k, j) * a ** (k-j) * h**j
                   for k in range(j, m+1)) for j in range(m+1)]
    return [sum(shifted[k] * F(comb(i, k), comb(m, k))
                for k in range(i+1)) for i in range(m+1)]


def symadd(*terms):
    out = Counter()
    for coefficient, p in terms:
        for key, value in p.items():
            out[key] += coefficient * value
    return {key: value for key, value in out.items() if value}


def symmul(p, q):
    out = Counter()
    for key, value in p.items():
        for other, weight in q.items():
            out[tuple(sorted(key + other))] += value * weight
    return dict(out)


def moments():
    """Expand all subset monomials before using sum x_i=0, sum x_i^2=1/4."""
    s = {k: {(i,) * k: F(1) for i in range(12)} for k in range(1, 5)}
    prob = {j: F(comb(12-j, 3-j), comb(12, 3)) for j in range(1, 4)}
    a, b, c = prob[1], prob[2], prob[3]
    s11 = symmul(s[1], s[1])
    expected = {
        1: symadd((-a, s[1])),
        2: symadd((a-b, s[2]), (b, s11)),
        3: symadd((-(a-3*b+2*c), s[3]),
                  (-3*(b-c), symmul(s[1], s[2])), (-c, symmul(s11, s[1]))),
        4: symadd((a-7*b+12*c, s[4]),
                  (4*b-12*c, symmul(s[1], s[3])),
                  (3*b-6*c, symmul(s[2], s[2])),
                  (6*c, symmul(s11, s[2]))),
    }
    counts = {}
    for degree in range(1, 5):
        expanded = Counter()
        for subset in combinations(range(12), 3):
            for ordered in product(subset, repeat=degree):
                expanded[tuple(sorted(ordered))] += F((-1)**degree, 220)
        assert dict(expanded) == expected[degree]
        counts[str(degree)] = len(expanded)
    m2 = (a-b) / 4
    neg_m3 = (a-3*b+2*c) / 8
    m4 = (3*b-6*c) / 16
    assert (m2, neg_m3, m4) == (F(9,176), F(27,1760), F(3,440))
    return {"inclusion_probabilities": prob, "expanded_monomial_counts": counts,
            "second": m2, "negative_third_upper": neg_m3, "fourth_upper": m4}


def run():
    pins = {name: digest(ROOT / name) for name in PINS}
    assert pins == PINS, "Changed source dependency"
    p = [F(19,125), -F(979,1000), F(768,625), -F(971,3125), F(53186,9375)]
    r = [F(5,4), -F(1)]
    p_neg = [value * (-1)**i for i, value in enumerate(p)]
    n = add(mul(mul(r, r), add(p_neg, [F(1,2)])), [-F(1)])
    assert n == [F(3,160), -F(321,3200), F(249,2000), -F(643,400),
                 F(27949,3000), -F(130052,9375), F(53186,9375)]
    p0 = [F(7,50), -F(128,125), F(768,625), -F(4096,3125), F(13312,3125)]
    assert add(p, scale(p0, -1)) == [F(3,250), F(9,200), F(0), F(1), F(106,75)]
    r_plus = [F(5,4), F(1)]
    lhs = add(mul(mul(r_plus, r_plus), add(p0, [F(1,2)])), [-F(1)])
    rhs = mul([F(0)]*4 + [F(1)/(F(5,4)**5)],
              add(scale(mul(r_plus, r_plus), 13), scale(r_plus, -4), [-F(5,4)]))
    assert lhs == rhs
    p_mid, p_slope = at(p, F(1,4)), at(deriv(p), F(1,4))
    assert p_mid == F(1627,1200000) > F(13,10000)
    assert p_slope == -F(10243,150000) and abs(p_slope) < F(7,100)
    assert F(36,625)-F(1,60) > 0 and F(13,10000)-F(49,48000) == F(67,240000) > 0
    q2 = add(deriv(p, 2), [-F(12,5)])
    assert q2 == [F(36,625), -F(5826,3125), F(212744,3125)]
    assert q2[1] >= -2 and q2[2] >= 60
    tables = []
    specs = [(1,F(0),F(1,5),[-F(1,10),-F(9,100),-F(99,1000),-F(98,1000),-F(79,1000),-F(9,200)], "upper"),
             (1,F(3,10),F(1,2),[F(3,50),F(11,100),F(17,100),F(21,100),F(23,100),F(1,5)], "lower"),
             (2,F(1,2),F(3,4),[-F(3,5),-F(19,10),-F(37,10),-F(29,5),-F(36,5)], "upper")]
    for order, a, b, bound, side in specs:
        values = bernstein(deriv(n, order), a, b)
        assert all(x < y if side == "upper" else x > y for x, y in zip(values, bound))
        tables.append({"derivative": order, "interval": [a,b], "coefficients": values,
                       "strict_bounds": bound, "side": side})
    center = {"N(1/5)": at(n,F(1,5)), "N4(1/5)": at(deriv(n,4),F(1,5)),
              "N4(3/10)": at(deriv(n,4),F(3,10)), "N3(3/10)": at(deriv(n,3),F(3,10)),
              "N2(1/5)": at(deriv(n,2),F(1,5)), "N(1/4)": at(n,F(1,4)),
              "N1(1/4)": at(deriv(n),F(1,4)), "N(3/4)": at(n,F(3,4))}
    assert center["N(1/5)"] > F(1,1000) and center["N4(1/5)"] < -27
    assert center["N4(3/10)"] < -91 and center["N3(3/10)"] > F(9,10)
    assert center["N2(1/5)"] > F(4,5) and deriv(n,6)[0] > 0
    assert center["N(1/4)"] == F(679,1200000) > F(1,2000)
    assert center["N1(1/4)"] == F(3061,600000) and abs(center["N1(1/4)"]) < F(1,100)
    assert F(1,2000)-F(1,100)**2 > 0 and center["N(3/4)"] == F(57,64000) > 0
    m = moments()
    average = p[0]+p[2]*m["second"]-p[3]*m["negative_third_upper"]+p[4]*m["fourth_upper"]
    assert p[2] > 0 > p[3] and p[4] > 0 and average == F(1420561,5500000)
    cap = 440*average
    lower = F(2*(7870*880+50*267300), 101*880+267300)
    gap = lower-cap
    assert cap == F(1420561,12500) and lower == F(184460,1619) and gap > 0
    assert 184460*12500-1420561*1619 == 5861741
    return {"status":"PASS", "scope":"exact polynomial/moment coefficients; written geometry and transport inherited",
            "owns_lock":False, "inputs":pins, "checker_sha256":digest(Path(__file__)),
            "P":p, "N_independent_expansion":n, "bernstein":tables, "central_values":center,
            "P_midpoint":[p_mid,p_slope], "moments":m, "average_upper":average,
            "triangle_cap":cap, "strict_triangle_lower":lower, "gap":gap}


def canonical(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key):canonical(item) for key,item in value.items()}
    if isinstance(value, (list,tuple)):
        return [canonical(item) for item in value]
    return value


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT)
    args = parser.parse_args()
    signal.alarm(45)
    data = canonical(run())
    if args.write:
        with args.output.open("x", encoding="utf-8") as stream:
            json.dump(data, stream, indent=2, sort_keys=True)
            stream.write("\n")
    else:
        assert json.loads(args.output.read_text()) == data, "Saved certificate mismatch"
    print("PASS positive-part scalar certificate; cap="+data["triangle_cap"]+"; lower="+data["strict_triangle_lower"])
