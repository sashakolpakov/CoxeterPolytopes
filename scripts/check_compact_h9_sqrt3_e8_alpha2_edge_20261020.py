#!/usr/bin/env python3
"""Complete finite interval on one E8 alpha2 edge; no global volume claim.

All eighteen spherical profiles are included. Two separate two-dimensional
integer boxes and a square dictionary exhaust the interval proved in the
pinned note. Cap compatibility is reported, NEVER used as a census filter.
Generation and replay share this implementation. The root coordinator must
supply the global scientific lock, memory/disk gates and external RSS cap.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import signal


ROOT = Path(__file__).resolve().parents[1]
INPUTS = {
    "results/compact_h9_sqrt3_scalar_facet_next_edge_20261020.md":
        "f6be0faf51bdc8188a3f8363fe2c5e3ee35e16c356ba49dd1546574b5195567d",
    "results/compact_h9_sqrt3_e8_integral_splitting_20261020.md":
        "0e8f62f9354fb3be65738db2bb6b8ffcbfdd864d535a70c10292f71d4210bacb",
    "results/compact_h9_sqrt3_e8_first_vinberg_batch_20261020.md":
        "9d71d939945ca043851ffe53b7eea64c7f7abd53067b94ece6f600c537c1b36b",
    "results/compact_h9_sqrt3_twelve_wall_escape_and_completion_20261020.md":
        "a0633db09d09f5de6bcd9d7b149ca05ed2773e801299e1bec5695e483a178f5f",
    "results/compact_h9_sqrt3_e8_alpha2_interval_20261020.md":
        "c19b41db49c82c66e3a2a4b0bc2c14cac3a246af09071e2c0431a106a979edf1",
}
DEFAULT = ROOT / "results/compact_h9_sqrt3_e8_alpha2_edge_20261020.json"
MAX_JSON = 256 * 1024
ZERO, ONE, S = (0, 0), (1, 0), (0, 1)
EPS, ETA, LAMBDA = (2, 1), (1, 1), (-1, 1)
D = (-14, 14)
UPPER_A, UPPER_Z = (23, 13), (73, 42)
PATH = (7, 2, 3, 4, 5, 6)  # beta, alpha3,...,alpha7; zero-based.


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x, y):
    return add(x, neg(y))


def scale(k, x):
    return k * x[0], k * x[1]


def mul(x, y):
    return x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def sq(x):
    return mul(x, x)


def conj(x):
    return x[0], -x[1]


def sign(x):
    a, b = x
    if not a:
        return (b > 0) - (b < 0)
    if not b or (a > 0) == (b > 0):
        return (a > 0) - (a < 0)
    delta = a * a - 3 * b * b
    assert delta, "A nonzero rational pair cannot represent zero in Q(sqrt3)"
    return ((delta > 0) - (delta < 0)) * ((a > 0) - (a < 0))


def divide_integer(x, k):
    if x[0] % k or x[1] % k:
        return None
    return x[0] // k, x[1] // k


def integral_quotient(x, y):
    numerator = mul(x, conj(y))
    denominator = mul(y, conj(y))
    assert denominator[1] == 0 and denominator[0]
    return divide_integer(numerator, denominator[0])


def summation(values):
    result = ZERO
    for value in values:
        result = add(result, value)
    return result


def weight_gram():
    # A7 inverse plus its beta-node Schur update, all exactly integral.
    p = (5, 10, 15, 12, 9, 6, 3)
    answer = [[0] * 8 for _ in range(8)]
    for i in range(7):
        for j in range(7):
            numerator = min(i + 1, j + 1) * (8 - max(i + 1, j + 1)) + p[i] * p[j]
            assert numerator % 8 == 0
            answer[i][j] = numerator // 8
        answer[i][7] = answer[7][i] = p[i]
    answer[7][7] = 8
    return answer


W = weight_gram()


def weight_dot(left, right):
    return summation(scale(W[i][j], mul(left[i], right[j]))
                     for i in range(8) for j in range(8))


def dot(left, right):
    # Coordinates are (a; eight fundamental-weight coefficients).
    return sub(weight_dot(left[1:], right[1:]), mul(LAMBDA, mul(left[0], right[0])))


def root(a, entries):
    coefficients = [ZERO] * 8
    for i, value in entries:
        coefficients[i] = add(coefficients[i], value)
    return [a] + coefficients


def caps():
    return [
        root((5, 3), [(5, neg(EPS))]),
        root((5, 3), [(6, (-3, -2))]),
        root((4, 2), [(0, neg(ETA)), (6, (-1, 0))]),
        root((9, 5), [(0, neg(EPS)), (6, (-3, -2))]),
    ]


def file_sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(16384), b""):
            digest.update(block)
    return digest.hexdigest()


def digest(value):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def profiles():
    # (name,q,A1 gap,path gap,one-based path position,expected R).
    rows = [
        ("none_norm2", (2, 0), ZERO, ZERO, 0, (28, 0)),
        ("none_norm2epsilon", scale(2, EPS), ZERO, ZERO, 0, scale(28, EPS)),
        ("A1_label3", (2, 0), ONE, ZERO, 0, (21, 0)),
        ("A1_label6", (2, 0), S, ZERO, 0, (7, 0)),
        ("A1_label4", scale(2, EPS), ETA, ZERO, 0, scale(14, EPS)),
        ("A1_label12", scale(2, EPS), EPS, ZERO, 0, (7, 0)),
    ]
    for j, rhs in enumerate((16, 8, 4, 4, 8, 16), start=1):
        rows.append((f"path3_j{j}", (2, 0), ZERO, ONE, j, (rhs, 0)))
    for j in (1, 6):
        rows.append((f"path4_j{j}", scale(2, EPS), ZERO, ETA, j, scale(4, EPS)))
    for j, rhs in ((1, 9), (2, 1), (5, 1), (6, 9)):
        rows.append((f"bridge3_j{j}", (2, 0), ONE, ONE, j, (rhs, 0)))
    assert len(rows) == 18
    result = []
    for name, q, a1_gap, path_gap, j, expected in rows:
        u = [ZERO] * 8
        u[0] = a1_gap
        if j:
            u[PATH[j - 1]] = path_gap
        kappa = summation(scale(W[1][i], u[i]) for i in range(8))
        c = weight_dot(u, u)
        rhs = add(sq(kappa), scale(14, sub(q, c)))
        assert rhs == expected
        assert sign(rhs) > 0 and sign(conj(rhs)) > 0
        assert sign(sub(rhs, (105, 0))) < 0
        assert sign(sub(conj(rhs), (28, 0))) <= 0
        result.append({"name": name, "q": q, "U": u, "kappa": kappa, "R": rhs})
    return result


def algebra_checks():
    cartan = [[2 if i == j else 0 for j in range(8)] for i in range(8)]
    for i in range(6):
        cartan[i][i + 1] = cartan[i + 1][i] = -1
    cartan[2][7] = cartan[7][2] = -1
    assert all(sum(cartan[i][k] * W[k][j] for k in range(8)) == (i == j)
               for i in range(8) for j in range(8))
    assert [W[i][i] for i in range(8)] == [4, 14, 30, 20, 12, 6, 2, 8]
    assert W[1] == [7, 14, 20, 16, 12, 8, 4, 10]
    upper = root(UPPER_A, [(1, (-4, -3)), (0, (-1, 0)), (7, (-1, 0))])
    assert dot(upper, upper) == (2, 0)
    assert integral_quotient(UPPER_A, ETA) == (8, 5)
    assert sub(sq(UPPER_Z), mul(D, sq(UPPER_A))) == (9, 0)
    assert [dot(upper, cap) for cap in caps()] == [(-14, -8), (-35, -20), (-3, -2), (-22, -13)]
    for q, cap in zip((scale(2, EPS), scale(2, EPS), (2, 0), scale(2, EPS)), caps()):
        assert dot(cap, cap) == q
        assert sign(sub(sq(dot(upper, cap)), scale(2, q))) > 0
    return upper


def boxes():
    zs = [(z0, z1) for z0 in range(-2, 253) for z1 in range(-1, 149)
          if sign((z0, z1)) > 0 and sign(sub((z0, z1), (500, 0))) < 0
          and sign(sub(sq(conj((z0, z1))), (36, 0))) < 0]
    aa = [(a0, a1) for a0 in range(80) for a1 in range(47)
          if (a0 - a1) % 2 == 0 and sign((a0, a1)) > 0
          and sign(sub((a0, a1), (158, 0))) < 0
          and sign(sub(sq(conj((a0, a1))), ONE)) < 0]
    index = {}
    for a in aa:
        index.setdefault(sq(a), []).append(a)
    return zs, aa, index


def interval_solutions(rhs, zs, square_index):
    result = []
    for z in zs:
        # 1/(14lambda)=eta/28, without field division or approximation.
        target = divide_integer(mul(ETA, sub(sq(z), rhs)), 28)
        if target is None:
            continue
        for a in square_index.get(target, []):
            assert sub(sq(z), mul(D, sq(a))) == rhs
            if sign(sub(scale(9, sq(z)), mul(rhs, sq(UPPER_Z)))) > 0:
                continue
            if sign(sub(scale(9, sq(a)), mul(rhs, sq(UPPER_A)))) > 0:
                continue
            if sign(sub(sq(conj(z)), conj(rhs))) > 0:
                continue
            if sign(sub(scale(14, mul(ETA, sq(conj(a)))), conj(rhs))) > 0:
                continue
            if sign(sub(mul(a, UPPER_Z), mul(UPPER_A, z))) > 0:
                continue
            result.append((z, a))
    return sorted(result)


def candidate(profile, z, a):
    h = divide_integer(sub(z, profile["kappa"]), 14)
    if h is None or sign(h) < 0:
        return None
    assert integral_quotient(a, ETA) is not None
    coefficients = [neg(x) for x in profile["U"]]
    coefficients[1] = sub(coefficients[1], h)
    vector = [a] + coefficients
    assert all(sign(x) <= 0 for x in coefficients)
    assert dot(vector, vector) == profile["q"]
    assert integral_quotient((2, 0), profile["q"]) is not None
    products = [dot(vector, cap) for cap in caps()]
    # These products are evidence only. DO NOT reject positive cap products.
    return {"profile": profile["name"], "q": profile["q"], "R": profile["R"],
            "Z": z, "a": a, "h": h, "root_v_weight_coordinates": vector,
            "all_eight_old_products": coefficients, "four_cap_products": products,
            "compatible_with_all_four_caps": all(sign(x) <= 0 for x in products)}


def collision_cmp(left, right):
    return sign(sub(mul(left["a"], right["Z"]), mul(right["a"], left["Z"])))


def simple_dominates(left, right):
    # left+whites implies right iff a_left*right-a_right*left is white-positive.
    delta = [sub(mul(left["a"], r), mul(right["a"], l))
             for l, r in zip(left["root_v_weight_coordinates"][1:],
                             right["root_v_weight_coordinates"][1:])]
    simple = [summation(scale(W[i][j], delta[j]) for j in range(8)) for i in range(8)]
    assert simple[1] == ZERO  # same collision on the released edge.
    return all(sign(value) >= 0 for i, value in enumerate(simple) if i != 1)


def same_projective_wall(left, right):
    return all(sub(mul(left["a"], r), mul(right["a"], l)) == ZERO
               for l, r in zip(left["root_v_weight_coordinates"][1:],
                               right["root_v_weight_coordinates"][1:]))


def census():
    upper = algebra_checks()
    zs, aa, square_index = boxes()
    summaries, records = [], []
    for profile in profiles():
        solutions = interval_solutions(profile["R"], zs, square_index)
        accepted = []
        for z, a in solutions:
            row = candidate(profile, z, a)
            if row is not None:
                accepted.append(row)
                records.append(row)
        summaries.append({"profile": profile["name"], "q": profile["q"],
                          "R": profile["R"], "kappa": profile["kappa"],
                          "interval_solutions": len(solutions),
                          "interval_solutions_sha256": digest(solutions),
                          "accepted_count": len(accepted), "accepted_sha256": digest(accepted)})
    assert any(row["root_v_weight_coordinates"] == upper for row in records), "Upper mirror omitted"
    best = records[0]
    for row in records[1:]:
        if collision_cmp(row, best) < 0:
            best = row
    winners = [row for row in records if collision_cmp(row, best) == 0]
    simple = [row for row in winners if all(simple_dominates(row, other) for other in winners)]
    unique = bool(simple) and all(same_projective_wall(simple[0], row) for row in simple)
    chosen = simple[0] if unique else None
    z, a = best["Z"], best["a"]
    point = root(z, [(1, neg(mul(LAMBDA, a)))])
    assert dot(point, point) == neg(mul(LAMBDA, best["R"]))
    assert sign(dot(point, point)) < 0 and sign(point[0]) > 0
    assert all(sign(dot(point, row["root_v_weight_coordinates"])) <= 0 for row in records)
    assert all(dot(point, row["root_v_weight_coordinates"]) == ZERO for row in winners)
    assert all(value == ZERO for i, value in enumerate(point[1:]) if i != 1)
    assert sign(point[2]) < 0
    point_caps = [dot(point, cap) for cap in caps()]
    assert all(sign(value) <= 0 for value in point_caps)
    status = ("PASS_COMPLETE_FIRST_E8_ALPHA2_EDGE_SQRT3_M" if unique
              else "UNKNOWN_SIMULTANEOUS_SIMPLE_WALL_NOT_UNIQUE")
    return {
        "status": status,
        "scope": "One complete finite first edge; no global compactness or nonreflectivity conclusion",
        "coverage": "All eighteen spherical/norm profiles; direct upper interval and two2D boxes",
        "cap_handling": "No cap-product census filter; compatibility reported after enumeration",
        "primitive_handling": "Nonprimitive norm2/2epsilon would give a forbidden unit-norm M vector",
        "box_counts": {"Z": len(zs), "a": len(aa)},
        "upper_root_v_weight_coordinates": upper, "upper_Z": UPPER_Z, "upper_a": UPPER_A,
        "rows": summaries, "all_candidates_count": len(records),
        "all_candidates_sha256": digest(records), "all_candidates": records,
        "simultaneous_candidates": winners, "simple_dominating_candidates": simple,
        "selected_simple_wall": chosen,
        "explicit_upper_is_first_simple_wall": bool(chosen and chosen["root_v_weight_coordinates"] == upper),
        "endpoint_v_weight_coordinates_homogeneous": point,
        "endpoint_square": dot(point, point), "endpoint_four_cap_products": point_caps,
        "written_inputs": ["Full M norm pool and a in eta O", "Complete finite A1+A6 extensions",
                           "All-arrangement first-wall and finite local positive-root interpretation",
                           "Explicit upper mirror and direct coefficient bounds"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT)
    args = parser.parse_args()
    if args.write and args.output.exists():
        raise SystemExit("Refusing to overwrite existing evidence")
    def expired(_signum, _frame):
        raise TimeoutError("45-second internal wall cap")
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(45)
    for relative, expected in INPUTS.items():
        assert file_sha(ROOT / relative) == expected, f"Input changed: {relative}"
    paths = [Path(__file__).resolve()] + [ROOT / p for p in INPUTS]
    before = {str(path.relative_to(ROOT)): file_sha(path) for path in paths}
    result = census()
    assert before == {str(path.relative_to(ROOT)): file_sha(path) for path in paths}
    result["sources"] = before
    payload = json.dumps(result, sort_keys=True, indent=2) + "\n"
    assert len(payload.encode()) <= MAX_JSON
    if args.write:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(payload)
    else:
        assert args.output.stat().st_size <= MAX_JSON
        with args.output.open(encoding="utf-8") as stream:
            saved = json.load(stream)
        assert saved == json.loads(payload), "Saved evidence differs from complete exact replay"
    signal.alarm(0)
    print(json.dumps({"status": result["status"], "mode": "generation" if args.write else "replay",
                      "candidates": result["all_candidates_count"],
                      "first_profile": (result["selected_simple_wall"] or {}).get("profile"),
                      "upper_is_first": result["explicit_upper_is_first_simple_wall"],
                      "certificate": str(args.output.relative_to(ROOT))}, sort_keys=True))
    if not result["status"].startswith("PASS_"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
