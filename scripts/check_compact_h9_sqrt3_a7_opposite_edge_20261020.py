#!/usr/bin/env python3
"""Complete bounded opposite A7-edge census, with no global volume claim.

Two2D integer boxes and a square dictionary cover all six profile equations.
The finite upper mirror gives DIRECT physical/conjugate bounds; no unit
classification, orbit truncation, floating-point predicate, or height sweep.
This standalone replay uses the same implementation as generation, not a
second independent arithmetic implementation. The root orchestrator supplies
the scientific lock, resource gates, and RSS cap. Never run concurrently.
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
    "results/compact_h9_sqrt3_a7_opposite_interval_20261020.md":
        "6e287ec6d16f9c1700ac717ce605a3ea956afefce1b0768182d711b2d87dcd58",
}
DEFAULT = ROOT / "results/compact_h9_sqrt3_a7_opposite_edge_20261020.json"
MAX_JSON = 256 * 1024
ZERO, ONE, S = (0, 0), (1, 0), (0, 1)
EPS, ETA, M = (2, 1), (1, 1), (4, 1)
D = (-2, 2)
UPPER_Z, UPPER_T, UPPER_X = (11, 8), (5, 3), (10, 6)


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
    assert delta
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
    answer = ZERO
    for value in values:
        answer = add(answer, value)
    return answer


def dot(left, right):
    return sub(summation([mul(a, b) for a, b in zip(left[1:], right[1:])]),
               mul(EPS, mul(left[0], right[0])))


def ambient(root):
    t, zs = root[0], root[1:]
    return [mul((423, 244), t)] + [add(z, mul((185, 107), t)) for z in zs]


def ambient_square(root):
    return sub(summation([sq(y) for y in root[1:]]), mul(S, sq(root[0])))


def file_sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(16384), b""):
            h.update(block)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def upper_check():
    root = [UPPER_T] + [(4, 2)] * 3 + [(3, 2)] * 5 + [EPS]
    assert summation(root[1:]) == mul(M, UPPER_T)
    assert dot(root, root) == (2, 0)
    assert sub(sq(UPPER_Z), mul(D, sq(UPPER_X))) == (9, 0)
    products = [sub(root[i + 2], root[i + 1]) for i in range(8)]
    assert products == [ZERO, ZERO, (-1, 0), ZERO, ZERO, ZERO, ZERO, neg(ETA)]
    return root


def profiles():
    return [
        ("A7+A1_norm2", (2, 0), ZERO, (0,), (144, 0)),
        ("A7+A1_norm2epsilon", scale(2, EPS), ZERO, (0,), scale(144, EPS)),
        ("A8", (2, 0), ONE, (1, 7), (81, 0)),
        ("D8", (2, 0), ONE, (2, 6), (36, 0)),
        ("E8", (2, 0), ONE, (3, 5), (9, 0)),
        ("B8", scale(2, EPS), ETA, (1, 7), scale(18, EPS)),
    ]


def boxes():
    # Bounds are proved in the pinned interval note, not inferred numerically.
    zs = [(a, b) for a in range(-5, 106) for b in range(-3, 71)
          if sign((a, b)) > 0 and sign(sub((a, b), (200, 0))) < 0
          and sign(sub(sq(conj((a, b))), (144, 0))) <= 0]
    xs = [(a, b) for a in range(-2, 85) for b in range(-1, 57)
          if sign((a, b)) > 0 and sign(sub((a, b), (164, 0))) < 0
          and sign(sub(sq(conj((a, b))), (36, 0))) < 0]
    index = {}
    for x in xs:
        index.setdefault(sq(x), []).append(x)
    return zs, xs, index


def interval_solutions(r, zs, square_index):
    answer = []
    for z in zs:
        #1/D=(1+s)/4, exact in K.
        target = divide_integer(mul(ETA, sub(sq(z), r)), 4)
        if target is None:
            continue
        for x in square_index.get(target, []):
            assert sub(sq(z), mul(D, sq(x))) == r
            if sign(sub(scale(9, sq(z)), mul(r, sq(UPPER_Z)))) > 0:
                continue
            if sign(sub(scale(9, sq(x)), mul(r, sq(UPPER_X)))) > 0:
                continue
            if sign(sub(sq(conj(z)), conj(r))) > 0:
                continue
            if sign(sub(mul(neg(conj(D)), sq(conj(x))), conj(r))) > 0:
                continue
            if sign(sub(mul(x, UPPER_Z), mul(UPPER_X, z))) > 0:
                continue
            answer.append((z, x))
    return sorted(answer)


def candidate(z, x, q, g, k, name):
    t = divide_integer(x, 2)
    if t is None:
        return None
    s8 = divide_integer(add(z, scale(8, mul(M, t))), 9)
    if s8 is None:
        return None
    p = divide_integer(sub(s8, scale(k, g)), 8)
    if p is None:
        return None
    z9 = sub(mul(M, t), s8)
    first_eight = [add(p, g)] * k + [p] * (8 - k)
    root = [t] + first_eight + [z9]
    assert summation(root[1:]) == mul(M, t)
    assert dot(root, root) == q
    old_products = [sub(root[i + 2], root[i + 1]) for i in range(8)]
    assert old_products[:7] == [neg(g) if i + 1 == k else ZERO for i in range(7)]
    if any(sign(v) > 0 for v in old_products):
        return None
    #All M/ambient reflection entries are integral because2/q is an O-unit.
    assert integral_quotient((2, 0), q) is not None
    lifted = ambient(root)
    assert ambient_square(lifted) == q
    assert sub(neg(mul(S, mul((1669, 964), lifted[0]))),
               neg(mul((732, 423), summation(lifted[1:])))) == ZERO
    return {"profile": name, "q": q, "gap": g, "k": k,
            "Z": z, "X": x, "T": t, "root_M": root,
            "root_L": lifted, "all_eight_old_products": old_products}


def collision_cmp(left, right):
    return sign(sub(mul(left["T"], right["Z"]), mul(right["T"], left["Z"])))


def simple_dominates(left, right):
    lp, rp = ZERO, ZERO
    for ly, ry in zip(left["root_M"][1:8], right["root_M"][1:8]):
        lp, rp = add(lp, ly), add(rp, ry)
        if sign(sub(mul(right["T"], lp), mul(left["T"], rp))) < 0:
            return False
    return True


def census():
    upper_root = upper_check()
    zs, xs, square_index = boxes()
    summaries, records = [], []
    for name, q, g, placements, r in profiles():
        assert sign(r) > 0 and sign(conj(r)) > 0
        assert sign(sub(r, (540, 0))) < 0 and sign(sub(conj(r), (144, 0))) <= 0
        solutions = interval_solutions(r, zs, square_index)
        accepted = []
        for z, x in solutions:
            for k in placements:
                row = candidate(z, x, q, g, k, name)
                if row is not None:
                    accepted.append(row)
                    records.append(row)
        summaries.append({"profile": name, "R": r, "interval_solutions": len(solutions),
                          "interval_solutions_sha256": digest(solutions),
                          "accepted_count": len(accepted), "accepted_sha256": digest(accepted)})
    assert records
    assert any(row["root_M"] == upper_root for row in records), "Explicit upper root omitted"
    best = records[0]
    for row in records[1:]:
        if collision_cmp(row, best) < 0:
            best = row
    winners = [row for row in records if collision_cmp(row, best) == 0]
    simple = [row for row in winners if all(simple_dominates(row, other) for other in winners)]
    assert simple
    chosen = simple[0]
    assert all(simple_dominates(other, chosen) for other in simple)
    z, t = chosen["Z"], chosen["T"]
    drift = mul((-1, 1), t)
    point = [scale(9, z)] + [add(mul(M, z), drift)] * 8 + [sub(mul(M, z), scale(8, drift))]
    r = sub(sq(z), scale(8, mul((-1, 1), sq(t))))
    assert summation(point[1:]) == mul(M, point[0])
    assert dot(point, point) == scale(-9, mul((-1, 1), r))
    assert sign(dot(point, point)) < 0 and sign(point[0]) > 0
    assert dot(point, chosen["root_M"]) == ZERO
    assert all(sign(dot(point, row["root_M"])) <= 0 for row in records)
    assert all(point[i] == point[1] for i in range(1, 9))
    assert sign(sub(point[9], point[8])) < 0
    return {
        "status": "PASS_COMPLETE_FIRST_A7_OPPOSITE_EDGE_SQRT3_M",
        "scope": "One complete finite first edge; no compactness or nonreflectivity conclusion",
        "coverage": "six primitive-norm/profile rows; direct finite interval; complete two2D boxes",
        "primitive_handling": "a common nonunit factor would give a forbidden unit-norm M root",
        "box_counts": {"Z": len(zs), "X": len(xs)},
        "upper_root_M": upper_root,
        "upper_Z": UPPER_Z, "upper_T": UPPER_T,
        "rows": summaries,
        "all_candidates_count": len(records), "all_candidates_sha256": digest(records),
        "simultaneous_candidates": winners,
        "selected_simple_wall": chosen,
        "simple_wall_projective_representatives": len(simple),
        "endpoint_M_homogeneous": point,
        "endpoint_M_square": dot(point, point),
        "endpoint_L_homogeneous": ambient(point),
        "written_inputs": ["full M norm pool2/2epsilon", "finite A7 spherical extensions",
                           "first simple sector and simultaneous positive-root cone",
                           "explicit upper root and direct interval coefficient bounds"],
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
    before = {str(p.relative_to(ROOT)): file_sha(p) for p in paths}
    result = census()
    assert before == {str(p.relative_to(ROOT)): file_sha(p) for p in paths}
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
                      "first_profile": result["selected_simple_wall"]["profile"],
                      "certificate": str(args.output.relative_to(ROOT))}, sort_keys=True))


if __name__ == "__main__":
    main()
