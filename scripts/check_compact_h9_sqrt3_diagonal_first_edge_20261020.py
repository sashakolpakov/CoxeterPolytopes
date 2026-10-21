#!/usr/bin/env python3
"""One exact first A8-edge census; no global volume conclusion.

All six norms are representatives modulo SQUARES OF O-UNITS, not six square
classes in K*. In particular2epsilon=(1+sqrt3)^2 is a square in K.

Completeness: the accompanying note proves the full norm/profile inventory,
an integral upper mirror, and reduction by the explicit relative unit U.
Each reduced z=A+B*s lies in -2<=A<=26,-1<=B<=17; each X=C+D*s lies in
-1<=C<=19,-1<=D<=12. These enlarged bounds follow from0<z<48,|z'|<5,
0<=X<35,|X'|<4 and3/2<s<7/4. They cover ALL ten R, since0<R<160 and
0<R'<=18. We index the TWO-dimensional X-square box, rather than taking a
four-dimensional product. Every retained representative needs only U^0..U^4.

No producer arithmetic, PARI, floating-point decisions, maximal-order basis,
fundamental-unit assertion, or exhaustive Vinberg height search is used.
--verify reconstructs the entire finite census; it is a separate replay of
this same implementation, not a second independent arithmetic implementation.
The root orchestrator must supply the global lock/resource gates and RSS cap.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import signal


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "results/compact_h9_sqrt3_diagonal_construction_route_20261020.md"
NOTE_SHA = "16732ee8454f20051c524cbaeb8a56ec8b1b9dbb916f7644d6b08adfd1f60297"
DEFAULT = ROOT / "results/compact_h9_sqrt3_diagonal_first_edge_20261020.json"
MAX_JSON = 256 * 1024
ZERO = (0, 0)
ONE = (1, 0)
S = (0, 1)
EPS = (2, 1)
ETA = (1, 1)


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    return x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def scale(k, x):
    return k * x[0], k * x[1]


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
    assert delta != 0
    return ((delta > 0) - (delta < 0)) * ((a > 0) - (a < 0))


def divide_integer(x, k):
    if x[0] % k or x[1] % k:
        return None
    return x[0] // k, x[1] // k


def integral_quotient(x, y):
    numerator = mul(x, conj(y))
    denominator = mul(y, conj(y))
    assert denominator[1] == 0 and denominator[0] != 0
    return divide_integer(numerator, denominator[0])


def power_step(z, x):
    return add(mul(EPS, z), mul(mul(S, ETA), x)), add(mul(ETA, z), mul(EPS, x))


def file_sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(16384), b""):
            h.update(block)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def reflection_is_integral(root, q):
    dual = [neg(mul(S, root[0]))] + root[1:]
    return all(integral_quotient(scale(2, mul(a, b)), q) is not None
               for a in root for b in dual)


def upper_witness():
    z, x = ONE, ZERO
    for _ in range(9):
        z, x = power_step(z, x)
    assert sub(sq(z), mul(S, sq(x))) == ONE
    assert divide_integer(add(z, ONE), 9) is not None
    a = divide_integer(x, 3)
    b = divide_integer(add(z, ONE), 9)
    assert a is not None and b is not None and sign(a) > 0 and sign(b) > 0
    root = [a] + [b] * 9
    q = scale(2, b)
    assert sub(scale(9, sq(b)), mul(S, sq(a))) == q
    assert reflection_is_integral(root, q)
    return z, x, root, q


def raw_boxes():
    z_box = [(a, b) for a in range(-2, 27) for b in range(-1, 18)
             if sign((a, b)) > 0 and sign(sub((a, b), (48, 0))) < 0
             and sign(sub(sq(conj((a, b))), (25, 0))) < 0]
    x_box = [(a, b) for a in range(-1, 20) for b in range(-1, 13)
             if sign((a, b)) >= 0 and sign(sub((a, b), (35, 0))) < 0
             and sign(sub(sq(conj((a, b))), (16, 0))) < 0]
    index = {}
    for x in x_box:
        index.setdefault(sq(x), []).append(x)
    return z_box, x_box, index


def representatives(r, z_box, square_index):
    answer = []
    for z in z_box:
        residual = sub(sq(z), r)
        if residual[0] % 3:
            continue
        # (A+B*s)/s=B+(A/3)*s.
        target = residual[1], residual[0] // 3
        for x in square_index.get(target, []):
            assert sub(sq(z), mul(S, sq(x))) == r
            if sign(sub(sq(z), mul(sq(EPS), r))) >= 0:
                continue
            if sign(sub(sq(x), mul(sq(ETA), r))) >= 0:
                continue
            if sign(sub(sq(conj(z)), conj(r))) > 0:
                continue
            if sign(sub(mul(S, sq(conj(x))), conj(r))) > 0:
                continue
            answer.append((z, x))
    return sorted(answer)


def profiles():
    norms = [ONE, EPS, (2, 0), scale(2, EPS), (3, 1), mul((3, 1), EPS)]
    rows = []
    for i, q in enumerate(norms):
        rows.append((f"A8+A1_norm_{i}", q, ZERO, (0,), scale(9, q)))
    rows.extend([
        ("A9", (2, 0), ONE, (1, 8), (10, 0)),
        ("D9", (2, 0), ONE, (2, 7), (4, 0)),
        ("B9_norm_1", ONE, ONE, (1, 8), ONE),
        ("B9_norm_2epsilon", scale(2, EPS), ETA, (1, 8), scale(2, EPS)),
    ])
    assert len(rows) == 10
    return rows


def build_root(z, x, q, g, k):
    a = divide_integer(x, 3)
    b = divide_integer(sub(z, scale(k, g)), 9)
    if a is None or b is None or sign(a) <= 0 or sign(b) < 0:
        return None
    ys = [add(b, g)] * k + [b] * (9 - k)
    root = [a] + ys
    total = neg(mul(S, sq(a)))
    for y in ys:
        total = add(total, sq(y))
    assert total == q
    # EVERY original zero wall, including the released ninth wall.
    old_products = [sub(ys[i + 1], ys[i]) for i in range(8)] + [neg(ys[8])]
    assert all(sign(v) <= 0 for v in old_products)
    assert old_products[:8] == [neg(g) if i + 1 == k else ZERO for i in range(8)]
    if not reflection_is_integral(root, q):
        return None
    return root


def collision_cmp(left, right):
    # lambda=s*a/z; all a,z strictly positive.
    return sign(sub(mul(left["a"], right["z"]), mul(right["a"], left["z"])))


def simple_dominates(left, right):
    # a_left*r_right-a_right*r_left=sum d_i*n_i, d_i>=0.
    lp, rp = ZERO, ZERO
    for ly, ry in zip(left["root"][1:9], right["root"][1:9]):
        lp, rp = add(lp, ly), add(rp, ry)
        if sign(sub(mul(right["a"], lp), mul(left["a"], rp))) < 0:
            return False
    return True


def census():
    upper_z, upper_x, upper_root, upper_q = upper_witness()
    z_box, x_box, square_index = raw_boxes()
    records, summaries = [], []
    for name, q, g, placements, r in profiles():
        reps = representatives(r, z_box, square_index)
        accepted = []
        powers_in_interval = 0
        for z0, x0 in reps:
            z, x = z0, x0
            for j in range(5):
                assert sub(sq(z), mul(S, sq(x))) == r
                if sign(sub(mul(x, add(upper_z, ONE)), mul(upper_x, z))) <= 0:
                    powers_in_interval += 1
                    for k in placements:
                        root = build_root(z, x, q, g, k)
                        if root is not None:
                            item = {"profile": name, "q": q, "gap": g, "k": k,
                                    "representative": (z0, x0), "power": j,
                                    "z": z, "X": x, "a": root[0], "root": root}
                            accepted.append(item)
                            records.append(item)
                z, x = power_step(z, x)
        summaries.append({"profile": name, "R": r, "representative_count": len(reps),
                          "representatives_sha256": digest(reps),
                          "powers_in_interval": powers_in_interval,
                          "accepted_root_count": len(accepted),
                          "accepted_roots_sha256": digest(accepted)})
    assert records, "No candidate: mathematical coverage or implementation error"
    best = records[0]
    for item in records[1:]:
        if collision_cmp(item, best) < 0:
            best = item
    winners = [item for item in records if collision_cmp(item, best) == 0]
    simple = [item for item in winners if all(simple_dominates(item, other) for other in winners)]
    assert simple, "No common simple sector among simultaneous candidates"
    chosen = simple[0]
    assert all(simple_dominates(other, chosen) for other in simple)
    assert sign(sub(sq(chosen["z"]), scale(9, mul(S, sq(chosen["a"]))))) > 0
    return {
        "status": "PASS_COMPLETE_FIRST_A8_EDGE_SQRT3_H9",
        "scope": "One complete first finite edge; no finite-volume or compact-polytope claim",
        "norm_equivalence": "squares of O-units for primitive integral normals; not square classes of K*",
        "primitive_handling": "all primitive norm representatives covered; additional integral nonprimitive candidates harmless",
        "box_counts": {"z": len(z_box), "X": len(x_box)},
        "upper_mirror": {"Z9": upper_z, "X9": upper_x, "root": upper_root, "square": upper_q},
        "rows": summaries,
        "all_candidates_count": len(records),
        "all_candidates_sha256": digest(records),
        "simultaneous_candidates": winners,
        "selected_simple_wall": chosen,
        "simple_wall_projective_representatives": len(simple),
        "written_inputs": ["finite rank9 extensions of A8", "full primitive norm pool",
                           "finite mirror arrangement and first simple sector", "U-strip coverage and five-power cutoff"],
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
    assert file_sha(NOTE) == NOTE_SHA, "Route source changed"
    source_paths = [Path(__file__).resolve(), NOTE]
    before = {str(path.relative_to(ROOT)): file_sha(path) for path in source_paths}
    result = census()
    after = {str(path.relative_to(ROOT)): file_sha(path) for path in source_paths}
    assert before == after, "Sources changed during census"
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
        assert saved == json.loads(payload), "Saved certificate differs from complete exact replay"
    signal.alarm(0)
    print(json.dumps({"status": result["status"], "mode": "generation" if args.write else "replay",
                      "candidates": result["all_candidates_count"],
                      "first_profile": result["selected_simple_wall"]["profile"],
                      "certificate": str(args.output.relative_to(ROOT))}, sort_keys=True))


if __name__ == "__main__":
    main()
