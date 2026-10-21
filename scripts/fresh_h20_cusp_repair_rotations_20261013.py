#!/usr/bin/env python3
"""Fresh supporting-plane rotations of the H20 mandatory bad ridge.

The new planes are positive combinations of old walls 0 and 74. All pair
decisions and every saved escape are exact; LP and greedy selection only
propose candidates. This is a bounded construction search, not a no-go for
arbitrary new mirrors.
"""

import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import ceil, sqrt
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import h20_gluing_independent_escape_20261013 as original

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/fresh_h20_cusp_repair_rotations_20261013.json"


def compare_rational_to_t(value, label):
    """Sign(value-t), with t=sqrt(2)*cot(pi/label), labels2,3,4,5."""
    if label == 2:
        return (value > 0) - (value < 0)
    if value <= 0:
        return -1
    if label in (3, 4):
        delta = value*value - (F(2, 3) if label == 3 else F(2))
        return (delta > 0) - (delta < 0)
    delta = value*value - 2
    if delta <= 0:
        return -1
    result = delta*delta - F(16, 5)
    return (result > 0) - (result < 0)


def linear_sign(a, b, label):
    a, b = F(a), F(b)
    if b == 0:
        return (a > 0) - (a < 0)
    return -compare_rational_to_t(-a/b, label) * (1 if b > 0 else -1)


INTERVALS = {}


def quadratic_sign(a, b, d, label):
    a, b, d = F(a), F(b), F(d)
    if label in (2, 3, 4):
        square = {2: F(0), 3: F(2, 3), 4: F(2)}[label]
        return linear_sign(a + d*square, b, label)
    if a == b == d == 0:
        return 0
    if label not in INTERVALS:
        low, high = F(1), F(2)
        for _ in range(160):
            mid = (low + high)/2
            if compare_rational_to_t(mid, label) < 0:
                low = mid
            else:
                high = mid
        INTERVALS[label] = low, high
    low, high = INTERVALS[label]
    lower = a + min(b*low, b*high) + min(d*low*low, d*high*high)
    upper = a + max(b*low, b*high) + max(d*low*low, d*high*high)
    assert lower > 0 or upper < 0, "Increase exact isolating interval precision."
    return 1 if lower > 0 else -1


def finite_angle(A, B, norm, N0, N2, label):
    # Test a declared finite alphabet. The target label5 is handled exactly
    # in its quartic field; label8 adds the relevant quadratic-field possibility.
    for m, ratio in ((2, F(0)), (3, F(1, 4)), (4, F(1, 2)), (6, F(3, 4))):
        if quadratic_sign(A*A - norm*N0*ratio, 2*A*B,
                          B*B - norm*N2*ratio, label) == 0:
            return m
    if label == 5 and A == 0 and B*B == norm*N2:
        return 5
    if label == 4:
        # cos^2(pi/8)=1/2+sqrt(2)/4.
        den = norm*(N0 + 2*N2)
        if F(A*A + 2*B*B, den) == F(1, 2) and F(2*A*B, den) == F(1, 4):
            return 8
    return None


def compatible(A, B, norm, N0, N2, label):
    if linear_sign(A, B, label) > 0:
        return False
    if quadratic_sign(A*A - norm*N0, 2*A*B, B*B - norm*N2, label) >= 0:
        return True
    return finite_angle(A, B, norm, N0, N2, label) is not None


def find_escape(ids, normal_a, normal_b, label, Q, f, center, normals, spatial):
    t = {2: 0.0, 3: sqrt(2/3), 4: sqrt(2), 5: sqrt(2 + 4*sqrt(5)/5)}[label]
    old = normals[ids]
    new_float = np.array(normal_a, dtype=float) + t*np.array(normal_b, dtype=float)
    all_float = np.vstack([np.array(old, dtype=float), new_float])
    A = all_float @ np.array(Q, dtype=float) @ spatial
    b = -(all_float @ np.array(Q @ center, dtype=float))
    scale = np.linalg.norm(A, axis=1)
    A, b = A/scale[:, None], b/scale
    rng = np.random.default_rng(2718)
    directions = list(np.eye(20)) + list(-np.eye(20)) + list(rng.normal(size=(20, 20)))
    ca, cb = int(normal_a @ Q @ center), int(normal_b @ Q @ center)
    assert linear_sign(ca, cb, label) < 0
    for objective, direction in enumerate(directions):
        result = linprog(-direction, A_ub=A, b_ub=b, bounds=(None, None), method="highs")
        if result.status == 3:
            result = linprog(-direction, A_ub=A, b_ub=b, bounds=(-10000, 10000), method="highs")
        if not result.success:
            continue
        proposed = np.array(center, dtype=float) + spatial @ result.x
        witness = original.corrected_integer_escape(proposed, Q, f, center, old)
        if witness is None:
            continue
        x = np.array(witness["vector"], dtype=object)
        da, db = int(normal_a @ Q @ x), int(normal_b @ Q @ x)
        shift = max(0, ceil((da + t*db)/(-(ca + t*cb))) + 1)
        while linear_sign(da + shift*ca, db + shift*cb, label) >= 0:
            shift += 1
        x += shift*center
        if x @ Q @ x <= 0:
            continue
        old_values = old @ Q @ x
        assert max(old_values) < 0 and x @ Q @ f == 0 and x @ Q @ center < 0
        da, db = int(normal_a @ Q @ x), int(normal_b @ Q @ x)
        assert linear_sign(da, db, label) < 0
        return {"vector": list(map(int, x)), "square": int(x @ Q @ x),
                "center_pairing": int(x @ Q @ center), "maximum_old_wall_pairing": int(max(old_values)),
                "new_wall_pairing_coefficients_1_t": [da, db], "objective_index": objective}
    return None


def verify(path):
    """Replay every assertion with exact arithmetic; never call an LP solver."""
    data = json.loads(path.read_text())
    Q, f, center, normals, _ = original.setup()
    assert data["common_strict_center"] == list(center)
    assert center @ Q @ center < 0 and center @ Q @ f == 0
    assert max(normals @ Q @ center) < 0
    gram = normals @ Q @ normals.T
    for row in data["results"]:
        label = row["target_angle_label"]
        ids = row["indices_zero_based"]
        assert len(ids) == len(set(ids))
        assert row["wall_count_including_new"] == len(ids) + 1
        aa, bb = (np.array(r, dtype=object) for r in row["new_normal_coefficients_1_t"])
        N0, N2 = row["new_norm_coefficients_1_t_squared"]
        if row["direction"] == "replace_0":
            assert np.array_equal(aa, 3*normals[0] + 2*normals[74])
            assert np.array_equal(bb, -normals[74])
            anchor, excluded = 74, 0
        else:
            assert row["direction"] == "replace_74"
            assert np.array_equal(aa, 2*normals[74] + 2*normals[0])
            assert np.array_equal(bb, -normals[0])
            anchor, excluded = 0, 74
        assert anchor in ids and excluded not in ids
        assert aa @ Q @ f == bb @ Q @ f == aa @ Q @ bb == 0
        assert aa @ Q @ aa == N0 and bb @ Q @ bb == N2
        # 0 <= t < 2: the new wall is a positive combination of the old
        # incompatible pair, hence contains the known six-tile polytope.
        assert compare_rational_to_t(F(2), label) > 0
        assert linear_sign(int(aa @ Q @ center), int(bb @ Q @ center), label) < 0
        pa, pb = normals @ Q @ aa, normals @ Q @ bb
        assert finite_angle(int(pa[anchor]), int(pb[anchor]), int(gram[anchor, anchor]), N0, N2, label) == label
        for i in ids:
            assert compatible(int(pa[i]), int(pb[i]), int(gram[i, i]), N0, N2, label)
        for i, j in combinations(ids, 2):
            ip, ni, nj = int(gram[i, j]), int(gram[i, i]), int(gram[j, j])
            assert ip <= 0
            assert ip*ip >= ni*nj or F(ip*ip, ni*nj) in (F(0), F(1, 4), F(1, 2), F(3, 4))
        witness = row["witness"]
        if witness is None:
            assert row["status"] == "NO_ESCAPE_FOUND_REQUIRES_VOLUME_AUDIT"
            continue
        assert row["status"] == "EXACT_STRICT_SPACELIKE_ESCAPE"
        x = np.array(witness["vector"], dtype=object)
        assert x @ Q @ f == 0
        assert int(x @ Q @ x) == witness["square"] > 0
        assert int(x @ Q @ center) == witness["center_pairing"] < 0
        assert int(max(normals[ids] @ Q @ x)) == witness["maximum_old_wall_pairing"] < 0
        wa, wb = int(aa @ Q @ x), int(bb @ Q @ x)
        assert [wa, wb] == witness["new_wall_pairing_coefficients_1_t"]
        assert linear_sign(wa, wb, label) < 0
    print(f"PASS: {len(data['results'])} exact compatible selections and their saved volume obstructions")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--certificate", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.verify:
        verify(args.certificate)
        return
    Q, f, center, normals, spatial = original.setup()
    source = json.loads(original.CANDIDATES.read_text())
    conflicts = {tuple(sorted(edge)) for edge in source["complete_conflict_edges"]}
    adjacent = [set() for _ in range(580)]
    for i, j in conflicts:
        adjacent[i].add(j)
        adjacent[j].add(i)
    old_gram = normals @ Q @ normals.T
    results = []
    for direction in ("replace_0", "replace_74"):
        if direction == "replace_0":
            aa, bb, N0, N2, anchor, excluded = 3*normals[0] + 2*normals[74], -normals[74], 12, 6, 74, 0
        else:
            aa, bb, N0, N2, anchor, excluded = 2*normals[74] + 2*normals[0], -normals[0], 8, 4, 0, 74
        assert aa @ Q @ aa == N0 and aa @ Q @ bb == 0 and bb @ Q @ bb == N2
        pair_a, pair_b = normals @ Q @ aa, normals @ Q @ bb
        for label in (5, 4, 3, 2):
            available = [i for i in range(580) if i != excluded and
                         compatible(int(pair_a[i]), int(pair_b[i]), int(old_gram[i, i]), N0, N2, label)]
            assert anchor in available
            seen = set()
            for seed in range(6):
                rng = np.random.default_rng(31415 + seed)
                order = list(available)
                rng.shuffle(order)
                order.sort(key=lambda i: len(adjacent[i].intersection(available)))
                chosen = [anchor]
                for i in order:
                    if i not in chosen and all(j not in adjacent[i] for j in chosen):
                        chosen.append(i)
                chosen = sorted(chosen)
                key = tuple(chosen)
                if key in seen:
                    continue
                seen.add(key)
                assert all(tuple(sorted((i, j))) not in conflicts for i, j in combinations(chosen, 2))
                witness = find_escape(chosen, aa, bb, label, Q, f, center, normals, spatial)
                row = {"direction": direction, "target_angle_label": label, "seed": seed,
                    "t_definition": {2: "0", 3: "sqrt(2/3)", 4: "sqrt(2)", 5: "sqrt(2+4*sqrt(5)/5)"}[label],
                    "new_normal_coefficients_1_t": [list(map(int, aa)), list(map(int, bb))],
                    "new_norm_coefficients_1_t_squared": [N0, N2], "available_old_walls": len(available),
                    "indices_zero_based": chosen, "wall_count_including_new": len(chosen) + 1,
                    "all_pair_compatibility": "EXACT_WITH_DECLARED_FINITE_ALPHABET_2_3_4_5_6_8",
                    "status": "EXACT_STRICT_SPACELIKE_ESCAPE" if witness else "NO_ESCAPE_FOUND_REQUIRES_VOLUME_AUDIT",
                    "witness": witness}
                results.append(row)
                print(json.dumps({k: v for k, v in row.items() if k in
                    ("direction", "target_angle_label", "seed", "wall_count_including_new", "available_old_walls", "status")}), flush=True)
    output = {"schema": "coxeter.fresh_h20_targeted_rotations.v1", "source_boundary": str(original.BOUNDARY.relative_to(ROOT)),
              "common_strict_center": list(map(int, center)), "results": results,
              "scope": "Eight genuinely new rotations of one bad ridge, with bounded greedy compatible selections. Exact escape witnesses reject only the recorded candidates. Other combinations, iterative rotations and new planes remain open."}
    args.certificate.write_text(json.dumps(output, indent=2) + "\n")
    verify(args.certificate)


if __name__ == "__main__":
    main()
