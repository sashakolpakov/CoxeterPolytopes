#!/usr/bin/env python3
"""A finite-volume-preserving quartic-field inward block cut in H20.

The new plane makes pi/5 with wall74 and pi/2 with walls2,3. An untouched
actual ridge with cos^2=1/9 proves the resulting polytope is not Coxeter.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np
import sympy as sp

import h20_gluing_independent_escape_20261013 as old

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/fresh_h20_cusp_repair_golden_block_20261013.json"
T = sp.Symbol("t")
MINIMAL = sp.Poly(1475*T**4 - 7080*T**3 + 10476*T**2 - 4752*T - 324, T, domain=sp.QQ)


def evaluate(poly, t):
    answer = F(0)
    for c in poly.all_coeffs():
        answer = answer*t + F(c)
    return answer


def interval():
    low, high = F(-3, 50), F(-59, 1000)
    assert evaluate(MINIMAL, low) > 0 > evaluate(MINIMAL, high)
    assert MINIMAL.count_roots(sp.Rational(low), sp.Rational(high)) == 1
    for _ in range(160):
        mid = (low + high)/2
        if evaluate(MINIMAL, mid) > 0:
            low = mid
        else:
            high = mid
    return low, high


LOW, HIGH = interval()


def sign(expr):
    p = sp.Poly(expr, T, domain=sp.QQ).rem(MINIMAL)
    if p.is_zero:
        return 0
    lower = upper = F(0)
    for (power,), coeff in p.terms():
        c = F(coeff)
        bounds = (c*LOW**power, c*HIGH**power)
        lower += min(bounds)
        upper += max(bounds)
    assert lower > 0 or upper < 0, "Increase rational isolating precision"
    return 1 if lower > 0 else -1


def build():
    Q, f, center, R, _ = old.setup()
    ids = [0, 74, 2, 3]
    gram = R[ids] @ Q @ R[ids].T
    assert gram.tolist() == [[4, -4, 0, 0], [-4, 6, -2, 0], [0, -2, 2, -1], [0, 0, -1, 2]]
    assert int(sp.Matrix(gram).det()) == -8
    assert MINIMAL.is_irreducible
    A, B = 3*R[0], 3*R[74] + 4*R[2] + 2*R[3]
    N = 36 - 72*T + 30*T*T
    assert [int(A @ Q @ A), 2*int(A @ Q @ B), int(B @ Q @ B)] == [36, -72, 30]
    assert A @ Q @ f == B @ Q @ f == 0
    assert sign(T) < 0 and sign(N) > 0
    center_pair = int(A @ Q @ center) + int(B @ Q @ center)*T
    assert center_pair == -48 - 42*T and sign(center_pair) < 0
    # sqrt(5) = numerator/denominator at the chosen root.
    numerator, denominator = 65*T*T - 156*T + 126, 45*T*T - 108*T + 54
    assert sign(numerator) > 0 and sign(denominator) > 0
    assert sign(numerator*numerator - 5*denominator*denominator) == 0
    target_ratios = {2: (0, 1), 3: (1, 4), 4: (1, 2), 6: (3, 4),
                     5: (3*denominator + numerator, 8*denominator),
                     10: (5*denominator + numerator, 8*denominator)}
    pairings_a, pairings_b = R @ Q @ A, R @ Q @ B
    labels, positive, wrong_angles = {}, [], []
    for i in range(580):
        if i == 0:
            continue
        ip = int(pairings_a[i]) + int(pairings_b[i])*T
        d = int(R[i] @ Q @ R[i])*N
        if sign(ip) > 0:
            positive.append(i)
            continue
        if sign(ip*ip-d) >= 0:
            labels[str(i)] = "disjoint_or_parallel"
            continue
        matched = None
        for m, (u, v) in target_ratios.items():
            if sign(v*ip*ip-u*d) == 0:
                matched = m
                break
        if matched is None:
            wrong_angles.append(i)
        else:
            labels[str(i)] = matched
    assert labels["74"] == 5 and labels["2"] == labels["3"] == 2
    old_audit = json.loads((ROOT / "data/enumerations/h20_gluing_independent_audit_20261013.json").read_text())
    x = np.array(old_audit["six_tile_patch"]["non_coxeter_ridge"]["relative_interior_timelike_point"], dtype=object)
    values = R @ Q @ x
    ridge_ids = [i for i, val in enumerate(values) if val == 0]
    assert len(ridge_ids) == 2 and 0 not in ridge_ids and max(values) == 0
    assert x @ Q @ x < 0 and x @ Q @ center < 0 and x @ Q @ f == 0
    wa, wb = int(A @ Q @ x), int(B @ Q @ x)
    assert sign(wa+wb*T) < 0
    H = R[ridge_ids] @ Q @ R[ridge_ids].T
    assert H.tolist() == [[6, -2], [-2, 6]]
    return {"schema": "coxeter.fresh_h20_golden_block.v1", "old_block_indices": ids,
            "old_block_gram": gram.tolist(), "old_block_determinant": -8,
            "t_minimal_polynomial_descending_coefficients": [1475, -7080, 10476, -4752, -324],
            "t_rational_isolating_interval": [str(LOW), str(HIGH)],
            "t_short_certified_interval": ["-3/50", "-59/1000"],
            "new_normal_coefficients_1_t": [list(map(int, A)), list(map(int, B))],
            "new_norm_coefficients_1_t_t_squared": [36, -72, 30],
            "center_pairing_coefficients_1_t": [-48, -42],
            "containment_identity": "r0 = (1/3)V + (-t)r74 + (-4t/3)r2 + (-2t/3)r3; all coefficients positive",
            "pair_gate_declared_finite_alphabet": [2, 3, 4, 5, 6, 10],
            "new_wall_compatible_old_indices_and_labels": labels,
            "new_wall_positive_pairing_indices": positive,
            "new_wall_other_intersecting_angle_indices": wrong_angles,
            "untouched_actual_non_coxeter_ridge": {"boundary_indices": ridge_ids, "gram": H.tolist(),
                "relative_interior_vector": list(map(int, x)), "square": int(x @ Q @ x),
                "new_wall_pairing_coefficients_1_t": [wa, wb], "cosine_squared": "1/9"},
            "conclusion": "Nonempty finite-volume inward-cut polytope with a genuinely new quartic-field plane; exact untouched genuine ridge certifies it is not Coxeter. Pair-gate failures are not a general obstruction to further cuts or redundancy removal."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    data = build()
    if args.verify:
        assert data == json.loads(OUTPUT.read_text())
    else:
        OUTPUT.write_text(json.dumps(data, indent=2) + "\n")
    print("PASS golden inward block: exact pi/5,pi/2,pi/2; finite-volume containment; untouched non-Coxeter ridge")
    print("New-wall pair gate:", len(data["new_wall_compatible_old_indices_and_labels"]), "compatible;",
          len(data["new_wall_positive_pairing_indices"]), "positive pairings;",
          len(data["new_wall_other_intersecting_angle_indices"]), "other intersecting angles")


if __name__ == "__main__":
    main()
