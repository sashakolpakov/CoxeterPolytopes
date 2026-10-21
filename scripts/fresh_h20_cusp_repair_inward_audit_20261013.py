#!/usr/bin/env python3
"""Independent exact source-anchored replay of the inward-cut archive.

No geometry or arithmetic routine is imported from the construction scripts.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/fresh_h20_mixed_inward_20261013.json"
OUTPUT = ROOT / "data/enumerations/fresh_h20_cusp_repair_inward_audit_20261013.json"


def sign(a, b=0):
    a, b = F(a), F(b)
    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    d = a*a - 2*b*b
    return ((d > 0) - (d < 0)) * (1 if a > 0 else -1)


def multiply(a, b):
    return a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0]


def inner(a, b, Q):
    return (int(a[0] @ Q @ b[0] + 2*a[1] @ Q @ b[1]),
            int(a[0] @ Q @ b[1] + a[1] @ Q @ b[0]))


def recombine(terms, R):
    answer = np.zeros((2, 22), dtype=object)
    for index, sa, sb in terms:
        a, b = F(sa), F(sb)
        assert sign(a, b) > 0
        answer[0] += a*R[index, 0] + 2*b*R[index, 1]
        answer[1] += a*R[index, 1] + b*R[index, 0]
    return answer


def pair_label(p, u, v):
    if sign(*p) > 0:
        return None
    numerator, denominator = multiply(p, p), multiply(u, v)
    if sign(numerator[0]-denominator[0], numerator[1]-denominator[1]) >= 0:
        return "disjoint_or_parallel"
    for label, ratio in ((2, F(0)), (3, F(1, 4)), (4, F(1, 2)), (6, F(3, 4))):
        if numerator == tuple(ratio*x for x in denominator):
            return label
    if tuple(4*x for x in numerator) == multiply((2, 1), denominator):
        return 8
    return None


def build():
    d = json.loads(SOURCE.read_text())
    raw = json.loads((ROOT / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json").read_text())["coxeter_input"]
    Q = -np.array([[int(x) for x in row] for row in raw["basis_gram"]], dtype=object)
    f = np.array(list(map(int, raw["wall_coordinates_as_columns_transposed"][0])), dtype=object)
    authoritative = json.loads((ROOT / "data/enumerations/new_dimension_h20_six_patch_boundary_20261013.json").read_text())["six_tile"]["primitive_outward_normals"]
    R = np.array(d["normal_pairs"], dtype=object)
    assert R[:580, 0].tolist() == authoritative and not any(R[:580, 1].flatten())
    center = np.array(d["common_strict_timelike_center"], dtype=object)
    source_center = json.loads((ROOT / "data/enumerations/h20_gluing_independent_escape_20261013.json").read_text())["common_strict_timelike_center"]
    assert center.tolist() == source_center
    assert center @ Q @ center < 0 and center @ Q @ f == 0
    assert max(R[:580, 0] @ Q @ center) < 0
    C = np.array([center, np.zeros(22, dtype=object)], dtype=object)
    active = list(range(580))
    steps = []
    for number, row in enumerate(d["steps"]):
        i, j, k = row["replaced_index"], row["retained_index"], row["new_index"]
        assert i in active and j in active and i != j and k == 580 + number
        assert np.array_equal(R[i], recombine([[k, *row["new_coefficient"]], [j, *row["retained_coefficient"]]], R))
        norm, ip = inner(R[k], R[k], Q), inner(R[k], R[j], Q)
        assert sign(*norm) > 0 and sign(*inner(R[k], C, Q)) < 0
        assert all(v == 0 for v in R[k] @ Q @ f)
        assert pair_label(ip, norm, inner(R[j], R[j], Q)) == row["target_angle_m"]
        active.remove(i)
        active.append(k)
        for deletion in row["redundancy_proofs"]:
            removed = deletion["removed_index"]
            assert removed in active
            assert all(t[0] in active and t[0] != removed for t in deletion["terms"])
            assert np.array_equal(R[removed], recombine(deletion["terms"], R))
            active.remove(removed)
        assert len(active) == row["retained_inequality_count"]
        steps.append({"replaced": i, "retained": j, "new": k, "target_angle_m": row["target_angle_m"],
                      "new_norm": list(norm), "center_pairing": list(inner(R[k], C, Q)),
                      "exact_redundancy_deletions": len(row["redundancy_proofs"]), "inequalities": len(active)})
    assert active == d["final_active_indices"]
    A, B = R[active, 0], R[active, 1]
    H0, H1 = A @ Q @ A.T + 2*B @ Q @ B.T, A @ Q @ B.T + B @ Q @ A.T
    conflicts = 0
    for i in range(len(active)):
        for j in range(i):
            if pair_label((int(H0[i,j]), int(H1[i,j])), (int(H0[i,i]), int(H1[i,i])), (int(H0[j,j]), int(H1[j,j]))) is None:
                conflicts += 1
    assert conflicts == d["final_conflict_count"]
    rational_ids = [i for i in active if not any(R[i, 1])]
    rational_rank = int(sp.Matrix(R[rational_ids, 0]).rank())
    return {"source": str(SOURCE.relative_to(ROOT)), "source_status": d["status"], "steps": steps,
            "center_square": int(center @ Q @ center), "final_inequalities": len(active),
            "final_conflicts_in_declared_alphabet": conflicts, "final_rational_normal_rank": rational_rank,
            "conclusion": "Every archived region is full-dimensional and contained in the original finite-volume patch. All removals satisfy exact positive cone identities. The final inequalities still fail the complete pair-compatibility gate; no Coxeter claim is made."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    result = build()
    if args.verify:
        assert result == json.loads(OUTPUT.read_text())
    else:
        OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS independent inward replay:", len(result["steps"]), "cuts;", result["final_conflicts_in_declared_alphabet"], "remaining pair conflicts")


if __name__ == "__main__":
    main()
