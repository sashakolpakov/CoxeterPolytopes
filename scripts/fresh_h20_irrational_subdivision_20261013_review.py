#!/usr/bin/env python3
"""Independent Fraction check of the subdivision obstruction's local geometry.

No producer or repository arithmetic helpers are imported. All saved face
points, parameter-dependent Grams, and exceptional ridge formulas are rebuilt
from the original integer form. The tiling/conjugate arguments require the
mathematical proof in the report and were reviewed separately.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = "fresh_h20_irrational_subdivision_20261013"


def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p


def multiply(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return trim(c)


def main():
    path = ROOT / "results" / (TAG + ".json")
    cert = json.loads(path.read_text())
    for p, digest in cert["sources"].items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == digest
    raw = json.loads((ROOT / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json").read_text())["coxeter_input"]
    Q = [[-F(z) for z in row] for row in raw["basis_gram"]]
    W = [[F(z) for z in row] for row in raw["wall_coordinates_as_columns_transposed"]]

    def dual(v):
        return [sum(x * y for x, y in zip(row, v)) for row in Q]

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    def pair(a, b):
        return dot(a, dual(b))

    C = [dual(w) for w in W]
    G = [[dot(w, c) for c in C] for w in W]
    old = json.loads((ROOT / "results/fresh_h20_cusp_sections_20261013.json").read_text())
    h = list(map(F, old["strict_center"]))
    point_cache = {}

    def face(row):
        x = list(map(F, row["point"]))
        key = tuple(x)
        if key not in point_cache:
            vals = [dot(x, c) for c in C]
            assert pair(x, x) < 0 and pair(x, h) < 0 and max(vals) == 0
            point_cache[key] = [i for i, v in enumerate(vals) if not v]
        assert point_cache[key] == row["walls"]
        return x

    counts = Counter()
    signed_planes = set()
    assert sorted(row["source_row"] for row in cert["rows"]) == list(range(118))
    for row in cert["rows"]:
        a, b = row["retained_leaf"], row["omitted_leaf"]
        e, d = list(map(F, row["base_normal"])), list(map(F, row["cusp"]))
        original = old["pencils"][row["source_row"]]
        assert a == original["retained_leaf"] and b == original["omitted_leaf"]
        assert e == list(map(F, original["base_normal"]))
        assert d == list(map(F, original["cusp"]))
        assert pair(e, e) == 1 and pair(e, d) == pair(d, d) == 0
        assert e == [(x - y) / 2 for x, y in zip(W[a], W[b])]
        ep, dp = [dot(e, c) for c in C], [dot(d, c) for c in C]
        assert ep[a] == 1 and dp[a] == 0
        signed_planes.add((tuple(e), tuple(d)))

        def gram(i, j):
            return trim([G[i][j] - ep[i] * ep[j], -ep[i] * dp[j] - dp[i] * ep[j], -dp[i] * dp[j]])

        def actual(proof, wanted):
            left, right = face(proof["left_face"]), face(proof["right_face"])
            A, B = -pair(d, left), -pair(d, right)
            assert A == F(proof["A"]) > 0 and B == F(proof["B"]) > 0
            assert pair(e, left) == 0 and pair(e, right) == B / 2
            zeros = set(proof["left_face"]["walls"]) & set(proof["right_face"]["walls"])
            assert zeros == set(wanted) == set(proof["actual_section_zero_walls_for_0_lt_t_lt_half"])
            assert set(proof["left_face"]["walls"]) == set(wanted) | {b}
            # These identities prove the interpolation for every real t.
            assert a in wanted and b not in wanted
            counts["actual_face_interpolations"] += 1

        if row["method"].startswith("SPHERICAL"):
            c, j = row["common_neighbor"], row["auxiliary_wall"]
            assert row["projected_Gram_order"] == [a, c, j]
            actual(row["actual_codimension3_face"], [a, c, j])
            assert [[gram(i, k) for k in [a, c, j]] for i in [a, c, j]] == [[trim(p) for p in r] for r in row["projected_Gram_polynomials"]]
            assert gram(a, a) == [1] and gram(c, c) == [2]
            assert gram(a, c) == [-1] and gram(c, j) == [0]
            N, D = multiply(gram(a, j), gram(a, j)), gram(j, j)
            assert N == trim(row["theta_cosine_squared"]["numerator"])
            assert D == trim(row["theta_cosine_squared"]["denominator"])
            assert N == [F(1, 4), -1, 1] and D[0] == F(7, 4) and D[2] == -1
            assert D[1] in [-1, 1]
            # D-4N is a concave quadratic positive at both endpoints.
            assert D[0] - 4 * N[0] > 0
            assert sum((x - 4 * y) * F(1, 2) ** i for i, (x, y) in enumerate(zip(D, N))) > 0
            counts["spherical_orientations"] += 1
        else:
            for label in ["A", "D"]:
                j = row[label + "_auxiliary_wall"]
                actual(row[label + "_ridge"], [a, j])
                N = [4 * x for x in multiply(gram(a, j), gram(a, j))]
                D = multiply(gram(a, a), gram(j, j))
                formula = row[label + "_4cosine_squared"]
                assert trim(N) == trim(formula["numerator"])
                assert D == trim(formula["denominator"])
            assert row["A_4cosine_squared"] == {"numerator": [1, -4, 4], "denominator": ["7/4", 1, -1]}
            assert row["D_4cosine_squared"] == {"numerator": [4, -32, 64], "denominator": [3, 8, -16]}
            counts["exceptional_orientations"] += 1
    assert counts == {"actual_face_interpolations": 122, "spherical_orientations": 114, "exceptional_orientations": 4}
    assert len(signed_planes) == 118
    assert all((tuple(-x for x in e), d) in signed_planes for e, d in signed_planes)
    # Strict rational bounds used in the exceptional conjugate argument.
    assert F(2) < F(3, 2) ** 2
    assert [F(3, 4) + 5 * t - 5 * t * t for t in [F(-1, 8), F(5, 8)]] == [F(3, 64), F(123, 64)]
    result = {"status": "PASS_INDEPENDENT_LOCAL_GEOMETRY_AND_FORMULAS",
              "certificate_sha256": sha256(path.read_bytes()).hexdigest(),
              "reviewer_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "counts": dict(counts), "distinct_face_points": len(point_cache),
              "scope": "Direct Fraction reconstruction of all actual-face interpolation inputs and every projected Gram/ridge polynomial. Spherical tiling and cyclotomic arguments separately reviewed; this is not a formal proof assistant."}
    (ROOT / "results" / (TAG + "_review.json")).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
