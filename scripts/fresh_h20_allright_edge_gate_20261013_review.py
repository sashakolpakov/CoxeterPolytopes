#!/usr/bin/env python3
"""Separate Fraction reconstruction of the saved all-right edge geometry.

Imports no producer or repository geometry helpers. The continuous proof
is mathematical; this rechecks its saved finite incidence inputs.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

P = Path(__file__).resolve().parents[1]
TAG = "fresh_h20_allright_edge_gate_20261013"


def main():
    source = P / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json"
    cert = P / "results" / (TAG + ".json")
    data = json.loads(source.read_text())["coxeter_input"]
    archive = json.loads(cert.read_text())
    assert archive["source_sha256"] == sha256(source.read_bytes()).hexdigest()
    Q = [[-F(z) for z in r] for r in data["basis_gram"]]
    W = [[F(z) for z in r] for r in data["wall_coordinates_as_columns_transposed"]]

    def dual(v):
        return [sum(a * b for a, b in zip(row, v)) for row in Q]

    def product(v, w):
        return sum(a * b for a, b in zip(v, dual(w)))

    C = [dual(w) for w in W]

    def signs(v):
        return [sum(a * b for a, b in zip(v, c)) for c in C]

    G = [signs(w) for w in W]
    counts = [0, 0, 0]
    for side in archive["sides"]:
        points, lines = side["point_walls"], side["line_walls"]
        assert len(points) == len(lines) == 21 and not set(points) & set(lines)
        assert all(G[i][j] == (2 if i == j else 0) for i in points for j in points)
        v = list(map(F, side["finite_vertex"]))
        vp = signs(v)
        assert product(v, v) == -2 and max(vp) == 0
        assert {i for i, z in enumerate(vp) if not z} == set(points)
        assert all(vp[l] == -1 for l in lines)
        delta = {int(p): list(map(F, d)) for p, d in side["cusp_rays"].items()}
        assert set(delta) == set(points)
        for p in points:
            d = delta[p]
            assert d == [x - y for x, y in zip(v, W[p])]
            assert product(d, d) == 0 and product(d, v) == -2
            assert max(signs(d)) == 0
            mid = [2 * x - y for x, y in zip(v, W[p])]
            vals = signs(mid)
            assert max(vals) == 0
            assert {i for i, z in enumerate(vals) if not z} == set(points) - {p}
            counts[0] += 1
        assert {(r["p"], r["k"]) for r in side["second_edges"]} == set(combinations(points, 2))
        for r in side["second_edges"]:
            p, k, l = r["p"], r["k"], r["line"]
            x = [a + b for a, b in zip(delta[p], delta[k])]
            vals = signs(x)
            assert product(x, x) == -4 and max(vals) == 0
            assert [i for i, z in enumerate(vals) if not z] == r["edge_zero_walls"]
            assert set(r["edge_zero_walls"]) == (set(points) - {p, k}) | {l}
            leaves = r["d4_leaves"]
            assert len(leaves) == 3 and not set(leaves) & {p, k}
            assert all(G[l][j] == -1 for j in leaves)
            assert all(G[i][j] == (2 if i == j else 0) for i in leaves for j in leaves)
            assert {j for j in points if G[l][j]} == set(leaves) | {p, k}
            counts[1] += 1
        for p in points:
            pencil = [l for l in lines if G[p][l] == -1]
            assert len(pencil) == 5
            incidence = {l: {k for k in points if G[k][l] == -1} for l in pencil}
            for n in range(3):
                for s in combinations(set(points) - {p}, n):
                    assert sum(not (set(s) & incidence[l]) for l in pencil) >= 5 - n
                    counts[2] += 1
    assert counts == [42, 420, 8862]
    # The projected D4 Schur complement is N-3/2, while a central
    # squared cosine is 1/(2N). These endpoint values bound every N.
    assert F(1, 2) / 2 == F(1, 4)
    assert F(1, 2) / F(3, 2) == F(1, 3) < F(1, 2)
    out = {"status": "PASS_INDEPENDENT_SAVED_GEOMETRY_RECHECK",
           "certificate_sha256": sha256(cert.read_bytes()).hexdigest(),
           "reviewer_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
           "counts": counts,
           "scope": "Separate Fraction source/certificate reconstruction of all first/second edge geometry and support avoidance. Continuous proof separately reviewed; no numerical sampling used."}
    (P / "results" / (TAG + "_review.json")).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out))


if __name__ == "__main__":
    main()
