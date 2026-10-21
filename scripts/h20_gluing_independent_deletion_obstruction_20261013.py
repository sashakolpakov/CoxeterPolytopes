#!/usr/bin/env python3
"""Standalone integer/rational replay; requires no LP, NumPy, or network.

Checks two mandatory incompatible facets in a fixed 580-wall H20 dictionary.
With --all-deletions, additionally checks the optional 580 escape archive.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/enumerations/h20_gluing_independent_deletion_obstruction_20261013.json"
ALL_DELETIONS = ROOT / "data/enumerations/h20_gluing_independent_all_deletions_20261013.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-deletions", action="store_true")
    args = parser.parse_args()
    proof = json.loads(CERTIFICATE.read_text())
    Q = proof["basis_gram"]
    R = proof["boundary_normals"]
    f, c = proof["ambient_facet_normal"], proof["strict_timelike_center"]
    assert len(Q) == len(f) == len(c) == 22 and len(R) == 580
    assert all(len(row) == 22 for row in Q + R)
    assert all(Q[i][j] == Q[j][i] for i in range(22) for j in range(22))

    def covector(x):
        return [sum(x[i]*Q[i][j] for i in range(22)) for j in range(22)]

    def dot(x, y):
        return sum(a*b for a, b in zip(covector(x), y))

    def pairings(x):
        qx = covector(x)
        return [sum(a*b for a, b in zip(qx, r)) for r in R]

    # Exact LDL elimination: 21 positive pivots followed by one negative.
    M = [[F(x) for x in row] for row in Q]
    signs = []
    for k in range(22):
        pivot = M[k][k]
        assert pivot != 0
        signs.append(1 if pivot > 0 else -1)
        for i in range(k + 1, 22):
            for j in range(k + 1, 22):
                M[i][j] -= M[i][k]*M[k][j]/pivot
    assert signs == [1]*21 + [-1]
    assert dot(f, f) == 2 and dot(c, f) == 0
    assert dot(c, c) == proof["center_square"] == -58
    assert max(pairings(c)) < 0
    assert all(dot(r, f) == 0 and dot(r, r) > 0 for r in R)
    i, j = proof["mandatory_indices"]
    assert (i, j) == (0, 74)
    gram = [[dot(R[a], R[b]) for b in (i, j)] for a in (i, j)]
    assert gram == proof["mandatory_pair_gram"] == [[4, -4], [-4, 6]]
    assert F(gram[0][1]**2, gram[0][0]*gram[1][1]) == F(2, 3)
    # The explicit point proves they are genuine facets meeting in a ridge,
    # with strict inequalities on every other member of the full dictionary.
    p = proof["strict_ridge_point"]
    values = pairings(p)
    assert dot(p, f) == 0 and dot(p, c) < 0
    assert dot(p, p) == proof["ridge_point_square"] == -364
    assert [a for a, value in enumerate(values) if value == 0] == [i, j]
    assert max(value for a, value in enumerate(values) if a not in (i, j)) < 0

    def check_escape(index, vector, recorded=None):
        values = pairings(vector)
        square, center_pairing = dot(vector, vector), dot(vector, c)
        maximum = max(value for a, value in enumerate(values) if a != index)
        assert dot(vector, f) == 0
        assert square > 0 and center_pairing < 0 and maximum < 0 and values[index] > 0
        if recorded is not None:
            assert square == recorded["square"]
            assert maximum == recorded["maximum_other_wall_pairing"]
            assert center_pairing == recorded["center_pairing"]
            assert values[index] == recorded["excluded_wall_pairing"]

    assert {row["deleted_index"] for row in proof["deletion_witnesses"]} == {i, j}
    for row in proof["deletion_witnesses"]:
        check_escape(row["deleted_index"], row["vector"], row)
    if args.all_deletions:
        more = json.loads(ALL_DELETIONS.read_text())["mandatory"]
        assert {int(k) for k in more} == set(range(580))
        for key, row in more.items():
            check_escape(int(key), row["vector"])
        print("PASS_ALL_580_SINGLE_DELETION_ESCAPES")
    print("PASS_EXACT_TWO_WITNESS_580_DICTIONARY_OBSTRUCTION")


if __name__ == "__main__":
    main()
