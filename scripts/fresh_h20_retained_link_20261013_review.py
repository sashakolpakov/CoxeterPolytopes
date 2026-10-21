#!/usr/bin/env python3
"""Independent SymPy reconstruction of the mixed-factor inverse census.

The producer uses its own Fraction elimination. This imports no producer,
reconstructs the declared connected deletions from the source census, and
checks their inverse diagonals and all attachment inequalities independently.
The actual-facet reduction is reviewed in the companion report, not formalized.
"""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from sympy import Matrix, Rational

P = Path(__file__).resolve().parents[1]
TAG = "fresh_h20_retained_link_20261013"


def connected(A):
    reached, todo = {0}, [0]
    while todo:
        i = todo.pop()
        for j in range(A.rows):
            if j not in reached and A[i, j]:
                reached.add(j)
                todo.append(j)
    return len(reached) == A.rows


def main():
    source = P / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json"
    certificate = P / "results" / (TAG + ".json")
    data = json.loads(source.read_text())["coxeter_input"]
    saved = json.loads(certificate.read_text())
    for path, digest in saved["sources"].items():
        assert sha256((P / path).read_bytes()).hexdigest() == digest
    Q, W = -Matrix(data["basis_gram"]), Matrix(data["wall_coordinates_as_columns_transposed"])
    G = W * Q * W.T
    walk = json.loads((P / "results/fresh_h21_source_corner_walk_20261013.json").read_text())
    expected = {}
    cusps = [v for v in walk["corner_orbits"] if v["type"] == "ideal"]
    assert len(cusps) == 18
    counts = Counter()
    for v in cusps:
        for factor in v["affine_components"]:
            walls = factor["walls"]
            counts["factors"] += 1
            if len(walls) == 2:
                assert v["id"] == 0 and walls == [162, 164]
                continue
            for omitted in walls:
                keep = [w for w in walls if w != omitted]
                B = G.extract(keep, keep)
                assert B.is_positive_definite
                counts["positive_deletions"] += 1
                if connected(B):
                    expected[(v["id"], tuple(walls), omitted)] = (keep, B)
    assert len(expected) == 247 and counts["positive_deletions"] == 411
    seen = set()
    minima = {}
    for row in saved["connected_deletion_records"]:
        key = (row["cusp"], tuple(row["factor"]), row["omitted_wall"])
        assert key not in seen
        seen.add(key)
        keep, B = expected[key]
        assert keep == row["retained"]
        inv = B.inv()
        assert all(x > 0 for x in inv)
        assert [a["wall"] for a in row["attachments"]] == keep
        for at in row["attachments"]:
            i, j = at["wall"], at["neighbor"]
            assert j in keep and i != j and G[i, j] < 0
            c = G[i, j] ** 2 / (G[i, i] * G[j, j])
            A = inv[keep.index(i), keep.index(i)]
            D = A * G[i, i]
            assert A == Rational(at["inverse_diagonal"])
            assert c == Rational(at["old_cosine_squared"])
            assert D == Rational(at["normalized_inverse_diagonal"])
            assert D >= 1 / (1 - c) == Rational(at["two_node_lower_bound"])
            assert c in [Rational(1, 4), Rational(1, 2)]
            counts["attachments"] += 1
            if row["cusp"] == 0:
                f = tuple(row["factor"])
                minima[f] = min(D, minima.get(f, D))
    assert seen == set(expected) and counts["attachments"] == 2067
    assert sorted((len(k) - 1, str(v)) for k, v in minima.items()) == [(5, "2"), (7, "7/4"), (7, "7/4")]
    out = {"status": "PASS_INDEPENDENT_INVERSE_CENSUS_RECONSTRUCTION",
           "certificate_sha256": sha256(certificate.read_bytes()).hexdigest(),
           "reviewer_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
           "counts": dict(counts),
           "scope": "All declared connected deletions and inverse attachment bounds reconstructed independently with SymPy exact matrices; no producer imported. Root separately reviewed the actual-facet, simplex/triangle and all-label arguments in the report. This is not a formal proof assistant verification."}
    (P / "results" / (TAG + "_review.json")).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out))


if __name__ == "__main__":
    main()
