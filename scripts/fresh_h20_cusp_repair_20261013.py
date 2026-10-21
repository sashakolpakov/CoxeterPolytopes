#!/usr/bin/env python3
"""Exact local affine-E7 and affine-D4 replacements at two H20 cusps.

These are supporting-plane proposals, not finite-volume H20 examples.
Both old mandatory-wall escapes survive all four genuinely new planes.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

import h20_gluing_independent_escape_20261013 as old

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/fresh_h20_cusp_repair_20261013.json"


def rows(A):
    return [[int(x) for x in r] for r in A]


def root_closure(K):
    roots = {tuple(row) for row in np.eye(len(K), dtype=int)}
    pending = list(roots)
    while pending:
        v = np.array(pending.pop(), dtype=object)
        for i in range(len(K)):
            numerator, denominator = 2*int(v @ K[:, i]), int(K[i, i])
            assert numerator % denominator == 0
            w = v.copy()
            w[i] -= numerator // denominator
            key = tuple(map(int, w))
            if key not in roots:
                roots.add(key)
                pending.append(key)
    return sorted(roots)


def cone_coefficients(generators, vectors):
    A = sp.Matrix(generators.T)
    assert A.rank() == len(generators)
    answer = []
    for r in vectors:
        x, parameters = A.gauss_jordan_solve(sp.Matrix(r))
        assert not parameters and all(t >= 0 for t in x)
        answer.append([str(t) for t in x])
    return answer


def choose_alcove(roots, delta, Q, center):
    assert delta @ Q @ delta == 0 and delta @ Q @ center < 0
    shifted = []
    for v in roots:
        assert v @ Q @ delta == 0
        k = int(v @ Q @ center) // int(-(delta @ Q @ center)) + 1
        r = v + k*delta
        assert r @ Q @ center < 0
        shifted.append(r)
    shifted.sort(key=lambda r: (int(-r @ Q @ center), tuple(r)))
    chosen = []
    for r in shifted:
        if all(r @ Q @ s <= 0 for s in chosen):
            chosen.append(r)
    chosen = np.array(chosen, dtype=object)
    return chosen, cone_coefficients(chosen, shifted)


def decompose_cusps(Q, f, center, R):
    source = json.loads((ROOT / "data/enumerations/h20_gluing_independent_cusps_20261013.json").read_text())
    result = []
    for cusp in source["cusps"]:
        u = np.array(cusp["vector"], dtype=object)
        values = R @ Q @ u
        assert u @ Q @ u == u @ Q @ f == 0 and u @ Q @ center < 0
        assert max(values) <= 0
        ids = [i for i, value in enumerate(values) if value == 0]
        assert ids == cusp["incident_boundary_indices"]
        G = R[ids] @ Q @ R[ids].T
        remaining, components = set(range(len(ids))), []
        while remaining:
            todo, component = [min(remaining)], []
            while todo:
                i = todo.pop()
                if i not in remaining:
                    continue
                remaining.remove(i)
                component.append(i)
                todo.extend(j for j in remaining if G[i, j] != 0)
            component.sort()
            H = G[np.ix_(component, component)]
            components.append({"boundary_indices": [ids[i] for i in component],
                               "gram": rows(H), "rank": int(sp.Matrix(H).rank())})
        assert sum(item["rank"] for item in components) == 19
        result.append({"vector": cusp["vector"], "components": components})
    return result


def build():
    Q, f, center, R, _ = old.setup()
    _, W, G, _ = old.audit.geometry()
    cusp_decomposition = decompose_cusps(Q, f, center, R)
    existing = {tuple(r): i for i, r in enumerate(R)}

    parent = [0, 25, 28, 29, 30, 32, 38, 48, 49]
    marks = [3, 6, 5, 1, 4, 2, 3, 2, 4]
    delta = sum((m*W[i] for m, i in zip(marks, parent)), np.zeros(22, dtype=object))
    assert list(delta) == cusp_decomposition[0]["vector"]
    finite = [i for i in parent if i != 29]
    roots = root_closure(G[np.ix_(finite, finite)])
    assert len(roots) == 240
    perpendicular = [np.array(v, dtype=object) @ W[finite] for v in roots]
    perpendicular = [v for v in perpendicular if v @ Q @ f == 0]
    assert len(perpendicular) == 126 and all(v @ Q @ v == 2 for v in perpendicular)
    e7, root_coefficients = choose_alcove(perpendicular, delta, Q, center)
    assert len(e7) == 8
    e7_ids = [26, 29, 30, 31, 32, 38, 100, 101]
    e7_gram = e7 @ Q @ e7.T
    assert sp.Matrix(e7_gram).rank() == 7
    assert sorted(sum(x < 0 for x in row) for row in e7_gram) == [1, 1, 1, 2, 2, 2, 2, 3]
    assert set(e7_gram.flatten()) <= {2, 0, -1}
    e7_result = {"source_cusp_index": 0, "parent_affine_component_one_based_labels": [i+1 for i in parent],
                 "parent_affine_marks": marks, "delta": list(map(int, delta)),
                 "finite_E8_root_count": len(roots), "orthogonal_E7_root_count": len(perpendicular),
                 "old_factor_indices": e7_ids, "new_normals": rows(e7), "new_gram": rows(e7_gram),
                 "existing_boundary_indices": [existing.get(tuple(r)) for r in e7],
                 "new_normals_as_nonnegative_old_combinations": cone_coefficients(R[e7_ids], e7),
                 "all_shifted_E7_roots_as_nonnegative_new_combinations": root_coefficients}

    d4_ids, scales, marks = [0, 1, 2, 3, 41], [3, 3, 3, 3, 1], [1, 2, 3, 2, 1]
    B = np.array([s*R[i] for s, i in zip(scales, d4_ids)], dtype=object)
    H = B @ Q @ B.T
    delta = sum((m*r for m, r in zip(marks, B)), np.zeros(22, dtype=object))
    assert sp.Matrix(np.array([delta, cusp_decomposition[6]["vector"]]).T).rank() == 1
    roots = root_closure(H[:4, :4])
    assert len(roots) == 48
    short = [np.array(v, dtype=object) @ B[:4] for v in roots]
    short = [v for v in short if v @ Q @ v == 18]
    assert len(short) == 24
    d4, root_coefficients = choose_alcove(short, delta, Q, center)
    assert len(d4) == 5
    d4_gram = d4 @ Q @ d4.T
    assert set(d4_gram.flatten()) <= {18, 0, -9} and sp.Matrix(d4_gram).rank() == 4
    assert sorted(sum(x < 0 for x in row) for row in d4_gram) == [1, 1, 1, 1, 4]
    primitive_d4 = np.array([old.audit.primitive(r) for r in d4], dtype=object)
    d4_result = {"source_cusp_index": 6, "old_factor_indices": d4_ids,
                 "old_roots_scale": scales, "parent_affine_marks": marks, "delta": list(map(int, delta)),
                 "finite_F4_root_count": len(roots), "short_D4_root_count": len(short),
                 "scaled_new_normals": rows(d4), "scaled_new_gram": rows(d4_gram),
                 "primitive_new_normals": rows(primitive_d4),
                 "existing_boundary_indices": [existing.get(tuple(r)) for r in primitive_d4],
                 "scaled_new_normals_as_nonnegative_old_combinations": cone_coefficients(R[d4_ids], d4),
                 "all_shifted_D4_roots_as_nonnegative_scaled_new_combinations": root_coefficients}

    new = [r for r in list(e7) + list(primitive_d4) if tuple(r) not in existing]
    assert len(new) == 4
    cert = json.loads((ROOT / "data/enumerations/h20_gluing_independent_deletion_obstruction_20261013.json").read_text())
    old_escape_checks = []
    for witness in cert["deletion_witnesses"]:
        x = np.array(witness["vector"], dtype=object)
        pairings = [int(r @ Q @ x) for r in new]
        assert max(pairings) < 0
        old_escape_checks.append({"deleted_old_wall": witness["deleted_index"], "new_wall_pairings": pairings})
    return {"schema": "coxeter.fresh_h20_cusp_repair.v1", "cusp_decomposition": cusp_decomposition,
            "affine_E7_replacement": e7_result, "affine_D4_replacement": d4_result,
            "four_new_primitive_normals": rows(new), "old_minimal_escapes_still_satisfy_new_walls": old_escape_checks,
            "scope": "Exact Euclidean cusp-link replacements yield four globally supporting planes. This is local geometry only: both old mandatory-wall escapes survive, so these four additions do not bypass the fixed-dictionary obstruction."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    data = build()
    if args.verify:
        assert data == json.loads(OUTPUT.read_text())
    else:
        OUTPUT.write_text(json.dumps(data, indent=2) + "\n")
    print("PASS: seven cusp decompositions; affine E7/D4 local repairs; four new supporting planes; both old escapes survive")


if __name__ == "__main__":
    main()
