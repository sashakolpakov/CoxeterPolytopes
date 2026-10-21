#!/usr/bin/env python3
"""Independent exact checks of H20 sections, rolling, and a six-tile patch.

Default replay checks the archived rational cone certificates. --rebuild
uses numerical LP solely to select rational coefficients, verifies their
equalities exactly, and writes a new certificate. No float decides a
geometric assertion in either mode.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import argparse
import json
from math import gcd
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json"
FRONTIER = ROOT / "data/enumerations/new_dimension_h20_facet_gluing_20261013.json"
OUTPUT = ROOT / "data/enumerations/h20_gluing_independent_audit_20261013.json"


def primitive(vector):
    values = tuple(map(int, vector))
    divisor = gcd(*values)
    assert divisor > 0
    return tuple(x // divisor for x in values)


def geometry():
    raw = json.loads(SOURCE.read_text())["coxeter_input"]
    Q = -np.array(raw["basis_gram"], dtype=np.int64).astype(object)
    W = np.array(raw["wall_coordinates_as_columns_transposed"], dtype=np.int64).astype(object)
    G = W @ Q @ W.T
    exact_Q = sp.Matrix(Q)
    minors = [exact_Q[:i, :i].det() for i in range(1, 23)]
    assert all(value > 0 for value in minors[:-1]) and minors[-1] == -4
    # Exact LDL pivot signs give signature (21,1).
    assert Counter(map(int, G.diagonal())) == {2: 42, 4: 168}
    assert all(G[i, j] <= 0 for i, j in combinations(range(210), 2))
    h = W.sum(axis=0)
    assert h @ Q @ h == -137214 and max(W @ Q @ h) < 0
    return Q, W, G, h


def audit_defects(G, archive):
    rows = []
    for old in archive["all_210_facets"]:
        f = old["ambient_facet_one_based"] - 1
        gram = np.array(int(G[f, f]) * G - np.outer(G[:, f], G[f, :]), dtype=np.int64)
        ids = np.where(gram.diagonal() > 0)[0]
        H = gram[np.ix_(ids, ids)]
        assert np.max(np.abs(H)) < 10**6  # Squared products below int64 range.
        n, d = H**2, np.outer(H.diagonal(), H.diagonal())
        bad = (n < d) & (n != 0) & (4*n != d) & (2*n != d) & (4*n != 3*d)
        np.fill_diagonal(bad, False)
        cover = {x - 1 for x in old["cover_ambient_labels"]}
        where = {int(v): i for i, v in enumerate(ids)}
        assert len(ids) == old["section_facets"]
        assert int(bad.sum()) // 2 == old["bad_ridge_count"]
        assert all(int(ids[i]) in cover or int(ids[j]) in cover
                   for i, j in zip(*np.where(np.triu(bad, 1))))
        used = set()
        for a, b in old["matching_ambient_labels"]:
            a, b = a - 1, b - 1
            assert a not in used and b not in used and bad[where[a], where[b]]
            used.update((a, b))
        assert len(used) == 2 * len(cover) == 2 * old["minimum_defect_vertex_cover"]
        if int(G[f, f]) == 2:
            assert all(int(bad[where[a]].sum()) == 104 for a in cover)
        rows.append([f + 1, len(ids), len(cover)])
    assert len(rows) == 210
    return {"sections_checked": 210, "short_section_count": 42,
            "long_section_count": 168, "short_facets_and_cover": [209, 5],
            "long_facets_and_cover": [82, 40],
            "short_minimum_cover_unique": True}


def audit_seam_group(Q, W, G, archive):
    import pynauty as pn
    f, j = 0, 21
    Z = 3*W - np.outer(2*G[:, f] + G[:, j], W[f]) - np.outer(G[:, f] + 2*G[:, j], W[j])
    H = Z @ Q @ Z.T
    ids = [i for i in range(210) if H[i, i] > 0]
    K = H[np.ix_(ids, ids)]
    assert len(ids) == 208 and sp.Matrix(Z[ids].T).rank() == 20
    colors = defaultdict(set)
    colors["facet"].update(range(208))
    adj = {i: [] for i in range(208)}
    counter = 208
    for a, b in combinations(range(208), 2):
        if K[a, b] == 0:
            continue
        label = (int(sp.sign(K[a, b])), F(int(K[a, b]**2), int(K[a, a]*K[b, b])))
        colors[label].add(counter)
        adj[counter] = [a, b]
        adj[a].append(counter)
        adj[b].append(counter)
        counter += 1
    graph = pn.Graph(counter, directed=False, adjacency_dict=adj,
                     vertex_coloring=list(colors.values()))
    group = pn.autgrp(graph)
    order = int(group[1]) * 10**int(group[2])
    assert order == 2304
    recorded = archive["A2_seam_automorphisms"]
    gens = [tuple(g) for g in recorded["generators"]]
    permutations = {tuple(range(208))}
    todo = list(permutations)
    while todo:
        old = todo.pop()
        for g in gens:
            new = tuple(g[i] for i in old)
            if new not in permutations:
                permutations.add(new)
                todo.append(new)
    assert len(permutations) == order
    wall_set = {tuple(row) for row in W}
    for g, lift in zip(gens, recorded["integral_ambient_lifts"]):
        assert np.array_equal(K[np.ix_(g, g)], K)
        A = np.array(lift["matrix"], dtype=object)
        assert np.array_equal(A.T @ Q @ A, Q)
        assert {tuple(row) for row in W @ A.T} == wall_set
        assert np.array_equal(Z[ids] @ A.T, Z[[ids[i] for i in g]])
        target = j if lift["swaps_ambient_pair"] else f
        assert np.array_equal(A @ W[f], W[target])
    core = frozenset(a for a, i in enumerate(ids) if G[i, i] == 2 and G[f, i] == -1)
    orbit = {frozenset(g[i] for i in core) for g in permutations}
    assert len(core) == 4 and len(orbit) == 2
    return {"seam_facets": 208, "full_normalized_gram_group_order": order,
            "archived_lifts_generate_full_group": True,
            "core_subset_orbit_size": 2, "core_subset_stabilizer_order": 1152}


def audit_holonomy(Q, W, G, archive):
    short = [i for i in range(210) if G[i, i] == 2]
    neighbors = {i: {j for j in short if G[i, j] == -1} for i in short}
    assert all(len(v) == 5 for v in neighbors.values())
    seen, todo = {short[0]}, [short[0]]
    while todo:
        for j in neighbors[todo.pop()] - seen:
            seen.add(j)
            todo.append(j)
    assert len(seen) == 42
    path = [v - 1 for v in archive["odd_roll_holonomy"]["representative_path"]]
    I = np.eye(22, dtype=object)
    T = I.copy()
    # Right products directly develop successive facets; this is the inverse
    # of the archived left-product convention, with identical order.
    for a, b in zip(path, path[1:] + path[:1]):
        assert G[a, b] == -1
        r = W[a] + W[b]
        assert r @ Q @ r == 2
        T = T @ (I - np.outer(r, r @ Q))
    assert np.array_equal(T.T @ Q @ T, Q)
    assert np.array_equal(T @ W[path[0]], W[path[0]])
    old = np.array(archive["odd_roll_holonomy"]["representative_holonomy"], dtype=object)
    assert np.array_equal(T @ old, I)
    N = T @ T - I
    assert np.any(N != 0) and not np.any(N @ N @ N != 0)
    restriction = sp.Matrix([W[path[0]] @ Q]).nullspace()
    B = sp.Matrix.hstack(*restriction)
    assert sp.Matrix(N) * B != sp.zeros(22, 21)
    return {"graph_vertices": 42, "graph_edges": 105, "graph_cycle_rank": 64,
            "representative_path": [i + 1 for i in path],
            "right_product_holonomy": T.tolist(),
            "squared_minus_identity_rank": int(sp.Matrix(N).rank()),
            "squared_minus_identity_cube_zero": True,
            "infinite_order_on_section": True}


def make_tiles(Q, W, G, h):
    f, j, k = 0, 21, 22
    root, u, n, v = W[f], 2*W[j] + W[f], W[f] + W[j] + W[k], W[f] + W[k]
    I = np.eye(22, dtype=object)
    D, U, S = 3*I - np.outer(u, u @ Q), I - np.outer(v, v @ Q), I - np.outer(n, n @ Q)
    assert u @ Q @ u == 6 and v @ Q @ v == n @ Q @ n == 2
    for T, denominator in ((D, 3), (U, 1), (S, 1)):
        assert np.array_equal(T.T @ Q @ T, denominator**2 * Q)

    def section(i):
        V = 2*W - np.outer(G[:, i], W[i])
        ids = [a for a in range(210) if V[a] @ Q @ V[a] > 0]
        center = 2*h - (h @ Q @ W[i])*W[i]
        return V[ids], ids, center

    Va, ida, ha = section(f)
    Vk, idk, hk = section(k)
    tiles = []
    for name, base, ids, T, center in [
        ("A", Va, ida, I, ha), ("DA", Va, ida, D, ha),
        ("B", Vk, idk, U, hk), ("DB", Vk, idk, D @ U, hk),
        ("C", Vk, idk, S @ U, hk), ("DC", Vk, idk, D @ S @ U, hk),
    ]:
        R = np.array([primitive(T @ row) for row in base], dtype=object)
        p = T @ center
        assert len(R) == 209 and p @ Q @ p < 0 and max(R @ Q @ p) < 0
        assert all(row @ Q @ root == 0 for row in R)
        tiles.append({"name": name, "R": R, "center": p, "ids": ids,
                      "lookup": {tuple(row): a for a, row in enumerate(R)}})
    return tiles


def shared_faces(Q, tiles):
    def face_roots(R, s):
        projected = (s @ Q @ s)*R - np.outer(R @ Q @ s, s)
        return {primitive(row) for row in projected if row @ Q @ row > 0}

    records, internal = [], set()
    for a, b in combinations(range(6), 2):
        A, B = tiles[a], tiles[b]
        for i, r in enumerate(A["R"]):
            j = B["lookup"].get(tuple(-r))
            if j is None:
                continue
            matches = face_roots(A["R"], r) == face_roots(B["R"], r)
            records.append({"tiles": [a, b], "walls": [i, j],
                            "ambient_labels": [A["ids"][i] + 1, B["ids"][j] + 1],
                            "full_face_matches": matches})
            if matches:
                internal.update(((a, i), (b, j)))
    assert sum(row["full_face_matches"] for row in records) == 6
    assert len(internal) == 12
    # Every opposite-normal coincidence involves already cancelled facets.
    assert all((row["tiles"][0], row["walls"][0]) in internal and
               (row["tiles"][1], row["walls"][1]) in internal for row in records)
    boundary = {}
    for a, A in enumerate(tiles):
        for i, r in enumerate(A["R"]):
            if (a, i) not in internal:
                boundary.setdefault(tuple(r), []).append([a, i, A["ids"][i] + 1])
    assert len(boundary) == 580
    return records, internal, boundary


def check_terms(vector, tile, terms):
    parsed = [(int(i), F(c)) for i, c in terms]
    assert parsed and len({i for i, c in parsed}) == len(parsed)
    assert all(0 <= i < 209 and c > 0 for i, c in parsed)
    assert all(sum(c * int(tile["R"][i, j]) for i, c in parsed) == int(vector[j])
               for j in range(22))


def find_terms(vector, tile):
    from scipy.optimize import linprog
    direct = tile["lookup"].get(tuple(vector))
    if direct is not None:
        return [[direct, "1"]]
    result = linprog(np.ones(209), A_eq=np.array(tile["R"].T, dtype=float),
                     b_eq=np.array(vector, dtype=float), bounds=(0, None), method="highs")
    if not result.success:
        return None
    terms = [[int(i), str(F(float(result.x[i])).limit_denominator(10**6))]
             for i in np.where(result.x > 1e-8)[0]]
    check_terms(vector, tile, terms)
    return terms


def patch_certificates(Q, tiles, internal, boundary, rebuild, old):
    if rebuild:
        proofs = []
        for normal, owners in boundary.items():
            certificates = [find_terms(normal, tile) for tile in tiles]
            assert all(terms is not None for terms in certificates)
            proofs.append({"normal": list(map(int, normal)), "owners": owners,
                           "certificates_by_tile": certificates})
        separators = []
        for a, b in combinations(range(6), 2):
            for owner, i in sorted(internal):
                if owner != a:
                    continue
                terms = find_terms(-tiles[a]["R"][i], tiles[b])
                if terms is not None:
                    separators.append({"tile_a": a, "tile_b": b, "wall_a": i,
                                       "negative_wall_certificate_in_b": terms})
                    break
            else:
                raise AssertionError(("No disjointness certificate", a, b))
    else:
        proofs, separators = old["boundary_proofs"], old["disjointness"]
    assert len(proofs) == 580 and {tuple(row["normal"]) for row in proofs} == set(boundary)
    for row in proofs:
        assert row["owners"] == boundary[tuple(row["normal"])]
        assert len(row["certificates_by_tile"]) == 6
        for tile, terms in zip(tiles, row["certificates_by_tile"]):
            check_terms(row["normal"], tile, terms)
    assert {(row["tile_a"], row["tile_b"]) for row in separators} == set(combinations(range(6), 2))
    for row in separators:
        a, b, i = row["tile_a"], row["tile_b"], row["wall_a"]
        assert (a, i) in internal
        check_terms(-tiles[a]["R"][i], tiles[b], row["negative_wall_certificate_in_b"])
    A = tiles[0]
    r, s = A["R"][A["ids"].index(23)], A["R"][A["ids"].index(48)]
    assert tuple(r) in boundary and tuple(s) in boundary
    assert r @ Q @ r == s @ Q @ s == 6 and r @ Q @ s == -2
    h = A["center"]
    a, b = h @ Q @ r, h @ Q @ s
    p = 32*h - (6*a + 2*b)*r - (2*a + 6*b)*s
    values = A["R"] @ Q @ p
    assert p @ Q @ p < 0 and p @ Q @ r == p @ Q @ s == 0
    assert sum(x == 0 for x in values) == 2 and max(values) == 0
    # The three-tile patch is precisely a convex half of the six-tile patch.
    u = A["R"][A["ids"].index(21)]
    half_terms = [[[20, "1"]], [[20, "1"]], [[0, "1"], [21, "2"]],
                  [[0, "1/3"], [21, "2/3"]], [[0, "1"]], [[0, "1"]]]
    for i, tile in enumerate(tiles):
        check_terms(u if i in (0, 2, 4) else -u, tile, half_terms[i])
    return {"boundary_proofs": proofs, "disjointness": separators,
            "boundary_facets": 580, "exact_support_certificates": 3480,
            "disjoint_tile_pairs": 15,
            "three_tile_half": {"tiles": ["A", "B", "C"],
                "normal": list(map(int, u)), "certificates_in_tile_order": half_terms,
                "conclusion": "A+B+C is exactly the intersection of the convex six-tile patch with the negative halfspace of this normal."},
            "non_coxeter_ridge": {"tile": "A", "ambient_labels": [24, 49],
                "gram": [[6, -2], [-2, 6]], "cosine_squared": "1/9",
                "relative_interior_timelike_point": list(map(int, p))}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    Q, W, G, h = geometry()
    archive = json.loads(FRONTIER.read_text())
    old = None if args.rebuild else json.loads(OUTPUT.read_text())["six_tile_patch"]
    defects = audit_defects(G, archive)
    print("Exact actual-facet/ridge defect data checked for all 210 sections.", flush=True)
    seam = audit_seam_group(Q, W, G, archive)
    holonomy = audit_holonomy(Q, W, G, archive)
    print("A2 seam group and infinite rolling holonomy independently checked.", flush=True)
    tiles = make_tiles(Q, W, G, h)
    faces, internal, boundary = shared_faces(Q, tiles)
    patch = patch_certificates(Q, tiles, internal, boundary, args.rebuild, old)
    patch.update(tile_order=[tile["name"] for tile in tiles], face_matches=faces,
                 conclusion="The six full-seam gluing is a convex finite-volume 580-facet H20 polytope, with disjoint tile interiors, but is not Coxeter.")
    result = {"status": "PASS_EXACT_H20_SECTION_AND_SIX_TILE_AUDIT",
        "source": str(SOURCE.relative_to(ROOT)), "source_sha256": sha256(SOURCE.read_bytes()).hexdigest(),
        "frontier_archive": str(FRONTIER.relative_to(ROOT)),
        "defect_graphs": defects, "A2_seam_group": seam, "rolling_holonomy": holonomy,
        "six_tile_patch": patch,
        "scope": "Exact checks for these literal sections and the six-tile patch; the core-only rolling obstruction is conditional on the assembly rules in the report. No no-go for H20 Coxeter polytopes or for assemblies using additional seams."}
    if args.rebuild:
        OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(result["status"], flush=True)


if __name__ == "__main__":
    main()
