#!/usr/bin/env python3
"""Exploratory exact reconstruction of the T1 rank-18 unfolding package.

This probe is intentionally separate from the frozen D34 census.  It exposes
the internal matrices needed to test whether the weak T1 carrier can be
promoted to a complete invariant-nef K3 realization.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import itertools
import json
import math

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

import scripts.certify_h9_14_t1_t2_d34_unfolding_census as census
import scripts.certify_h9_14_t2_surface_unfolding_pilot as pilot


def reconstruct() -> dict[str, object]:
    corpus, _ = census._corpus()
    representative = next(row for row in corpus if row["id"] == "A14_8_1_C01")
    raw = representative["gram"]
    edges, graph, _ = census._relation_data(raw)
    scaling_rows, scalings, _, _ = census._minimal_scalings(edges, graph)
    e8_gram, e8_classes = pilot._e8_data()

    survivor = None
    folded_gram = None
    for scaling_row, scaling in zip(scaling_rows, scalings, strict=True):
        if not scaling_row["one_two_sheet_compatible"]:
            continue
        folded = pilot._gram_for_scaling(scaling, edges)
        _, internal = census._fixed_candidates(
            scaling_row["number"], folded, e8_gram, e8_classes
        )
        for candidate in internal:
            if candidate["survives_minimal_two_sheet_weak_K3_gate"]:
                if survivor is not None:
                    raise ArithmeticError("T1 survivor is no longer unique")
                survivor = candidate
                folded_gram = folded
    if survivor is None or folded_gram is None:
        raise ArithmeticError("T1 survivor disappeared")

    fixed_gram = survivor["gram"]
    fixed_roots = survivor["roots"]
    long_indices = [
        index for index in range(14) if folded_gram[index, index] == -4
    ]
    long_basis_indices = [
        index - 1 for index in survivor["long_class_basis_facets_one_based"]
    ]
    fixed_basis_classes = [
        pilot._vector_bit(fixed_roots[:, index]) for index in long_basis_indices
    ]
    target_basis_classes = survivor["e8_embedding_basis_classes"]
    anti_roots: dict[int, sp.Matrix] = {}
    for index in long_indices:
        relation = pilot._relation(
            pilot._vector_bit(fixed_roots[:, index]), fixed_basis_classes
        )
        target = 0
        for position, basis_class in enumerate(target_basis_classes):
            if (relation >> position) & 1:
                target ^= basis_class
        anti_roots[index] = sp.Matrix(max(e8_classes[target]))

    anti_gram = -2 * e8_gram
    direct_gram = sp.diag(fixed_gram, anti_gram)
    direct_roots: list[sp.Matrix] = []
    root_orbits: list[list[int]] = []
    glue: list[sp.Matrix] = []
    for facet in range(14):
        fixed = sp.Matrix.vstack(fixed_roots[:, facet], sp.zeros(8, 1))
        if folded_gram[facet, facet] == -2:
            direct_roots.append(fixed)
            root_orbits.append([len(direct_roots) - 1])
        else:
            anti = sp.Matrix.vstack(sp.zeros(10, 1), anti_roots[facet])
            plus, minus = (fixed + anti) / 2, (fixed - anti) / 2
            direct_roots.extend((plus, minus))
            root_orbits.append([len(direct_roots) - 2, len(direct_roots) - 1])
            glue.append(plus)
    basis = hermite_normal_form(
        (2 * sp.eye(18)).row_join(
            sp.Matrix.hstack(*[(2 * vector).applyfunc(int) for vector in glue])
        )
    ) / 2
    n_gram = pilot._integral(
        sp.simplify(basis.T * direct_gram * basis), "N Gram is not integral"
    )
    roots = pilot._integral(
        sp.simplify(basis.inv() * sp.Matrix.hstack(*direct_roots)),
        "root coordinates are not integral",
    )
    involution = pilot._integral(
        sp.simplify(basis.inv() * sp.diag(sp.eye(10), -sp.eye(8)) * basis),
        "involution is not integral",
    )
    fixed_embedding = pilot._integral(
        basis.inv() * sp.Matrix.vstack(sp.eye(10), sp.zeros(8, 10)),
        "fixed embedding is not integral",
    )
    anti_embedding = pilot._integral(
        basis.inv() * sp.Matrix.vstack(sp.zeros(10, 8), sp.eye(8)),
        "anti embedding is not integral",
    )
    return {
        "raw": raw,
        "folded_gram": folded_gram,
        "fixed_gram": fixed_gram,
        "fixed_roots": fixed_roots,
        "anti_gram": anti_gram,
        "direct_gram": direct_gram,
        "basis": basis,
        "n_gram": n_gram,
        "roots": roots,
        "involution": involution,
        "fixed_embedding": fixed_embedding,
        "anti_embedding": anti_embedding,
        "root_orbits": root_orbits,
        "survivor": survivor,
    }


def _primitive_column(vector: sp.Matrix) -> sp.Matrix:
    denominator = math.lcm(*(int(value.q) for value in vector))
    integral = [int(denominator * value) for value in vector]
    divisor = math.gcd(*map(abs, integral))
    return sp.Matrix([value // divisor for value in integral])


def cone_rays(fixed_gram: sp.Matrix, fixed_roots: sp.Matrix) -> list[sp.Matrix]:
    inequality = fixed_roots.T * fixed_gram
    rays: dict[tuple[int, ...], sp.Matrix] = {}
    for active in itertools.combinations(range(inequality.rows), 9):
        block = inequality[list(active), :]
        if block.rank() != 9:
            continue
        vector = _primitive_column(block.nullspace()[0])
        slacks = inequality * vector
        if all(value <= 0 for value in slacks):
            vector = -vector
            slacks = -slacks
        if not all(value >= 0 for value in slacks):
            continue
        square = int((vector.T * fixed_gram * vector)[0])
        if square < 0:
            continue
        key = tuple(map(int, vector))
        rays[key] = vector
    return [rays[key] for key in sorted(rays)]


def small_fixed_vectors(
    fixed_gram: sp.Matrix, bound: int, target_square: int
) -> list[sp.Matrix]:
    """Enumerate fixed vectors of a target square with nine coordinates bounded.

    Coordinate 6 (zero based 5) is solved linearly from the norm equation.
    When its coefficient vanishes, its value is also restricted to the same
    exploratory box.
    """

    free = [index for index in range(10) if index != 5]
    roots: dict[tuple[int, ...], sp.Matrix] = {}
    for values in itertools.product(range(-bound, bound + 1), repeat=9):
        vector = sp.zeros(10, 1)
        for index, value in zip(free, values, strict=True):
            vector[index] = value
        coefficient = 2 * sum(
            int(fixed_gram[5, index]) * int(vector[index]) for index in free
        )
        constant = int((vector.T * fixed_gram * vector)[0])
        candidates: list[int]
        if coefficient:
            numerator = target_square - constant
            if numerator % coefficient:
                continue
            candidates = [numerator // coefficient]
        elif constant == target_square:
            candidates = list(range(-bound, bound + 1))
        else:
            continue
        for value in candidates:
            vector[5] = value
            if (vector.T * fixed_gram * vector)[0] != target_square:
                raise ArithmeticError("linear root solver failed")
            first = next((int(x) for x in vector if x), 0)
            oriented = vector if first > 0 else -vector
            roots[tuple(map(int, oriented))] = oriented.copy()
    return [roots[key] for key in sorted(roots)]


def hidden_fixed_root_search(package: dict[str, object], bound: int) -> dict[str, object]:
    fixed_gram = package["fixed_gram"]
    fixed_roots = package["fixed_roots"]
    rays = cone_rays(fixed_gram, fixed_roots)
    short_roots = small_fixed_vectors(fixed_gram, bound, -2)
    long_vectors = small_fixed_vectors(fixed_gram, bound, -4)
    _, e8_classes = pilot._e8_data()
    anti_roots = {
        tuple(vector)
        for vectors in e8_classes.values()
        for vector in vectors
    }
    anti_embedding = package["anti_embedding"]
    fixed_embedding = package["fixed_embedding"]
    allowed_anti_parities = {
        tuple(int(value) % 2 for value in anti_embedding * sp.Matrix(vector))
        for vector in anti_roots
    }
    eligible_long = [
        vector
        for vector in long_vectors
        if tuple(int(value) % 2 for value in fixed_embedding * vector)
        in allowed_anti_parities
    ]
    roots = [*short_roots, *eligible_long]
    pairing_matrix = (
        sp.Matrix.hstack(*roots).T
        * fixed_gram
        * sp.Matrix.hstack(*rays)
    )
    crossing: list[dict[str, object]] = []
    supporting: list[dict[str, object]] = []
    for root_index, root in enumerate(roots):
        pairings = [int(value) for value in pairing_matrix.row(root_index)]
        if min(pairings) < 0 < max(pairings):
            crossing.append(
                {
                    "root": list(map(int, root)),
                    "square": int((root.T * fixed_gram * root)[0]),
                    "minimum_ray_pairing": min(pairings),
                    "maximum_ray_pairing": max(pairings),
                    "negative_ray": list(map(int, rays[pairings.index(min(pairings))])),
                    "positive_ray": list(map(int, rays[pairings.index(max(pairings))])),
                }
            )
        elif min(pairings) == 0 or max(pairings) == 0:
            supporting.append(
                {
                    "root": list(map(int, root)),
                    "square": int((root.T * fixed_gram * root)[0]),
                    "minimum_ray_pairing": min(pairings),
                    "maximum_ray_pairing": max(pairings),
                }
            )
    return {
        "bound": bound,
        "cone_ray_count": len(rays),
        "ray_square_histogram": dict(
            sorted(
                Counter(
                    int((ray.T * fixed_gram * ray)[0]) for ray in rays
                ).items()
            )
        ),
        "small_fixed_root_count_up_to_sign": len(roots),
        "small_square_minus_two_count_up_to_sign": len(short_roots),
        "small_square_minus_four_count_up_to_sign": len(long_vectors),
        "eligible_square_minus_four_count_up_to_sign": len(eligible_long),
        "crossing_root_count": len(crossing),
        "first_crossing_roots": crossing[:20],
        "supporting_root_count": len(supporting),
        "first_supporting_roots": supporting[:20],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrices", action="store_true")
    parser.add_argument("--hidden-search-bound", type=int, default=0)
    arguments = parser.parse_args()
    package = reconstruct()
    roots = package["roots"]
    gram = package["n_gram"]
    action = package["involution"]
    root_gram = roots.T * gram * roots
    graph = nx.Graph()
    graph.add_nodes_from(range(roots.cols))
    for left in range(roots.cols):
        for right in range(left):
            if root_gram[left, right]:
                graph.add_edge(left, right, pairing=int(root_gram[left, right]))
    automorphisms = list(nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter())
    action_lookup = {tuple(roots[:, index]): index for index in range(roots.cols)}
    action_permutation = [
        action_lookup[tuple(action * roots[:, index])] for index in range(roots.cols)
    ]
    payload = {
                "N_det": int(gram.det()),
                "fixed_det": int(package["fixed_gram"].det()),
                "root_count": roots.cols,
                "root_span_index": abs(int(hermite_normal_form(roots).det())),
                "degree_histogram": dict(sorted(Counter(dict(graph.degree()).values()).items())),
                "edge_count": graph.number_of_edges(),
                "graph_automorphism_order": len(automorphisms),
                "action_permutation_one_based": [value + 1 for value in action_permutation],
                "action_orbit_sizes": sorted(len(orbit) for orbit in nx.connected_components(nx.Graph([(i, action_permutation[i]) for i in range(roots.cols)]))),
    }
    if arguments.matrices:
        payload.update(
            {
                "fixed_gram": census._rows(package["fixed_gram"]),
                "fixed_roots": census._rows(package["fixed_roots"]),
                "N_gram": census._rows(gram),
                "roots": census._rows(roots),
                "involution": census._rows(action),
            }
        )
    if arguments.hidden_search_bound:
        payload["hidden_fixed_root_search"] = hidden_fixed_root_search(
            package, arguments.hidden_search_bound
        )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
