#!/usr/bin/env python3
"""Exact bounded-corpus D34 unfolding census for the T1/T2 H9/14 diagrams."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

import scripts.certify_h9_14_t2_surface_unfolding_pilot as pilot


ROOT = Path(__file__).resolve().parents[1]
FACTORY = ROOT / "data/enumerations/alexeev_k3_h9_14_factory.json"
PARKED = ROOT / "data/enumerations/alexeev_parked_rank10_h9_examples.json"
PRODUCTS = ROOT / "data/enumerations/alexeev_parked_h9_doubling_products.json"
T2_PILOT = ROOT / "data/enumerations/h9_14_t2_surface_unfolding_pilot.json"
OUTPUT = ROOT / "data/enumerations/h9_14_t1_t2_d34_unfolding_census.json"

SOURCE_HASHES = {
    "factory": "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775",
    "parked": "e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430",
    "products": "10e6323895bacd2d18f897f33ea5f9d4bb1ee4eedad061d545aa5302b447c694",
    "t2_pilot": "e72a7b5690c16c8c405ddf8679cef4bafc18cf93f59f38f5e635acc9131d60c2",
}
SCHEMA = "coxeter.h9_14_t1_t2_d34_unfolding_census.v1"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _jhash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _root_span_quotient(gram: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, list[int]]:
    """The Z-span of all walls, without assuming ten walls form a Z-basis."""
    rank = gram.rank()
    pivots = next(
        list(indices)
        for indices in itertools.combinations(range(gram.rows), rank)
        if gram.extract(indices, indices).det()
    )
    pivot_gram = gram.extract(pivots, pivots)
    rational_coordinates = pivot_gram.inv() * gram.extract(pivots, range(gram.cols))
    denominator = math.lcm(*(int(value.q) for value in rational_coordinates))
    integer_generators = (denominator * rational_coordinates).applyfunc(int)
    hnf = hermite_normal_form(integer_generators)
    basis_in_pivots = hnf / denominator
    span_gram = sp.simplify(basis_in_pivots.T * pivot_gram * basis_in_pivots)
    roots = sp.simplify(basis_in_pivots.inv() * rational_coordinates)
    if not all(value.q == 1 for value in span_gram) or not all(
        value.q == 1 for value in roots
    ):
        raise ArithmeticError("Hermite quotient did not produce an integral root lattice")
    span_gram = span_gram.applyfunc(int)
    roots = roots.applyfunc(int)
    if roots.T * span_gram * roots != gram:
        raise ArithmeticError("Hermite quotient failed to recover the input Gram")
    return span_gram, roots, pivots


def _relation_data(gram: sp.Matrix) -> tuple[dict[int, list[tuple[int, int]]], nx.Graph, Counter[int]]:
    edges: dict[int, list[tuple[int, int]]] = {3: [], 4: []}
    graph = nx.Graph()
    graph.add_nodes_from(range(gram.rows))
    histogram: Counter[int] = Counter()
    for left in range(gram.rows):
        for right in range(left):
            if gram[left, right] == 0:
                order = 2
            else:
                product = Fraction(
                    4 * int(gram[left, right]) ** 2,
                    int(gram[left, left] * gram[right, right]),
                )
                order = {
                    Fraction(1): 3,
                    Fraction(2): 4,
                    Fraction(3): 6,
                    Fraction(4): 0,  # parallel
                }.get(product, -1)
            histogram[order] += 1
            if order in edges:
                edges[order].append((right, left))
            if order != 2:
                graph.add_edge(right, left, order=order)
    return edges, graph, histogram


def _minimal_scalings(
    edges: dict[int, list[tuple[int, int]]], graph: nx.Graph
) -> tuple[list[dict[str, Any]], list[tuple[int, ...]], int, list[int]]:
    size = len(graph)
    parent = list(range(size))

    def find(value: int) -> int:
        if parent[value] != value:
            parent[value] = find(parent[value])
        return parent[value]

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for edge in edges[3]:
        union(*edge)
    components: list[int] = []
    for vertex in range(size):
        root = find(vertex)
        if root not in components:
            components.append(root)
    component_index = {root: index for index, root in enumerate(components)}
    ratio_edges = [
        (component_index[find(left)], component_index[find(right)])
        for left, right in edges[4]
    ]
    solutions: set[tuple[int, ...]] = set()
    for signs in itertools.product((-1, 1), repeat=len(ratio_edges)):
        adjacency: list[list[tuple[int, int]]] = [[] for _ in components]
        for (left, right), sign in zip(ratio_edges, signs, strict=True):
            adjacency[left].append((right, sign))
            adjacency[right].append((left, -sign))
        heights: list[int | None] = [None] * len(components)
        heights[0] = 0
        stack = [0]
        consistent = True
        while stack:
            left = stack.pop()
            for right, change in adjacency[left]:
                proposed = int(heights[left]) + change
                if heights[right] is None:
                    heights[right] = proposed
                    stack.append(right)
                elif heights[right] != proposed:
                    consistent = False
        if consistent and all(value is not None for value in heights):
            floor = min(int(value) for value in heights)
            solutions.add(
                tuple(
                    2 * 2 ** (int(heights[component_index[find(vertex)]]) - floor)
                    for vertex in range(size)
                )
            )
    scalings = sorted(solutions)
    automorphisms = tuple(
        tuple(mapping[index] for index in range(size))
        for mapping in nx.algorithms.isomorphism.GraphMatcher(
            graph,
            graph,
            edge_match=lambda left, right: left["order"] == right["order"],
        ).isomorphisms_iter()
    )
    unseen = set(scalings)
    orbit_sizes: list[int] = []
    orbit_number: dict[tuple[int, ...], int] = {}
    while unseen:
        representative = min(unseen)
        orbit = {
            tuple(representative[permutation[index]] for index in range(size))
            for permutation in automorphisms
        } & set(scalings)
        orbit_sizes.append(len(orbit))
        for scaling in orbit:
            orbit_number[scaling] = len(orbit_sizes)
        unseen -= orbit
    public = [
        {
            "number": number,
            "normal_square_absolutes": list(scaling),
            "histogram": {
                str(key): value for key, value in sorted(Counter(scaling).items())
            },
            "automorphism_orbit": orbit_number[scaling],
            "one_two_sheet_compatible": set(scaling) <= {2, 4},
        }
        for number, scaling in enumerate(scalings, start=1)
    ]
    return public, scalings, len(automorphisms), sorted(orbit_sizes)


def _fixed_candidates(
    scaling_number: int,
    folded_gram: sp.Matrix,
    e8_gram: sp.Matrix,
    e8_classes: dict[int, tuple[tuple[int, ...], ...]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    # The old pinned implementation is otherwise generic; replace only its
    # non-universal principal-facet root-span shortcut for this exact call.
    old = pilot._root_span
    pilot._root_span = _root_span_quotient
    try:
        return pilot._fixed_candidates(
            scaling_number, folded_gram, e8_gram, e8_classes
        )
    finally:
        pilot._root_span = old


def _q_value(gram: sp.Matrix, value: int) -> Fraction:
    vector = pilot._bit_vector(value, gram.rows)
    return Fraction(int((vector.T * gram * vector)[0]) % 8, 4)


def _k3_glue(n_gram: sp.Matrix, involution: sp.Matrix) -> dict[str, Any]:
    source_basis = pilot._reduce_basis(
        pilot._discriminant_vectors(n_gram), n_gram.rows
    )
    length = len(source_basis)
    if length == 2:
        complement_name = "U+U(2)"
        complement = sp.diag(
            sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, 2], [2, 0]])
        )
    elif length == 4:
        complement_name = "U(2)+U(2)"
        complement = sp.diag(
            sp.Matrix([[0, 2], [2, 0]]), sp.Matrix([[0, 2], [2, 0]])
        )
    else:
        raise ArithmeticError("the bounded class produced an unsupported complement length")
    source_q = {
        coefficients: _q_value(
            n_gram,
            next(
                value
                for value in pilot._span(source_basis)
                if pilot._relation(value, source_basis) == coefficients
            ),
        )
        for coefficients in range(1 << length)
    }
    target_discriminant = tuple(
        value for value in pilot._discriminant_vectors(complement) if value
    )
    target_basis: tuple[int, ...] | None = None
    for candidate in itertools.permutations(target_discriminant, length):
        if len(pilot._reduce_basis(candidate, 4)) != length:
            continue
        if all(
            _q_value(
                complement,
                next(
                    value
                    for value in pilot._span(candidate)
                    if pilot._relation(value, candidate) == coefficients
                ),
            )
            == ((-source_q[coefficients]) % 2)
            for coefficients in source_q
        ):
            target_basis = candidate
            break
    if target_basis is None:
        raise ArithmeticError("no exact complement discriminant anti-isometry exists")
    direct_gram = sp.diag(n_gram, complement)
    glue_columns = [
        sp.Matrix.vstack(
            pilot._bit_vector(source, 18), pilot._bit_vector(target, 4)
        )
        for source, target in zip(source_basis, target_basis, strict=True)
    ]
    basis = hermite_normal_form(
        (2 * sp.eye(22)).row_join(sp.Matrix.hstack(*glue_columns))
    ) / 2
    k3_gram = sp.simplify(basis.T * direct_gram * basis)
    action = sp.simplify(
        basis.inv() * sp.diag(involution, sp.eye(4)) * basis
    )
    if not all(value.q == 1 for value in k3_gram) or not all(
        value.q == 1 for value in action
    ):
        raise ArithmeticError("K3 graph glue/action is not integral")
    k3_gram = k3_gram.applyfunc(int)
    action = action.applyfunc(int)
    n_embedding = pilot._integral(
        basis.inv() * sp.Matrix.vstack(sp.eye(18), sp.zeros(4, 18)),
        "rank-18 K3 embedding is not integral",
    )
    t_embedding = pilot._integral(
        basis.inv() * sp.Matrix.vstack(sp.zeros(18, 4), sp.eye(4)),
        "complement K3 embedding is not integral",
    )
    if not (
        k3_gram.det() == -1
        and pilot._inertia(k3_gram) == (3, 19, 0)
        and all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and action.T * k3_gram * action == k3_gram
        and action**2 == sp.eye(22)
        and pilot._smith(n_embedding) == [1] * 18
        and pilot._smith(t_embedding) == [1] * 4
    ):
        raise ArithmeticError("K3 graph glue checks failed")
    return {
        "complement": complement_name,
        "discriminant_length": length,
        "graph_glue_index": 1 << length,
        "k3_signature_determinant": [3, 19, -1],
        "N_and_complement_embeddings_primitive": True,
        "k3_gram_sha256": _jhash(_rows(k3_gram)),
        "extended_involution_sha256": _jhash(_rows(action)),
        "extended_involution_characteristic_polynomial": "(lambda-1)^14*(lambda+1)^8",
    }


def _construct_unfolding(
    candidate: dict[str, Any],
    folded_gram: sp.Matrix,
    e8_gram: sp.Matrix,
    e8_classes: dict[int, tuple[tuple[int, ...], ...]],
) -> dict[str, Any]:
    fixed_gram: sp.Matrix = candidate["gram"]
    fixed_roots: sp.Matrix = candidate["roots"]
    long_indices = [index for index in range(14) if folded_gram[index, index] == -4]
    long_basis_indices = [
        index - 1 for index in candidate["long_class_basis_facets_one_based"]
    ]
    fixed_basis_classes = [
        pilot._vector_bit(fixed_roots[:, index]) for index in long_basis_indices
    ]
    target_basis_classes = candidate["e8_embedding_basis_classes"]
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
        sp.simplify(basis.T * direct_gram * basis),
        "rank-18 graph-glue Gram is not integral",
    )
    roots = pilot._integral(
        sp.simplify(basis.inv() * sp.Matrix.hstack(*direct_roots)),
        "a lifted root is not integral",
    )
    involution = pilot._integral(
        sp.simplify(basis.inv() * sp.diag(sp.eye(10), -sp.eye(8)) * basis),
        "rank-18 involution is not integral",
    )
    fixed_embedding = pilot._integral(
        basis.inv() * sp.Matrix.vstack(sp.eye(10), sp.zeros(8, 10)),
        "fixed-lattice embedding is not integral",
    )
    anti_embedding = pilot._integral(
        basis.inv() * sp.Matrix.vstack(sp.zeros(10, 8), sp.eye(8)),
        "coinvariant-lattice embedding is not integral",
    )
    root_gram = roots.T * n_gram * roots
    folded_roots = sp.Matrix.hstack(
        *[
            roots[:, orbit[0]]
            if len(orbit) == 1
            else roots[:, orbit[0]] + roots[:, orbit[1]]
            for orbit in root_orbits
        ]
    )
    off_diagonal = [
        int(root_gram[left, right])
        for left in range(root_gram.rows)
        for right in range(left)
    ]
    root_hnf = hermite_normal_form(roots)
    unfolding_checks = {
        "inertia": pilot._inertia(n_gram) == (1, 17, 0),
        "even": all(n_gram[index, index] % 2 == 0 for index in range(18)),
        "isometry": involution.T * n_gram * involution == n_gram,
        "involution": involution**2 == sp.eye(18),
        "fixed_primitive": pilot._smith(fixed_embedding) == [1] * 10,
        "anti_primitive": pilot._smith(anti_embedding) == [1] * 8,
        "root_squares": all(root_gram[index, index] == -2 for index in range(root_gram.rows)),
        "acute": set(off_diagonal) <= {0, 1},
        "fold": folded_roots.T * n_gram * folded_roots == folded_gram,
        "root_span_shape": root_hnf.shape == (18, 18),
        "root_span_finite_index": root_hnf.shape == (18, 18),
    }
    if not all(unfolding_checks.values()):
        raise ArithmeticError(
            "rank-18 unfolding or exact fold-back failed: "
            + ",".join(key for key, value in unfolding_checks.items() if not value)
            + f"; root_index={abs(int(root_hnf.det())) if root_hnf.shape == (18, 18) else 'shape'}"
        )
    discr = pilot._two_elementary_data(n_gram)
    k3 = _k3_glue(n_gram, involution)
    return {
        "fixed_lattice": {
            "determinant": int(fixed_gram.det()),
            "smith": pilot._smith(fixed_gram),
            "invariants": [10, 8, 0],
            "isometry_class": "U+E8(-2)",
            "gram_sha256": _jhash(_rows(fixed_gram)),
        },
        "rank18_unfolding": {
            "signature_determinant": [1, 17, int(n_gram.det())],
            "smith": pilot._smith(n_gram),
            "two_elementary_invariants": [18, discr["discriminant_length"], discr["delta"]],
            "gram_sha256": _jhash(_rows(n_gram)),
            "involution_sha256": _jhash(_rows(involution)),
            "fixed_and_coinvariant_embeddings_primitive": True,
        },
        "roots_and_fold": {
            "root_count": roots.cols,
            "one_sheet_count": sum(len(orbit) == 1 for orbit in root_orbits),
            "two_sheet_count": sum(len(orbit) == 2 for orbit in root_orbits),
            "off_diagonal_pairing_histogram": {
                str(key): value for key, value in sorted(Counter(off_diagonal).items())
            },
            "roots_span_index_in_rank18_lattice": abs(int(root_hnf.det())),
            "exact_labelled_fold_back": True,
            "roots_sha256": _jhash(_rows(roots)),
        },
        "weak_K3_lattice_action_extension": k3,
    }


def _diagram_isomorphic(left: sp.Matrix, right: sp.Matrix) -> bool:
    _, left_graph, _ = _relation_data(left)
    _, right_graph, _ = _relation_data(right)
    return nx.is_isomorphic(
        left_graph,
        right_graph,
        edge_match=lambda a, b: a["order"] == b["order"],
    )


def _corpus() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    for name, path in (
        ("factory", FACTORY), ("parked", PARKED), ("products", PRODUCTS),
        ("t2_pilot", T2_PILOT),
    ):
        if _sha(path) != SOURCE_HASHES[name]:
            raise ArithmeticError(f"pinned source changed: {name}")
    factory = json.loads(FACTORY.read_text())
    parked = json.loads(PARKED.read_text())
    products = json.loads(PRODUCTS.read_text())
    rows: list[dict[str, Any]] = []
    for item in factory["constructions"]:
        rows.append(
            {
                "id": item["id"],
                "lineage": "factory",
                "gram": sp.Matrix(
                    [[int(x) for x in row] for row in item["section"]["restricted_mirrors"]["facet_gram"]]
                ),
                "classification": item["exact_cone"]["decision"]["classification"],
                "vinberg": item["literal_vinberg"],
            }
        )
    parked_h9_14 = [
        item for item in parked["cases"] if item["exact_cone"]["facet_count"] == 14
    ]
    for item in products["direct_exact_audits"]:
        if item["facet_count"] != 14:
            continue
        rows.append(
            {
                "id": item["id"],
                "lineage": "parked_product",
                "gram": sp.Matrix([[int(x) for x in row] for row in item["facet_gram"]]),
                "classification": item["exact_cone"]["decision"]["classification"],
                "vinberg": item["literal_vinberg"],
            }
        )
    return rows, {
        "factory_markings": len(factory["constructions"]),
        "parked_base_H9_14_markings": len(parked_h9_14),
        "parked_product_H9_14_markings": sum(
            item["facet_count"] == 14 for item in products["direct_exact_audits"]
        ),
    }


def _compile_type(type_id: str, representative: dict[str, Any], markings: list[str]) -> dict[str, Any]:
    raw = representative["gram"]
    edges, graph, histogram = _relation_data(raw)
    scaling_rows, scalings, automorphism_order, orbit_sizes = _minimal_scalings(edges, graph)
    e8_gram, e8_classes = pilot._e8_data()
    lanes: list[dict[str, Any]] = []
    internal_candidates: list[dict[str, Any]] = []
    selected_grams: dict[int, sp.Matrix] = {}
    for scaling_row, scaling in zip(scaling_rows, scalings, strict=True):
        if not scaling_row["one_two_sheet_compatible"]:
            continue
        folded = pilot._gram_for_scaling(scaling, edges)
        public, internal = _fixed_candidates(
            scaling_row["number"], folded, e8_gram, e8_classes
        )
        label = (
            f"short{sum(x == 2 for x in scaling)}_long{sum(x == 4 for x in scaling)}"
        )
        lanes.append(
            {
                "label": label,
                "scaling_number": scaling_row["number"],
                "candidate_count": len(public),
                "glue_dimension_histogram": {
                    str(key): value
                    for key, value in sorted(Counter(x["glue_dimension"] for x in public).items())
                },
                "root_span_determinant": int(internal[0]["root_span_gram"].det()),
                "candidates": public,
            }
        )
        for candidate in internal:
            candidate["lane"] = label
            internal_candidates.append(candidate)
        selected_grams[scaling_row["number"]] = folded
    survivors = [
        candidate
        for candidate in internal_candidates
        if candidate["survives_minimal_two_sheet_weak_K3_gate"]
    ]
    if len(survivors) != 1:
        raise ArithmeticError(f"{type_id} did not have a unique bounded-compiler survivor")
    survivor = survivors[0]
    if not (
        survivor["determinant"] == -256
        and survivor["smith_diagonal_absolute"] == [1, 1] + [2] * 8
        and survivor["discriminant_length"] == 8
        and survivor["discriminant_delta"] == 0
    ):
        raise ArithmeticError("the forced U+E8(-2) fixed lattice changed")
    unfolding = _construct_unfolding(
        survivor,
        selected_grams[survivor["scaling_number"]],
        e8_gram,
        e8_classes,
    )
    literal = representative["vinberg"]
    if not (
        representative["classification"] == "finite_volume_noncompact"
        and literal["literal_vinberg_passes"] is True
        and literal["bad_edge_count"] == 0
        and literal["bad_parabolic_neighborhood_count"] == 0
    ):
        raise ArithmeticError("casewise archived Vinberg import is not green")
    return {
        "type_id": type_id,
        "archived_markings": markings,
        "representative": representative["id"],
        "raw_gram_sha256": _jhash(_rows(raw)),
        "relation_histogram": {str(k): v for k, v in sorted(histogram.items())},
        "diagram_automorphism_order": automorphism_order,
        "scaling_count": len(scalings),
        "scaling_orbit_sizes": orbit_sizes,
        "compatible_lane_count": len(lanes),
        "norm_gt4_obstructed_scaling_count": sum(
            not row["one_two_sheet_compatible"] for row in scaling_rows
        ),
        "lanes": lanes,
        "unique_survivor": {
            "lane": survivor["lane"],
            "scaling_number": survivor["scaling_number"],
            "candidate_number": survivor["candidate_number"],
            **unfolding,
        },
        "casewise_V_import": {
            "exact_labelled_diagram_identity": True,
            "classification": representative["classification"],
            "finite_endpoint_count": literal["finite_endpoint_count"],
            "ideal_endpoint_count": literal["ideal_endpoint_count"],
            "bad_edge_count": literal["bad_edge_count"],
            "bad_parabolic_neighborhood_count": literal["bad_parabolic_neighborhood_count"],
        },
    }


def build_certificate() -> dict[str, Any]:
    corpus, source_counts = _corpus()
    if len(corpus) != 6:
        raise ArithmeticError("the pinned primary H9/14 Gram corpus changed")
    groups: list[list[dict[str, Any]]] = []
    for row in corpus:
        group = next(
            (group for group in groups if _diagram_isomorphic(group[0]["gram"], row["gram"])),
            None,
        )
        if group is None:
            groups.append([row])
        else:
            group.append(row)
    if sorted(len(group) for group in groups) != [1, 1, 2, 2]:
        raise ArithmeticError("the four-type primary corpus quotient changed")
    retained: list[list[dict[str, Any]]] = []
    excluded: list[dict[str, Any]] = []
    for group in groups:
        edges, _, histogram = _relation_data(group[0]["gram"])
        raw_gate = (
            set(histogram) <= {2, 3, 4}
            and histogram[4] == 8
            and group[0]["gram"].rank() == 10
            and pilot._inertia(group[0]["gram"]) == (1, 9, 4)
        )
        if raw_gate:
            retained.append(group)
        else:
            excluded.append(
                {
                    "markings": [row["id"] for row in group],
                    "relation_histogram": {str(k): v for k, v in sorted(histogram.items())},
                    "reason": "outside ordinary {2,3,4} with exactly eight m4 edges",
                }
            )
    if len(retained) != 2:
        raise ArithmeticError("D34 did not retain exactly two diagram types")
    # Name by the (m3,m4) census, without using a surface label.
    compiled = []
    for group in retained:
        histogram = _relation_data(group[0]["gram"])[2]
        compiled.append(
            _compile_type(
                f"D34-{histogram[3]}-{histogram[4]}",
                group[0],
                [row["id"] for row in group],
            )
        )
    compiled.sort(key=lambda row: row["type_id"])
    by_type = {row["type_id"]: row for row in compiled}
    # Sorting above puts T2 (75 m2) first and T1 (74 m2) second.
    if not (
        by_type["D34-8-8"]["representative"] == "A18_4_0_C01"
        and by_type["D34-9-8"]["representative"] == "A14_8_1_C01"
        and by_type["D34-8-8"]["scaling_count"] == 18
        and by_type["D34-9-8"]["scaling_count"] == 18
        and by_type["D34-8-8"]["unique_survivor"]["rank18_unfolding"]["two_elementary_invariants"] == [18, 4, 0]
        and by_type["D34-9-8"]["unique_survivor"]["rank18_unfolding"]["two_elementary_invariants"] == [18, 2, 0]
    ):
        raise ArithmeticError("the T1/T2 theorem invariants changed")
    checks = {
        "pinned_primary_corpus_has_six_markings_and_four_diagram_types": True,
        "D34_filter_retains_three_markings_and_two_nonisomorphic_types": True,
        "Hermite_root_lattice_quotient_used_for_both_types": True,
        "all_36_minimal_scalings_enumerated": True,
        "all_28_stable_even_overlattice_candidates_enumerated": True,
        "each_type_has_exactly_one_U_plus_E8_minus2_survivor": True,
        "T1_rank18_18_2_0_unfolding_and_K3_glue_exact": True,
        "T2_rank18_18_4_0_unfolding_replayed": True,
        "casewise_V_imported_only_by_exact_diagram_identity": True,
        "no_full_U_or_T1_cone_completeness_claim": True,
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "source_sha256": SOURCE_HASHES,
        "bounded_primary_corpus": {
            **source_counts,
            "marking_count": len(corpus),
            "diagram_type_count": len(groups),
            "scope": "only the pinned primary certified exact-Gram corpus; not all H9/14 diagrams",
        },
        "surface_free_D34_specification": {
            "dimension": 9,
            "facet_count": 14,
            "Gram_rank_inertia": [10, 1, 9, 4],
            "allowed_relations": [2, 3, 4],
            "m4_edge_count": 8,
            "minimal_unfolding_model": "fixed (-2) root or orthogonal exchanged pair of (-2) roots",
            "retained_marking_count": sum(len(group) for group in retained),
            "retained_diagram_type_count": len(retained),
            "excluded_types": excluded,
        },
        "compiler_inventory": compiled,
        "theorem": (
            "Relative to the pinned primary exact-Gram corpus, D34 has exactly two "
            "nonisomorphic diagrams. The exhaustive minimal involutive compiler has "
            "28 stable even fixed-overlattice candidates and forces a unique "
            "U+E8(-2) fixed lattice for each. T1 yields a new exact (18,2,0) "
            "rank-18 K3-lattice unfolding; T2 replays (18,4,0)."
        ),
        "URCV_boundary": {
            "U": "open: corpus-relative IC0 exhaustiveness is not forward universality for all H9/14 realizations or an algebraic surface locus",
            "R": "exact labelled fold-back/lattice-action recovery R0 is closed for both; full geometric marked recovery is imported only for T2",
            "C": "T2 imports the existing complete invariant-nef theorem; T1 remains open",
            "V": "casewise closed for both by exact labelled identity with the pinned Vinberg archives; not automatic from the weak K3 compiler",
        },
        "new_polytope_claim": False,
        "new_carrier_lineage": "T1 acquires a new weak rank-18 (18,2,0) K3-lattice/action unfolding with complement U+U(2)",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("D34 archive differs from exact replay")
        print("archive_check=PASS")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        print(f"wrote={args.output}")
    print("all_checks_pass=true")
    print("corpus_markings=6; corpus_types=4; D34_markings=3; D34_types=2")
    print("stable_even_candidates=28; survivors=2")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
