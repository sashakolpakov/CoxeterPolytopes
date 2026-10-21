#!/usr/bin/env python3
"""Exact lineage test for the promising Bredon N9 and N11 triality rows.

This does not search for a new K3 wall stock.  It tests whether the two
Lorentzian simply-laced unfoldings already occur among the classified
complete high-Picard-rank K3 nef-wall stocks, and, when they do, verifies the
integral lattice and triality action rather than stopping at an unweighted
graph match.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import networkx as nx
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
BREDON = ROOT / "data/enumerations/bredon_adeg24_triality_census.json"
CLASSIFICATION = ROOT / "data/enumerations/k3_high_rank_full_root_classification.json"
VOLUME = ROOT / "data/enumerations/k3_rank11_14_literal_vinberg_volume.json"
INVENTORY = ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
R12 = ROOT / "data/enumerations/roulleau_non2elementary_involution_folds.json"
R14 = ROOT / "data/enumerations/roulleau_ue8d4_all_subgroup_folds.json"
ORDER6 = ROOT / "data/enumerations/ue8d4_order6_weierstrass_family.json"


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _integer_rows(matrix: sp.Matrix) -> list[list[int]]:
    if not all(value.is_Integer for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _target_orientation(
    payload: dict[str, Any], row_id: str, orientation: int
) -> tuple[dict[str, Any], dict[str, Any]]:
    row = next(record for record in payload["rows"] if record["id"] == row_id)
    oriented = next(
        record for record in row["orientations"] if record["orientation"] == orientation
    )
    return row, oriented


def _gram_graph(gram: sp.Matrix) -> nx.Graph:
    if gram != gram.T or any(gram[i, i] <= 0 for i in range(gram.rows)):
        raise ArithmeticError("expected a symmetric positive-normal Gram matrix")
    graph = nx.Graph()
    for vertex in range(gram.rows):
        graph.add_node(vertex, norm=int(gram[vertex, vertex]))
    for left in range(gram.rows):
        for right in range(left + 1, gram.rows):
            value = int(gram[left, right])
            if value:
                graph.add_edge(left, right, pairing=value)
    return graph


def _isomorphisms(left: sp.Matrix, right: sp.Matrix) -> list[tuple[int, ...]]:
    left_graph, right_graph = _gram_graph(left), _gram_graph(right)
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        left_graph,
        right_graph,
        node_match=lambda first, second: first["norm"] == second["norm"],
        edge_match=lambda first, second: first["pairing"] == second["pairing"],
    )
    mappings = {
        tuple(mapping[index] for index in range(left.rows))
        for mapping in matcher.isomorphisms_iter()
    }
    return sorted(mappings)


def _graph_invariants(gram: sp.Matrix) -> dict[str, Any]:
    graph = _gram_graph(gram)
    component_count = nx.number_connected_components(graph)
    automorphism_order = sum(
        1
        for _ in nx.algorithms.isomorphism.GraphMatcher(
            graph,
            graph,
            node_match=lambda first, second: first["norm"] == second["norm"],
            edge_match=lambda first, second: first["pairing"] == second["pairing"],
        ).isomorphisms_iter()
    )
    return {
        "vertices": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "connected_components": component_count,
        "cycle_rank": graph.number_of_edges() - graph.number_of_nodes() + component_count,
        "triangle_count": sum(nx.triangles(graph).values()) // 3,
        "degree_histogram": {
            str(degree): count
            for degree, count in sorted(Counter(dict(graph.degree()).values()).items())
        },
        "weighted_graph_automorphism_order": automorphism_order,
    }


def _wall_coordinates(gram: sp.Matrix, basis_indices: Iterable[int]) -> sp.Matrix:
    indices = tuple(basis_indices)
    basis_gram = gram.extract(indices, indices)
    if not basis_gram.det():
        raise ArithmeticError("selected wall basis is singular")
    coordinates = gram[:, indices] * basis_gram.inv()
    if not all(value.is_Integer for value in coordinates):
        raise ArithmeticError("wall coordinates are not integral")
    return coordinates


def _unfolded_orbits(vertex_count: int, long_vertices: set[int]) -> list[list[int]]:
    result: list[list[int]] = []
    cursor = 0
    for vertex in range(vertex_count):
        size = 3 if vertex in long_vertices else 1
        result.append(list(range(cursor, cursor + size)))
        cursor += size
    return result


def _relation_matrix(gram: sp.Matrix) -> list[list[str]]:
    result = [["1" if i == j else "2" for j in range(gram.cols)] for i in range(gram.rows)]
    labels = {1: "3", 2: "4", 3: "6"}
    for left in range(gram.rows):
        for right in range(left + 1, gram.cols):
            numerator = 4 * gram[left, right] ** 2
            denominator = gram[left, left] * gram[right, right]
            product = sp.Rational(numerator, denominator)
            if product == 0:
                label = "2"
            elif product in (1, 2, 3):
                label = labels[int(product)]
            else:
                raise ArithmeticError(f"non-Coxeter pair in folded Gram: {product}")
            result[left][right] = result[right][left] = label
    return result


def _unique_stock_candidates(
    classification: dict[str, Any], rank: int, facets: int
) -> list[dict[str, Any]]:
    return [
        {
            "key": row["key"],
            "lattice": row["lattice"],
            "picard_rank": row["picard_rank"],
            "facet_count": row["facet_count"],
            "complete_root_chamber": row["complete_root_chamber"],
        }
        for row in classification["cases"]
        if row["picard_rank"] == rank and row["facet_count"] == facets
    ]


def _source_cross_checks(
    key: str,
    source_gram: list[list[int]],
    volume: dict[str, Any],
    inventory: dict[str, Any],
) -> dict[str, Any]:
    volume_row = next(row for row in volume["cases"] if row["key"] == key)
    inventory_row = next(row for row in inventory["stocks"] if row["key"] == key)
    gram_hash = _canonical_sha256(source_gram)
    if not (
        volume_row["curve_gram_sha256"] == gram_hash
        and volume_row["simple_roots_complete"]
        and inventory_row["current_machine_ledger_status"]
        == "represented_exact_full_nef_stock"
        and inventory_row["wall_count"] == len(source_gram)
    ):
        raise ArithmeticError(f"{key}: complete-stock provenance did not cross-check")
    return {
        "source_gram_sha256": gram_hash,
        "literal_vinberg_decision": volume_row["decision"],
        "literal_vinberg_finite_vertices": volume_row["finite_vertex_count"],
        "literal_vinberg_cusps": volume_row["ideal_vertex_count"],
        "literal_vinberg_elliptic_edges": volume_row["elliptic_edge_count"],
        "simple_roots_complete": volume_row["simple_roots_complete"],
        "inventory_status": inventory_row["current_machine_ledger_status"],
        "matrix_provenance": volume_row["matrix_provenance"],
    }


def _analyze_n9(
    bredon: dict[str, Any],
    classification: dict[str, Any],
    volume: dict[str, Any],
    inventory: dict[str, Any],
    r12: dict[str, Any],
) -> dict[str, Any]:
    row, oriented = _target_orientation(bredon, "ADEG-N9-01-branched-path", 0)
    target = sp.Matrix(oriented["unfolded_gram"])
    source_case = next(case for case in r12["cases"] if case["id"] == "R12_UA2E8_C2_FIXED")
    source_negative = sp.Matrix(source_case["source"]["curve_gram"])
    source = -source_negative
    candidates = _unique_stock_candidates(classification, target.rank(), target.rows)
    mappings = _isomorphisms(target, source)

    target_basis_indices = oriented["root_generated_lattice"]["independent_root_indices"]
    target_lattice_determinant = int(target.extract(target_basis_indices, target_basis_indices).det())
    source_basis_labels = [label - 1 for label in source_case["source"]["basis_wall_labels"]]
    source_lattice_determinant = int(
        source_negative.extract(source_basis_labels, source_basis_labels).det()
    )
    source_coordinates = _wall_coordinates(source_negative, source_basis_labels)
    if [int(value) for value in source_coordinates.row(source_negative.rows - 1)] != source_case[
        "source"
    ]["dependent_last_wall_coordinates"]:
        raise ArithmeticError("R12 source dependent wall stopped reconstructing")

    target_invariants = _graph_invariants(target)
    source_invariants = _graph_invariants(source)
    if not (
        candidates == [
            {
                "key": "R12_UA2E8",
                "lattice": "U+A2(-1)+E8(-1)",
                "picard_rank": 12,
                "facet_count": 13,
                "complete_root_chamber": True,
            }
        ]
        and not mappings
        and target_invariants["degree_histogram"] != source_invariants["degree_histogram"]
        and target_lattice_determinant == -12
        and source_lattice_determinant == -3
    ):
        raise ArithmeticError("N9 lineage separation changed")

    return {
        "bredon_row": row["id"],
        "orientation": oriented["orientation"],
        "folded_dimension_facets": [row["dimension"], row["facet_count"]],
        "unfolded_rank_root_count": [target.rank(), target.rows],
        "classified_stock_candidates_with_same_rank_and_wall_count": candidates,
        "tested_stock": "R12_UA2E8",
        "target_graph_invariants": target_invariants,
        "stock_graph_invariants": source_invariants,
        "exact_weighted_gram_isomorphism_count": len(mappings),
        "target_root_lattice_determinant": target_lattice_determinant,
        "stock_full_NS_determinant": source_lattice_determinant,
        "integral_lattice_isometry": False,
        "triality_compatible_stock_automorphism": False,
        "failure_certificate": [
            "degree histogram differs (the target has a degree-four vertex)",
            "target graph is a tree while the stock graph contains one triangle",
            "wall-lattice determinants are -12 and -3",
            "target graph automorphism group has order divisible by 3; stock graph group has order 2",
        ],
        "complete_stock_provenance": _source_cross_checks(
            "R12_UA2E8", source_case["source"]["curve_gram"], volume, inventory
        ),
        "verdict": (
            "no match to the unique classified complete rank-12/13-wall K3 nef stock; "
            "the N9 same-P carrier remains open"
        ),
    }


def _analyze_n11(
    bredon: dict[str, Any],
    classification: dict[str, Any],
    volume: dict[str, Any],
    inventory: dict[str, Any],
    r14: dict[str, Any],
    order6: dict[str, Any],
) -> dict[str, Any]:
    row, oriented = _target_orientation(bredon, "ADEG-N11-01-branched-path", 1)
    target = sp.Matrix(oriented["unfolded_gram"])
    source_negative = sp.Matrix(r14["source"]["wall_gram"])
    source = -source_negative
    candidates = _unique_stock_candidates(classification, target.rank(), target.rows)
    mappings = _isomorphisms(target, source)
    if not mappings:
        raise ArithmeticError("N11 unfolded graph stopped matching R14_UE8D4")
    mapping = mappings[0]
    expected_mapping = (1, 2, 4, 5, 6, 3, 7, 9, 8, 0, 10, 11, 12, 13, 14)
    if mapping != expected_mapping:
        raise ArithmeticError(f"canonical N11 wall mapping changed: {mapping}")
    if any(target[i, j] != source[mapping[i], mapping[j]] for i in range(15) for j in range(15)):
        raise ArithmeticError("N11 mapping is not entrywise Gram preserving")

    target_basis_indices = tuple(oriented["root_generated_lattice"]["independent_root_indices"])
    source_basis_gram = sp.Matrix(r14["source"]["basis_gram"])
    source_wall_rows = sp.Matrix(r14["source"]["wall_coordinate_rows"])
    basis_change = source_wall_rows[[mapping[index] for index in target_basis_indices], :]
    target_basis_gram = target.extract(target_basis_indices, target_basis_indices)
    target_wall_rows = _wall_coordinates(target, target_basis_indices)
    if not (
        abs(int(basis_change.det())) == 1
        and basis_change * source_basis_gram * basis_change.T == -target_basis_gram
        and target_wall_rows * basis_change == source_wall_rows[list(mapping), :]
    ):
        raise ArithmeticError("N11 integral wall-lattice isometry failed")

    # The Bredon long wall unfolds into the last three target roots.  Under
    # the wall isometry these are exactly the three outer D4 roots.
    target_orbits = _unfolded_orbits(row["facet_count"], set(oriented["long_vertices"]))
    c3_class = next(
        record for record in r14["subgroup_census"]["classes"] if record["order"] == 3
    )
    source_orbits = [
        [wall - 1 for wall in orbit] for orbit in c3_class["wall_orbits_one_based"]
    ]
    folded_mapping: list[int] = []
    for orbit in target_orbits:
        image = {mapping[index] for index in orbit}
        matches = [index for index, source_orbit in enumerate(source_orbits) if image == set(source_orbit)]
        if len(matches) != 1:
            raise ArithmeticError("N11 unfolded isomorphism does not respect triality orbits")
        folded_mapping.append(matches[0])

    target_folded = sp.Matrix(oriented["folded_integral_gram"])
    type_row = next(row for row in r14["types"] if row["type_id"] == "R14_UE8D4_C3_S3")
    source_orbit_gram = -sp.Matrix(type_row["geometry"]["orbit_gram"])
    if any(
        target_folded[i, j]
        != source_orbit_gram[folded_mapping[i], folded_mapping[j]]
        for i in range(target_folded.rows)
        for j in range(target_folded.cols)
    ):
        raise ArithmeticError("N11 folded Gram comparison failed")
    source_relations = type_row["normalized_relation_matrix"]
    target_relations = _relation_matrix(target_folded)
    if any(
        target_relations[i][j]
        != source_relations[folded_mapping[i]][folded_mapping[j]]
        for i in range(target_folded.rows)
        for j in range(target_folded.cols)
    ):
        raise ArithmeticError("N11 folded Coxeter relation comparison failed")

    # Literal integral C3 equivariance on the rank-fourteen source lattice.
    target_generator = tuple(range(12)) + (13, 14, 12)
    source_generator = tuple(
        value - 1 for value in type_row["ns_action"]["generator_permutations_one_based"][0]
    )
    if any(mapping[target_generator[i]] != source_generator[mapping[i]] for i in range(15)):
        raise ArithmeticError("N11 graph isomorphism is not C3-equivariant")
    target_action = target_wall_rows[
        [target_generator[index] for index in target_basis_indices], :
    ]
    source_action = sp.Matrix(type_row["ns_action"]["generator_action_matrices"][0])
    if not (
        target_action**3 == sp.eye(14)
        and target_action * basis_change == basis_change * source_action
        and type_row["ns_action"]["all_integral_lattice_checks_pass"]
    ):
        raise ArithmeticError("N11 integral C3 action comparison failed")

    # Also exhibit a determinant-one isometry of the rank-twelve folded
    # lattices, not just the ambient rank-fourteen wall spans.
    target_fixed_basis = tuple(range(12))
    target_fixed_rows = _wall_coordinates(target_folded, target_fixed_basis)
    source_fixed_rows = sp.Matrix(
        [record["fixed_coordinates"] for record in type_row["reflection_audit"]["checks"]]
    )
    fixed_basis_change = source_fixed_rows[
        [folded_mapping[index] for index in target_fixed_basis], :
    ]
    source_fixed_gram = sp.Matrix(type_row["fixed_lattice"]["gram"])
    if not (
        abs(int(fixed_basis_change.det())) == 1
        and fixed_basis_change * source_fixed_gram * fixed_basis_change.T
        == -target_folded[:12, :12]
        and target_fixed_rows * fixed_basis_change
        == source_fixed_rows[folded_mapping, :]
    ):
        raise ArithmeticError("N11 fixed-lattice change of basis failed")

    handoff = order6["hyperbolic_handoff"]
    if not (
        order6["all_checks_pass"]
        and order6["family"]["global_order"] == 6
        and order6["brandhorst_hofmann_handoff"]["archived_type_id"]
        == "R14_UE8D4_C3_S3"
        and order6["brandhorst_hofmann_handoff"]["same_d4_action_entrywise"]
        and order6["picard_action"]["fixed_lattice"]["class"]
        == "U+E8(-1)+A2(-1)"
        and handoff == {
            "abstract_polytope_novelty": "known; this is a new explicit realization in the project",
            "bad_edges": 0,
            "classification": "finite_volume_noncompact",
            "cusps": 1,
            "dimension": 11,
            "elliptic_edges": 162,
            "facets": 13,
            "finite_vertices": 27,
        }
    ):
        raise ArithmeticError("N11 equation-level K3 handoff changed")

    if candidates != [
        {
            "key": "R14_UE8D4",
            "lattice": "U+E8(-1)+D4(-1)",
            "picard_rank": 14,
            "facet_count": 15,
            "complete_root_chamber": True,
        }
    ]:
        raise ArithmeticError("N11 classified-stock denominator changed")

    return {
        "bredon_row": row["id"],
        "orientation": oriented["orientation"],
        "folded_dimension_facets": [row["dimension"], row["facet_count"]],
        "unfolded_rank_root_count": [target.rank(), target.rows],
        "classified_stock_candidates_with_same_rank_and_wall_count": candidates,
        "tested_stock": "R14_UE8D4",
        "target_graph_invariants": _graph_invariants(target),
        "stock_graph_invariants": _graph_invariants(source),
        "exact_weighted_gram_isomorphism_count": len(mappings),
        "canonical_target_to_source_wall_permutation_one_based": [value + 1 for value in mapping],
        "long_target_wall_unfolds_to_target_roots_one_based": [value + 1 for value in target_orbits[12]],
        "their_source_images_one_based": [mapping[value] + 1 for value in target_orbits[12]],
        "source_c3_orbit_one_based": [value + 1 for value in source_orbits[12]],
        "triality_orbits_respected": True,
        "folded_target_to_source_facet_permutation_one_based": [
            value + 1 for value in folded_mapping
        ],
        "folded_integral_gram_entrywise_equal": True,
        "folded_normalized_relation_matrix_entrywise_equal": True,
        "unfolded_basis_change_determinant": int(basis_change.det()),
        "unfolded_basis_change_target_to_source": _integer_rows(basis_change),
        "unfolded_integral_lattice_isometry": True,
        "c3_integral_action_entrywise_intertwined": True,
        "fixed_lattice_basis_change_determinant": int(fixed_basis_change.det()),
        "fixed_lattice_basis_change_target_to_source": _integer_rows(fixed_basis_change),
        "fixed_lattice_determinant": int(source_fixed_gram.det()),
        "fixed_lattice": "U+E8(-1)+A2(-1)",
        "equation_level_surface": {
            "family": order6["family"]["affine_equation"],
            "automorphism": order6["family"]["automorphism"],
            "global_order": order6["family"]["global_order"],
            "order_on_NS": 3,
            "same_d4_action_entrywise": order6["brandhorst_hofmann_handoff"][
                "same_d4_action_entrywise"
            ],
        },
        "same_polytope_vinberg_certificate": handoff,
        "complete_stock_provenance": _source_cross_checks(
            "R14_UE8D4", r14["source"]["wall_gram"], volume, inventory
        ),
        "known_lineage": {
            "roulleau_type": "R14_UE8D4_C3_S3",
            "scharlau_kirschmer": "List11 entry 4",
            "coxiter": "graphs/Tumarkin pyramids n+2/11-tum04_12_05.coxiter",
            "reflection_subgroup_index": 1,
        },
        "verdict": (
            "exact Bredon-N11 -> complete K3 nef stock -> C3 fixed section -> "
            "same Bredon polytope return; known abstract H11/13 chamber"
        ),
    }


def analyze() -> dict[str, Any]:
    bredon = _load(BREDON)
    classification = _load(CLASSIFICATION)
    volume = _load(VOLUME)
    inventory = _load(INVENTORY)
    r12 = _load(R12)
    r14 = _load(R14)
    order6 = _load(ORDER6)
    result = {
        "certificate": "bredon_n9_n11_unfolded_k3_stock_lineage",
        "scope": (
            "exact comparison with the classified complete finite-volume K3 root "
            "chambers of Picard ranks 11 through 20"
        ),
        "n9": _analyze_n9(bredon, classification, volume, inventory, r12),
        "n11": _analyze_n11(
            bredon, classification, volume, inventory, r14, order6
        ),
        "summary": {
            "exact_stock_matches": 1,
            "exact_stock_failures": 1,
            "new_abstract_polytope_count": 0,
            "new_named_same_P_polytope_to_surface_identifications": 1,
        },
    }
    result["canonical_sha256"] = _canonical_sha256(result)
    return result


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
