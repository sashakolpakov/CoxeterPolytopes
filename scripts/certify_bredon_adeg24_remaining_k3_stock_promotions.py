#!/usr/bin/env python3
"""Compare the remaining stable ADEG triality unfoldings with K3 nef stocks.

The input side is polytope-first: all stable-range packages in the frozen
ADEG24 triality census, except the already closed N11 promotion and the
separately audited P_star row.  The comparison denominator contains the 17
Roulleau rank-5--8 complete stocks selected by the target headers and every
complete K3 nef-root stock currently frozen locally in ranks 9--19.

Rank/root-count, Smith invariants, and exact weighted-Gram isomorphism are
successive gates.  An exact match is promoted only after checking an integral
wall-lattice isometry, the transported C3 action, an invariant ample class,
the complete-nef/no-hidden-wall source certificate, and exact orbit-sum
return to the original folded Gram.

The transported C3 is a divisor-lattice/nef-cone symmetry.  This script does
not assert that it is induced by a K3 automorphism or construct a quotient
surface; that stronger Hodge/Torelli lift remains a separate gate.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/bredon_adeg24_remaining_k3_stock_promotions.json"
SCHEMA = "coxeter.bredon_adeg24_remaining_k3_stock_promotions.v2"

BREDON = ROOT / "data/enumerations/bredon_adeg24_triality_census.json"
LOWER = ROOT / "data/enumerations/roulleau_rank5_8_matching_full_nef_stocks.json"
RANK9 = ROOT / "data/enumerations/roulleau_rank9_all_full_nef.json"
RANK10_SIMPLEX = ROOT / "data/enumerations/roulleau_rank10_two_full_nef_simplices.json"
RANK10_REMAINING = ROOT / "data/enumerations/roulleau_rank10_remaining_full_nef.json"
HIGH = ROOT / "data/enumerations/k3_high_rank_full_root_classification.json"
R12_UA2E8 = ROOT / "data/enumerations/roulleau_non2elementary_involution_folds.json"
R12_GRAPH = ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"

EXCLUDED_ROWS = {
    "ADEG-N11-01-branched-path": "already closed by the N9/N11 stock-lineage theorem",
    "ADEG-N9-02-P-star": "separately closed by the P_star K3/triality audit",
}

SOURCE_AMPLE_COORDINATES = {
    "R9-UD4A3-NEF-H8-10": (7, 7, 7, 15, 10, 6, 3, 1, 0),
    "R10-UD4D4-NEF-H9-11": (7, 15, 7, 7, 10, 6, 3, 1, 0, 0),
}


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _integer_rows(matrix: sp.Matrix) -> list[list[int]]:
    if not all(value.is_Integer for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _nonzero_smith(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.shape))
        if smith[index, index] != 0
    ]


def _gram_graph(gram: sp.Matrix) -> nx.Graph:
    if gram != gram.T or any(gram[index, index] != 2 for index in range(gram.rows)):
        raise ArithmeticError("expected a symmetric square-two positive-normal Gram")
    graph = nx.Graph()
    graph.add_nodes_from(range(gram.rows))
    for left in range(gram.rows):
        for right in range(left + 1, gram.rows):
            value = int(gram[left, right])
            if value:
                graph.add_edge(left, right, pairing=value)
    return graph


def _graph_invariants(gram: sp.Matrix) -> dict[str, Any]:
    graph = _gram_graph(gram)
    components = nx.number_connected_components(graph)
    return {
        "vertices": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "connected_components": components,
        "cycle_rank": graph.number_of_edges() - graph.number_of_nodes() + components,
        "triangle_count": sum(nx.triangles(graph).values()) // 3,
        "degree_histogram": {
            str(degree): count
            for degree, count in sorted(Counter(dict(graph.degree()).values()).items())
        },
    }


def _isomorphisms(left: sp.Matrix, right: sp.Matrix) -> list[tuple[int, ...]]:
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _gram_graph(left),
        _gram_graph(right),
        edge_match=lambda first, second: first["pairing"] == second["pairing"],
    )
    return sorted(
        {
            tuple(mapping[index] for index in range(left.rows))
            for mapping in matcher.isomorphisms_iter()
        }
    )


def _wall_coordinates(gram: sp.Matrix, basis_indices: Iterable[int]) -> sp.Matrix:
    indices = tuple(basis_indices)
    basis = gram.extract(indices, indices)
    if not basis.det():
        raise ArithmeticError("selected wall basis is singular")
    coordinates = gram[:, indices] * basis.inv()
    if not all(value.is_Integer for value in coordinates):
        raise ArithmeticError("wall coordinates are not integral")
    return coordinates


def _independent_principal_indices(gram: sp.Matrix) -> tuple[int, ...]:
    rank = gram.rank()
    pivots = tuple(int(value) for value in gram.rref()[1])
    if len(pivots) == rank and gram.extract(pivots, pivots).det() != 0:
        return pivots
    # All exact-promotion headers are small.  The fallback is retained for
    # denominator rows whose RREF column pivots do not form a principal minor.
    import itertools

    for indices in itertools.combinations(range(gram.rows), rank):
        if gram.extract(indices, indices).det() != 0:
            return tuple(indices)
    raise ArithmeticError("symmetric wall Gram has no full-rank principal minor")


def _unfolded_orbits(facet_count: int, long_vertices: set[int]) -> list[list[int]]:
    result: list[list[int]] = []
    cursor = 0
    for vertex in range(facet_count):
        size = 3 if vertex in long_vertices else 1
        result.append(list(range(cursor, cursor + size)))
        cursor += size
    return result


def _orbit_sum_matrix(orbits: list[list[int]], root_count: int) -> sp.Matrix:
    matrix = sp.zeros(len(orbits), root_count)
    for orbit_index, orbit in enumerate(orbits):
        for root in orbit:
            matrix[orbit_index, root] = 1
    return matrix


def _triality_permutation(orbits: list[list[int]], root_count: int) -> tuple[int, ...]:
    permutation = list(range(root_count))
    for orbit in orbits:
        if len(orbit) == 3:
            permutation[orbit[0]], permutation[orbit[1]], permutation[orbit[2]] = (
                orbit[1],
                orbit[2],
                orbit[0],
            )
        elif len(orbit) != 1:
            raise ArithmeticError("triality orbit has size other than one or three")
    return tuple(permutation)


def _low_stock_records(payload: dict[str, Any], source: Path) -> list[dict[str, Any]]:
    if not payload["all_checks_pass"]:
        raise ArithmeticError(f"{source.name}: source archive does not pass")
    records = []
    for case in payload["cases"]:
        if not case["all_checks_pass"] or "full_nef_K3" not in case["classification"]:
            raise ArithmeticError(f"{case['id']}: source is not a complete K3 nef stock")
        certificate = case["matrix_certificate"]
        negative = sp.Matrix(certificate["wall_gram"])
        rank = int(case["picard_rank"])
        if not (
            negative.rank() == rank
            and exact_inertia(negative) == (1, rank - 1, negative.rows - rank)
            and all(negative[index, index] == -2 for index in range(negative.rows))
        ):
            raise ArithmeticError(f"{case['id']}: source wall Gram failed")
        smith = _nonzero_smith(negative)
        stored_smith = certificate.get(
            "basis_smith_invariants", certificate.get("smith_invariants")
        )
        if smith != stored_smith:
            raise ArithmeticError(f"{case['id']}: source Smith data changed")
        basis_raw = certificate.get("basis_indices_one_based")
        if basis_raw is None:
            possible = certificate.get("integral_basis")
            basis_raw = possible if isinstance(possible, list) else None
        basis_indices = (
            tuple(int(value) - 1 for value in basis_raw)
            if basis_raw is not None
            else (tuple(range(rank)) if negative.rows == rank else None)
        )
        if basis_indices is not None and len(basis_indices) != rank:
            raise ArithmeticError(f"{case['id']}: source basis size changed")
        records.append(
            {
                "stock_id": case["id"],
                "rank": rank,
                "root_count": int(case["facet_count"]),
                "lattice": case["lattice"],
                "smith": smith,
                "negative_gram": negative,
                "basis_indices": basis_indices,
                "archive": source.name,
                "archive_wall_gram_sha256": certificate["wall_gram_sha256"],
                "full_nef_reason": case["complete_full_nef_reason"],
                "literal_vinberg": case["literal_vinberg"],
                "source_ample": certificate.get(
                    "ample_coordinates_in_basis",
                    SOURCE_AMPLE_COORDINATES.get(case["id"]),
                ),
            }
        )
    return records


def _high_exact_grams() -> dict[str, tuple[sp.Matrix, tuple[int, ...]]]:
    ua2e8 = _load(R12_UA2E8)
    case = next(row for row in ua2e8["cases"] if row["id"] == "R12_UA2E8_C2_FIXED")
    ua2e8_gram = sp.Matrix(case["source"]["curve_gram"])
    ua2e8_basis = tuple(int(value) - 1 for value in case["source"]["basis_wall_labels"])

    graph = _load(R12_GRAPH)
    ue8a1 = next(row for row in graph["stocks"] if row["stock_key"] == "R12_UE8A1x2")
    ue8a1_gram = sp.Matrix(ue8a1["wall_gram"])
    # The archived wall-coordinate rows use the first twelve walls as basis.
    ue8a1_basis = tuple(range(12))
    if not (ua2e8["all_checks_pass"] and graph["all_checks_pass"] and ue8a1["all_checks_pass"]):
        raise ArithmeticError("rank-twelve exact source archives failed")
    return {
        "R12_UA2E8": (ua2e8_gram, ua2e8_basis),
        "R12_UE8A1x2": (ue8a1_gram, ue8a1_basis),
    }


def _stock_denominator() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    records.extend(_low_stock_records(_load(LOWER), LOWER))
    records.extend(_low_stock_records(_load(RANK9), RANK9))
    records.extend(_low_stock_records(_load(RANK10_SIMPLEX), RANK10_SIMPLEX))
    records.extend(_low_stock_records(_load(RANK10_REMAINING), RANK10_REMAINING))

    classification = _load(HIGH)
    exact = _high_exact_grams()
    if not classification["all_checks_pass"] or len(classification["cases"]) != 19:
        raise ArithmeticError("high-rank complete-root classification changed")
    for case in classification["cases"]:
        if not case["complete_root_chamber"]:
            raise ArithmeticError(f"{case['key']}: high-rank stock is not complete")
        negative = basis_indices = None
        smith = None
        if case["key"] in exact:
            negative, basis_indices = exact[case["key"]]
            smith = _nonzero_smith(negative)
            if not (
                negative.rank() == case["picard_rank"]
                and len(smith) == case["picard_rank"]
                and abs(int(sp.prod(smith))) == case["absolute_discriminant"]
            ):
                raise ArithmeticError(f"{case['key']}: exact high-rank Gram changed")
        records.append(
            {
                "stock_id": case["key"],
                "rank": int(case["picard_rank"]),
                "root_count": int(case["facet_count"]),
                "lattice": case["lattice"],
                "smith": smith,
                "negative_gram": negative,
                "basis_indices": basis_indices,
                "archive": HIGH.name,
                "archive_wall_gram_sha256": (
                    _canonical_sha256(_integer_rows(negative)) if negative is not None else None
                ),
                "full_nef_reason": case["complete_root_evidence"],
                "literal_vinberg": {"decision": case["literal_vinberg_status"]},
                "source_ample": None,
            }
        )
    if len(records) != 55 or len({record["stock_id"] for record in records}) != 55:
        raise ArithmeticError("local complete-stock denominator is not 55 distinct rows")
    return records


def _public_stock(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "stock_id": record["stock_id"],
        "rank": record["rank"],
        "root_count": record["root_count"],
        "lattice": record["lattice"],
        "smith": record["smith"],
        "exact_gram_locally_loaded": record["negative_gram"] is not None,
        "archive": record["archive"],
        "archive_wall_gram_sha256": record["archive_wall_gram_sha256"],
    }


def _promotion(
    row: dict[str, Any],
    orientation: dict[str, Any],
    target: sp.Matrix,
    stock: dict[str, Any],
    isomorphisms: list[tuple[int, ...]],
) -> dict[str, Any]:
    if not isomorphisms:
        raise ArithmeticError("promotion called without a Gram isomorphism")
    mapping = isomorphisms[0]
    source_positive = -stock["negative_gram"]
    if any(
        target[left, right] != source_positive[mapping[left], mapping[right]]
        for left in range(target.rows)
        for right in range(target.cols)
    ):
        raise ArithmeticError("weighted-Gram mapping is not entrywise exact")

    target_basis_indices = tuple(
        int(value)
        for value in orientation["root_generated_lattice"]["independent_root_indices"]
    )
    source_basis_indices = tuple(stock["basis_indices"])
    target_coordinates = _wall_coordinates(target, target_basis_indices)
    source_coordinates = _wall_coordinates(source_positive, source_basis_indices)
    target_basis_gram = target.extract(target_basis_indices, target_basis_indices)
    source_basis_gram = source_positive.extract(source_basis_indices, source_basis_indices)
    basis_change = source_coordinates[[mapping[index] for index in target_basis_indices], :]
    if not (
        abs(int(basis_change.det())) == 1
        and basis_change * source_basis_gram * basis_change.T == target_basis_gram
        and target_coordinates * basis_change == source_coordinates[list(mapping), :]
    ):
        raise ArithmeticError("exact graph match did not extend to a lattice isometry")

    orbits = _unfolded_orbits(row["facet_count"], set(orientation["long_vertices"]))
    orbit_sum = _orbit_sum_matrix(orbits, target.rows)
    folded = sp.Matrix(orientation["folded_integral_gram"])
    if orbit_sum * target * orbit_sum.T != folded:
        raise ArithmeticError("target orbit sums stopped recovering the folded Gram")
    source_orbit_sum = sp.zeros(len(orbits), source_positive.rows)
    for orbit_index, orbit in enumerate(orbits):
        for root in orbit:
            source_orbit_sum[orbit_index, mapping[root]] = 1
    if source_orbit_sum * source_positive * source_orbit_sum.T != folded:
        raise ArithmeticError("source orbit sums do not return the original folded Gram")

    target_permutation = _triality_permutation(orbits, target.rows)
    source_permutation = list(range(source_positive.rows))
    for target_root, source_root in enumerate(mapping):
        source_permutation[source_root] = mapping[target_permutation[target_root]]
    source_permutation = tuple(source_permutation)
    if any(
        source_positive[left, right]
        != source_positive[source_permutation[left], source_permutation[right]]
        for left in range(source_positive.rows)
        for right in range(source_positive.cols)
    ):
        raise ArithmeticError("transported triality is not a source wall symmetry")

    source_wall_coordinates_negative = _wall_coordinates(
        stock["negative_gram"], source_basis_indices
    )
    action = source_wall_coordinates_negative[
        [source_permutation[index] for index in source_basis_indices], :
    ]
    negative_basis_gram = stock["negative_gram"].extract(
        source_basis_indices, source_basis_indices
    )
    fixed_rank = len((action - sp.eye(action.rows)).nullspace())
    if not (
        action**3 == sp.eye(action.rows)
        and action * negative_basis_gram * action.T == negative_basis_gram
        and fixed_rank == row["dimension"] + 1
    ):
        raise ArithmeticError("transported triality failed its integral lattice checks")

    discriminant_difference = negative_basis_gram.inv() * (
        action - sp.eye(action.rows)
    )
    discriminant_action_trivial = all(value.q == 1 for value in discriminant_difference)

    ample = stock["source_ample"]
    if ample is None:
        raise ArithmeticError("an exact promotion lacks a frozen source ample class")
    ample_row = sp.Matrix(1, len(ample), ample)
    source_pairings = ample_row * negative_basis_gram * source_wall_coordinates_negative.T
    if not (
        (ample_row * negative_basis_gram * ample_row.T)[0] > 0
        and min(source_pairings) > 0
    ):
        raise ArithmeticError("frozen source ample class is not strict")
    invariant_ample = ample_row + ample_row * action + ample_row * action**2
    invariant_pairings = (
        invariant_ample * negative_basis_gram * source_wall_coordinates_negative.T
    )
    invariant_ample_square = int(
        (invariant_ample * negative_basis_gram * invariant_ample.T)[0]
    )
    if not (
        invariant_ample * action == invariant_ample
        and invariant_ample_square > 0
        and min(invariant_pairings) > 0
    ):
        raise ArithmeticError("triality average did not give a strict invariant ample class")

    source_orbits: list[list[int]] = []
    seen: set[int] = set()
    for start in range(source_positive.rows):
        if start in seen:
            continue
        orbit = []
        current = start
        while current not in orbit:
            orbit.append(current)
            seen.add(current)
            current = source_permutation[current]
        source_orbits.append(sorted(orbit))
    source_orbits.sort()

    return {
        "stock_id": stock["stock_id"],
        "stock_lattice": stock["lattice"],
        "exact_weighted_gram_isomorphism_count": len(isomorphisms),
        "canonical_target_to_stock_wall_permutation_one_based": [
            value + 1 for value in mapping
        ],
        "target_unfolded_orbits_one_based": [
            [value + 1 for value in orbit] for orbit in orbits
        ],
        "stock_triality_orbits_one_based": [
            [value + 1 for value in orbit] for orbit in source_orbits
        ],
        "wall_lattice_basis_change_determinant": int(basis_change.det()),
        "wall_lattice_basis_change_target_to_stock": _integer_rows(basis_change),
        "integral_wall_lattice_isometry": True,
        "triality_action_on_stock_basis": _integer_rows(action),
        "triality_order": 3,
        "triality_fixed_rank": fixed_rank,
        "triality_preserves_full_nef_wall_set": True,
        "triality_action_on_discriminant_is_trivial": discriminant_action_trivial,
        "invariant_ample_coordinates_in_stock_basis": [
            int(value) for value in invariant_ample
        ],
        "invariant_ample_square": invariant_ample_square,
        "invariant_ample_minimum_wall_pairing": min(
            int(value) for value in invariant_pairings
        ),
        "coinvariant_minus_two_rootless_by_invariant_ample": True,
        "complete_full_nef_source": {
            "archive": stock["archive"],
            "wall_gram_sha256": stock["archive_wall_gram_sha256"],
            "reason": stock["full_nef_reason"],
            "literal_vinberg_decision": stock["literal_vinberg"]["decision"],
        },
        "same_P_return": {
            "source_orbit_sum_gram_equals_original_folded_gram": True,
            "fixed_nef_section_has_exactly_original_facet_orbits": True,
            "hidden_nef_wall_count": 0,
            "level": "complete K3 nef cone plus divisor-lattice triality fixed section",
        },
        "geometric_order_three_K3_automorphism": {
            "certified": False,
            "reason": (
                "the transported nef-preserving NS isometry acts nontrivially on the "
                "discriminant; a matching transcendental Hodge action and Torelli lift "
                "have not been constructed"
            ),
        },
    }


def _target_rows(bredon: dict[str, Any]) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    targets = []
    for row in bredon["rows"]:
        if row["id"] in EXCLUDED_ROWS:
            continue
        for orientation in row["orientations"]:
            if orientation["nikulin_stable_range_primitive_k3_embedding"]:
                targets.append((row, orientation))
    if len(targets) != 30:
        raise ArithmeticError(f"expected 30 remaining stable packages, got {len(targets)}")
    return targets


def _compare(
    row: dict[str, Any], orientation: dict[str, Any], stocks: list[dict[str, Any]]
) -> dict[str, Any]:
    target = sp.Matrix(orientation["unfolded_gram"])
    target_smith = _nonzero_smith(target)
    if target_smith != orientation["root_generated_lattice"]["smith_diagonal"]:
        raise ArithmeticError("target Smith data stopped replaying")
    header = [
        stock
        for stock in stocks
        if stock["rank"] == target.rank() and stock["root_count"] == target.rows
    ]
    if any(stock["negative_gram"] is None or stock["smith"] is None for stock in header):
        raise ArithmeticError("a matching local header lacks its exact source Gram")
    lattice = [stock for stock in header if stock["smith"] == target_smith]
    exact: list[tuple[dict[str, Any], list[tuple[int, ...]]]] = []
    for stock in lattice:
        mappings = _isomorphisms(target, -stock["negative_gram"])
        if mappings:
            exact.append((stock, mappings))
    promotions = [
        _promotion(row, orientation, target, stock, mappings)
        for stock, mappings in exact
    ]
    if promotions:
        status = "exact_complete_stock_match_and_same_P_lattice_nef_return"
    elif header:
        status = "matching_headers_all_rejected_by_lattice_or_exact_gram"
    else:
        status = "no_matching_stock_in_local_complete_denominator"
    return {
        "package": f"{row['id']}::orientation-{orientation['orientation']}",
        "bredon_row": row["id"],
        "orientation": orientation["orientation"],
        "folded_dimension": row["dimension"],
        "folded_facet_count": row["facet_count"],
        "unfolded_rank": target.rank(),
        "unfolded_root_count": target.rows,
        "unfolded_smith": target_smith,
        "unfolded_graph_invariants": _graph_invariants(target),
        "matching_header_stock_ids": [stock["stock_id"] for stock in header],
        "matching_lattice_stock_ids": [stock["stock_id"] for stock in lattice],
        "exact_weighted_gram_stock_ids": [stock["stock_id"] for stock, _ in exact],
        "promotions": promotions,
        "status": status,
        "gap": (
            None
            if promotions
            else (
                "all locally frozen stocks with this rank/root-count header fail exact invariants"
                if header
                else "no complete K3 nef-root stock with this header is frozen locally"
            )
        ),
    }


def analyze() -> dict[str, Any]:
    bredon = _load(BREDON)
    stocks = _stock_denominator()
    comparisons = [
        _compare(row, orientation, stocks)
        for row, orientation in _target_rows(bredon)
    ]
    status_counts = Counter(row["status"] for row in comparisons)
    promoted = [
        (row["package"], promotion["stock_id"])
        for row in comparisons
        for promotion in row["promotions"]
    ]
    expected_promoted = [
        ("ADEG-N2-01-path-66::orientation-1", "R5-P4D4-NEF-H4-5"),
        ("ADEG-N3-01-path-633::orientation-0", "R6-UD4-NEF-H5-6"),
        ("ADEG-N3-04-star-6::orientation-0", "R6-U2D4-NEF-H5-6"),
        ("ADEG-N3-05-triangle-leaf-6::orientation-1", "R8-UA2x3-NEF-H7-10"),
        ("ADEG-N5-01-path-633336::orientation-0", "R10-UD4D4-NEF-H9-11"),
        ("ADEG-N5-02-path-triangle::orientation-0", "R8-UA2D4-NEF-H7-9"),
        ("ADEG-N6-01-diamond-tail::orientation-1", "R9-UD4A3-NEF-H8-10"),
        ("ADEG-N7-02-path-star::orientation-0", "R10-UD4D4-NEF-H9-11"),
    ]
    if sorted(promoted) != sorted(expected_promoted):
        raise ArithmeticError(f"promotion set changed: {promoted}")
    expected_status = {
        "no_matching_stock_in_local_complete_denominator": 7,
        "matching_headers_all_rejected_by_lattice_or_exact_gram": 15,
        "exact_complete_stock_match_and_same_P_lattice_nef_return": 8,
    }
    if dict(status_counts) != expected_status:
        raise ArithmeticError(f"comparison status census changed: {status_counts}")

    payload = {
        "schema": SCHEMA,
        "scope": {
            "polytope_first_targets": (
                "all 30 stable-range ADEG24 triality packages left after excluding "
                "the already closed N11 row and separately audited P_star"
            ),
            "stock_denominator": (
                "17 target-header-selected Roulleau complete stocks in ranks 5--8 "
                "plus all 38 complete K3 nef-root stocks frozen locally in ranks "
                "9--19; not a worldwide or all-rank K3 stock census"
            ),
            "return_level": (
                "exact complete-nef divisor-lattice fixed section; no geometric "
                "order-three K3 automorphism or quotient-surface claim"
            ),
        },
        "source_archives": {
            path.name: _file_sha256(path)
            for path in (
                BREDON,
                LOWER,
                RANK9,
                RANK10_SIMPLEX,
                RANK10_REMAINING,
                HIGH,
                R12_UA2E8,
                R12_GRAPH,
            )
        },
        "excluded_rows": EXCLUDED_ROWS,
        "stock_denominator": [_public_stock(stock) for stock in stocks],
        "comparisons": comparisons,
        "summary": {
            "remaining_stable_packages": len(comparisons),
            "local_complete_k3_nef_stocks": len(stocks),
            "local_stock_rank_distribution": {
                str(rank): count
                for rank, count in sorted(Counter(stock["rank"] for stock in stocks).items())
            },
            "packages_without_matching_local_header": status_counts[
                "no_matching_stock_in_local_complete_denominator"
            ],
            "packages_rejected_after_header_match": status_counts[
                "matching_headers_all_rejected_by_lattice_or_exact_gram"
            ],
            "exact_complete_stock_promotions": status_counts[
                "exact_complete_stock_match_and_same_P_lattice_nef_return"
            ],
            "distinct_promoted_stocks": len({stock for _, stock in promoted}),
            "promoted_package_stock_pairs": [list(pair) for pair in promoted],
            "geometric_order_three_K3_automorphism_promotions": 0,
            "new_abstract_polytopes": 0,
        },
        "status": (
            "eight exact polytope-first complete-nef lattice-carrier returns; "
            "geometric triality lifts and all-rank stock coverage remain open"
        ),
    }
    payload["canonical_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    arguments = parser.parse_args()
    payload = analyze()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(encoded, encoding="utf-8")
    if arguments.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != encoded:
            raise SystemExit("archive does not match exact replay")
    if not arguments.quiet:
        print(encoded, end="")


if __name__ == "__main__":
    main()
