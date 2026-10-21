#!/usr/bin/env python3
"""Identify the three non-C007 D0004 Coxeter survivors exactly.

The complete archived-stabilizer screen for Shimada's discriminant-16 atlas
leaves four Coxeter cases on chamber D0004.  C007 is the 27-wall H10 chamber.
This certificate compares the other three normalized wall Gram matrices with
the already certified D5 chambers.  Equality of the signed normalized Gram
matrix after a facet permutation proves equality of the hyperbolic chambers
up to isometry; facet count alone is never used as an identification.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import networkx as nx


ROOT = Path(__file__).resolve().parents[1]
CASE_DIRECTORY = (
    ROOT / "data" / "enumerations" / "shimada_disc16_atlas_archived_high_sections"
)
D5_SCREEN = (
    ROOT / "data" / "enumerations" / "shimada_disc16_d5_full_subgroup_fishing.json"
)
H10_CONE = (
    ROOT / "data" / "enumerations" / "shimada_disc16_d5_h10_30_exact_cone_certificate.json"
)
H10_VINBERG = (
    ROOT / "data" / "enumerations" / "shimada_disc16_d5_h10_30_vinberg_face_certificate.json"
)
H14_CONE = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h14_independent_exact_cone_certificate.json"
)
H14_ALEXEEV = (
    ROOT / "data" / "enumerations" / "shimada_disc16_d5_h14_alexeev_exact_match.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d4_duplicate_survivor_identifications.json"
)

EXPECTED_SHA256 = {
    "D0004_C001.json": "1b65b47f7b3b13b9fee830515056cfd60bdb39f0342f14e1d19a9dfa88adaf1f",
    "D0004_C003.json": "b39a9f071c92f546f79f37a61ba9bfb09abd513129436bad7e8a737921934c01",
    "D0004_C004.json": "8a32c3ebd06faff25225e685b393a6bce310532294834936e72446f82a9430ba",
    "shimada_disc16_d5_full_subgroup_fishing.json": "58ac74b11b59e6216c0aac2a34b8ce91c91ef3d891c41ea805ed3c3e9155b566",
    "shimada_disc16_d5_h10_30_exact_cone_certificate.json": "1160d8e5e4425a20a8a12bc9df6c6dd7fbd5031c5ebe8965b8d4357c82dad083",
    "shimada_disc16_d5_h10_30_vinberg_face_certificate.json": "4e0fc3d33da785537527df982d422312c0596ca8e93ced73f812934133255f8b",
    "shimada_disc16_d5_h14_independent_exact_cone_certificate.json": "c6574488aa1eb78da4cdfcbd2e9037c9099d6789ba3fb33724feccd9a256bae5",
    "shimada_disc16_d5_h14_alexeev_exact_match.json": "b8cb43fecc81c9036f72d5d8325a43f1085a2561e6eefce7f343f88d473b12ef",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_pinned(path: Path) -> dict[str, Any]:
    digest = _sha256(path)
    expected = EXPECTED_SHA256[path.name]
    if digest != expected:
        raise ArithmeticError(f"source changed: {path}: {digest} != {expected}")
    return json.loads(path.read_text(encoding="utf-8"))


def _fraction_matrix(payload: Sequence[Sequence[str | int]]) -> list[list[Fraction]]:
    matrix = [[Fraction(entry) for entry in row] for row in payload]
    size = len(matrix)
    if not size or any(len(row) != size for row in matrix):
        raise ArithmeticError("wall Gram matrix is not square")
    if any(matrix[i][j] != matrix[j][i] for i in range(size) for j in range(size)):
        raise ArithmeticError("wall Gram matrix is not symmetric")
    if any(matrix[i][i] <= 0 for i in range(size)):
        raise ArithmeticError("wall Gram matrix has a nonpositive diagonal")
    if any(matrix[i][j] > 0 for i in range(size) for j in range(i + 1, size)):
        raise ArithmeticError("wall Gram matrix is not in the acute sign convention")
    return matrix


def _label(matrix: Sequence[Sequence[Fraction]], left: int, right: int) -> str:
    if left == right:
        return "1"
    product = matrix[left][right]
    if product == 0:
        return "0"
    ratio = product * product / (matrix[left][left] * matrix[right][right])
    # All nonzero off-diagonal entries are negative.  Retaining this explicit
    # sign makes the certificate valid independently of that global check.
    return f"-sqrt({ratio})"


def _weighted_graph(matrix: Sequence[Sequence[Fraction]]) -> nx.Graph:
    size = len(matrix)
    graph = nx.Graph()
    for index in range(size):
        signature = Counter(
            _label(matrix, index, other) for other in range(size) if other != index
        )
        graph.add_node(index, signature=tuple(sorted(signature.items())))
    for left in range(size):
        for right in range(left + 1, size):
            graph.add_edge(left, right, normalized_product=_label(matrix, left, right))
    return graph


def _match(
    source_payload: Sequence[Sequence[str | int]],
    target_payload: Sequence[Sequence[str | int]],
) -> dict[str, Any]:
    source = _fraction_matrix(source_payload)
    target = _fraction_matrix(target_payload)
    if len(source) != len(target):
        raise ArithmeticError("facet counts differ")
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _weighted_graph(source),
        _weighted_graph(target),
        node_match=nx.algorithms.isomorphism.categorical_node_match("signature", None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "normalized_product", None
        ),
    )
    if not matcher.is_isomorphic():
        raise ArithmeticError("signed normalized Gram matrices are not isomorphic")
    mapping = {int(left): int(right) for left, right in matcher.mapping.items()}
    size = len(source)
    if set(mapping) != set(range(size)) or set(mapping.values()) != set(range(size)):
        raise ArithmeticError("graph isomorphism did not return a facet permutation")
    if any(
        _label(source, left, right) != _label(target, mapping[left], mapping[right])
        for left in range(size)
        for right in range(size)
    ):
        raise ArithmeticError("facet permutation failed entrywise replay")
    label_histogram = Counter(
        _label(source, left, right)
        for left in range(size)
        for right in range(left + 1, size)
    )
    return {
        "facet_count": size,
        "source_to_target_facet_permutation_one_based": [
            mapping[index] + 1 for index in range(size)
        ],
        "signed_normalized_off_diagonal_histogram": dict(
            sorted(label_histogram.items())
        ),
        "complete_edge_coloured_graph_isomorphic": True,
        "all_normalized_gram_entries_equal_after_permutation": True,
        "same_hyperbolic_chamber_up_to_isometry": True,
    }


def build_certificate() -> dict[str, Any]:
    cases = {
        index: _load_pinned(CASE_DIRECTORY / f"D0004_C{index:03d}.json")
        for index in (1, 3, 4)
    }
    d5 = _load_pinned(D5_SCREEN)
    h10_cone = _load_pinned(H10_CONE)
    h10_vinberg = _load_pinned(H10_VINBERG)
    h14_cone = _load_pinned(H14_CONE)
    h14_alexeev = _load_pinned(H14_ALEXEEV)

    targets = {
        1: ("D5-H14-60", d5["above_H9"]["complete_section_screens"]["class_02_H14"]),
        3: ("D5-H10-30", d5["above_H9"]["complete_section_screens"]["class_03_H10"]),
        4: ("D5-H10-30", d5["above_H9"]["complete_section_screens"]["class_03_H10"]),
    }
    identifications: list[dict[str, Any]] = []
    for index in (1, 3, 4):
        case = cases[index]
        target_name, target = targets[index]
        comparison = _match(
            case["coxeter_gate"]["positive_wall_gram"], target["positive_wall_gram"]
        )
        dimension = int(case["case"]["hyperbolic_dimension"])
        if dimension != int(target["hyperbolic_dimension"]):
            raise ArithmeticError("matched chambers have different dimensions")
        identifications.append(
            {
                "source": f"D0004/C{index:03d}",
                "target": target_name,
                "hyperbolic_dimension": dimension,
                "comparison": comparison,
                "finite_volume_noncompact_inherited_from_exact_target_isometry": True,
            }
        )

    checks = {
        "all_three_source_cases_pass_complete_exact_coxeter_gate": all(
            case["coxeter_gate"]["all_pairs_coxeter"] for case in cases.values()
        ),
        "h10_target_has_both_exact_volume_certificates": (
            h10_cone["decision"]["classification"] == "finite_volume_noncompact"
            and h10_cone["exact_double_description"]["hyperideal_vertex_count"] == 0
            and h10_vinberg["decision"]["classification"] == "finite_volume_noncompact"
            and h10_vinberg["literal_endpoint_audit"]["bad_edge_count"] == 0
            and h10_vinberg["literal_endpoint_audit"]["bad_parabolic_neighborhood_count"] == 0
        ),
        "h14_target_has_exact_volume_certificate": (
            h14_cone["decision"]["classification"] == "finite_volume_noncompact"
            and h14_cone["exact_double_description"]["hyperideal_vertex_count"] == 0
        ),
        "h14_target_is_exactly_alexeev_15_7_1": h14_alexeev["decision"]
        ["same_hyperbolic_chamber_up_to_isometry"],
        "all_three_normalized_gram_permutations_replay": all(
            record["comparison"]["all_normalized_gram_entries_equal_after_permutation"]
            for record in identifications
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"duplicate-survivor checks failed: {checks}")

    return {
        "schema": "coxeter.shimada_disc16_d4_duplicate_survivors.v1",
        "sources": {
            path.name: {
                "path": str(path.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256[path.name],
            }
            for path in (
                *(CASE_DIRECTORY / f"D0004_C{index:03d}.json" for index in (1, 3, 4)),
                D5_SCREEN,
                H10_CONE,
                H10_VINBERG,
                H14_CONE,
                H14_ALEXEEV,
            )
        },
        "identifications": identifications,
        "target_volume_facts": {
            "D5-H10-30": {
                "finite_vertices": h10_cone["exact_double_description"]["finite_vertex_count"],
                "cusps": h10_cone["exact_double_description"]["ideal_vertex_count"],
                "hyperideal_rays": h10_cone["exact_double_description"]["hyperideal_vertex_count"],
                "vinberg_edges": h10_vinberg["literal_endpoint_audit"]["edge_count"],
            },
            "D5-H14-60": {
                "finite_vertices": h14_cone["exact_double_description"]["finite_vertex_count"],
                "cusps": h14_cone["exact_double_description"]["ideal_vertex_count"],
                "hyperideal_rays": h14_cone["exact_double_description"]["hyperideal_vertex_count"],
            },
        },
        "decision": {
            "D0004_C001_is_D5_H14_60": True,
            "D0004_C003_is_D5_H10_30": True,
            "D0004_C004_is_D5_H10_30": True,
            "three_cast_instances_add_zero_abstract_chamber_types": True,
            "D0004_C007_is_the_only_new_literal_type_among_the_four_D0004_survivors": True,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)


if __name__ == "__main__":
    main()
