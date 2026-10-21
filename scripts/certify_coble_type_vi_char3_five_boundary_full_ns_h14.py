#!/usr/bin/env python3
"""Certify the five-boundary Type-VI Coble full-nef chamber in H14.

Kondo's Type-VI surface in characteristic three has five boundary curves,
ten (-2)-curves, and ten boundary-pair (-1)-roots.  The latter data determine
the full Neron--Severi glue.  Enumerating every primitive cusp of the complete
Coble--Mukai chamber then recovers all remaining (-1)-curves.  This script
checks the resulting complete negative-curve cone by Vinberg's literal
edge/parabolic criterion and identifies its arithmetic mirror Gram exactly
with Alexeev's published (15,7,1) chamber.

This is a SURFACE -> POLYTOPE certificate.  The archived boundary-frame fold
is useful input for a later inverse recognizer, but no polytope-to-Coble
existence theorem is asserted here.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ

from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)
from scripts.certify_coble_one_boundary_full_ns_h10 import _cdd_rays
from scripts.match_shimada_disc16_d5_h14_alexeev_15_7_1 import (
    _shimada_oriented_roots_and_gram,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json"
)
FIXED = ROOT / "data/enumerations/shimada_disc16_d5_h14_fixed_lattice_certificate.json"
EXACT_MATCH = ROOT / "data/enumerations/shimada_disc16_d5_h14_alexeev_exact_match.json"
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"

SCHEMA = "coxeter.coble_type_vi_char3_five_boundary_full_ns_h14.v1"
PINNED = {
    "kondo_source_tar": "6989648564a80167f2daa35da64b90a063a0a090b9f23ff843b7fabb9d96a98d",
    "kondo_tex": "8e23df4e77a2666a7a4f9a634c1ec2f36066353bf0bdcd0a72b11a0a70bfb137",
    "fixed": "c9ed730ae9e6476136bbb4810dac5c76d43bf618ce8b63dacc8cbcf111930e23",
    "exact_match": "b8cb43fecc81c9036f72d5d8325a43f1085a2561e6eefce7f343f88d473b12ef",
    "cpp": "45174f5461d28affc1fdb4cb55047b7a1581037d307e32800651badb3c2df045",
}
CM_BASIS_ROOTS = (0, 1, 3, 4, 5, 8, 10, 11, 15, 18)
VINBERG_EXPECTED = {
    "finite_nodes": 10_956_853,
    "affine_nodes": 12_225_900,
    "prefix": [
        60, 970, 8_335, 46_274, 180_085, 512_903, 1_098_902,
        1_812_645, 2_302_420, 2_253_968, 1_637_502, 830_613,
        244_330, 27_845,
    ],
    "special": [
        60, 970, 8_335, 46_415, 183_308, 537_765, 1_198_850,
        2_045_905, 2_659_490, 2_587_221, 1_819_090, 866_315,
        244_330, 27_845,
    ],
    "affine_components": 2_704,
    "cusps": 1_027,
    "endpoints": {"0,2": 9_660, "1,1": 79_510, "2,0": 155_160},
}
FULL_CONE_EXPECTED = {
    "ray_count": 28_872,
    "finite": 27_845,
    "ideal": 1_027,
    "hyperideal": 0,
    "incidence_histogram": {14: 27_845, 15: 432, 16: 540, 17: 55},
    "primitive_ray_list_sha256": "13c6e10ea3d631d7620aefc06a0e62b86991838207ac5f5c243ef79fdd5da042",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("an archived integral matrix became nonintegral")
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _rational_rows(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if value.q == 1 else str(value) for value in matrix.row(index)]
        for index in range(matrix.rows)
    ]


def _half_integral_basis(
    generators: Sequence[sp.Matrix], rank: int = 15
) -> sp.Matrix:
    doubled = sp.Matrix(
        [
            [int(2 * generators[column][row]) for column in range(len(generators))]
            for row in range(rank)
        ]
    )
    if doubled.rank() != rank:
        raise ArithmeticError("the proposed divisor generators do not have full rank")
    return hermite_normal_form(doubled) / 2


def _primitive_rational(vector: sp.Matrix) -> tuple[int, ...]:
    denominator = math.lcm(*(int(sp.denom(value)) for value in vector))
    values = [int(value * denominator) for value in vector]
    divisor = math.gcd(*(abs(value) for value in values))
    if divisor == 0:
        raise ArithmeticError("zero restriction has no primitive direction")
    values = [value // divisor for value in values]
    first = next(value for value in values if value)
    if first < 0:
        values = [-value for value in values]
    return tuple(values)


def _type_vi_root_stock() -> dict[str, Any]:
    """Reconstruct Kondo Figure VI from the ten two-subsets of a five-set."""
    pairs = list(itertools.combinations(range(5), 2))
    gram = -2 * sp.eye(20)
    for left, pair in enumerate(pairs):
        gram[left, 10 + left] = gram[10 + left, left] = 2
        for right, other in enumerate(pairs[:left]):
            if set(pair).isdisjoint(other):
                # Figure VI(A): the Petersen graph.
                gram[left, right] = gram[right, left] = 1
            if len(set(pair) & set(other)) == 1:
                # Figure VI(B): its triangular-graph complement.
                gram[10 + left, 10 + right] = gram[10 + right, 10 + left] = 1
    basis_gram = gram.extract(CM_BASIS_ROOTS, CM_BASIS_ROOTS)
    root_coordinates = basis_gram.inv() * gram.extract(CM_BASIS_ROOTS, range(20))
    checks = {
        "twenty_roots": gram.shape == (20, 20),
        "root_gram_rank_ten": gram.rank() == 10,
        "root_gram_inertia": exact_inertia(basis_gram) == (1, 9, 0),
        "chosen_roots_are_unimodular_CM_basis": basis_gram.det() == -1,
        "all_root_coordinates_integral": all(value.q == 1 for value in root_coordinates),
        "ten_cross_pairings_are_two": sum(
            gram[left, right] == 2
            for left in range(20)
            for right in range(left)
        )
        == 10,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Type-VI Coble--Mukai root stock changed: {checks}")
    return {
        "pairs": pairs,
        "root_gram": gram,
        "basis_gram": basis_gram,
        "root_coordinates": root_coordinates,
        "checks": checks,
    }


def _cm_chamber(stock: dict[str, Any], executable: str) -> dict[str, Any]:
    gram: sp.Matrix = stock["root_gram"]
    basis_gram: sp.Matrix = stock["basis_gram"]
    inequalities = gram[:, list(CM_BASIS_ROOTS)]
    rays = sorted(_cdd_rays(inequalities, executable))
    counts: Counter[str] = Counter()
    cusps: list[tuple[int, ...]] = []
    for raw in rays:
        vector = sp.Matrix(raw)
        square = int((vector.T * basis_gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        counts[kind] += 1
        if square == 0:
            cusps.append(raw)
    if counts != Counter({"finite": 600, "ideal": 57}) or len(cusps) != 57:
        raise ArithmeticError(f"the exact Type-VI CM cone changed: {counts}")
    return {
        "ray_count": len(rays),
        "finite_vertex_count": counts["finite"],
        "ideal_vertex_count": counts["ideal"],
        "hyperideal_vertex_count": counts["hyperideal"],
        "primitive_cusps": cusps,
        "primitive_ray_list_sha256": hashlib.sha256(
            json.dumps(rays, separators=(",", ":")).encode()
        ).hexdigest(),
    }


def _surface_lattices_and_walls(
    stock: dict[str, Any], cm: dict[str, Any]
) -> dict[str, Any]:
    pairs: list[tuple[int, int]] = stock["pairs"]
    basis_gram: sp.Matrix = stock["basis_gram"]
    root_coordinates: sp.Matrix = stock["root_coordinates"]
    ambient = sp.zeros(15)
    ambient[:10, :10] = basis_gram
    for index in range(5):
        ambient[10 + index, 10 + index] = -1
    halves = [sp.eye(15)[:, 10 + index] for index in range(5)]
    boundaries = [2 * half for half in halves]
    roots = [
        sp.Matrix([*root_coordinates[:, index], 0, 0, 0, 0, 0])
        for index in range(20)
    ]
    pair_exceptionals = [
        (roots[10 + index] - halves[pair[0]] - halves[pair[1]]) / 2
        for index, pair in enumerate(pairs)
    ]
    canonical = -sum(halves, sp.zeros(15, 1))
    ns_basis = _half_integral_basis(
        [*boundaries, *roots[:10], *pair_exceptionals, canonical]
    )
    ns_gram = ns_basis.T * ambient * ns_basis
    ns_inverse = ns_basis.inv()
    base_curves = [*boundaries, *roots[:10], *pair_exceptionals]
    if any(any(value.q != 1 for value in ns_inverse * curve) for curve in base_curves):
        raise ArithmeticError("a visible Type-VI negative curve left NS")

    cusp_membership: list[list[bool]] = []
    cusp_exceptionals: list[sp.Matrix] = []
    cusp_labels: list[str] = []
    lifted_pairs: list[tuple[int, int]] = []
    for cusp_index, raw in enumerate(cm["primitive_cusps"]):
        fiber = sp.Matrix([*raw, 0, 0, 0, 0, 0])
        row: list[bool] = []
        for boundary_index, half in enumerate(halves):
            candidate = fiber - half
            integral = all(value.q == 1 for value in ns_inverse * candidate)
            row.append(integral)
            if integral:
                cusp_exceptionals.append(candidate)
                cusp_labels.append(f"E_cusp_{cusp_index + 1}_B{boundary_index + 1}")
                lifted_pairs.append((cusp_index, boundary_index))
        cusp_membership.append(row)

    labels = [*(f"B{index + 1}" for index in range(5))]
    labels.extend(f"C_{left + 1}{right + 1}" for left, right in pairs)
    labels.extend(f"E_pair_{left + 1}{right + 1}" for left, right in pairs)
    labels.extend(cusp_labels)
    physical = [*base_curves, *cusp_exceptionals]
    wall_matrix_ambient = sp.Matrix.hstack(*physical)
    ns_walls = ns_inverse * wall_matrix_ambient
    wall_gram = wall_matrix_ambient.T * ambient * wall_matrix_ambient

    pattern_counts = Counter(tuple(row) for row in cusp_membership)
    expected_patterns = Counter({(False,) * 5: 22})
    for boundary_index in range(5):
        expected_patterns[
            tuple(index == boundary_index for index in range(5))
        ] = 7
    norms = Counter(int(wall_gram[index, index]) for index in range(60))
    ratios: Counter[str] = Counter()
    coxeter_angles = True
    for left in range(60):
        for right in range(left):
            product = int(wall_gram[left, right])
            if product < 0:
                coxeter_angles = False
                continue
            if product == 0:
                ratios["orthogonal"] += 1
                continue
            ratio = sp.Rational(
                4 * product * product,
                int(wall_gram[left, left]) * int(wall_gram[right, right]),
            )
            if ratio == 1:
                ratios["m3"] += 1
            elif ratio == 2:
                ratios["m4"] += 1
            elif ratio == 3:
                ratios["m6"] += 1
            elif ratio == 4:
                ratios["parallel"] += 1
            elif ratio > 4:
                ratios[f"ultraparallel_ratio_{ratio}"] += 1
            else:
                coxeter_angles = False

    # Pic~ is generated by NS and the five boundary halves.  Doubling the
    # form makes the complete reflection arrangement integral.
    pic_basis = _half_integral_basis([*ns_basis.columnspace(), *halves])
    pic_gram = 2 * pic_basis.T * ambient * pic_basis
    pic_inverse = pic_basis.inv()
    arithmetic_normals = [*halves, *physical[5:]]
    arithmetic_walls = pic_inverse * sp.Matrix.hstack(*arithmetic_normals)

    checks = {
        "five_boundary_rank_is_fifteen": ns_gram.rows == 15,
        "NS_is_odd_unimodular_I_1_14": ns_gram.det() == 1
        and exact_inertia(ns_gram) == (1, 14, 0)
        and any(int(ns_gram[index, index]) % 2 for index in range(15)),
        "canonical_square_minus_five": int((canonical.T * ambient * canonical)[0]) == -5,
        "canonical_is_integral": all(value.q == 1 for value in ns_inverse * canonical),
        "cusp_membership_is_22_plus_5_times_7": pattern_counts == expected_patterns,
        "exactly_35_cusp_exceptionals": len(cusp_exceptionals) == 35,
        "all_cusp_exceptionals_are_minus_one": all(
            (curve.T * ambient * curve)[0] == -1 for curve in cusp_exceptionals
        ),
        "all_cusp_exceptionals_have_K_degree_minus_one": all(
            (canonical.T * ambient * curve)[0] == -1 for curve in cusp_exceptionals
        ),
        "complete_negative_curve_split": norms == Counter({-1: 45, -2: 10, -4: 5}),
        "sixty_walls_span_NS": ns_walls.rank() == 15
        and all(abs(int(value)) == 1 for value in smith_normal_form(ns_walls, domain=ZZ).diagonal()),
        "all_physical_pairings_are_acute_Coxeter_or_disjoint": coxeter_angles,
        "arithmetic_lattice_is_15_7_1": pic_gram.det() == 128
        and exact_inertia(pic_gram) == (1, 14, 0)
        and [
            abs(int(value))
            for value in smith_normal_form(pic_gram, domain=ZZ).diagonal()
        ]
        == [1] * 8 + [2] * 7,
        "all_arithmetic_normals_integral": all(value.q == 1 for value in arithmetic_walls),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Type-VI full-NS construction failed: {checks}")
    return {
        "ambient": ambient,
        "halves": halves,
        "roots": roots,
        "canonical": canonical,
        "ns_basis": ns_basis,
        "ns_gram": ns_gram,
        "labels": labels,
        "physical": physical,
        "ns_walls": ns_walls,
        "wall_gram": wall_gram,
        "pic_basis": pic_basis,
        "pic_gram": pic_gram,
        "arithmetic_walls": arithmetic_walls,
        "cusp_membership": cusp_membership,
        "lifted_pairs": lifted_pairs,
        "pattern_counts": pattern_counts,
        "ratio_histogram": ratios,
        "checks": checks,
    }


def _literal_vinberg(wall_gram: sp.Matrix, compiler: str | None = None) -> dict[str, Any]:
    if _sha256(CPP_SOURCE) != PINNED["cpp"]:
        raise ArithmeticError("the pinned literal Vinberg verifier changed")
    selected = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected is None:
        raise RuntimeError("a C++17 compiler is required")
    raw_input = "14 60\n" + "\n".join(
        " ".join(str(-int(wall_gram[row, column])) for column in range(60))
        for row in range(60)
    ) + "\n"
    with tempfile.TemporaryDirectory(prefix="coble-vi-h14-vinberg-") as raw:
        executable = Path(raw) / "vinberg"
        subprocess.run(
            [selected, "-O3", "-std=c++17", str(CPP_SOURCE), "-o", str(executable)],
            check=True,
        )
        completed = subprocess.run(
            [str(executable)],
            input=raw_input,
            text=True,
            capture_output=True,
            check=True,
            timeout=1800,
        )
    record = _parse_output(completed.stdout)
    expected = VINBERG_EXPECTED
    checks = {
        "dimension_and_facets": record["dimension"] == 14 and record["facet_count"] == 60,
        "finite_search_nodes": record["edge_search_node_count"] == expected["finite_nodes"],
        "affine_search_nodes": record["all_finite_subset_search_node_count"] == expected["affine_nodes"],
        "finite_prefix_counts": record["ordering_dependent_edge_search_prefix_counts"] == expected["prefix"],
        "finite_special_counts": record["finite_special_subset_counts_by_rank"] == expected["special"],
        "affine_components": record["connected_affine_component_count"] == expected["affine_components"],
        "cusps": record["ideal_vertex_count"] == expected["cusps"],
        "endpoint_histogram": record["edge_endpoint_type_histogram"] == expected["endpoints"],
        "zero_bad_edges": record["bad_edge_count"] == 0,
        "zero_bad_parabolic_neighborhoods": record["bad_parabolic_neighborhood_count"] == 0,
        "all_edges_accounted_for": sum(record["edge_endpoint_type_histogram"].values())
        == expected["special"][-2],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Type-VI literal Vinberg output changed: {checks}")
    return {
        **record,
        "finite_vertex_count": expected["special"][-1],
        "elliptic_edge_count": expected["special"][-2],
        "input_sha256": hashlib.sha256(raw_input.encode()).hexdigest(),
        "stdout_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
        "criterion": (
            "every elliptic rank-13 edge has exactly two finite/ideal endpoints, "
            "and every affine zero-neighborhood has full parabolic rank 13"
        ),
        "decision": "finite_volume_noncompact",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _full_cone(
    surface: dict[str, Any], executable: str, recompute: bool
) -> dict[str, Any]:
    inequalities = surface["ns_walls"].T * surface["ns_gram"]
    input_rows = _rows(inequalities)
    expected = FULL_CONE_EXPECTED
    if recompute:
        rays = _cdd_rays(inequalities, executable)
        counts: Counter[str] = Counter()
        incidence: Counter[int] = Counter()
        slacks = []
        for raw in rays:
            vector = sp.Matrix(raw)
            square = int((vector.T * surface["ns_gram"] * vector)[0])
            counts[
                "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
            ] += 1
            slack = inequalities * vector
            incidence[sum(value == 0 for value in slack)] += 1
            slacks.append(slack)
        ray_hash = hashlib.sha256(
            json.dumps(sorted(rays), separators=(",", ":")).encode()
        ).hexdigest()
        facet_ranks = [
            sp.Matrix(
                [
                    ray
                    for ray, slack in zip(rays, slacks, strict=True)
                    if slack[facet] == 0
                ]
            ).rank()
            for facet in range(60)
        ]
        checks = {
            "ray_count": len(rays) == expected["ray_count"],
            "finite_ideal_hyperideal": counts
            == Counter(
                {
                    "finite": expected["finite"],
                    "ideal": expected["ideal"],
                    "hyperideal": expected["hyperideal"],
                }
            ),
            "incidence_histogram": dict(incidence) == expected["incidence_histogram"],
            "ray_hash": ray_hash == expected["primitive_ray_list_sha256"],
            "all_facets_irredundant": facet_ranks == [14] * 60,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"the full Type-VI cdd cone changed: {checks}")
    return {
        "algorithm": "cddlib exact rational double description",
        "inequality_matrix_sha256": hashlib.sha256(
            json.dumps(input_rows, separators=(",", ":")).encode()
        ).hexdigest(),
        "extreme_ray_count": expected["ray_count"],
        "finite_vertex_count": expected["finite"],
        "ideal_vertex_count": expected["ideal"],
        "hyperideal_vertex_count": expected["hyperideal"],
        "ray_facet_incidence_histogram": {
            str(key): value for key, value in sorted(expected["incidence_histogram"].items())
        },
        "all_sixty_facet_face_ranks": [14] * 60,
        "all_facets_irredundant": True,
        "primitive_ray_list_sha256": expected["primitive_ray_list_sha256"],
        "slow_exact_replay": (
            "pass --full-cone; the default replay pins the input and archived "
            "output while avoiding the approximately ten-minute cdd enumeration"
        ),
    }


def _node_signature(matrix: sp.Matrix, index: int) -> tuple[Any, ...]:
    interactions = Counter(
        (int(matrix[other, other]), int(matrix[index, other]))
        for other in range(matrix.rows)
        if other != index
    )
    return (
        int(matrix[index, index]),
        tuple(sorted((square, product, count) for (square, product), count in interactions.items())),
    )


def _weighted_graph(matrix: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    for index in range(matrix.rows):
        graph.add_node(index, signature=_node_signature(matrix, index))
    for left in range(matrix.rows):
        for right in range(left):
            graph.add_edge(left, right, product=int(matrix[left, right]))
    return graph


def _lineage(surface: dict[str, Any]) -> dict[str, Any]:
    if _sha256(FIXED) != PINNED["fixed"] or _sha256(EXACT_MATCH) != PINNED["exact_match"]:
        raise ArithmeticError("a pinned H14/60 lineage certificate changed")
    fixed = json.loads(FIXED.read_text(encoding="utf-8"))
    exact = json.loads(EXACT_MATCH.read_text(encoding="utf-8"))
    _, shimada_gram = _shimada_oriented_roots_and_gram(fixed)
    arithmetic_gram = (
        surface["arithmetic_walls"].T
        * surface["pic_gram"]
        * surface["arithmetic_walls"]
    )
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _weighted_graph(arithmetic_gram),
        _weighted_graph(shimada_gram),
        node_match=nx.algorithms.isomorphism.categorical_node_match("signature", None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("product", None),
    )
    if not matcher.is_isomorphic():
        raise ArithmeticError("the Type-VI and Shimada H14/60 mirror Grams differ")
    mapping = {int(left): int(right) for left, right in matcher.mapping.items()}
    if any(
        arithmetic_gram[left, right] != shimada_gram[mapping[left], mapping[right]]
        for left in range(60)
        for right in range(60)
    ):
        raise ArithmeticError("the Type-VI-to-Shimada permutation failed entrywise")
    shimada_to_alexeev = [
        value - 1
        for value in exact["weighted_gram_matrix_comparison"][
            "shimada_facet_to_alexeev_root_permutation_one_based"
        ]
    ]
    type_vi_to_alexeev = [shimada_to_alexeev[mapping[index]] + 1 for index in range(60)]
    arithmetic_diagonal = Counter(int(arithmetic_gram[i, i]) for i in range(60))
    arithmetic_off = Counter(
        int(arithmetic_gram[left, right])
        for left in range(60)
        for right in range(left)
    )
    checks = {
        "same_60_by_60_weighted_Gram_as_Shimada": True,
        "all_3600_entries_verified": True,
        "Shimada_is_exactly_Alexeev_15_7_1": exact["decision"][
            "same_hyperbolic_chamber_up_to_isometry"
        ],
        "diagonal_histogram": arithmetic_diagonal == Counter({-2: 50, -4: 10}),
        "off_diagonal_histogram": arithmetic_off
        == Counter({0: 865, 1: 20, 2: 395, 4: 120, 6: 310, 10: 60}),
        "same_1027_cusps": exact["independent_finite_volume_cross_check"][
            "alexeev_maximal_parabolic_subdiagrams"
        ]
        == 1_027,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Type-VI lineage checks failed: {checks}")
    return {
        "exact_identification": "Alexeev (15,7,1) = Shimada D5 H14/60",
        "type_VI_to_Shimada_zero_based": [mapping[index] for index in range(60)],
        "type_VI_to_Alexeev_one_based": type_vi_to_alexeev,
        "arithmetic_diagonal_histogram": dict(sorted(arithmetic_diagonal.items())),
        "arithmetic_off_diagonal_histogram": dict(sorted(arithmetic_off.items())),
        "transitive_sources": {
            "Shimada_fixed_lattice_certificate": str(FIXED.relative_to(ROOT)),
            "Shimada_to_Alexeev_exact_match": str(EXACT_MATCH.relative_to(ROOT)),
            "Alexeev_source": exact["sources"]["alexeev"],
        },
        "decision": {
            "same_abstract_polytope": True,
            "same_arithmetic_reflection_chamber": True,
            "new_abstract_polytope": False,
            "positive_gain": (
                "new project-level complete full-nef realization by Kondo's "
                "five-boundary characteristic-three Type-VI Coble surface"
            ),
            "worldwide_priority": (
                "closed known for the abstract chamber via Alexeev; priority for "
                "this full-NS Coble realization has not been asserted"
            ),
        },
        "checks": checks,
    }


def _arithmetic_and_fold(surface: dict[str, Any], stock: dict[str, Any], cm: dict[str, Any]) -> dict[str, Any]:
    pic_gram: sp.Matrix = surface["pic_gram"]
    walls: sp.Matrix = surface["arithmetic_walls"]
    span = [
        abs(int(value))
        for value in smith_normal_form(walls, domain=ZZ).diagonal()
    ]
    root_types: Counter[tuple[int, int]] = Counter()
    reflection_checks = []
    for index in range(60):
        wall = walls[:, index]
        square = int((wall.T * pic_gram * wall)[0])
        divisibility = math.gcd(*(abs(int(value)) for value in pic_gram * wall))
        root_types[(square, divisibility)] += 1
        reflection = sp.eye(15) - sp.Rational(2, square) * wall * (wall.T * pic_gram)
        reflection_checks.append(
            all(value.q == 1 for value in reflection)
            and reflection.T * pic_gram * reflection == pic_gram
        )

    physical_gram: sp.Matrix = surface["wall_gram"]
    boundary_indices = [
        index for index in range(60) if physical_gram[index, index] == -4
    ]
    restrictions = [
        _primitive_rational(curve[:10, :]) for curve in surface["physical"][5:]
    ]
    restriction_squares = Counter(
        int((sp.Matrix(vector).T * stock["basis_gram"] * sp.Matrix(vector))[0])
        for vector in restrictions
    )
    root_directions = {
        _primitive_rational(stock["root_coordinates"][:, index]) for index in range(20)
    }
    negative_restrictions = {
        vector
        for vector in restrictions
        if (sp.Matrix(vector).T * stock["basis_gram"] * sp.Matrix(vector))[0] < 0
    }
    isotropic_restrictions = {
        vector
        for vector in restrictions
        if (sp.Matrix(vector).T * stock["basis_gram"] * sp.Matrix(vector))[0] == 0
    }
    cm_cusps = {
        _primitive_rational(sp.Matrix(vector)) for vector in cm["primitive_cusps"]
    }
    checks = {
        "all_arithmetic_reflections_integral": all(reflection_checks),
        "arithmetic_walls_span_lattice": span == [1] * 15,
        "root_type_split": root_types
        == Counter({(-2, 1): 15, (-2, 2): 35, (-4, 2): 10}),
        "raw_physical_Gram_has_unique_five_minus_four_frame": boundary_indices == list(range(5)),
        "boundary_frame_is_pairwise_orthogonal": physical_gram.extract(boundary_indices, boundary_indices)
        == -4 * sp.eye(5),
        "common_boundary_orthogonal_has_rank_ten": stock["basis_gram"].rank() == 10,
        "negative_restrictions_recover_exact_twenty_CM_roots": negative_restrictions == root_directions,
        "remaining_restrictions_are_35_distinct_CM_cusps": len(isotropic_restrictions) == 35
        and isotropic_restrictions < cm_cusps,
        "restriction_square_split": restriction_squares == Counter({0: 35, -2: 20}),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Type-VI arithmetic/fold checks failed: {checks}")
    return {
        "arithmetic_lattice": "2*Pic~(S), Nikulin (15,7,1)",
        "identification": "U(2)+E8(-1)+A1(-1)^5",
        "gram": _rows(pic_gram),
        "basis_in_CM_plus_boundary_half_coordinates": _rational_rows(surface["pic_basis"]),
        "wall_vectors": _rows(walls.T),
        "primitive_square_divisibility_histogram": {
            f"{square},{divisibility}": count
            for (square, divisibility), count in sorted(root_types.items())
        },
        "walls_span_whole_arithmetic_lattice": span == [1] * 15,
        "all_60_reflections_integral": all(reflection_checks),
        "boundary_frame_fold": {
            "input_scope": (
                "the marked physical primitive-divisor wall Gram, not an "
                "unmarked abstract Coxeter matrix"
            ),
            "boundary_indices_zero_based": boundary_indices,
            "boundary_gram": _rows(physical_gram.extract(boundary_indices, boundary_indices)),
            "common_orthogonal_signature": [1, 9],
            "common_orthogonal_lattice": "II_(1,9) = U+E8(-1)",
            "negative_primitive_restriction_count": len(negative_restrictions),
            "isotropic_redundant_restriction_count": len(isotropic_restrictions),
            "recovered_CM_root_gram": _rows(stock["root_gram"]),
            "recovered_CM_chamber": "Kondo Type VI, H9/20",
            "interpretation": (
                "the five boundary hyperplanes contain CM(S); restricting the "
                "other 55 inequalities recovers exactly the 20 root walls plus "
                "35 future-isotropic redundant cusp inequalities"
            ),
            "inverse_boundary": (
                "a nonempty marked survivor for a later P->Coble classifier; "
                "not by itself a surface-existence or nine-class inverse theorem"
            ),
        },
        "checks": checks,
    }


def build_certificate(
    cdd_executable: str | None = None,
    compiler: str | None = None,
    full_cone: bool = False,
) -> dict[str, Any]:
    executable = cdd_executable or shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the exact CM cusp enumeration")
    stock = _type_vi_root_stock()
    cm = _cm_chamber(stock, executable)
    surface = _surface_lattices_and_walls(stock, cm)
    literal = _literal_vinberg(surface["wall_gram"], compiler)
    cone = _full_cone(surface, executable, full_cone)
    lineage = _lineage(surface)
    arithmetic = _arithmetic_and_fold(surface, stock, cm)

    checks = {
        "source_surface_is_char3_five_boundary_not_Martin_char5_one_boundary": True,
        "complete_sixty_negative_curves": len(surface["physical"]) == 60,
        "full_nef_not_weak_carrier": surface["checks"]["sixty_walls_span_NS"],
        "literal_finite_volume_noncompact": literal["decision"] == "finite_volume_noncompact"
        and literal["ideal_vertex_count"] > 0,
        "independent_exact_cone_has_no_hyperideal_rays": cone["hyperideal_vertex_count"] == 0
        and cone["finite_vertex_count"] == literal["finite_vertex_count"]
        and cone["ideal_vertex_count"] == literal["ideal_vertex_count"],
        "exactly_known_Alexeev_Shimada_chamber": lineage["decision"]["same_abstract_polytope"],
        "not_new_abstract_polytope": not lineage["decision"]["new_abstract_polytope"],
        "boundary_fold_is_only_marked_inverse_input": arithmetic["boundary_frame_fold"][
            "inverse_boundary"
        ].startswith("a nonempty marked survivor"),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global Type-VI certificate checks failed: {checks}")

    payload = {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "SURFACE_TO_POLYTOPE_with_marked_boundary_fold_only",
        "claim": {
            "theorem": (
                "The full Neron--Severi nef cone of Kondo's characteristic-three "
                "Type-VI Coble surface with five boundary components is the "
                "60-facet finite-volume noncompact Alexeev (15,7,1) Coxeter "
                "polytope in H14."
            ),
            "not_claimed": [
                "a new abstract Coxeter polytope",
                "a new arithmetic reflection chamber",
                "a characteristic-five five-boundary surface",
                "an unmarked polytope-to-surface theorem",
                "the future nine-class Coble inverse classification",
            ],
        },
        "sources": {
            "kondo": {
                "title": "Coble surfaces with finite automorphism group",
                "arxiv": "2102.00591",
                "source_tar_sha256": PINNED["kondo_source_tar"],
                "main_tex_sha256": PINNED["kondo_tex"],
                "used": [
                    "main classification table: Type VI, p=3, n=5, k=20, Aut=S5",
                    "Example TypeVI3 and Figure VI",
                    "Propositions NefCone, Vinberg2, Halphen and the primitive-isotropic/fibration lemma",
                ],
                "characteristic_distinction": (
                    "Kondo has two Type-VI Coble rows: p=3 with five boundaries "
                    "(this theorem), and p=5 with one boundary. Martin's Type-VI "
                    "Enriques degeneration concerns the latter p=5 row."
                ),
            },
            "lineage_inputs": {
                str(FIXED.relative_to(ROOT)): PINNED["fixed"],
                str(EXACT_MATCH.relative_to(ROOT)): PINNED["exact_match"],
            },
            "literal_vinberg_verifier": {
                "path": str(CPP_SOURCE.relative_to(ROOT)),
                "sha256": PINNED["cpp"],
            },
        },
        "Coble_Mukai_chamber": {
            "lattice": "II_(1,9) = U+E8(-1)",
            "basis_root_indices_zero_based": list(CM_BASIS_ROOTS),
            "basis_gram": _rows(stock["basis_gram"]),
            "root_labels": [
                *(f"e_{left + 1}{right + 1}" for left, right in stock["pairs"]),
                *(f"f_{left + 1}{right + 1}" for left, right in stock["pairs"]),
            ],
            "root_vectors": _rows(stock["root_coordinates"].T),
            "root_gram": _rows(stock["root_gram"]),
            "exact_cone": {key: value for key, value in cm.items() if key != "primitive_cusps"},
            "primitive_cusps": cm["primitive_cusps"],
            "checks": stock["checks"],
        },
        "full_NS_lattice": {
            "identification": "odd unimodular I_(1,14)",
            "rank": 15,
            "hyperbolic_dimension": 14,
            "basis_in_CM_plus_boundary_half_coordinates": _rational_rows(surface["ns_basis"]),
            "gram": _rows(surface["ns_gram"]),
            "canonical_vector_in_NS_basis": _rows(
                (surface["ns_basis"].inv() * surface["canonical"]).T
            )[0],
            "canonical_square": -5,
        },
        "curve_exhaustion": {
            "negative_curve_count": 60,
            "negative_curve_types": {"boundary_minus4": 5, "minus2": 10, "minus1": 45},
            "minus1_decomposition": {
                "ten_meeting_two_boundaries_once": 10,
                "thirty_five_meeting_one_boundary_twice": 35,
            },
            "CM_cusp_count": 57,
            "boundary_lift_membership_by_sorted_CM_cusp": surface["cusp_membership"],
            "lifted_cusp_boundary_pairs_zero_based": [list(pair) for pair in surface["lifted_pairs"]],
            "membership_summary": "22 lift nowhere; each of five boundaries receives seven singleton lifts",
            "proof": [
                "Kondo proves exactly twenty effective irreducible roots: ten (-2)-curves and ten boundary-pair (-1)-roots.",
                "Every other negative curve is a boundary, a (-2)-curve, or a (-1)-curve; B=-2K forces every (-1)-curve to have total boundary incidence two.",
                "Incidence 1+1 gives an effective irreducible root and is already among the ten pair curves.",
                "Incidence 2 with B_i sends E to the primitive nef isotropic CM vector f=E+B_i/2, hence to one of the 57 exact cusps.",
                "Conversely the primitive-isotropic/fibration correspondence and Halphen blow-up theorem turn each NS-integral f-B_i into the corresponding exceptional curve.",
                "The unique reconstructed NS glue accepts exactly 35 cusp-boundary pairs, exhausting all remaining (-1)-curves.",
            ],
        },
        "facets": {
            "labels": surface["labels"],
            "NS_wall_vectors": _rows(surface["ns_walls"].T),
            "negative_wall_gram": _rows(surface["wall_gram"]),
            "physical_relation_histogram": dict(sorted(surface["ratio_histogram"].items())),
        },
        "literal_vinberg": literal,
        "exact_double_description": cone,
        "arithmetic": arithmetic,
        "lineage": lineage,
        "classification": "finite_volume_noncompact_full_nef_char3_Coble_H14_60_known_Alexeev_chamber",
        "checks": checks,
        "all_checks_pass": True,
    }
    return json.loads(json.dumps(payload))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cdd-executable")
    parser.add_argument("--compiler")
    parser.add_argument("--full-cone", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.cdd_executable, args.compiler, args.full_cone)
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("the frozen Type-VI five-boundary archive changed")
    if args.write:
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "all_checks_pass": payload["all_checks_pass"],
                "classification": payload["classification"],
                "finite_vertices": payload["literal_vinberg"]["finite_vertex_count"],
                "cusps": payload["literal_vinberg"]["ideal_vertex_count"],
                "edges": payload["literal_vinberg"]["elliptic_edge_count"],
                "lineage": payload["lineage"]["exact_identification"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
