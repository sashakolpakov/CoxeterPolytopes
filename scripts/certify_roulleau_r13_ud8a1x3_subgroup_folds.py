#!/usr/bin/env python3
"""Certify every high fixed section of Roulleau's 39-wall rank-13 stock.

The complete K3 nef chamber has
``NS(X)=U+D8(-1)+A1(-1)^3``.  Roulleau gives the first nineteen curve
intersections, an integral basis ``A1,...,A13``, the class ``D2``, and the
closed formula for ten complementary pairs ``A_J,B_J``.  Those data are
enough to reconstruct all thirty-nine divisor walls without guessing a
missing incidence picture.

This certificate computes the full weighted-wall group, all subgroup
conjugacy classes, the exact Coxeter gate, Vinberg's literal finite-volume
criterion, reflection-stable integral overclosures, and bounded diagram
lineage.  Diagram/lattice actions are kept separate from automorphisms of a
general K3 surface.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Sequence

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ as POLY_ZZ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from scripts.certify_alexeev_remaining_fixed_sections_lineage import _solomon
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _overlattice_basis,
    _q_histogram,
    _reflection,
    _stable_subgroups,
)
from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)
from scripts.certify_jacobian_kummer_h10_maximal_reflection_lineage import (
    _permutation as _relation_permutation,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import _relation_matrix
from scripts.certify_roulleau_non2elementary_involution_folds import _scan_coxiter
from scripts.certify_roulleau_ud8d4_all_subgroup_folds import (
    _gap_subgroup_census,
    _greedy_generators,
)
from scripts.search_k3_high_rank_involution_folds import (
    _compose,
    _group_checks,
    _weighted_graph_automorphisms,
)
from scripts.search_k3_nef_subgroup_folds import (
    Subgroup,
    _orbit_gram,
    _orbits,
    _order,
    _pairwise_disjoint,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_r13_ud8a1x3_all_subgroup_folds.json"
SOURCE_INVENTORY = ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
SOURCE_MAGMA = Path("/private/tmp/roulleau_source/anc/ExamplesMagma.txt")
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
SK_ROOT = Path("/private/tmp/Scharlau-Kirschmer-Lorentz")
EVEN_DET6_FULL = ROOT / "data/enumerations/roulleau_ue8a1x3_all_subgroup_folds.json"
SCHEMA = "coxeter.roulleau_r13_ud8a1x3.all_subgroup_folds.v1"
RANK = 13
WALL_COUNT = 39


D2 = sp.Matrix([-1, 1, 1, 1, 1, 3, 3, 3, 1, 2, 1, -1, -2])
D86 = sp.Matrix([2, 1, 2, 5, 5, 11, 8, 6, 2, 3, 1, 0, 1])

EXPECTED_HIGH = {
    4: {"order": 2, "class_size": 10, "fixed_rank": 12, "status": "positive"},
    6: {"order": 2, "class_size": 15, "fixed_rank": 11, "status": "isotropic_reject"},
    7: {"order": 3, "class_size": 10, "fixed_rank": 11, "status": "positive"},
    12: {"order": 4, "class_size": 15, "fixed_rank": 11, "status": "isotropic_reject"},
    22: {"order": 6, "class_size": 10, "fixed_rank": 11, "status": "positive"},
}

EXPECTED_TYPES = {
    4: {
        "id": "R13-UD8A1x3-C2-H11-31",
        "classes": [4],
        "dimension": 11,
        "facets": 31,
        "finite_counts": [31, 318, 1726, 5856, 13436, 21684, 25003, 20491, 11510, 3948, 552],
        "cusps": 41,
        "affine_components": 121,
        "endpoint_histogram": {"0,2": 186, "1,1": 1452, "2,0": 2310},
        "fixed_determinant": -32,
        "fixed_smith": [1] * 8 + [2, 2, 2, 4],
        "stable_orders": [1, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4],
        "stable_even_orders": [1, 2, 2, 2],
        "maximal_index": 4,
        "maximal_count": 3,
        "maximal_determinant": -2,
        "maximal_parities": [False, False, False],
        "sk_targets": [("List11", 2, "I_{11,1}([1^{11}2^{1}]_{2})")],
        "section_euler": "0",
        "relation_sha256": "ef3ad1896e862f6a1d8cdb57f31f0519d73754ab91578b2bddb5203d8eab2df4",
        "root_types": {"-2,1,1": 15, "-2,2,1": 8, "-4,2,2": 2, "-4,4,2": 6},
    },
    7: {
        "id": "R13-UD8A1x3-C3S3-H10-23",
        "classes": [7, 22],
        "dimension": 10,
        "facets": 23,
        "finite_counts": [23, 188, 826, 2270, 4186, 5348, 4765, 2893, 1086, 170],
        "cusps": 17,
        "affine_components": 42,
        "endpoint_histogram": {"0,2": 42, "1,1": 388, "2,0": 656},
        "fixed_determinant": 24,
        "fixed_smith": [1] * 8 + [2, 2, 6],
        "stable_orders": [1, 2, 2, 2],
        "stable_even_orders": [1, 2, 2],
        "maximal_index": 2,
        "maximal_count": 3,
        "maximal_determinant": 6,
        "maximal_parities": [False, True, True],
        "sk_targets": [
            ("List10", 8, "I_{10,1}([1^{-10}2^{1}]_{7} 3^{-1})"),
            ("List10", 9, "II_{10,1}([2^{1}]_{7} 3^{-1})"),
        ],
        "section_euler": "-1891/1393459200",
        "relation_sha256": "a8ef1e076e386d5378ae7ca193b6a5058e649f2637bae43cbd5b5a4159b21c66",
        "root_types": {"-2,1,1": 13, "-2,2,1": 2, "-6,3,3": 2, "-6,6,3": 6},
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(matrix[row, col]) for col in range(matrix.cols)] for row in range(matrix.rows)]


def _smith(matrix: sp.Matrix) -> list[int]:
    return [abs(int(value)) for value in smith_normal_form(matrix, domain=ZZ).diagonal()]


def _edge(gram: sp.Matrix, left: int, right: int, weight: int = 1) -> None:
    gram[left - 1, right - 1] = gram[right - 1, left - 1] = weight


def _first_nineteen_gram() -> sp.Matrix:
    gram = -2 * sp.eye(19)
    for left, right in (
        (12, 13), (12, 14), (18, 12), (12, 15), (15, 3),
        (3, 6), (6, 7), (12, 16), (12, 17), (7, 8),
        (8, 9), (8, 10), (10, 11), (18, 11), (11, 19),
        (6, 1), (6, 2), (6, 4), (6, 5),
    ):
        _edge(gram, left, right)
    for left, right in ((1, 13), (2, 14), (3, 15), (4, 16), (5, 17)):
        _edge(gram, left, right, 2)
    return gram


def _coordinate_model() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, list[tuple[int, int]]]:
    first = _first_nineteen_gram()
    basis_gram = first[:RANK, :RANK]
    coordinates = [basis_gram.inv() * first[:RANK, wall] for wall in range(19)]
    labels: list[tuple[int, int]] = []
    for left, right in combinations(range(5), 2):
        curve = (
            -coordinates[left]
            - coordinates[right]
            + sum((coordinates[index] for index in range(5)), sp.zeros(RANK, 1))
            + 2 * coordinates[5]
            + coordinates[6]
            - coordinates[8]
        )
        coordinates.extend((curve, 2 * D2 - curve))
        labels.append((left + 1, right + 1))
    walls = sp.Matrix.hstack(*coordinates)
    wall_gram = walls.T * basis_gram * walls
    return basis_gram, walls, wall_gram, labels


def _source_audit(
    inventory: dict[str, Any],
    basis_gram: sp.Matrix,
    walls: sp.Matrix,
    wall_gram: sp.Matrix,
    labels: list[tuple[int, int]],
) -> dict[str, Any]:
    stock = next(row for row in inventory["stocks"] if row["key"] == "R13_UD8A1x3")
    d2_pairings = [int(value) for value in list(D2.T * basis_gram * walls)]
    d86_pairings = [int(value) for value in list(D86.T * basis_gram * walls)]
    petersen_checks = []
    for left, label_left in enumerate(labels):
        a_left = 19 + 2 * left
        b_left = a_left + 1
        for right, label_right in enumerate(labels):
            a_right = 19 + 2 * right
            b_right = a_right + 1
            intersection = len(set(label_left) & set(label_right))
            expected_a = -2 if left == right else (2 if intersection == 0 else 0)
            petersen_checks.append(
                wall_gram[a_left, a_right] == expected_a
                and wall_gram[b_left, b_right] == expected_a
                and wall_gram[a_left, b_right] == 4 - expected_a
            )
    tex_lines = SOURCE_TEX.read_text().splitlines()
    tex_excerpt = "\n".join(tex_lines[9249:9490])
    magma_lines = SOURCE_MAGMA.read_text().splitlines()
    magma_excerpt = "\n".join(magma_lines[144:214])
    checks = {
        "inventory_certified": inventory["all_checks_pass"] is True,
        "inventory_stock": (
            stock["picard_rank"] == 13
            and stock["wall_count"] == 39
            and stock["lattice"] == "U+D8(-1)+A1(-1)^3"
        ),
        "source_hashes_match_inventory": (
            _sha256(SOURCE_TEX) == inventory["source"]["sha256"]["tex"]
            and _sha256(SOURCE_MAGMA) == inventory["source"]["sha256"]["examples_magma"]
        ),
        "tex_formula_lane_present": (
            "contains $39$ $(-2)$-curves" in tex_excerpt
            and "$A_{1},\\dots,A_{13}$ generate" in tex_excerpt
            and "A_{J}=-A_{i}-A_{j}" in tex_excerpt
            and "B_{J}=2D_{2}-A_{J}" in tex_excerpt
        ),
        "magma_39_wall_ordering_present": (
            "Lattice \"U\",\"D8\",\"A1\",\"A1\",\"A1\"" in magma_excerpt
            and "SE:=[17,16,15,6,7,8" in magma_excerpt
        ),
        "basis_gram": basis_gram.det() == 32 and exact_inertia(basis_gram) == (1, 12, 0),
        "all_coordinates_integral": all(value.is_Integer for value in walls),
        "complete_gram_rank_inertia": (
            wall_gram.rank() == 13 and exact_inertia(wall_gram) == (1, 12, 26)
        ),
        "all_39_classes_are_minus_two": all(wall_gram[index, index] == -2 for index in range(39)),
        "first_19_printed_graph_replayed": wall_gram[:19, :19] == _first_nineteen_gram(),
        "D2_square_and_pairings": (
            (D2.T * basis_gram * D2)[0] == 2
            and d2_pairings == [1] * 5 + [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1] + [2] * 20
        ),
        "D86_ample_square_and_strict_pairings": (
            (D86.T * basis_gram * D86)[0] == 86 and min(d86_pairings) > 0
        ),
        "ten_complementary_pair_relations": all(
            walls[:, 19 + 2 * index] + walls[:, 20 + 2 * index] == 2 * D2
            for index in range(10)
        ),
        "weighted_Petersen_and_cross_formulas": all(petersen_checks),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"R13/39 source reconstruction changed: {checks}")
    return {
        "stock_inventory_record": stock,
        "authoritative_coordinate_formula_lane": {
            "first_19_graph_tex_lines": [9258, 9340],
            "basis_D86_D2_and_pair_relations_tex_lines": [9342, 9382],
            "AJ_BJ_formula_tex_lines": [9404, 9421],
            "canonical_lexicographic_J_labels_for_A20_to_A39": [list(label) for label in labels],
            "basis_gram": _payload(basis_gram),
            "wall_coordinate_columns": _payload(walls),
            "wall_gram": _payload(wall_gram),
            "rank": 13,
            "inertia": [1, 12, 26],
        },
        "printed_graph_status": (
            "the printed graph is used only for A1--A19 and is exactly consistent; "
            "all twenty extra walls come from the displayed AJ/BJ formulas"
        ),
        "source_completeness": (
            "published induced-chamber completion inside II_(1,17), not merely a degree cutoff"
        ),
        "checks": checks,
    }


def _action(permutation: Sequence[int], walls: sp.Matrix) -> sp.Matrix:
    return sp.Matrix.hstack(*(walls[:, permutation[index]] for index in range(RANK)))


def _fixed_basis(subgroup: Subgroup, walls: sp.Matrix) -> sp.Matrix:
    equations = sp.Matrix.vstack(*(_action(element, walls) - sp.eye(RANK) for element in subgroup))
    diagonal, _, right = smith_normal_decomp(
        DomainMatrix.from_Matrix(equations).convert_to(POLY_ZZ)
    )
    normal = diagonal.to_Matrix()
    rank = sum(normal[index, index] != 0 for index in range(min(normal.rows, normal.cols)))
    transform = right.to_Matrix()
    basis = transform[:, rank:]
    if not (abs(int(transform.det())) == 1 and equations * basis == sp.zeros(equations.rows, basis.cols)):
        raise ArithmeticError("fixed Smith kernel is not saturated")
    return basis


def _group_audit(
    group: tuple[tuple[int, ...], ...],
    walls: sp.Matrix,
    basis_gram: sp.Matrix,
) -> tuple[dict[str, Any], dict[tuple[int, ...], sp.Matrix]]:
    actions: dict[tuple[int, ...], sp.Matrix] = {}
    for permutation in group:
        action = _action(permutation, walls)
        if not (
            all(value.is_Integer for value in action)
            and abs(int(action.det())) == 1
            and action.T * basis_gram * action == basis_gram
            and all(action * walls[:, index] == walls[:, permutation[index]] for index in range(WALL_COUNT))
        ):
            raise ArithmeticError("weighted graph symmetry did not lift integrally")
        actions[permutation] = action
    blocks = [{index, index + 12} for index in range(5)]
    induced = set()
    kernel = []
    for permutation in group:
        image = []
        for block in blocks:
            target = {permutation[index] for index in block}
            image.append(blocks.index(target))
        induced.add(tuple(image))
        if image == list(range(5)):
            kernel.append(permutation)
    center = [
        element
        for element in group
        if all(_compose(element, other) == _compose(other, element) for other in group)
    ]
    central = next(element for element in center if _order(element) == 2)
    central_fixed_rank = RANK - (actions[central] - sp.eye(RANK)).rank()
    checks = {
        **_group_checks(group),
        "order_240": len(group) == 240,
        "five_paired_line_blocks_give_S5": len(induced) == math.factorial(5),
        "block_kernel_is_central_C2": len(kernel) == len(center) == 2 and set(kernel) == set(center),
        "central_geometric_candidate_has_only_H4": central_fixed_rank == 5,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"R13/39 group recognition changed: {checks}")
    return (
        {
            "order": len(group),
            "recognition": "S5 x C2",
            "id_group": [240, 189],
            "element_order_histogram": {
                str(key): value for key, value in sorted(Counter(_order(row) for row in group).items())
            },
            "central_C2_permutation_one_based": [value + 1 for value in central],
            "central_C2_fixed_rank": central_fixed_rank,
            "checks": checks,
        },
        actions,
    )


def _strict_facet_witnesses(
    subgroup: Subgroup,
    fixed_basis: sp.Matrix,
    fixed_gram: sp.Matrix,
    roots: sp.Matrix,
    walls: sp.Matrix,
) -> dict[str, Any]:
    invariant_ample = sum((_action(element, walls) * D86 for element in subgroup), sp.zeros(RANK, 1))
    ample = fixed_basis.gauss_jordan_solve(invariant_ample)[0]
    rows = []
    for index in range(roots.cols):
        root = roots[:, index]
        norm = (root.T * fixed_gram * root)[0]
        pairing = (ample.T * fixed_gram * root)[0]
        witness = ample - pairing / norm * root
        pairings = list(witness.T * fixed_gram * roots)
        square = (witness.T * fixed_gram * witness)[0]
        rows.append(
            {
                "facet_one_based": index + 1,
                "witness_square": str(square),
                "target_pairing": str(pairings[index]),
                "minimum_other_pairing": str(min(value for pos, value in enumerate(pairings) if pos != index)),
                "strict_timelike_facet_witness": (
                    square > 0
                    and pairings[index] == 0
                    and all(value > 0 for pos, value in enumerate(pairings) if pos != index)
                ),
            }
        )
    return {
        "method": "orthogonal projection of the subgroup-averaged ample D86 class",
        "records": rows,
        "all_orbit_inequalities_are_genuine_facets": all(row["strict_timelike_facet_witness"] for row in rows),
    }


def _positive_internal(
    class_index: int,
    subgroup: Subgroup,
    orbits: tuple[tuple[int, ...], ...],
    basis_gram: sp.Matrix,
    walls: sp.Matrix,
    wall_gram: sp.Matrix,
) -> dict[str, Any]:
    fixed_basis = _fixed_basis(subgroup, walls)
    fixed_gram = fixed_basis.T * basis_gram * fixed_basis
    sums = sp.zeros(WALL_COUNT, len(orbits))
    for column, orbit in enumerate(orbits):
        for row in orbit:
            sums[row, column] = 1
    roots = fixed_basis.gauss_jordan_solve(walls * sums)[0]
    orbit_gram = roots.T * fixed_gram * roots
    if orbit_gram != _orbit_gram(wall_gram, orbits):
        raise ArithmeticError("fixed-lattice and direct orbit Grams disagree")
    reflection_checks = []
    reflections = []
    root_types: Counter[tuple[int, int, int]] = Counter()
    for index, orbit in enumerate(orbits):
        root = roots[:, index]
        norm = int((root.T * fixed_gram * root)[0])
        divisibility = math.gcd(*(abs(int(value)) for value in fixed_gram * root))
        reflection = _reflection(fixed_gram, root)
        disjoint = _pairwise_disjoint(wall_gram, orbit)
        reflection_checks.append(
            {
                "wall_orbit_one_based": [value + 1 for value in orbit],
                "norm": norm,
                "divisibility": divisibility,
                "pairwise_disjoint": disjoint,
                "integral_involutive_isometry": (
                    all(value.is_Integer for value in reflection)
                    and reflection * reflection == sp.eye(fixed_gram.rows)
                    and reflection.T * fixed_gram * reflection == fixed_gram
                ),
            }
        )
        reflections.append(reflection)
        root_types[(norm, divisibility, len(orbit))] += 1
    positive_gram = -orbit_gram
    relations = _relation_matrix(_payload(positive_gram))
    return {
        "class_index": class_index,
        "subgroup": subgroup,
        "orbits": orbits,
        "fixed_basis": fixed_basis,
        "fixed_gram": fixed_gram,
        "roots": roots,
        "orbit_gram": orbit_gram,
        "positive_gram": positive_gram,
        "relations": relations,
        "reflections": reflections,
        "reflection_checks": reflection_checks,
        "root_types": {
            f"{norm},{divisibility},{size}": count
            for (norm, divisibility, size), count in sorted(root_types.items())
        },
        "strict_facet_witnesses": _strict_facet_witnesses(
            subgroup, fixed_basis, fixed_gram, roots, walls
        ),
    }


def _vinberg(executable: Path, gram: sp.Matrix, expected: dict[str, Any]) -> dict[str, Any]:
    rows = [f"{expected['dimension']} {gram.rows}"]
    rows.extend(" ".join(str(int(gram[left, right])) for right in range(gram.cols)) for left in range(gram.rows))
    completed = subprocess.run(
        [str(executable)], input="\n".join(rows) + "\n", text=True, capture_output=True, check=True
    )
    result = _parse_output(completed.stdout)
    counts = result["finite_special_subset_counts_by_rank"]
    endpoints = result["edge_endpoint_type_histogram"]
    checks = {
        "dimension_facets": result["dimension"] == expected["dimension"] and result["facet_count"] == expected["facets"],
        "finite_special_counts": counts == expected["finite_counts"],
        "finite_subset_search_complete": result["all_finite_subset_search_node_count"] == 1 + sum(counts),
        "cusps_and_affine_components": (
            result["ideal_vertex_count"] == expected["cusps"]
            and result["connected_affine_component_count"] == expected["affine_components"]
        ),
        "zero_bad_edges": result["bad_edge_count"] == 0,
        "zero_bad_parabolic_neighborhoods": result["bad_parabolic_neighborhood_count"] == 0,
        "endpoint_histogram": endpoints == expected["endpoint_histogram"],
        "two_endpoints_per_edge": sum((int(key[0]) + int(key[2])) * value for key, value in endpoints.items()) == 2 * counts[-2],
        "finite_vertex_edge_incidence": 2 * endpoints.get("2,0", 0) + endpoints.get("1,1", 0) == expected["dimension"] * counts[-1],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"R13/39 Vinberg audit changed: {checks}")
    return result | {
        "elliptic_edge_count": counts[-2],
        "finite_vertex_count": counts[-1],
        "finite_volume": True,
        "compact": False,
        "checks": checks,
        "all_checks_pass": True,
    }


def _bilinear_diagonal_histogram(gram: sp.Matrix) -> dict[str, int]:
    moduli, _, discriminant_gram = _discriminant_data(gram)
    histogram: Counter[str] = Counter()
    for element in product(*(range(modulus) for modulus in moduli)):
        vector = sp.Matrix(element)
        histogram[str(sp.factor((vector.T * discriminant_gram * vector)[0]) % 1)] += 1
    return dict(sorted(histogram.items()))


def _sk_row(list_name: str, entry: int, symbol: str) -> dict[str, Any]:
    path = SK_ROOT / list_name
    rows = [line.strip().rstrip(",") for line in path.read_text().splitlines() if line.strip().startswith("<")]
    row = rows[entry - 1]
    if symbol not in row:
        raise ArithmeticError(f"SK row changed: {list_name}/{entry}: {row}")
    bracket_groups = __import__("re").findall(r"\[([^]]*)\]", row)
    counts = [int(value) for value in __import__("re").findall(r"-?\d+", bracket_groups[-1])]
    cusp = int(__import__("re").search(r",\s*(\d+)\s*>$", row).group(1))
    return {
        "list": list_name,
        "list_sha256": _sha256(path),
        "entry_one_based": entry,
        "genus_symbol": symbol,
        "row": row,
        "full_chamber_facets": counts[0],
        "full_chamber_finite_vertices": counts[-1],
        "full_chamber_cusps": cusp,
    }


def _stable_lineage(internal: dict[str, Any], expected: dict[str, Any]) -> dict[str, Any]:
    fixed_gram = internal["fixed_gram"]
    stable = _stable_subgroups(fixed_gram, internal["reflections"])
    stable.pop("groups_internal")
    left_inverse = stable.pop("left_smith_inverse")
    maxima = []
    for position in stable["maximal_stable_subgroup_positions_zero_based"]:
        basis = _overlattice_basis(
            fixed_gram,
            left_inverse,
            stable["stable_subgroups"][position]["cosets_in_smith_coordinates"],
        )
        gram = sp.simplify(basis.T * fixed_gram * basis)
        actions = [sp.simplify(basis.inv() * reflection * basis) for reflection in internal["reflections"]]
        even = all(gram[index, index] % 2 == 0 for index in range(gram.rows))
        checks = {
            "index": abs(sp.Rational(1, basis.det())) == expected["maximal_index"],
            "determinant": int(gram.det()) == expected["maximal_determinant"],
            "smith": _smith(gram) == [1] * (gram.rows - 1) + [abs(expected["maximal_determinant"])],
            "reflections_integral": all(value.is_Integer for action in actions for value in action),
            "reflections_are_isometries": all(action.T * gram * action == gram for action in actions),
        }
        if not all(checks.values()):
            raise ArithmeticError(f"stable maximum changed: {checks}")
        target = next(
            _sk_row(list_name, entry, symbol)
            for list_name, entry, symbol in expected["sk_targets"]
            if (list_name == "List10" and ((entry == 8) == (not even))) or list_name == "List11"
        )
        exact_full_bridge = None
        numerical_index = None
        if expected["dimension"] == 10 and even:
            full_source = json.loads(EVEN_DET6_FULL.read_text())
            full_type = next(
                row for row in full_source["types"]
                if row["type_id"] == "R13_UE8A1x3_C3_S3"
            )
            target_q = full_type["maximal_reflection_lineage"]["unique_maximal_lattice"][
                "discriminant_q_histogram_mod_two"
            ]
            full_euler = Fraction(full_type["solomon"]["rational_euler_characteristic"])
            section_euler = Fraction(expected["section_euler"])
            numerical_index = section_euler / full_euler
            bridge_checks = {
                "target_is_exact_full_index_one_chamber": (
                    full_type["maximal_reflection_lineage"]["full_reflection_group"][
                        "reflection_subgroup_index"
                    ] == 1
                ),
                "rank_signature_parity_determinant_smith": (
                    gram.rows == full_type["fixed_lattice"]["rank"] == 11
                    and exact_inertia(gram) == (1, 10, 0)
                    and int(gram.det()) == full_type["fixed_lattice"]["determinant"] == 6
                    and _smith(gram) == full_type["fixed_lattice"]["smith_invariants"] == [1] * 10 + [6]
                    and even and full_type["fixed_lattice"]["even"]
                ),
                "discriminant_quadratic_form": _q_histogram(gram) == target_q,
                "unique_indefinite_even_genus_range": gram.rows >= 1 + 2,
                "exact_Solomon_ratio": numerical_index == 2046,
            }
            if not all(bridge_checks.values()):
                raise ArithmeticError(f"even determinant-six full-chamber bridge changed: {bridge_checks}")
            exact_full_bridge = {
                "source_path": str(EVEN_DET6_FULL.relative_to(ROOT)),
                "source_sha256": _sha256(EVEN_DET6_FULL),
                "source_type_id": full_type["type_id"],
                "classification_theorem": (
                    "indefinite even lattice is unique in its genus because rank 11 "
                    "is at least discriminant length 1 plus 2"
                ),
                "full_chamber_facets": full_type["geometry"]["facet_count"],
                "full_chamber_finite_vertices": full_type["geometry"]["finite_vertex_count"],
                "full_chamber_cusps": full_type["geometry"]["ideal_vertex_count"],
                "section_euler": str(section_euler),
                "full_euler": str(full_euler),
                "section_to_full_index": int(numerical_index),
                "checks": bridge_checks,
            }
        maxima.append(
            {
                "stable_subgroup_position_zero_based": position,
                "index": expected["maximal_index"],
                "basis_columns_in_fixed_lattice": _payload(basis),
                "gram": _payload(gram),
                "determinant": int(gram.det()),
                "smith_invariants": _smith(gram),
                "even": even,
                "discriminant_bilinear_diagonal_mod_1_histogram": _bilinear_diagonal_histogram(gram),
                "discriminant_q_mod_2_histogram": _q_histogram(gram) if even else None,
                "Scharlau_Kirschmer_target": target,
                "exact_full_chamber_bridge": exact_full_bridge,
                "section_to_full_reflection_group_index": (
                    int(numerical_index) if numerical_index is not None else None
                ),
                "containment_index_status": (
                    "finite; numerical value open because both H11 Coxeter characteristics vanish"
                    if expected["dimension"] == 11
                    else (
                        "exact Solomon ratio through the known full even determinant-six chamber"
                        if even
                        else "finite; numerical value awaits an explicit full-chamber relation matrix for the pinned odd SK row"
                    )
                ),
                "checks": checks,
            }
        )
    parities = sorted((row["even"] for row in maxima))
    expected_parities = sorted(expected["maximal_parities"])
    checks = {
        "stable_orders": stable["stable_integral_subgroup_orders"] == expected["stable_orders"],
        "stable_even_orders": stable["stable_even_subgroup_orders"] == expected["stable_even_orders"],
        "maximal_count": len(maxima) == expected["maximal_count"],
        "maximal_parities": parities == expected_parities,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"stable-overlattice denominator changed: {checks}")
    return {
        "exhaustive_discriminant_module_audit": stable,
        "maximal_reflection_stable_overclosures": maxima,
        "new_arithmetic_commensurability_class": False,
        "checks": checks,
    }


def _type_record(
    internal: dict[str, Any],
    class_indices: list[int],
    executable: Path,
) -> dict[str, Any]:
    expected = EXPECTED_TYPES[internal["class_index"]]
    fixed_gram = internal["fixed_gram"]
    relation_sha = _canonical_sha256(internal["relations"])
    vinberg = _vinberg(executable, internal["positive_gram"], expected)
    stable = _stable_lineage(internal, expected)
    coxiter = _scan_coxiter(internal["relations"], expected["dimension"], expected["facets"])
    solomon = _solomon(_payload(internal["positive_gram"]))
    checks = {
        "classes": class_indices == expected["classes"],
        "dimension_facets": fixed_gram.rows - 1 == expected["dimension"] and internal["roots"].cols == expected["facets"],
        "fixed_lattice": (
            int(fixed_gram.det()) == expected["fixed_determinant"]
            and _smith(fixed_gram) == expected["fixed_smith"]
            and all(fixed_gram[index, index] % 2 == 0 for index in range(fixed_gram.rows))
        ),
        "root_types": internal["root_types"] == expected["root_types"],
        "all_reflections": all(
            row["pairwise_disjoint"] and row["integral_involutive_isometry"]
            for row in internal["reflection_checks"]
        ),
        "all_facets_irredundant": internal["strict_facet_witnesses"]["all_orbit_inequalities_are_genuine_facets"],
        "relation_sha": relation_sha == expected["relation_sha256"],
        "solomon": solomon["rational_euler_characteristic"] == expected["section_euler"],
        "coxiter_nonmatch": not coxiter["ordinary_exact_matches"] and not coxiter["doubling_exact_matches"],
        "vinberg": vinberg["all_checks_pass"],
        "stable_lineage": all(stable["checks"].values()),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"positive type changed for {expected['id']}: {checks}")
    return {
        "id": expected["id"],
        "representative_subgroup_class_one_based": internal["class_index"],
        "source_subgroup_classes_one_based": class_indices,
        "hyperbolic_dimension": expected["dimension"],
        "facet_count": expected["facets"],
        "volume_type": "finite_volume_noncompact",
        "fixed_lattice": {
            "saturated_basis_columns_in_NS_basis": _payload(internal["fixed_basis"]),
            "gram": _payload(fixed_gram),
            "rank": fixed_gram.rows,
            "determinant": int(fixed_gram.det()),
            "smith_invariants": _smith(fixed_gram),
            "even": True,
            "discriminant_q_mod_2_histogram": _q_histogram(fixed_gram),
        },
        "folded_walls": {
            "wall_orbits_one_based": [[value + 1 for value in orbit] for orbit in internal["orbits"]],
            "root_coordinate_columns": _payload(internal["roots"]),
            "negative_gram": _payload(internal["orbit_gram"]),
            "positive_gram": _payload(internal["positive_gram"]),
            "normalized_relation_matrix": internal["relations"],
            "relation_matrix_sha256": relation_sha,
            "root_type_norm_divisibility_orbit_size_histogram": internal["root_types"],
            "reflection_checks": internal["reflection_checks"],
            "strict_facet_witnesses": internal["strict_facet_witnesses"],
        },
        "geometry": {"vinberg": vinberg},
        "stable_overlattices": stable,
        "lineage": {
            "CoxIter": coxiter,
            "solomon": solomon,
            "bounded_literal_status": (
                "absent from pinned current CoxIter and index-two-doubling snapshots; worldwide priority open"
            ),
            "prior_positive_ledger_same_dimension_facet_pair_count": 0,
            "new_arithmetic_commensurability_class": False,
        },
        "surface_geometry": {
            "ambient_object": "complete nef cone of Roulleau's K3 surface with NS=U+D8(-1)+A1(-1)^3",
            "general_surface_automorphism_group": "C2 x C2",
            "published_action_tag": "tn",
            "selected_survivor_is_central_diagram_C2": False,
            "selected_nonidentity_subgroup_realized_on_general_surface": False,
            "surface_automorphism_not_needed_for_polytope_proof": True,
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def _same_h10_integral_chamber(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any]:
    permutation = _relation_permutation(left["relations"], right["relations"])
    if permutation is None:
        raise ArithmeticError("C3 and S3 H10 diagrams stopped being isomorphic")
    target_roots = sp.Matrix.hstack(*(right["roots"][:, permutation[index]] for index in range(len(permutation))))
    independent: list[int] = []
    for column in range(left["roots"].cols):
        if left["roots"][:, [*independent, column]].rank() > len(independent):
            independent.append(column)
    isometry = target_roots[:, independent] * left["roots"][:, independent].inv()
    checks = {
        "rank_spanning_roots": len(independent) == left["fixed_gram"].rows,
        "integral_unimodular": all(value.is_Integer for value in isometry) and abs(int(isometry.det())) == 1,
        "gram_isometry": isometry.T * right["fixed_gram"] * isometry == left["fixed_gram"],
        "all_23_roots_conjugated": isometry * left["roots"] == target_roots,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"C3/S3 integral chamber identification changed: {checks}")
    return {
        "source_class": 7,
        "target_class": 22,
        "facet_permutation_zero_based": permutation,
        "independent_source_root_columns_zero_based": independent,
        "integral_unimodular_isometry": _payload(isometry),
        "determinant": int(isometry.det()),
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    inventory = json.loads(SOURCE_INVENTORY.read_text())
    basis_gram, walls, wall_gram, labels = _coordinate_model()
    source = _source_audit(inventory, basis_gram, walls, wall_gram, labels)
    group, search_nodes = _weighted_graph_automorphisms(wall_gram)
    group_record, _ = _group_audit(group, walls, basis_gram)
    group_record["exhaustive_backtracking_node_count"] = search_nodes
    generators = _greedy_generators(group)
    census = _gap_subgroup_census(group, generators)
    if not (
        census["header"]["structure"] == "C2 x S5"
        and census["header"]["id_group"] == [240, 189]
        and census["header"]["subgroup_class_count"] == 57
    ):
        raise ArithmeticError(f"GAP subgroup denominator changed: {census['header']}")

    class_rows = []
    positive_internal: dict[int, dict[str, Any]] = {}
    actual_high: dict[int, dict[str, Any]] = {}
    for row in census["rows"]:
        subgroup = row["representative_subgroup"]
        orbits = _orbits(subgroup, WALL_COUNT)
        folded = _orbit_gram(wall_gram, orbits)
        fixed_rank = folded.rank()
        norms = [int(folded[index, index]) for index in range(folded.rows)]
        high = len(subgroup) > 1 and fixed_rank >= 11
        status = "below_H10_or_identity"
        if high:
            status = "isotropic_orbit_normal_reject" if 0 in norms else "exact_coxeter_candidate"
            actual_high[row["class_index"]] = {
                "order": row["order"],
                "class_size": row["conjugacy_class_size"],
                "fixed_rank": fixed_rank,
                "status": "isotropic_reject" if 0 in norms else "positive",
            }
            if 0 not in norms:
                positive_internal[row["class_index"]] = _positive_internal(
                    row["class_index"], subgroup, orbits, basis_gram, walls, wall_gram
                )
        class_rows.append(
            {
                "class_index_one_based": row["class_index"],
                "subgroup_order": row["order"],
                "conjugacy_class_size": row["conjugacy_class_size"],
                "structure": row["structure"],
                "fixed_rank": fixed_rank,
                "hyperbolic_dimension": fixed_rank - 1,
                "orbit_count": len(orbits),
                "isotropic_orbit_indices_one_based": [index + 1 for index, norm in enumerate(norms) if norm == 0],
                "status": status,
                "wall_orbits_one_based": [[value + 1 for value in orbit] for orbit in orbits] if high else None,
            }
        )
    if actual_high != EXPECTED_HIGH:
        raise ArithmeticError(f"R13/39 high subgroup denominator changed: {actual_high}")

    h10_identification = _same_h10_integral_chamber(positive_internal[7], positive_internal[22])
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="roulleau-r13-ud8a1x3-vinberg-") as raw:
        executable = Path(raw) / "vinberg"
        subprocess.run(
            [compiler, "-O3", "-std=c++17", "-Wall", "-Wextra", "-Werror", str(CPP_SOURCE), "-o", str(executable)],
            check=True,
        )
        types = [
            _type_record(positive_internal[4], [4], executable),
            _type_record(positive_internal[7], [7, 22], executable),
        ]
    checks = {
        "source": all(source["checks"].values()),
        "weighted_group": all(group_record["checks"].values()),
        "all_57_subgroup_classes": len(class_rows) == 57,
        "five_high_classes": len(actual_high) == 5,
        "two_isotropic_rejects": sum(row["status"] == "isotropic_orbit_normal_reject" for row in class_rows) == 2,
        "three_positive_classes_two_types": len(positive_internal) == 3 and len(types) == 2,
        "C3_S3_same_integral_chamber": all(h10_identification["checks"].values()),
        "all_types_finite_volume": all(row["all_checks_pass"] for row in types),
        "no_new_arithmetic_classes": all(not row["lineage"]["new_arithmetic_commensurability_class"] for row in types),
        "no_general_surface_action_claim": all(not row["surface_geometry"]["selected_nonidentity_subgroup_realized_on_general_surface"] for row in types),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"R13/39 global audit changed: {checks}")
    return {
        "schema": SCHEMA,
        "artifact_id": "roulleau_r13_ud8a1x3_all_subgroup_folds",
        "source": {
            "inventory_path": str(SOURCE_INVENTORY.relative_to(ROOT)),
            "inventory_sha256": _sha256(SOURCE_INVENTORY),
            "tex_path": str(SOURCE_TEX),
            "tex_sha256": _sha256(SOURCE_TEX),
            "ancillary_magma_path": str(SOURCE_MAGMA),
            "ancillary_magma_sha256": _sha256(SOURCE_MAGMA),
            "cpp_path": str(CPP_SOURCE.relative_to(ROOT)),
            "cpp_sha256": _sha256(CPP_SOURCE),
            "compiler": Path(compiler).name,
        },
        "source_reconstruction": source,
        "weighted_wall_group": group_record | {
            "gap_version": census["header"]["gap_version"],
            "gap_subgroup_class_count": census["header"]["subgroup_class_count"],
            "greedy_generators_zero_based": [list(row) for row in generators],
        },
        "subgroup_census": {
            "conjugacy_class_count": len(class_rows),
            "proper_fixed_rank_at_least_11_count": len(actual_high),
            "classes": class_rows,
        },
        "C3_S3_integral_chamber_identification": h10_identification,
        "positive_types": types,
        "summary": {
            "source_wall_count": 39,
            "weighted_group_order": 240,
            "subgroup_conjugacy_class_count": 57,
            "proper_high_class_count": 5,
            "isotropic_orbit_reject_count": 2,
            "positive_class_count": 3,
            "distinct_positive_type_count": 2,
            "dimension_facet_pairs": [[11, 31], [10, 23]],
            "finite_volume_noncompact_type_count": 2,
            "exact_CoxIter_match_count": 0,
            "new_arithmetic_commensurability_class_count": 0,
            "worldwide_literal_priority": "open",
            "general_surface_geometric_survivor_count": 0,
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.check_archive:
        if payload != json.loads(args.output.read_text()):
            raise ArithmeticError("rebuilt R13 UD8 A1^3 archive changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
