#!/usr/bin/env python3
"""Certify every high fixed section of Roulleau's rank-14 27-wall stock.

The stock is the complete K3 nef chamber on
``U + E8(-1) + A1(-1)^4``.  Roulleau gives a curve basis, four additional
curve coordinates, and four complementary-pair equivalences.  Those data
uniquely reconstruct all 27 divisor classes.  They are used here instead of
the final cross-incidence TikZ picture, whose literal labels are incompatible
with Picard rank 14; that source inconsistency is itself frozen below.

The exact weighted wall group, all its subgroups modulo conjugacy, every
fixed section above H9, literal Vinberg volume, stable overclosures, and
bounded CoxIter/Scharlau--Kirschmer lineage are then audited exhaustively.
No nonidentity graph action is asserted to be a surface automorphism.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
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
    _overlattice_basis,
    _q_histogram,
    _reflection,
    _stable_subgroups,
)
from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _relation_matrix,
)
from scripts.certify_roulleau_non2elementary_involution_folds import (
    _scan_coxiter,
)
from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    _subgroup_conjugacy_classes,
)
from scripts.search_k3_high_rank_involution_folds import (
    _compose,
    _group_checks,
    _weighted_graph_automorphisms,
)
from scripts.search_k3_nef_subgroup_folds import (
    Subgroup,
    _all_subgroups,
    _certify_fold,
    _orbit_gram,
    _orbits,
    _order,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/roulleau_r14_ue8a1x4_all_subgroup_folds.json"
)
SOURCE_INVENTORY = (
    ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
)
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
SK_ROOT = Path("/private/tmp/Scharlau-Kirschmer-Lorentz")
SCHEMA = "coxeter.roulleau_r14_ue8a1x4.all_subgroup_folds.v1"


D132 = sp.Matrix([2, 5, 9, 4, 10, 12, 15, 19, 9, 15, 12, 5, 5, 0])
D2 = sp.Matrix([1, 2, 3, 1, 3, 3, 3, 3, 1, 2, 1, 0, 0, 0])
EVEN_EXTRA_COORDINATES = {
    20: [0, 0, 0, -1, 1, 2, 3, 4, 2, 3, 2, 0, 0, 1],
    22: [2, 4, 6, 2, 6, 6, 6, 6, 3, 3, 0, -1, -1, -1],
    24: [0, 0, 0, -1, 1, 2, 3, 4, 2, 3, 2, 1, 0, 0],
    26: [0, 0, 0, -1, 1, 2, 3, 4, 2, 3, 2, 0, 1, 0],
}


EXPECTED = {
    2: {
        "id": "R14-UE8A1x4-C2T-H12-23",
        "kind": "transposition_C2",
        "dimension": 12,
        "facets": 23,
        "subgroup_order": 2,
        "class_size": 6,
        "fixed_determinant": 16,
        "fixed_smith": [1] * 10 + [2, 2, 4],
        "stable_orders": [1, 2, 2, 2, 4],
        "stable_even_orders": [1, 2],
        "maximal_index": 4,
        "maximal_determinant": 1,
        "maximal_smith": [1] * 13,
        "maximal_even": False,
        "standard_lattice": "odd unimodular I_(12,1)",
        "sk_list": "List12",
        "sk_entry": 1,
        "sk_symbol": "I_{12,1}()",
        "sk_counts": [14, 91, 364, 1000, 1992, 2958, 3312, 2793, 1749, 786, 234, 35],
        "sk_cusps": 2,
        "full_euler": "691/382588157952000",
        "section_euler": "691/45984153600",
        "section_index": 8320,
        "finite_counts": [23, 215, 1147, 3996, 9728, 17158, 22308, 21439, 15036, 7405, 2308, 320],
        "affine_components": 44,
        "cusps": 17,
        "endpoint_histogram": {"0,2": 52, "1,1": 672, "2,0": 1584},
        "relation_sha256": "b6cdd4a6942030715563084eddfc187d76cba60ae653ddefbe60d1444d1063f5",
        "coxiter_headers": [0, 0],
    },
    3: {
        "id": "R14-UE8A1x4-C2D-H11-19",
        "kind": "double_transposition_C2",
        "dimension": 11,
        "facets": 19,
        "subgroup_order": 2,
        "class_size": 3,
        "fixed_determinant": -16,
        "fixed_smith": [1] * 10 + [4, 4],
        "stable_orders": [1, 2, 2, 2, 4],
        "stable_even_orders": [1, 2],
        "maximal_index": 4,
        "maximal_determinant": -1,
        "maximal_smith": [1] * 12,
        "maximal_even": False,
        "standard_lattice": "odd unimodular I_(11,1)",
        "sk_list": "List11",
        "sk_entry": 1,
        "sk_symbol": "I_{11,1}()",
        "sk_counts": [13, 78, 285, 705, 1242, 1596, 1506, 1035, 504, 162, 26],
        "sk_cusps": 2,
        "full_euler": "0",
        "section_euler": "0",
        "section_index": None,
        "finite_counts": [19, 149, 667, 1940, 3912, 5654, 5930, 4485, 2370, 795, 114],
        "affine_components": 27,
        "cusps": 13,
        "endpoint_histogram": {"0,2": 31, "1,1": 274, "2,0": 490},
        "relation_sha256": "822ed1ce719a59521da7eca5546790bd356cd0b941a49f7dcecb9bbf8f70173a",
        "coxiter_headers": [0, 0],
    },
    7: {
        "id": "R14-UE8A1x4-C3-H11-19",
        "kind": "three_cycle_C3",
        "dimension": 11,
        "facets": 19,
        "subgroup_order": 3,
        "class_size": 4,
        "fixed_determinant": -12,
        "fixed_smith": [1] * 10 + [2, 6],
        "stable_orders": [1, 2],
        "stable_even_orders": [1, 2],
        "maximal_index": 2,
        "maximal_determinant": -3,
        "maximal_smith": [1] * 11 + [3],
        "maximal_even": True,
        "standard_lattice": "U + E8(-1) + A2(-1)",
        "sk_list": "List11",
        "sk_entry": 4,
        "sk_symbol": "II_{11,1}(3^{-1})",
        "sk_counts": [13, 78, 285, 705, 1242, 1596, 1506, 1035, 504, 162, 27],
        "sk_cusps": 1,
        "full_euler": "0",
        "section_euler": "0",
        "section_index": None,
        "finite_counts": [19, 151, 687, 2028, 4134, 6006, 6294, 4735, 2486, 833, 128],
        "affine_components": 19,
        "cusps": 7,
        "endpoint_histogram": {"0,2": 14, "1,1": 230, "2,0": 589},
        "relation_sha256": "25fdeb9504867a25e26ae75aa66d79e5b05e2cea77393b8ad266ff9530bbbd69",
        "coxiter_headers": [0, 0],
        "original_SK_collision_entry": 12,
    },
    11: {
        "id": "R14-UE8A1x4-V4-H10-15",
        "kind": "normal_Klein_four_V4",
        "dimension": 10,
        "facets": 15,
        "subgroup_order": 4,
        "class_size": 1,
        "fixed_determinant": 8,
        "fixed_smith": [1] * 10 + [8],
        "stable_orders": [1, 2],
        "stable_even_orders": [1, 2],
        "maximal_index": 2,
        "maximal_determinant": 2,
        "maximal_smith": [1] * 10 + [2],
        "maximal_even": True,
        "standard_lattice": "U + E8(-1) + A1(-1)",
        "sk_list": "List10",
        "sk_entry": 3,
        "sk_symbol": "II_{10,1}([2^{1}]_{1})",
        "sk_counts": [12, 65, 210, 450, 672, 714, 540, 285, 99, 18],
        "sk_cusps": 1,
        "full_euler": "-1/183936614400",
        "section_euler": "-1/174182400",
        "section_index": 1056,
        "finite_counts": [15, 95, 345, 810, 1302, 1470, 1170, 645, 230, 37],
        "affine_components": 10,
        "cusps": 5,
        "endpoint_histogram": {"0,2": 7, "1,1": 76, "2,0": 147},
        "relation_sha256": "caddc8bb3c5592b38c6a015bd6eeba06b84dd38799ac097d0c38c6bf3dcc1a0d",
        "coxiter_headers": [0, 1],
    },
}


EXPECTED_HIGH_CLASS_TO_REPRESENTATIVE = {
    2: 2,
    3: 3,
    7: 7,
    8: 3,
    11: 11,
    12: 11,
    17: 7,
    20: 11,
    28: 11,
    30: 11,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(matrix[row, col]) for col in range(matrix.cols)] for row in range(matrix.rows)]


def _smith(matrix: sp.Matrix) -> list[int]:
    return [abs(int(value)) for value in smith_normal_form(matrix, domain=ZZ).diagonal()]


def _edge(matrix: sp.Matrix, left: int, right: int, weight: int = 1) -> None:
    matrix[left - 1, right - 1] = matrix[right - 1, left - 1] = weight


def _first_nineteen_gram() -> sp.Matrix:
    gram = -2 * sp.eye(19)
    chain = [1, 2, 3, 5, 6, 7, 8, 10, 11]
    for left, right in zip(chain, chain[1:]):
        _edge(gram, left, right)
    _edge(gram, 3, 4)
    _edge(gram, 8, 9)
    for right in range(12, 16):
        _edge(gram, 11, right)
    for right in range(16, 20):
        _edge(gram, 1, right)
    for left, right in zip(range(12, 16), range(16, 20)):
        _edge(gram, left, right, 2)
    return gram


def _coordinate_model() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    first = _first_nineteen_gram()
    basis_gram = first[:14, :14]
    coordinates = [basis_gram.inv() * first[:14, wall] for wall in range(19)]
    for label in (20, 22, 24, 26):
        even = sp.Matrix(EVEN_EXTRA_COORDINATES[label])
        coordinates.extend((even, 2 * D2 - even))
    wall_coordinates = sp.Matrix.hstack(*coordinates)
    wall_gram = wall_coordinates.T * basis_gram * wall_coordinates
    return basis_gram, wall_coordinates, wall_gram


def _printed_tikz_gram(*, thick_precedence: bool) -> sp.Matrix:
    gram = -2 * sp.eye(27)
    first = _first_nineteen_gram()
    gram[:19, :19] = first
    for right in (20, 22, 24, 26):
        _edge(gram, 4, right, 2)
    for right in (21, 23, 25, 27):
        _edge(gram, 9, right, 2)
    for left, right in ((20, 21), (22, 23), (24, 25), (26, 27)):
        _edge(gram, left, right, 6)
    thin_four = (
        (20, 23), (20, 25), (20, 27),
        (22, 21), (22, 25), (22, 27),
        (24, 23), (24, 25), (24, 27),
        (26, 23), (26, 25), (26, 21), (24, 21),
    )
    for left, right in thin_four:
        if thick_precedence and {left, right} == {24, 25}:
            continue
        _edge(gram, left, right, 4)
    cross = {
        12: (20, 21, 23, 26),
        13: (22, 24, 25, 27),
        14: (20, 21, 22, 27),
        15: (23, 24, 25, 26),
        16: (21, 22, 23, 24),
        17: (20, 25, 26, 27),
        18: (20, 22, 23, 25),
        19: (21, 24, 26, 27),
    }
    for left, rights in cross.items():
        for right in rights:
            _edge(gram, left, right, 2)
    return gram


def _source_audit(
    inventory: dict[str, Any],
    basis_gram: sp.Matrix,
    walls: sp.Matrix,
    wall_gram: sp.Matrix,
) -> dict[str, Any]:
    stock = next(row for row in inventory["stocks"] if row["key"] == "R14_UE8A1x4")
    source_hash = _sha256(SOURCE_TEX)
    declared_hash = inventory["source"]["sha256"]["tex"]
    d2_pairings = list((D2.T * basis_gram * walls))
    d132_pairings = list((D132.T * basis_gram * walls))
    intended_d132 = [1] * 11 + [2, 2, 12, 12, 12, 12, 2, 2] + [28, 18, 28, 18, 18, 28, 18, 28]
    intended_d2 = [
        0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
        1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2,
    ]
    intended_extra = sp.zeros(8)
    for index in range(8):
        intended_extra[index, index] = -2
    for even in range(4):
        for odd in range(4):
            intended_extra[2 * even, 2 * odd + 1] = intended_extra[2 * odd + 1, 2 * even] = (
                6 if even == odd else 4
            )
    literal = _printed_tikz_gram(thick_precedence=False)
    repaired = _printed_tikz_gram(thick_precedence=True)
    literal_differences = [
        [left + 1, right + 1, int(wall_gram[left, right]), int(literal[left, right])]
        for left in range(27)
        for right in range(left)
        if wall_gram[left, right] != literal[left, right]
    ]
    repaired_differences = sum(
        wall_gram[left, right] != repaired[left, right]
        for left in range(27)
        for right in range(left)
    )
    checks = {
        "pinned_source_hash": source_hash == declared_hash,
        "inventory_row": (
            stock["picard_rank"] == 14
            and stock["wall_count"] == 27
            and stock["lattice"] == "U+E8(-1)+A1(-1)^4"
        ),
        "basis_determinant": basis_gram.det() == -16,
        "coordinate_lane_rank_and_inertia": (
            wall_gram.rank() == 14 and exact_inertia(wall_gram) == (1, 13, 13)
        ),
        "all_27_coordinate_classes_are_minus_two": all(
            wall_gram[index, index] == -2 for index in range(27)
        ),
        "all_coordinates_integral": all(value.is_Integer for value in walls),
        "first_19_graph_replayed": wall_gram[:19, :19] == _first_nineteen_gram(),
        "D2_square_and_pairings": (
            (D2.T * basis_gram * D2)[0] == 2
            and [int(value) for value in d2_pairings] == intended_d2
        ),
        "D132_square_and_coordinate_forced_pairings": (
            (D132.T * basis_gram * D132)[0] == 132
            and [int(value) for value in d132_pairings] == intended_d132
        ),
        "four_complementary_pair_equivalences": all(
            walls[:, even - 1] + walls[:, even] == 2 * D2
            for even in (20, 22, 24, 26)
        ),
        "extra_curve_weight_4_6_graph": wall_gram[19:27, 19:27] == intended_extra,
        "literal_cross_TikZ_is_impossible_in_rank_14": (
            literal.rank() == 24 and exact_inertia(literal) == (6, 18, 3)
        ),
        "thick_precedence_does_not_repair_TikZ": (
            repaired.rank() == 22 and exact_inertia(repaired) == (5, 17, 5)
        ),
        "coordinate_and_literal_lanes_differ": len(literal_differences) == 25,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Roulleau coordinate/source audit changed: {checks}")
    return {
        "stock_inventory_record": stock,
        "tex_sha256": source_hash,
        "authoritative_lane": {
            "method": (
                "A1--A19 graph at lines 9836--9943; NS basis statement at line 9944; "
                "D2 and four pair equivalences at lines 9954--9963; explicit A20, "
                "A22, A24, A26 coordinates at lines 9973--9977"
            ),
            "basis_gram": _payload(basis_gram),
            "wall_coordinate_columns": _payload(walls),
            "wall_gram": _payload(wall_gram),
            "rank": wall_gram.rank(),
            "inertia": list(exact_inertia(wall_gram)),
        },
        "source_inconsistency": {
            "printed_cross_TikZ_lines": [10027, 10114],
            "duplicate_A24_A25_segment_lines": (
                "one thick vertical segment plus a repeated thin segment in the same picture"
            ),
            "last_draw_literal_rank": literal.rank(),
            "last_draw_literal_inertia": list(exact_inertia(literal)),
            "thick_precedence_rank": repaired.rank(),
            "thick_precedence_inertia": list(exact_inertia(repaired)),
            "coordinate_vs_last_draw_lower_triangle_difference_count": len(literal_differences),
            "coordinate_vs_last_draw_differences": literal_differences,
            "coordinate_vs_thick_precedence_difference_count": repaired_differences,
            "D132_prose_index_typo": (
                "lines 9949--9952 overlap A18 in the 2 and 12 lists and omit A14; "
                "the displayed D132 coordinates force 12 on A14--A17"
            ),
            "conclusion": (
                "the printed cross-incidence labels cannot be the rho=14 divisor "
                "configuration; the exact coordinate/equivalence lane is authoritative"
            ),
        },
        "checks": checks,
    }


def _action(permutation: Sequence[int], walls: sp.Matrix) -> sp.Matrix:
    action = sp.Matrix.hstack(*(walls[:, permutation[index]] for index in range(14)))
    return action


def _fixed_basis(subgroup: Subgroup, walls: sp.Matrix) -> sp.Matrix:
    equations = sp.Matrix.vstack(
        *(_action(element, walls) - sp.eye(14) for element in subgroup)
    )
    diagonal, _, right = smith_normal_decomp(
        DomainMatrix.from_Matrix(equations).convert_to(POLY_ZZ)
    )
    normal = diagonal.to_Matrix()
    rank = sum(
        normal[index, index] != 0 for index in range(min(normal.rows, normal.cols))
    )
    transform = right.to_Matrix()
    basis = transform[:, rank:]
    if not (
        abs(int(transform.det())) == 1
        and equations * basis == sp.zeros(equations.rows, basis.cols)
    ):
        raise ArithmeticError("fixed Smith basis is not saturated")
    return basis


def _group_recognition(group: tuple[tuple[int, ...], ...], walls: sp.Matrix) -> dict[str, Any]:
    kernel = tuple(element for element in group if all(element[index] == index for index in range(11)))
    induced = {
        tuple(element[index] - 11 for index in range(11, 15)) for element in kernel
    }
    center = tuple(
        element
        for element in group
        if all(_compose(element, other) == _compose(other, element) for other in group)
    )
    central = next(element for element in center if element != tuple(range(27)))
    central_fixed_rank = 14 - (_action(central, walls) - sp.eye(14)).rank()
    checks = {
        **_group_checks(group),
        "order_48": len(group) == 48,
        "pointwise_A1_to_A11_kernel_is_S4": len(kernel) == len(induced) == math.factorial(4),
        "center_is_C2": len(center) == 2 and _order(central) == 2,
        "direct_product_order": len(kernel) * len(center) == len(group),
        "central_abstract_graph_action_has_only_H5": central_fixed_rank == 6,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"weighted group recognition changed: {checks}")
    return {
        "order": len(group),
        "recognition": "S4 x C2",
        "element_order_histogram": {
            str(key): value for key, value in sorted(Counter(_order(row) for row in group).items())
        },
        "wall_permutations_one_based": [
            [value + 1 for value in element] for element in group
        ],
        "S4_kernel_order": len(kernel),
        "central_C2_permutation_one_based": [value + 1 for value in central],
        "central_C2_fixed_rank": central_fixed_rank,
        "checks": checks,
    }


def _group_type(subgroup: Subgroup) -> str:
    histogram = Counter(_order(element) for element in subgroup)
    if len(subgroup) == 1:
        return "trivial"
    if len(subgroup) in (2, 3):
        return f"C{len(subgroup)}"
    if len(subgroup) == 4:
        return "V4" if histogram[2] == 3 else "C4"
    return f"order_{len(subgroup)}_orders_{dict(sorted(histogram.items()))}"


def _subgroup_census(
    group: tuple[tuple[int, ...], ...],
    walls: sp.Matrix,
    wall_gram: sp.Matrix,
) -> tuple[dict[str, Any], dict[int, tuple[Subgroup, tuple[tuple[int, ...], ...]]]]:
    subgroups = _all_subgroups(group)
    classes = _subgroup_conjugacy_classes(group, subgroups)
    representatives: dict[int, tuple[Subgroup, tuple[tuple[int, ...], ...]]] = {}
    class_rows = []
    high_map: dict[int, int] = {}
    representative_orbits: dict[int, tuple[tuple[int, ...], ...]] = {}
    for class_index, members in enumerate(classes, 1):
        subgroup = subgroups[members[0]]
        orbits = _orbits(subgroup, 27)
        pre = _certify_fold(
            wall_gram,
            subgroup,
            orbits,
            run_vinberg=False,
            full_postcheck=False,
        )
        fixed_basis = _fixed_basis(subgroup, walls)
        if fixed_basis.cols != pre["fixed_rank"]:
            raise ArithmeticError("orbit Gram and saturated fixed rank disagree")
        representative = None
        if (
            class_index > 1
            and pre["hyperbolic_dimension"] > 9
            and pre["status"] == "exact_coxeter_candidate"
        ):
            for known, known_orbits in representative_orbits.items():
                if orbits == known_orbits:
                    representative = known
                    break
            if representative is None:
                representative = class_index
                representative_orbits[class_index] = orbits
                representatives[class_index] = (subgroup, orbits)
            high_map[class_index] = representative
        class_rows.append(
            {
                "class_index_one_based": class_index,
                "subgroup_indices_one_based": [member + 1 for member in members],
                "class_size": len(members),
                "representative_subgroup_order": len(subgroup),
                "representative_group_type": _group_type(subgroup),
                "representative_element_order_histogram": {
                    str(key): value
                    for key, value in sorted(Counter(_order(row) for row in subgroup).items())
                },
                "wall_orbits_one_based": [[value + 1 for value in orbit] for orbit in orbits],
                "fixed_rank": fixed_basis.cols,
                "hyperbolic_dimension": pre["hyperbolic_dimension"],
                "screen_status": pre["status"],
                "screen_gates": pre["gates"],
                "collapsed_positive_representative_class": representative,
            }
        )
    checks = {
        "all_98_subgroups": len(subgroups) == 98,
        "all_33_conjugacy_classes": len(classes) == 33,
        "class_sizes_sum_to_98": sum(len(row) for row in classes) == 98,
        "ten_nontrivial_high_classes": len(high_map) == 10,
        "four_distinct_high_orbit_diagrams": sorted(representatives) == sorted(EXPECTED),
        "expected_high_class_collapse": high_map == EXPECTED_HIGH_CLASS_TO_REPRESENTATIVE,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"subgroup denominator changed: {checks}")
    return (
        {
            "subgroup_count": len(subgroups),
            "subgroup_order_histogram": {
                str(key): value for key, value in sorted(Counter(len(row) for row in subgroups).items())
            },
            "conjugacy_class_count": len(classes),
            "classes": class_rows,
            "high_class_to_distinct_diagram_representative": {
                str(key): value for key, value in sorted(high_map.items())
            },
            "checks": checks,
        },
        representatives,
    )


def _strict_facet_witnesses(
    subgroup: Subgroup,
    fixed_basis: sp.Matrix,
    fixed_gram: sp.Matrix,
    roots: sp.Matrix,
    walls: sp.Matrix,
) -> dict[str, Any]:
    invariant_ample = sum(
        (_action(element, walls) * D132 for element in subgroup),
        sp.zeros(14, 1),
    )
    ample_fixed = fixed_basis.gauss_jordan_solve(invariant_ample)[0]
    records = []
    for index in range(roots.cols):
        root = roots[:, index]
        norm = sp.Rational((root.T * fixed_gram * root)[0])
        pairing = sp.Rational((ample_fixed.T * fixed_gram * root)[0])
        witness = ample_fixed - pairing / norm * root
        pairings = [
            sp.Rational(value) for value in list(witness.T * fixed_gram * roots)
        ]
        square = sp.Rational((witness.T * fixed_gram * witness)[0])
        strict = [value for position, value in enumerate(pairings) if position != index]
        records.append(
            {
                "facet_index_one_based": index + 1,
                "witness_square": str(square),
                "minimum_other_wall_pairing": str(min(strict)),
                "target_pairing": str(pairings[index]),
                "strict_timelike_facet_witness": (
                    square > 0 and pairings[index] == 0 and all(value > 0 for value in strict)
                ),
            }
        )
    return {
        "method": "orthogonal projection of the subgroup-averaged ample D132 class",
        "invariant_ample_fixed_coordinates": [str(value) for value in ample_fixed],
        "witnesses": records,
        "all_orbit_inequalities_are_genuine_facets": all(
            row["strict_timelike_facet_witness"] for row in records
        ),
    }


def _vinberg_audit(executable: Path, gram: sp.Matrix, expected: dict[str, Any]) -> dict[str, Any]:
    dimension = expected["dimension"]
    rows = [f"{dimension} {gram.rows}"]
    rows.extend(" ".join(str(int(gram[left, right])) for right in range(gram.cols)) for left in range(gram.rows))
    completed = subprocess.run(
        [str(executable)],
        input="\n".join(rows) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    result = _parse_output(completed.stdout)
    counts = result["finite_special_subset_counts_by_rank"]
    endpoints = result["edge_endpoint_type_histogram"]
    checks = {
        "dimension_and_facets": (
            result["dimension"] == dimension and result["facet_count"] == expected["facets"]
        ),
        "all_finite_special_subsets": (
            counts == expected["finite_counts"]
            and result["all_finite_subset_search_node_count"] == 1 + sum(counts)
        ),
        "all_affine_components_and_cusps": (
            result["connected_affine_component_count"] == expected["affine_components"]
            and result["ideal_vertex_count"] == expected["cusps"]
            and result["bad_parabolic_neighborhood_count"] == 0
        ),
        "every_edge_has_two_endpoints": (
            endpoints == expected["endpoint_histogram"]
            and result["bad_edge_count"] == 0
            and sum(endpoints.values()) == counts[-2]
            and sum(
                (int(kind.split(",")[0]) + int(kind.split(",")[1])) * value
                for kind, value in endpoints.items()
            ) == 2 * counts[-2]
        ),
        "finite_vertex_edge_incidence": (
            2 * endpoints.get("2,0", 0) + endpoints.get("1,1", 0)
            == dimension * counts[-1]
        ),
        "noncompact": result["ideal_vertex_count"] > 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Vinberg audit changed for {expected['id']}: {checks}")
    return result | {
        "elliptic_edge_count": counts[-2],
        "finite_vertex_count": counts[-1],
        "finite_volume": True,
        "compact": False,
        "classification": "finite_volume_noncompact_Coxeter_polytope",
        "checks": checks,
        "all_checks_pass": True,
    }


def _sk_row(expected: dict[str, Any]) -> dict[str, Any]:
    path = SK_ROOT / expected["sk_list"]
    rows = [line.strip().rstrip(",") for line in path.read_text().splitlines() if line.strip().startswith("<")]
    row = rows[expected["sk_entry"] - 1]
    bracket_groups = re.findall(r"\[([^]]*)\]", row)
    row_counts = (
        [int(value) for value in re.findall(r"-?\d+", bracket_groups[-1])]
        if bracket_groups
        else []
    )
    checks = {
        "entry_exists": len(rows) >= expected["sk_entry"],
        "determinant_and_symbol": (
            row.startswith(f'<{abs(expected["maximal_determinant"])}, "{expected["sk_symbol"]}"')
        ),
        "full_special_counts": row_counts == expected["sk_counts"],
        "full_cusp_count": row.endswith(f" {expected['sk_cusps']}>") or row.endswith(f" {expected['sk_cusps']} >"),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"SK lineage row changed for {expected['id']}: {checks}; {row}")
    return {
        "path": str(path),
        "sha256": _sha256(path),
        "entry_one_based": expected["sk_entry"],
        "row": row,
        "genus_symbol": expected["sk_symbol"],
        "full_chamber_facets": expected["sk_counts"][0],
        "full_chamber_finite_vertices": expected["sk_counts"][-1],
        "full_chamber_cusps": expected["sk_cusps"],
        "checks": checks,
    }


def _analyze_diagram(
    class_index: int,
    subgroup: Subgroup,
    orbits: tuple[tuple[int, ...], ...],
    basis_gram: sp.Matrix,
    walls: sp.Matrix,
    wall_gram: sp.Matrix,
    vinberg_executable: Path,
) -> dict[str, Any]:
    expected = EXPECTED[class_index]
    fixed_basis = _fixed_basis(subgroup, walls)
    fixed_gram = fixed_basis.T * basis_gram * fixed_basis
    orbit_sums = sp.zeros(27, len(orbits))
    for column, orbit in enumerate(orbits):
        for row in orbit:
            orbit_sums[row, column] = 1
    roots = fixed_basis.gauss_jordan_solve(walls * orbit_sums)[0]
    orbit_gram = roots.T * fixed_gram * roots
    if orbit_gram != _orbit_gram(wall_gram, orbits):
        raise ArithmeticError("fixed-root and wall-orbit Grams disagree")
    reflections = [_reflection(fixed_gram, roots[:, index]) for index in range(roots.cols)]
    root_types: Counter[tuple[int, int, int]] = Counter()
    reflection_checks = []
    for index, orbit in enumerate(orbits):
        root = roots[:, index]
        norm = int((root.T * fixed_gram * root)[0])
        divisibility = math.gcd(*(abs(int(value)) for value in fixed_gram * root))
        action = reflections[index]
        disjoint = all(wall_gram[left, right] == 0 for left in orbit for right in orbit if left < right)
        root_types[(norm, divisibility, len(orbit))] += 1
        reflection_checks.append(
            {
                "wall_orbit_one_based": [value + 1 for value in orbit],
                "norm": norm,
                "divisibility": divisibility,
                "pairwise_disjoint_curve_orbit": disjoint,
                "reflection_integral": all(value.is_Integer for value in action),
                "reflection_involutive_isometry": (
                    action * action == sp.eye(fixed_gram.rows)
                    and action.T * fixed_gram * action == fixed_gram
                ),
            }
        )
    facet_witnesses = _strict_facet_witnesses(
        subgroup, fixed_basis, fixed_gram, roots, walls
    )
    positive_gram = -orbit_gram
    relations = _relation_matrix(_payload(positive_gram))
    relation_sha = _canonical_sha256(relations)
    solomon = _solomon(_payload(positive_gram))

    stable = _stable_subgroups(fixed_gram, reflections)
    groups = stable.pop("groups_internal")
    left_inverse = stable.pop("left_smith_inverse")
    maximal_positions = stable["maximal_stable_subgroup_positions_zero_based"]
    if len(maximal_positions) != 1:
        raise ArithmeticError("maximal reflection-stable overclosure is not unique")
    position = maximal_positions[0]
    maximal_basis = _overlattice_basis(
        fixed_gram,
        left_inverse,
        stable["stable_subgroups"][position]["cosets_in_smith_coordinates"],
    )
    maximal_gram = sp.simplify(maximal_basis.T * fixed_gram * maximal_basis)
    maximal_actions = [sp.simplify(maximal_basis.inv() * action * maximal_basis) for action in reflections]
    maximal_checks = {
        "index": abs(sp.Rational(1, maximal_basis.det())) == expected["maximal_index"],
        "determinant_square_law": (
            maximal_gram.det() == fixed_gram.det() / expected["maximal_index"] ** 2
        ),
        "expected_determinant": int(maximal_gram.det()) == expected["maximal_determinant"],
        "expected_smith": _smith(maximal_gram) == expected["maximal_smith"],
        "expected_parity": (
            all(maximal_gram[index, index] % 2 == 0 for index in range(maximal_gram.rows))
            == expected["maximal_even"]
        ),
        "all_reflections_integral_on_overclosure": all(
            value.is_Integer for action in maximal_actions for value in action
        ),
        "all_reflections_preserve_overclosure": all(
            action.T * maximal_gram * action == maximal_gram for action in maximal_actions
        ),
    }

    vinberg = _vinberg_audit(vinberg_executable, positive_gram, expected)
    coxiter = _scan_coxiter(relations, expected["dimension"], expected["facets"])
    sk = _sk_row(expected)
    section_euler = Fraction(solomon["rational_euler_characteristic"])
    full_euler = Fraction(expected["full_euler"])
    section_index = None if full_euler == 0 else section_euler / full_euler
    lineage_checks = {
        "relation_hash": relation_sha == expected["relation_sha256"],
        "no_exact_CoxIter_match": not coxiter["ordinary_exact_matches"] and not coxiter["doubling_exact_matches"],
        "CoxIter_header_counts": [
            coxiter["ordinary_header_candidate_count"],
            coxiter["doubling_header_candidate_count"],
        ] == expected["coxiter_headers"],
        "section_euler": str(section_euler) == expected["section_euler"],
        "full_euler": str(full_euler) == expected["full_euler"],
        "index_status": (
            (section_index is None and expected["section_index"] is None)
            or (section_index == expected["section_index"])
        ),
    }
    checks = {
        "dimension_facets": (
            fixed_gram.rows - 1 == expected["dimension"] and roots.cols == expected["facets"]
        ),
        "fixed_lattice": (
            int(fixed_gram.det()) == expected["fixed_determinant"]
            and _smith(fixed_gram) == expected["fixed_smith"]
            and all(fixed_gram[index, index] % 2 == 0 for index in range(fixed_gram.rows))
        ),
        "integral_fixed_roots": all(value.is_Integer for value in roots),
        "all_curve_orbits_disjoint_and_reflections_integral": all(
            row["pairwise_disjoint_curve_orbit"]
            and row["reflection_integral"]
            and row["reflection_involutive_isometry"]
            for row in reflection_checks
        ),
        "all_orbit_walls_genuine_facets": facet_witnesses["all_orbit_inequalities_are_genuine_facets"],
        "stable_orders": stable["stable_integral_subgroup_orders"] == expected["stable_orders"],
        "stable_even_orders": stable["stable_even_subgroup_orders"] == expected["stable_even_orders"],
        "maximal_overclosure": all(maximal_checks.values()),
        "literal_Vinberg_finite_volume": vinberg["all_checks_pass"],
        "lineage": all(lineage_checks.values()) and all(sk["checks"].values()),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"diagram audit changed for {expected['id']}: {checks}")

    original_collision = None
    if "original_SK_collision_entry" in expected:
        path = SK_ROOT / "List11"
        rows = [line.strip().rstrip(",") for line in path.read_text().splitlines() if line.strip().startswith("<")]
        row = rows[expected["original_SK_collision_entry"] - 1]
        collision_bracket_groups = re.findall(r"\[([^]]*)\]", row)
        collision_counts = (
            [int(value) for value in re.findall(r"-?\d+", collision_bracket_groups[-1])]
            if collision_bracket_groups
            else []
        )
        original_collision = {
            "entry_one_based": expected["original_SK_collision_entry"],
            "row": row,
            "finite_special_counts_and_cusps_match_exactly": (
                collision_counts == expected["finite_counts"]
                and (row.endswith(" 7>") or row.endswith(" 7 >"))
            ),
            "interpretation": (
                "exact full-census signature collision on the original determinant-12 fixed lattice; "
                "an explicit chamber conjugacy is not asserted"
            ),
        }
        if not original_collision["finite_special_counts_and_cusps_match_exactly"]:
            raise ArithmeticError("original determinant-12 SK collision changed")

    return {
        "id": expected["id"],
        "representative_subgroup_class_one_based": class_index,
        "subgroup": {
            "kind": expected["kind"],
            "order": len(subgroup),
            "elements_one_based": [
                [value + 1 for value in element] for element in sorted(subgroup)
            ],
            "wall_orbits_one_based": [[value + 1 for value in orbit] for orbit in orbits],
        },
        "fixed_lattice": {
            "saturated_basis_columns_in_NS_basis": _payload(fixed_basis),
            "gram": _payload(fixed_gram),
            "rank": fixed_gram.rows,
            "determinant": int(fixed_gram.det()),
            "smith_invariants": _smith(fixed_gram),
            "even": True,
            "discriminant_q_mod_2_histogram": _q_histogram(fixed_gram),
        },
        "folded_walls": {
            "root_coordinate_columns": _payload(roots),
            "negative_gram": _payload(orbit_gram),
            "positive_gram": _payload(positive_gram),
            "relation_matrix": relations,
            "relation_matrix_sha256": relation_sha,
            "root_type_norm_divisibility_orbit_size_histogram": {
                f"{norm},{divisibility},{size}": count
                for (norm, divisibility, size), count in sorted(root_types.items())
            },
            "reflection_checks": reflection_checks,
            "strict_facet_witnesses": facet_witnesses,
        },
        "geometry": {
            "hyperbolic_dimension": expected["dimension"],
            "facet_count": expected["facets"],
            "vinberg": vinberg,
        },
        "stable_overlattices": {
            "exhaustive_discriminant_module_audit": stable,
            "unique_maximal_reflection_stable_overclosure": {
                "index": expected["maximal_index"],
                "basis_columns_in_fixed_lattice": _payload(maximal_basis),
                "gram": _payload(maximal_gram),
                "determinant": int(maximal_gram.det()),
                "smith_invariants": _smith(maximal_gram),
                "even": expected["maximal_even"],
                "standard_lattice": expected["standard_lattice"],
                "discriminant_q_mod_2_histogram": _q_histogram(maximal_gram),
                "checks": maximal_checks,
            },
        },
        "lineage": {
            "Scharlau_Kirschmer": sk,
            "CoxIter": coxiter,
            "solomon": solomon,
            "section_to_full_reflection_group_index": (
                int(section_index) if section_index is not None else None
            ),
            "index_status": (
                "exact Solomon characteristic ratio"
                if section_index is not None
                else "finite; numerical value open because both H11 characteristics vanish"
            ),
            "new_arithmetic_commensurability_class": False,
            "original_fixed_lattice_SK_collision": original_collision,
            "bounded_literal_status": (
                "absent from pinned current CoxIter and index-two-doubling snapshots; "
                "worldwide priority open"
            ),
            "checks": lineage_checks,
        },
        "surface_geometry": {
            "ambient_object": "complete nef cone of a K3 surface with NS=U+E8(-1)+A1(-1)^4",
            "published_general_surface_automorphism_group": "C2",
            "published_action_tag": "t",
            "meaning_of_t": "the hyperelliptic involution acts trivially on NS",
            "selected_nonidentity_graph_subgroup_realized_on_general_surface": False,
            "special_period_realization": "not assessed",
            "surface_automorphism_not_needed_for_polytope_proof": True,
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def build_certificate() -> dict[str, Any]:
    inventory = json.loads(SOURCE_INVENTORY.read_text())
    if not inventory["all_checks_pass"]:
        raise ArithmeticError("Roulleau source inventory failed")
    basis_gram, walls, wall_gram = _coordinate_model()
    source = _source_audit(inventory, basis_gram, walls, wall_gram)
    group, search_nodes = _weighted_graph_automorphisms(wall_gram)
    group_record = _group_recognition(group, walls)
    group_record["exhaustive_backtracking_node_count"] = search_nodes
    subgroup_census, representatives = _subgroup_census(group, walls, wall_gram)

    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="roulleau-r14-ue8a1x4-vinberg-") as raw:
        executable = Path(raw) / "vinberg"
        subprocess.run(
            [compiler, "-O3", "-std=c++17", "-Wall", "-Wextra", "-Werror", str(CPP_SOURCE), "-o", str(executable)],
            check=True,
        )
        cases = [
            _analyze_diagram(
                class_index,
                representatives[class_index][0],
                representatives[class_index][1],
                basis_gram,
                walls,
                wall_gram,
                executable,
            )
            for class_index in sorted(representatives)
        ]
    checks = {
        "source_coordinate_lane": all(source["checks"].values()),
        "weighted_group": all(group_record["checks"].values()),
        "complete_subgroup_denominator": all(subgroup_census["checks"].values()),
        "four_positive_diagrams": len(cases) == 4 and all(case["all_checks_pass"] for case in cases),
        "dimensions_and_facets": [
            [case["geometry"]["hyperbolic_dimension"], case["geometry"]["facet_count"]]
            for case in cases
        ] == [[12, 23], [11, 19], [11, 19], [10, 15]],
        "no_new_arithmetic_classes": all(
            not case["lineage"]["new_arithmetic_commensurability_class"] for case in cases
        ),
        "no_general_surface_automorphism_claim": all(
            not case["surface_geometry"]["selected_nonidentity_graph_subgroup_realized_on_general_surface"]
            for case in cases
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global R14/27 certificate changed: {checks}")
    return {
        "schema": SCHEMA,
        "artifact_id": "roulleau_r14_ue8a1x4_all_subgroup_folds",
        "source": {
            "inventory_path": str(SOURCE_INVENTORY.relative_to(ROOT)),
            "inventory_sha256": _sha256(SOURCE_INVENTORY),
            "tex_path": str(SOURCE_TEX),
            "tex_sha256": _sha256(SOURCE_TEX),
            "cpp_path": str(CPP_SOURCE.relative_to(ROOT)),
            "cpp_sha256": _sha256(CPP_SOURCE),
            "compiler": Path(compiler).name,
        },
        "source_reconstruction": source,
        "weighted_wall_group": group_record,
        "subgroup_census": subgroup_census,
        "cases": cases,
        "summary": {
            "source_wall_count": 27,
            "weighted_group_order": 48,
            "subgroup_count": 98,
            "subgroup_conjugacy_class_count": 33,
            "nontrivial_above_H9_class_count": 10,
            "distinct_positive_diagram_count": 4,
            "finite_volume_noncompact_count": 4,
            "dimensions_and_facets": [[12, 23], [11, 19], [11, 19], [10, 15]],
            "exact_even_dimensional_indices": {"H12_23": 8320, "H10_15": 1056},
            "H11_numerical_indices_open": 2,
            "new_arithmetic_commensurability_class_count": 0,
            "exact_CoxIter_match_count": 0,
            "worldwide_literal_priority": "open",
            "general_surface_geometric_nonidentity_action_count": 0,
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
            raise ArithmeticError("rebuilt R14 UE8A1^4 archive changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
