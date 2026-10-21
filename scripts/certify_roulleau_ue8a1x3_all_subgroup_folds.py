#!/usr/bin/env python3
"""Certify every subgroup fold of Roulleau's 17-wall rank-13 K3 stock.

The source Neron--Severi lattice is ``U+E8(-1)+A1(-1)^3``.  Its complete
17-curve wall graph has weighted automorphism group ``S3``.  This certificate
checks the published coordinates, exhausts subgroup conjugacy classes, applies
Vinberg's finite-volume criterion to each distinct proper fixed cone, and then
identifies the full integral reflection group by independent all-root Vinberg
enumerations.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import reduce
import hashlib
import itertools
import json
from math import gcd
from pathlib import Path
from typing import Any, Sequence

import sympy as sp

from coxeter_polytopes.search.qf2_h11_vinberg import (
    lattice_gram as qf2_lattice_gram,
    pair as qf2_pair,
    run_search as run_qf2_vinberg,
)
from coxeter_polytopes.search.u_e8_m6_vinberg import (
    lattice_gram as ue8m6_lattice_gram,
    pair as ue8m6_pair,
    run_search as run_ue8m6_vinberg,
)
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _overlattice_basis,
    _payload,
    _q_histogram,
    _reflection,
    _smith,
    _stable_subgroups,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_direct,
    _euler_from_connected_components,
    _relation_matrix,
)
from scripts.certify_roulleau_non2elementary_involution_folds import (
    _exact_gram_permutation,
    _scan_coxiter,
)
from scripts.certify_roulleau_ue8d4_all_subgroup_folds import (
    _fixed_data,
    _lift_actions,
    _subgroup_classes,
)
from scripts.search_k3_high_rank_involution_folds import (
    _coxeter_audit,
    _geometry_certificate,
    _group_checks,
    _matrix_payload,
    _weighted_graph_automorphisms,
)
from scripts.search_k3_nef_subgroup_folds import _orbit_gram, _orbits, _order
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_ue8a1x3_all_subgroup_folds.json"
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
SOURCE_TEX_SHA256 = "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
SOURCE_INVENTORY = ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
SK_ROOT = Path("/private/tmp/Scharlau-Kirschmer-Lorentz")
SCHEMA = "coxeter.roulleau_ue8a1x3_all_subgroup_folds.v1"
RANK = 13
WALL_COUNT = 17


CHAIN = (1, 2, 4, 5, 6, 7, 8, 9, 10, 12)
SOURCE_EDGES = tuple((*edge, 1) for edge in zip(CHAIN, CHAIN[1:])) + (
    (12, 15, 2),
    (15, 17, 2),
    (17, 1, 1),
    (4, 3, 1),
    (10, 11, 1),
    (11, 14, 2),
    (14, 17, 2),
    (10, 13, 1),
    (13, 16, 2),
    (16, 17, 2),
)


EXPECTED_TYPES: dict[str, dict[str, Any]] = {
    "R13_UE8A1x3_C2": {
        "class_indices": [2],
        "dimension": 11,
        "facets": 15,
        "edges": 341,
        "finite_vertices": 54,
        "cusps": 3,
        "parabolic_audits": 7,
        "relations": {"0": 3, "1": 1, "2": 89, "3": 11, "4": 1},
        "finite_special_counts": [1, 15, 101, 406, 1089, 2058, 2814, 2808, 2030, 1033, 341, 54],
        "euler": "0",
        "fixed_determinant": -8,
        "fixed_smith": [1] * 10 + [2, 4],
        "stable_orders": [1, 2],
        "stable_even_orders": [1],
        "maximal_index": 2,
        "maximal_determinant": -2,
        "maximal_parity": "odd",
        "primitive_contents": [1] * 11 + [2, 1, 2, 1],
        "intrinsic_name": "<2>+<-1>^11",
        "intrinsic_bound": 36,
        "intrinsic_height": "36",
        "sk_list": "List11",
        "sk_entry": 2,
        "sk_symbol": "I_{11,1}([1^{11}2^{1}]_{2})",
        "coxiter_matches": [
            "graphs/group units qf/11-2,1,1,1,1,1,1,1,1,1,1,1.coxiter",
        ],
        "doubling_matches": [
            "testing_index2/11-2,1,1,1,1,1,1,1,1,1,1,1-11-15_0.coxiter",
            "testing_index2/11-2,1,1,1,1,1,1,1,1,1,1,1-15-11_0.coxiter",
        ],
    },
    "R13_UE8A1x3_C3_S3": {
        "class_indices": [3, 4],
        "dimension": 10,
        "facets": 13,
        "edges": 176,
        "finite_vertices": 31,
        "cusps": 2,
        "parabolic_audits": 4,
        "relations": {"0": 1, "1": 1, "2": 65, "3": 10, "6": 1},
        "finite_special_counts": [1, 13, 76, 264, 606, 966, 1092, 876, 488, 176, 31],
        "euler": "-61/91968307200",
        "fixed_determinant": 6,
        "fixed_smith": [1] * 10 + [6],
        "stable_orders": [1],
        "stable_even_orders": [1],
        "maximal_index": 1,
        "maximal_determinant": 6,
        "maximal_parity": "even",
        "primitive_contents": [1] * 13,
        "intrinsic_name": "U+E8(-1)+<-6>",
        "intrinsic_bound": 8,
        "intrinsic_height": "8",
        "sk_list": "List10",
        "sk_entry": 9,
        "sk_symbol": "II_{10,1}([2^{1}]_{7} 3^{-1})",
        "coxiter_matches": [],
        "doubling_matches": [],
    },
}


SURFACE_STATUS = (
    "not_realized_on_the_general_K3: Roulleau records Aut(X)=C2 and tag t; "
    "the nontrivial surface involution acts trivially on NS, whereas every "
    "nontrivial S3 diagram action used here acts nontrivially on NS"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source_gram() -> sp.Matrix:
    gram = -2 * sp.eye(WALL_COUNT)
    for left, right, product in SOURCE_EDGES:
        gram[left - 1, right - 1] = gram[right - 1, left - 1] = product
    return gram


def _ambient_reconstruction(gram: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix]:
    basis_gram = gram[:RANK, :RANK]
    wall_rows = gram[:, :RANK] * basis_gram.inv()
    expected_extra = [
        [2, 4, 3, 6, 5, 4, 3, 2, 1, 0, -1, 0, 0],
        [2, 4, 3, 6, 5, 4, 3, 2, 1, 0, 0, -1, 0],
        [2, 4, 3, 6, 5, 4, 3, 2, 1, 0, 0, 0, -1],
        [-2, -3, -2, -4, -3, -2, -1, 0, 1, 2, 1, 1, 1],
    ]
    ample = sp.Matrix([2, 5, 4, 9, 10, 12, 15, 19, 24, 30, 14, 14, 9])
    nef_two = sp.Matrix([0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1])
    ample_intersections = [int(value) for value in ample.T * basis_gram * wall_rows.T]
    nef_intersections = [int(value) for value in nef_two.T * basis_gram * wall_rows.T]
    checks = {
        "basis_determinant": int(basis_gram.det()) == 8,
        "basis_signature": exact_inertia(basis_gram) == (1, 12, 0),
        "first_thirteen_are_basis": wall_rows[:RANK, :] == sp.eye(RANK),
        "four_dependent_coordinates": [
            [int(value) for value in wall_rows.row(index)]
            for index in range(RANK, WALL_COUNT)
        ] == expected_extra,
        "all_wall_coordinates_integral": all(value.is_Integer for value in wall_rows),
        "wall_gram_reconstructed": wall_rows * basis_gram * wall_rows.T == gram,
        "published_ample_square_294": int((ample.T * basis_gram * ample)[0]) == 294,
        "published_ample_intersections": ample_intersections
        == [1] * 10 + [2, 2, 12, 28, 28, 18, 2],
        "published_nef_square_2": int((nef_two.T * basis_gram * nef_two)[0]) == 2,
        "published_nef_zero_set": [
            index + 1 for index, value in enumerate(nef_intersections) if value == 0
        ] == [2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 17],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"source-coordinate reconstruction failed: {checks}")
    return basis_gram, wall_rows


def _integral_overlattice_census(gram: sp.Matrix) -> dict[str, Any]:
    """Exhaust every integral overlattice via the finite discriminant group."""

    moduli, _, discriminant_gram = _discriminant_data(gram)
    elements = list(itertools.product(*(range(modulus) for modulus in moduli)))
    if len(elements) > 16:
        raise ArithmeticError("discriminant group too large for literal subgroup census")
    index_of = {element: index for index, element in enumerate(elements)}

    def add(left: int, right: int) -> int:
        return index_of[
            tuple(
                (elements[left][position] + elements[right][position])
                % moduli[position]
                for position in range(len(moduli))
            )
        ]

    good = []
    for mask in range(1 << len(elements)):
        subgroup = frozenset(
            index for index in range(len(elements)) if (mask >> index) & 1
        )
        if 0 not in subgroup or any(
            add(left, right) not in subgroup
            for left in subgroup
            for right in subgroup
        ):
            continue
        if any(
            (sp.Matrix(elements[left]).T * discriminant_gram * sp.Matrix(elements[right]))[0].q
            != 1
            for left in subgroup
            for right in subgroup
        ):
            continue
        good.append(subgroup)
    good.sort(key=lambda subgroup: (len(subgroup), sorted(subgroup)))
    maximal = [
        position
        for position, subgroup in enumerate(good)
        if not any(subgroup < other for other in good)
    ]
    return {
        "discriminant_smith_moduli": moduli,
        "discriminant_group_order": len(elements),
        "integral_overlattice_subgroup_count": len(good),
        "integral_overlattice_indices": [len(subgroup) for subgroup in good],
        "subgroups_in_smith_coordinates": [
            [list(elements[index]) for index in sorted(subgroup)] for subgroup in good
        ],
        "maximal_positions_zero_based": maximal,
        "unique_maximal_integral_overlattice": len(maximal) == 1,
    }


def _primitive_columns(matrix: sp.Matrix) -> tuple[sp.Matrix, list[int]]:
    roots = []
    contents = []
    for column in range(matrix.cols):
        root = matrix[:, column]
        content = reduce(gcd, (abs(int(value)) for value in root), 0)
        if content == 0:
            raise ArithmeticError("zero wall normal")
        roots.append(root / content)
        contents.append(content)
    return sp.Matrix.hstack(*roots), contents


def _intrinsic_search(identifier: str) -> tuple[Any, sp.Matrix, Any]:
    expected = EXPECTED_TYPES[identifier]
    if identifier == "R13_UE8A1x3_C2":
        return (
            run_qf2_vinberg(expected["intrinsic_bound"]),
            sp.Matrix(qf2_lattice_gram()),
            qf2_pair,
        )
    return (
        run_ue8m6_vinberg(expected["intrinsic_bound"]),
        sp.Matrix(ue8m6_lattice_gram()),
        ue8m6_pair,
    )


def _full_reflection_lineage(
    identifier: str, fixed_gram: sp.Matrix, roots: Sequence[sp.Matrix]
) -> dict[str, Any]:
    expected = EXPECTED_TYPES[identifier]
    reflections = [_reflection(fixed_gram, root) for root in roots]
    stable = _stable_subgroups(fixed_gram, reflections)
    groups = stable.pop("groups_internal")
    left_inverse = stable.pop("left_smith_inverse")
    if not (
        stable["stable_integral_subgroup_orders"] == expected["stable_orders"]
        and stable["stable_even_subgroup_orders"] == expected["stable_even_orders"]
        and len(stable["maximal_stable_subgroup_positions_zero_based"]) == 1
    ):
        raise ArithmeticError("reflection-stable overlattice census changed")
    position = stable["maximal_stable_subgroup_positions_zero_based"][0]
    maximal_basis = _overlattice_basis(
        fixed_gram,
        left_inverse,
        stable["stable_subgroups"][position]["cosets_in_smith_coordinates"],
    )
    maximal_gram = sp.simplify(maximal_basis.T * fixed_gram * maximal_basis)
    maximal_roots_raw = sp.Matrix.hstack(
        *(maximal_basis.inv() * root for root in roots)
    )
    maximal_roots, contents = _primitive_columns(maximal_roots_raw)
    maximal_reflections = [
        _reflection(maximal_gram, maximal_roots[:, column])
        for column in range(maximal_roots.cols)
    ]
    parity = (
        "even"
        if all(maximal_gram[index, index] % 2 == 0 for index in range(maximal_gram.rows))
        else "odd"
    )
    overclosures = _integral_overlattice_census(fixed_gram)

    search, intrinsic_gram, intrinsic_pair = _intrinsic_search(identifier)
    candidate_gram = [
        [
            int((maximal_roots[:, left].T * maximal_gram * maximal_roots[:, right])[0])
            for right in range(maximal_roots.cols)
        ]
        for left in range(maximal_roots.cols)
    ]
    intrinsic_root_gram = [
        [intrinsic_pair(left, right) for right in search.roots]
        for left in search.roots
    ]
    permutation = _exact_gram_permutation(candidate_gram, intrinsic_root_gram)
    if permutation is None:
        raise ArithmeticError("candidate primitive roots do not match intrinsic Vinberg roots")
    independent = []
    for column in range(maximal_roots.cols):
        if maximal_roots[:, [*independent, column]].rank() > len(independent):
            independent.append(column)
    if len(independent) != maximal_gram.rows:
        raise ArithmeticError("primitive roots do not span the maximal lattice")
    ordered_intrinsic = sp.Matrix.hstack(
        *(sp.Matrix(search.roots[permutation[column]]) for column in range(len(permutation)))
    )
    isometry = ordered_intrinsic[:, independent] * maximal_roots[:, independent].inv()
    intrinsic_checks = {
        "search_stopped_at_finite_volume": search.stopping.fired,
        "root_count": len(search.roots) == expected["facets"],
        "stopping_height": str(search.levels[-1].height) == expected["intrinsic_height"],
        "finite_vertices": search.stopping.finite_vertex_count
        == expected["finite_vertices"],
        "cusps": search.stopping.cusp_count == expected["cusps"],
        "no_hyperideal_rays": search.stopping.hyperideal_ray_count == 0,
        "no_opposite_sheet_rays": search.stopping.opposite_sheet_ray_count == 0,
        "isometry_integral_unimodular": all(value.is_Integer for value in isometry)
        and abs(int(isometry.det())) == 1,
        "isometry_preserves_gram": isometry.T * intrinsic_gram * isometry
        == maximal_gram,
        "wall_by_wall_vector_conjugacy": isometry * maximal_roots == ordered_intrinsic,
    }
    maximal_checks = {
        "index": len(groups[position]) == expected["maximal_index"],
        "determinant": int(maximal_gram.det()) == expected["maximal_determinant"],
        "parity": parity == expected["maximal_parity"],
        "primitive_contents": contents == expected["primitive_contents"],
        "reflections_integral": all(
            value.is_Integer for reflection in maximal_reflections for value in reflection
        ),
        "reflections_preserve_gram": all(
            reflection.T * maximal_gram * reflection == maximal_gram
            for reflection in maximal_reflections
        ),
        "unique_maximal_integral_overlattice": overclosures[
            "unique_maximal_integral_overlattice"
        ],
        "all_intrinsic_checks": all(intrinsic_checks.values()),
    }
    if not all(maximal_checks.values()):
        raise ArithmeticError(f"full reflection lineage failed: {maximal_checks}")
    return {
        "reflection_stable_overlattice_census": stable,
        "all_integral_overlattice_census": overclosures,
        "unique_maximal_lattice": {
            "index_over_fixed_lattice": len(groups[position]),
            "basis_columns_in_fixed_lattice": _payload(maximal_basis),
            "gram": _payload(maximal_gram),
            "determinant": int(maximal_gram.det()),
            "smith_invariants": _smith(maximal_gram),
            "parity": parity,
            "discriminant_q_histogram_mod_two": (
                _q_histogram(maximal_gram) if parity == "even" else None
            ),
            "primitive_facet_roots": _matrix_payload(maximal_roots),
            "raw_root_contents": contents,
        },
        "intrinsic_all_root_vinberg": {
            "lattice": expected["intrinsic_name"],
            "lattice_gram": _matrix_payload(intrinsic_gram),
            "height_bound": str(search.height_bound),
            "stopping_height": str(search.levels[-1].height),
            "root_count": len(search.roots),
            "finite_vertex_count": search.stopping.finite_vertex_count,
            "cusp_count": search.stopping.cusp_count,
            "levels": [
                {
                    "height": str(level.height),
                    "keys": [list(key) for key in level.keys],
                    "candidate_count": level.candidate_count,
                    "accepted_roots": [list(root) for root in level.accepted_roots],
                }
                for level in search.levels
            ],
            "candidate_to_intrinsic_facet_permutation_zero_based": permutation,
            "maximal_to_intrinsic_isometry": _matrix_payload(isometry),
            "isometry_determinant": int(isometry.det()),
            "checks": intrinsic_checks,
        },
        "full_reflection_group": {
            "status": "exact_intrinsic_Vinberg_full_chamber_conjugacy",
            "reflection_subgroup_index": 1,
            "fixed_lattice_scope_reason": (
                "the maximal integral overclosure is unique, hence characteristic "
                "under every fixed-lattice isometry; the folded generators preserve "
                "the fixed lattice and generate the full reflection group of that "
                "overclosure"
                if expected["maximal_index"] > 1
                else "the fixed lattice has no proper integral overclosure and the "
                "folded walls are the intrinsic all-root Vinberg chamber"
            ),
        },
        "checks": maximal_checks,
    }


def build_certificate() -> dict[str, Any]:
    if _sha256(SOURCE_TEX) != SOURCE_TEX_SHA256:
        raise ArithmeticError("pinned Roulleau TeX hash changed")
    inventory = json.loads(SOURCE_INVENTORY.read_text())
    stock = next(row for row in inventory["stocks"] if row["key"] == "R13_UE8A1x3")
    source_lines = SOURCE_TEX.read_text().splitlines()
    excerpt = "\n".join(source_lines[9128:9249])
    source_checks = {
        "inventory_certified": inventory["all_checks_pass"] is True,
        "stock_rank": stock["picard_rank"] == RANK,
        "stock_walls": stock["wall_count"] == WALL_COUNT,
        "stock_discriminant": stock["absolute_discriminant"] == 8,
        "section_token": r"U\oplus\mathbf{E}_{8}\oplus\mathbf{A}_{1}^{\oplus3}" in excerpt,
        "wall_count_token": "contains $17$ $(-2)$-curves" in excerpt,
        "basis_token": "$A_{1},\\dots,A_{13}$ generate" in excerpt,
        "ample_token": "D_{294}=(2,5,4,9,10,12,15,19,24,30,14,14,9)" in excerpt,
    }
    if not all(source_checks.values()):
        raise ArithmeticError(f"source checks failed: {source_checks}")

    wall_gram = _source_gram()
    basis_gram, wall_rows = _ambient_reconstruction(wall_gram)
    group, search_nodes = _weighted_graph_automorphisms(wall_gram)
    if not (
        len(group) == 6
        and Counter(_order(element) for element in group) == {1: 1, 2: 3, 3: 2}
        and all(_group_checks(group).values())
    ):
        raise ArithmeticError("unexpected weighted-diagram group")
    actions = _lift_actions(
        group, wall_rows, basis_gram, rank=RANK, wall_count=WALL_COUNT
    )
    classes_internal = _subgroup_classes(group, wall_count=WALL_COUNT)
    if [row["order"] for row in classes_internal] != [1, 2, 3, 6]:
        raise ArithmeticError("S3 subgroup conjugacy denominator changed")

    class_records = []
    proper_internal = []
    for row in classes_internal:
        subgroup = row.pop("representative_subgroup")
        generators = row.pop("representative_generators")
        orbits = _orbits(subgroup, WALL_COUNT)
        fixed_rank = _orbit_gram(wall_gram, orbits).rank()
        public = {
            **row,
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "wall_orbits_one_based": [[value + 1 for value in orbit] for orbit in orbits],
            "representative_generators_zero_based": [list(generator) for generator in generators],
        }
        class_records.append(public)
        if len(subgroup) > 1:
            data = _fixed_data(
                generators,
                subgroup,
                actions,
                basis_gram,
                wall_rows,
                wall_gram,
                rank=RANK,
                wall_count=WALL_COUNT,
                surface_realization_status=SURFACE_STATUS,
            )
            proper_internal.append({**public, "_data": data})

    types = []
    for identifier in ("R13_UE8A1x3_C2", "R13_UE8A1x3_C3_S3"):
        expected = EXPECTED_TYPES[identifier]
        members = [
            row for row in proper_internal
            if row["class_index"] in expected["class_indices"]
        ]
        representative = members[0]
        data = representative["_data"]
        orbit_gram = data["orbit_gram"]
        same_cone_checks = {
            "same_wall_orbits": all(
                member["_data"]["cycles_one"] == data["cycles_one"] for member in members
            ),
            "same_orbit_gram": all(
                member["_data"]["orbit_gram"] == orbit_gram for member in members
            ),
            "all_action_checks": all(
                member["_data"]["ns_action"]["all_integral_lattice_checks_pass"]
                for member in members
            ),
        }
        coxeter = _coxeter_audit(
            orbit_gram, list(range(orbit_gram.rows)), data["cycles_one"]
        )
        geometry = _geometry_certificate(
            RANK,
            data["cycles_one"],
            orbit_gram,
            coxeter,
            data["reflection_audit"],
            data["ns_action"],
        )
        relations = _relation_matrix(
            [
                [str(-orbit_gram[left, right]) for right in range(orbit_gram.rows)]
                for left in range(orbit_gram.rows)
            ]
        )
        relation_histogram = dict(
            sorted(
                Counter(
                    relations[left][right]
                    for left in range(len(relations))
                    for right in range(left + 1, len(relations))
                ).items()
            )
        )
        solomon = _euler_from_connected_components(
            relations, _connected_masks_direct(relations)
        )
        coxiter = _scan_coxiter(
            relations, expected["dimension"], expected["facets"]
        )
        ordinary_paths = sorted(row["path"] for row in coxiter["ordinary_exact_matches"])
        doubling_paths = sorted(row["path"] for row in coxiter["doubling_exact_matches"])
        lineage = _full_reflection_lineage(
            identifier, data["fixed_gram"], data["roots"]
        )
        sk_path = SK_ROOT / expected["sk_list"]
        sk_line = sk_path.read_text().splitlines()[expected["sk_entry"]]
        checks = {
            "class_partition": [row["class_index"] for row in members]
            == expected["class_indices"],
            "same_fixed_cone_for_merged_classes": all(same_cone_checks.values()),
            "fixed_lattice": int(data["fixed_gram"].det())
            == expected["fixed_determinant"]
            and _smith(data["fixed_gram"]) == expected["fixed_smith"]
            and all(
                data["fixed_gram"][index, index] % 2 == 0
                for index in range(data["fixed_gram"].rows)
            ),
            "geometry": geometry["all_checks_pass"],
            "vinberg_volume_counts": geometry["hyperbolic_dimension"]
            == expected["dimension"]
            and geometry["facet_count"] == expected["facets"]
            and geometry["vinberg_volume_criterion"]["elliptic_edge_count"]
            == expected["edges"]
            and geometry["finite_vertex_count"] == expected["finite_vertices"]
            and geometry["ideal_vertex_count"] == expected["cusps"]
            and geometry["vinberg_volume_criterion"][
                "parabolic_neighborhood_audit_count"
            ] == expected["parabolic_audits"],
            "relations": relation_histogram == expected["relations"],
            "solomon": solomon["finite_special_subset_count_by_rank"]
            == expected["finite_special_counts"]
            and solomon["rational_euler_characteristic"] == expected["euler"],
            "coxiter_matches": ordinary_paths == sorted(expected["coxiter_matches"])
            and doubling_paths == sorted(expected["doubling_matches"]),
            "scharlau_kirschmer_row": expected["sk_symbol"] in sk_line,
            "full_reflection_lineage": all(lineage["checks"].values())
            and lineage["full_reflection_group"]["reflection_subgroup_index"] == 1,
            "surface_nonrealization_separate": data["ns_action"][
                "surface_realization_status"
            ] == SURFACE_STATUS,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"type checks failed for {identifier}: {checks}")
        types.append(
            {
                "type_id": identifier,
                "subgroup_class_indices_one_based": expected["class_indices"],
                "subgroup_structures": [
                    row["element_order_histogram"] for row in members
                ],
                "fixed_cone_identification_checks": same_cone_checks,
                "fixed_lattice": {
                    "basis_rows_in_A1_to_A13_basis": _matrix_payload(data["fixed_basis"]),
                    "gram": _matrix_payload(data["fixed_gram"]),
                    "rank": data["fixed_gram"].rows,
                    "determinant": int(data["fixed_gram"].det()),
                    "smith_invariants": _smith(data["fixed_gram"]),
                    "even": True,
                    "smith_kernel_audit": data["fixed_kernel"],
                    "orbit_root_coordinate_columns": _matrix_payload(data["root_matrix"]),
                },
                "reflection_audit": data["reflection_audit"],
                "neron_severi_diagram_actions": [
                    {
                        "class_index_one_based": member["class_index"],
                        "action": member["_data"]["ns_action"],
                    }
                    for member in members
                ],
                "geometry": geometry,
                "normalized_relation_matrix": relations,
                "relation_histogram": relation_histogram,
                "solomon": solomon,
                "maximal_reflection_lineage": {
                    **lineage,
                    "scharlau_kirschmer": {
                        "list": expected["sk_list"],
                        "entry_one_based": expected["sk_entry"],
                        "genus_symbol": expected["sk_symbol"],
                        "list_row": sk_line.strip(),
                        "list_sha256": _sha256(sk_path),
                    },
                    "coxiter": coxiter,
                },
                "status_separation": {
                    "source_coordinate_validity": "proved_exactly",
                    "geometric_surface_action": "not_realized_on_the_general_K3",
                    "abstract_hyperbolic_chamber": "proved_finite_volume_noncompact",
                    "full_integral_reflection_group": "proved_exact_index_one",
                    "worldwide_abstract_priority": "known_not_new",
                    "algebro_geometric_realization_priority": "open_literature_check",
                },
                "checks": checks,
                "all_checks_pass": all(checks.values()),
            }
        )

    checks = {
        "source_metadata": all(source_checks.values()),
        "wall_gram_rank_signature": wall_gram.rank() == RANK
        and exact_inertia(wall_gram) == (1, 12, 4),
        "weighted_diagram_group_S3": len(group) == 6,
        "four_subgroup_classes": len(class_records) == 4,
        "three_proper_classes": len(proper_internal) == 3,
        "two_distinct_finite_volume_types": len(types) == 2
        and all(row["all_checks_pass"] for row in types),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global checks failed: {checks}")
    return {
        "schema": SCHEMA,
        "artifact_id": "roulleau_ue8a1x3_all_subgroup_folds",
        "source": {
            "stock_key": "R13_UE8A1x3",
            "lattice": "U+E8(-1)+A1(-1)^3",
            "tex_section_lines": [9129, 9249],
            "tex_sha256": SOURCE_TEX_SHA256,
            "source_edges_one_based": [list(edge) for edge in SOURCE_EDGES],
            "wall_gram": _matrix_payload(wall_gram),
            "basis_gram": _matrix_payload(basis_gram),
            "wall_coordinate_rows": _matrix_payload(wall_rows),
            "checks": source_checks,
        },
        "weighted_diagram_group": {
            "structure": "S3",
            "order": len(group),
            "element_order_histogram": {"1": 1, "2": 3, "3": 2},
            "exhaustive_backtracking_node_count": search_nodes,
            "group_checks": _group_checks(group),
            "general_K3_surface_automorphism_status": SURFACE_STATUS,
        },
        "subgroup_census": {
            "conjugacy_class_count": len(class_records),
            "proper_class_count": len(proper_internal),
            "classes": class_records,
        },
        "types": types,
        "summary": {
            "proper_subgroup_class_count": 3,
            "distinct_fixed_cone_count": 2,
            "finite_volume_noncompact_count": 2,
            "full_reflection_index_one_count": 2,
            "known_abstract_chamber_count": 2,
            "new_abstract_polytope_count": 0,
            "new_algebro_geometric_realization_priority": "open_literature_check",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.is_file() or OUTPUT.read_text() != rendered:
            raise SystemExit("recomputed certificate differs from archive")
    elif args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
