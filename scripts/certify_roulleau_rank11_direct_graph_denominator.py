#!/usr/bin/env python3
"""Close the two direct-graph rank-11 stocks in Roulleau's K3 atlas."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from coxeter_polytopes.search.u_e8_a1_vinberg import (
    lattice_gram as ue8a1_lattice_gram,
    pair as ue8a1_pair,
    run_search as run_ue8a1_vinberg,
)
from scripts.certify_alexeev_remaining_stable_overlattices import (
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
from scripts.certify_roulleau_ue8a1x3_all_subgroup_folds import (
    _integral_overlattice_census,
)
from scripts.certify_ss_fermat_char3_survivor_fixed_lattices import (
    _saturated_fixed_basis,
)
from scripts.search_k3_high_rank_involution_folds import (
    _coxeter_audit,
    _geometry_certificate,
    _group_checks,
    _matrix_payload,
    _weighted_graph_automorphisms,
)
from scripts.search_k3_nef_subgroup_folds import (
    _all_subgroups,
    _generated_subgroup,
    _identity,
    _orbit_gram,
    _orbits,
    _order,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_rank11_direct_graph_denominator.json"
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
SOURCE_TEX_SHA256 = "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
SOURCE_INVENTORY = ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"


CASES: tuple[dict[str, Any], ...] = (
    {
        "id": "R11_UE8A1",
        "lattice": "U+E8(-1)+A1(-1)",
        "walls": 12,
        "section_lines": [8053, 8119],
        "basis_labels": list(range(1, 12)),
        "edges": [
            *[(left, right, 1) for left, right in zip(
                (1, 2, 3, 5, 6, 7, 8, 9, 10, 11),
                (2, 3, 5, 6, 7, 8, 9, 10, 11),
            )],
            (3, 4, 1),
            (11, 12, 2),
        ],
        "dependent": {12: [2, 4, 6, 3, 5, 4, 3, 2, 1, 0, -1]},
        "ample": [48, 97, 147, 73, 125, 104, 84, 65, 47, 30, 14],
        "ample_square": 848,
        "ample_intersections": [1] * 10 + [2, 28],
        "basis_determinant": 2,
        "group_order": 1,
        "group_structure": "trivial",
        "subgroup_orders": [1],
        "proper_fixed_ranks": [],
        "geometry": {
            "facets": 12,
            "edges": 99,
            "finite": 18,
            "cusps": 1,
            "parabolic": 2,
            "relations": {"0": 1, "2": 55, "3": 10},
            "finite_special": [1, 12, 65, 210, 450, 672, 714, 540, 285, 99, 18],
            "euler": "-1/183936614400",
        },
        "coxiter": [
            "graphs/10-thierryP_1-10.coxiter",
            "graphs/Tumarkin pyramids n+2/10-tum04_12_01-k=2-l=3.coxiter",
        ],
    },
    {
        "id": "R11_UD8A1",
        "lattice": "U+D8(-1)+A1(-1)",
        "walls": 14,
        "section_lines": [8120, 8205],
        "basis_labels": [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12],
        "edges": [
            *[(left, right, 1) for left, right in zip(
                (1, 2, 3, 4, 5, 6, 13),
                (2, 3, 4, 5, 6, 13),
            )],
            (2, 12, 1),
            (12, 11, 1),
            (11, 10, 1),
            (6, 7, 1),
            (7, 8, 2),
            (8, 9, 2),
            (13, 14, 2),
            (10, 9, 2),
            (9, 14, 2),
        ],
        "dependent": {
            10: [-2, -4, -3, -2, -1, 0, 1, 1, 0, -2, -3],
            13: [-3, -6, -5, -4, -3, -2, 0, 1, 1, -2, -4],
            14: [3, 6, 5, 4, 3, 2, 1, 0, -1, 2, 4],
        },
        "ample": [-23, -45, -35, -24, -12, 1, 15, 15, 9, -16, -31],
        "ample_square": 208,
        "ample_intersections": [1, 1, 1, 1, 1, 1, 1, 18, 12, 2, 1, 1, 1, 18],
        "basis_determinant": 8,
        "group_order": 2,
        "group_structure": "C2",
        "subgroup_orders": [1, 2],
        "proper_fixed_ranks": [10],
        "geometry": {
            "facets": 14,
            "edges": 218,
            "finite": 35,
            "cusps": 4,
            "parabolic": 9,
            "relations": {"0": 5, "2": 76, "3": 10},
            "finite_special": [1, 14, 86, 309, 726, 1176, 1344, 1086, 608, 218, 35],
            "euler": "-31/5573836800",
        },
        "coxiter": [],
    },
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _curve_gram(case: dict[str, Any]) -> sp.Matrix:
    gram = -2 * sp.eye(case["walls"])
    for left, right, product in case["edges"]:
        gram[left - 1, right - 1] = gram[right - 1, left - 1] = product
    return gram


def _source_coordinates(
    case: dict[str, Any], gram: sp.Matrix
) -> tuple[sp.Matrix, sp.Matrix, dict[str, bool]]:
    basis_indices = [label - 1 for label in case["basis_labels"]]
    basis_gram = gram.extract(basis_indices, basis_indices)
    rows = gram[:, basis_indices] * basis_gram.inv()
    ample = sp.Matrix(case["ample"])
    dependent = {
        label: [int(value) for value in rows.row(label - 1)]
        for label in case["dependent"]
    }
    checks = {
        "basis_rank": basis_gram.rank() == 11,
        "basis_signature": exact_inertia(basis_gram) == (1, 10, 0),
        "basis_determinant": int(basis_gram.det()) == case["basis_determinant"],
        "basis_coordinate_rows": rows[basis_indices, :] == sp.eye(11),
        "dependent_coordinates": dependent == case["dependent"],
        "all_coordinates_integral": all(value.is_Integer for value in rows),
        "full_gram_reconstruction": rows * basis_gram * rows.T == gram,
        "published_ample_square": int((ample.T * basis_gram * ample)[0])
        == case["ample_square"],
        "published_ample_intersections": [
            int(value) for value in ample.T * basis_gram * rows.T
        ] == case["ample_intersections"],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"source reconstruction failed for {case['id']}: {checks}")
    return basis_gram, rows, checks


def _lift_group(
    group: tuple[tuple[int, ...], ...],
    basis_labels: list[int],
    rows: sp.Matrix,
    basis_gram: sp.Matrix,
) -> dict[tuple[int, ...], sp.Matrix]:
    basis_indices = [label - 1 for label in basis_labels]
    actions = {}
    for permutation in group:
        action = rows[[permutation[index] for index in basis_indices], :]
        expected = rows[[permutation[index] for index in range(rows.rows)], :]
        if not (
            all(value.is_Integer for value in action)
            and abs(int(action.det())) == 1
            and rows * action == expected
            and action * basis_gram * action.T == basis_gram
        ):
            raise ArithmeticError("weighted graph symmetry failed its integral NS lift")
        actions[permutation] = action
    return actions


def _subgroup_records(
    group: tuple[tuple[int, ...], ...],
    actions: dict[tuple[int, ...], sp.Matrix],
    gram: sp.Matrix,
) -> list[dict[str, Any]]:
    records = []
    for subgroup in sorted(_all_subgroups(group), key=lambda row: (len(row), sorted(row))):
        generators = []
        generated = frozenset({_identity(gram.rows)})
        for element in sorted(subgroup):
            if element not in generated:
                generators.append(element)
                generated = _generated_subgroup(generators)
        orbits = _orbits(subgroup, gram.rows)
        orbit_rank = _orbit_gram(gram, orbits).rank()
        if generators:
            fixed_basis, smith_audit = _saturated_fixed_basis(
                [actions[generator] for generator in generators], 11
            )
            action_fixed_rank = fixed_basis.rows
        else:
            smith_audit = None
            action_fixed_rank = 11
        if action_fixed_rank != orbit_rank:
            raise ArithmeticError("orbit and integral-action fixed ranks disagree")
        records.append(
            {
                "order": len(subgroup),
                "element_order_histogram": {
                    str(order): count
                    for order, count in sorted(
                        Counter(_order(element) for element in subgroup).items()
                    )
                },
                "wall_orbits_one_based": [
                    [index + 1 for index in orbit] for orbit in orbits
                ],
                "fixed_rank": orbit_rank,
                "hyperbolic_dimension": orbit_rank - 1,
                "representative_generators_zero_based": [
                    list(generator) for generator in generators
                ],
                "generator_action_matrices": [
                    _matrix_payload(actions[generator]) for generator in generators
                ],
                "saturated_fixed_basis_audit": smith_audit,
            }
        )
    return records


def _full_chamber_geometry(
    case: dict[str, Any],
    gram: sp.Matrix,
    basis_gram: sp.Matrix,
    rows: sp.Matrix,
) -> tuple[dict[str, Any], list[list[str]], dict[str, Any], dict[str, Any]]:
    roots = [sp.Matrix(rows.row(index)).T for index in range(rows.rows)]
    reflection_checks = []
    for index, root in enumerate(roots):
        reflection = _reflection(basis_gram, root)
        reflection_checks.append(
            {
                "wall_one_based": index + 1,
                "root_coordinates": [int(value) for value in root],
                "norm": int((root.T * basis_gram * root)[0]),
                "integral_involutive_isometry": all(value.is_Integer for value in reflection)
                and reflection * reflection == sp.eye(11)
                and reflection.T * basis_gram * reflection == basis_gram,
            }
        )
    reflection_audit = {
        "all_folded_reflections_are_integral": all(
            row["integral_involutive_isometry"] for row in reflection_checks
        ),
        "checks": reflection_checks,
    }
    singleton_orbits = [[index + 1] for index in range(rows.rows)]
    coxeter = _coxeter_audit(
        gram, list(range(rows.rows)), singleton_orbits
    )
    geometry = _geometry_certificate(
        11,
        singleton_orbits,
        gram,
        coxeter,
        reflection_audit,
        {"all_integral_lattice_checks_pass": True},
    )
    relations = _relation_matrix(
        [[str(-gram[left, right]) for right in range(gram.rows)] for left in range(gram.rows)]
    )
    solomon = _euler_from_connected_components(
        relations, _connected_masks_direct(relations)
    )
    coxiter = _scan_coxiter(relations, 10, rows.rows)
    expected = case["geometry"]
    histogram = dict(
        sorted(
            Counter(
                relations[left][right]
                for left in range(len(relations))
                for right in range(left + 1, len(relations))
            ).items()
        )
    )
    checks = {
        "geometry": geometry["all_checks_pass"],
        "counts": geometry["facet_count"] == expected["facets"]
        and geometry["vinberg_volume_criterion"]["elliptic_edge_count"]
        == expected["edges"]
        and geometry["finite_vertex_count"] == expected["finite"]
        and geometry["ideal_vertex_count"] == expected["cusps"]
        and geometry["vinberg_volume_criterion"]["parabolic_neighborhood_audit_count"]
        == expected["parabolic"],
        "relations": histogram == expected["relations"],
        "solomon": solomon["finite_special_subset_count_by_rank"]
        == expected["finite_special"]
        and solomon["rational_euler_characteristic"] == expected["euler"],
        "coxiter": sorted(row["path"] for row in coxiter["ordinary_exact_matches"])
        == sorted(case["coxiter"])
        and not coxiter["doubling_exact_matches"],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"full chamber checks failed for {case['id']}: {checks}")
    return geometry, relations, solomon, {"audit": coxiter, "checks": checks}


def _ue8a1_lineage(basis_gram: sp.Matrix, rows: sp.Matrix) -> dict[str, Any]:
    search = run_ue8a1_vinberg(40)
    intrinsic_gram = sp.Matrix(ue8a1_lattice_gram())
    candidate_roots = rows.T
    candidate_gram = [
        [int((candidate_roots[:, left].T * basis_gram * candidate_roots[:, right])[0])
         for right in range(candidate_roots.cols)]
        for left in range(candidate_roots.cols)
    ]
    intrinsic_root_gram = [
        [ue8a1_pair(left, right) for right in search.roots] for left in search.roots
    ]
    permutation = _exact_gram_permutation(candidate_gram, intrinsic_root_gram)
    if permutation is None:
        raise ArithmeticError("U+E8+A1 walls failed intrinsic Vinberg matching")
    independent = []
    for column in range(candidate_roots.cols):
        if candidate_roots[:, [*independent, column]].rank() > len(independent):
            independent.append(column)
    ordered = sp.Matrix.hstack(
        *(sp.Matrix(search.roots[permutation[column]]) for column in range(len(permutation)))
    )
    isometry = ordered[:, independent] * candidate_roots[:, independent].inv()
    checks = {
        "stopped": search.stopping.fired,
        "root_count": len(search.roots) == 12,
        "finite_vertices": search.stopping.finite_vertex_count == 18,
        "cusps": search.stopping.cusp_count == 1,
        "integral_unimodular_isometry": all(value.is_Integer for value in isometry)
        and abs(int(isometry.det())) == 1,
        "gram_isometry": isometry.T * intrinsic_gram * isometry == basis_gram,
        "wall_by_wall": isometry * candidate_roots == ordered,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"intrinsic U+E8+A1 lineage failed: {checks}")
    return {
        "status": "known_exact_full_reflection_chamber",
        "reflection_subgroup_index": 1,
        "intrinsic_lattice": "U+E8(-1)+A1(-1)",
        "vinberg_stopping_height": str(search.levels[-1].height),
        "facet_permutation_zero_based": permutation,
        "source_to_intrinsic_isometry": _matrix_payload(isometry),
        "checks": checks,
    }


def _ud8a1_lineage(
    case: dict[str, Any], basis_gram: sp.Matrix, rows: sp.Matrix
) -> dict[str, Any]:
    roots = [sp.Matrix(rows.row(index)).T for index in range(rows.rows)]
    reflections = [_reflection(basis_gram, root) for root in roots]
    stable = _stable_subgroups(basis_gram, reflections)
    groups = stable.pop("groups_internal")
    left_inverse = stable.pop("left_smith_inverse")
    maxima = []
    for position in stable["maximal_stable_subgroup_positions_zero_based"]:
        basis = _overlattice_basis(
            basis_gram,
            left_inverse,
            stable["stable_subgroups"][position]["cosets_in_smith_coordinates"],
        )
        gram = sp.simplify(basis.T * basis_gram * basis)
        even = all(gram[index, index] % 2 == 0 for index in range(gram.rows))
        actions = [sp.simplify(basis.inv() * reflection * basis) for reflection in reflections]
        if not (
            len(groups[position]) == 2
            and abs(int(gram.det())) == 2
            and all(value.is_Integer for action in actions for value in action)
            and all(action.T * gram * action == gram for action in actions)
        ):
            raise ArithmeticError("D8+A1 stable maximum failed")
        maxima.append(
            {
                "stable_position_zero_based": position,
                "index_over_source_lattice": 2,
                "basis_columns_in_source_lattice": _payload(basis),
                "gram": _payload(gram),
                "determinant": int(gram.det()),
                "smith_invariants": _smith(gram),
                "parity": "even" if even else "odd",
                "discriminant_q_histogram_mod_two": _q_histogram(gram) if even else None,
                "classification": (
                    "U+E8(-1)+A1(-1) by the indefinite even discriminant-form classification"
                    if even else "odd determinant-two maximum; exact full chamber not identified here"
                ),
                "containing_full_reflection_index": 1023 if even else None,
                "index_proof": (
                    "(-31/5573836800)/(-1/183936614400)=1023"
                    if even else "finite because both chambers have finite covolume; value open"
                ),
            }
        )
    all_overclosures = _integral_overlattice_census(basis_gram)
    checks = {
        "stable_orders": stable["stable_integral_subgroup_orders"] == [1, 2, 2, 2],
        "stable_even_orders": stable["stable_even_subgroup_orders"] == [1, 2, 2],
        "three_incomparable_maxima": len(maxima) == 3,
        "all_integral_overclosures_are_stable": all_overclosures[
            "integral_overlattice_indices"
        ] == [1, 2, 2, 2],
        "two_even_one_odd": Counter(row["parity"] for row in maxima)
        == {"even": 2, "odd": 1},
        "even_indices": [
            row["containing_full_reflection_index"]
            for row in maxima if row["parity"] == "even"
        ] == [1023, 1023],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"D8+A1 lineage changed: {checks}")
    return {
        "status": "published_K3_minus_two_Weyl_chamber; maximal_integral_lineage_partly_closed",
        "source_reflection_group_euler_characteristic": case["geometry"]["euler"],
        "reflection_stable_overlattice_census": stable,
        "all_integral_overlattice_census": all_overclosures,
        "maximal_stable_lattices": maxima,
        "bounded_literal_priority": (
            "no ordinary or doubling CoxIter match, but the chamber is explicitly "
            "published by Roulleau and therefore is not a new abstract example"
        ),
        "checks": checks,
    }


def _analyze(case: dict[str, Any], stock: dict[str, Any]) -> dict[str, Any]:
    gram = _curve_gram(case)
    basis_gram, rows, source_checks = _source_coordinates(case, gram)
    group, search_nodes = _weighted_graph_automorphisms(gram)
    if not (
        len(group) == case["group_order"]
        and all(_group_checks(group).values())
    ):
        raise ArithmeticError("weighted graph group changed")
    actions = _lift_group(group, case["basis_labels"], rows, basis_gram)
    subgroup_records = _subgroup_records(group, actions, gram)
    proper_fixed_ranks = [
        row["fixed_rank"] for row in subgroup_records if row["order"] > 1
    ]
    geometry, relations, solomon, coxiter = _full_chamber_geometry(
        case, gram, basis_gram, rows
    )
    lineage = (
        _ue8a1_lineage(basis_gram, rows)
        if case["id"] == "R11_UE8A1"
        else _ud8a1_lineage(case, basis_gram, rows)
    )
    checks = {
        "inventory_metadata": stock["picard_rank"] == 11
        and stock["wall_count"] == case["walls"]
        and stock["source"]["section_lines"] == case["section_lines"],
        "source_coordinates": all(source_checks.values()),
        "wall_gram_rank_signature": gram.rank() == 11
        and exact_inertia(gram) == (1, 10, case["walls"] - 11),
        "weighted_group": len(group) == case["group_order"],
        "subgroup_orders": [row["order"] for row in subgroup_records]
        == case["subgroup_orders"],
        "proper_fixed_ranks": proper_fixed_ranks == case["proper_fixed_ranks"],
        "zero_proper_sections_above_H9": all(rank <= 10 for rank in proper_fixed_ranks),
        "full_chamber": all(coxiter["checks"].values()),
        "lineage": all(lineage["checks"].values()),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"case checks failed for {case['id']}: {checks}")
    return {
        "id": case["id"],
        "source": {
            "lattice": case["lattice"],
            "picard_rank": 11,
            "wall_count": case["walls"],
            "section_lines": case["section_lines"],
            "edges_one_based": [list(edge) for edge in case["edges"]],
            "basis_wall_labels": case["basis_labels"],
            "wall_gram": _matrix_payload(gram),
            "basis_gram": _matrix_payload(basis_gram),
            "wall_coordinate_rows": _matrix_payload(rows),
            "checks": source_checks,
        },
        "weighted_wall_group": {
            "structure": case["group_structure"],
            "order": len(group),
            "element_order_histogram": {
                str(order): count
                for order, count in sorted(Counter(_order(element) for element in group).items())
            },
            "exhaustive_backtracking_node_count": search_nodes,
            "group_checks": _group_checks(group),
            "integral_NS_actions": [
                {
                    "permutation_zero_based": list(permutation),
                    "matrix": _matrix_payload(action),
                    "faithful_nonidentity": permutation == _identity(case["walls"])
                    or action != sp.eye(11),
                }
                for permutation, action in sorted(actions.items())
            ],
            "general_surface_action_status": (
                "Roulleau Aut(X)=C2 tag t acts trivially on NS; the weighted graph "
                "group is trivial"
                if len(group) == 1
                else "the nontrivial diagram involution acts nontrivially on NS and is "
                "not the tag-t involution of the general K3"
            ),
        },
        "subgroup_denominator": {
            "subgroup_count": len(subgroup_records),
            "records": subgroup_records,
            "proper_fixed_rank_at_least_11_count": sum(
                row["order"] > 1 and row["fixed_rank"] >= 11
                for row in subgroup_records
            ),
            "closure_reason": (
                "a faithful nonidentity action on rank-11 NS has fixed rank at most 10; "
                "the exhaustive records realize this bound"
            ),
        },
        "published_full_nef_chamber": {
            "geometry": geometry,
            "normalized_relation_matrix": relations,
            "solomon": solomon,
            "coxiter": coxiter["audit"],
            "lineage": lineage,
            "priority": "published_known_positive_example_not_a_new_abstract_polytope",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def build_certificate() -> dict[str, Any]:
    if _sha256(SOURCE_TEX) != SOURCE_TEX_SHA256:
        raise ArithmeticError("Roulleau TeX hash changed")
    inventory = json.loads(SOURCE_INVENTORY.read_text())
    stock_by_key = {row["key"]: row for row in inventory["stocks"]}
    lines = SOURCE_TEX.read_text().splitlines()
    for case in CASES:
        start, end = case["section_lines"]
        excerpt = "\n".join(lines[start - 1 : end])
        stock = stock_by_key[case["id"]]
        if not (
            f"contains ${case['walls']}$" in excerpt
            and stock["tex_lattice_token"] in excerpt
        ):
            raise ArithmeticError(f"source witness failed for {case['id']}")
    records = [_analyze(case, stock_by_key[case["id"]]) for case in CASES]
    checks = {
        "inventory_certified": inventory["all_checks_pass"] is True,
        "two_cases": len(records) == 2,
        "both_full_chambers_finite_volume_noncompact": all(
            row["published_full_nef_chamber"]["geometry"]["finite_volume_certified"]
            and row["published_full_nef_chamber"]["geometry"]["compact"] is False
            for row in records
        ),
        "zero_proper_folds_above_H9": sum(
            row["subgroup_denominator"]["proper_fixed_rank_at_least_11_count"]
            for row in records
        ) == 0,
        "all_case_checks": all(row["all_checks_pass"] for row in records),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global checks failed: {checks}")
    return {
        "schema": "coxeter.roulleau_rank11_direct_graph_denominator.v1",
        "artifact_id": "roulleau_rank11_direct_graph_denominator",
        "source": {
            "roulleau_tex_sha256": SOURCE_TEX_SHA256,
            "inventory_path": str(SOURCE_INVENTORY.relative_to(ROOT)),
        },
        "cases": records,
        "summary": {
            "source_stock_count": 2,
            "source_finite_volume_noncompact_chamber_count": 2,
            "weighted_group_orders": [1, 2],
            "total_subgroup_count": 3,
            "proper_fixed_rank_at_least_11_count": 0,
            "new_higher_dimensional_fold_count": 0,
            "new_abstract_polytope_count": 0,
            "rank11_direct_graph_denominator_status": "closed",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-archive", action="store_true")
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
