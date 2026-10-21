#!/usr/bin/env python3
"""Certify every high subgroup fold of Roulleau's graph-only rank-12 stocks."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Any, Sequence

import networkx as nx
import sympy as sp

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
from scripts.certify_roulleau_non2elementary_involution_folds import _scan_coxiter
from scripts.certify_ss_fermat_char3_survivor_fixed_lattices import (
    _saturated_fixed_basis,
)
from scripts.search_k3_high_rank_involution_folds import (
    _coxeter_audit,
    _fold_reflection_audit,
    _geometry_certificate,
    _group_checks,
    _matrix_payload,
    _weighted_graph_automorphisms,
)
from scripts.search_k3_nef_subgroup_folds import (
    Permutation,
    Subgroup,
    _all_subgroups,
    _compose,
    _generated_subgroup,
    _identity,
    _orbit_gram,
    _orbits,
    _order,
    _pairwise_disjoint,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"
R14_27_ARTIFACT = ROOT / "data/enumerations/roulleau_r14_ue8a1x4_all_subgroup_folds.json"
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
SOURCE_TEX_SHA256 = "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"


def _gram(size: int, edges: Sequence[tuple[int, int, int]]) -> sp.Matrix:
    gram = -2 * sp.eye(size)
    for left, right, product in edges:
        if gram[left - 1, right - 1] != 0:
            raise ArithmeticError(f"duplicate source edge {left}-{right}")
        gram[left - 1, right - 1] = product
        gram[right - 1, left - 1] = product
    return gram


R12_UE8A1X2_EDGES = tuple(
    (left, right, 1)
    for left, right in zip((1, 2, 3, 5, 6, 7, 8, 9, 10), (2, 3, 5, 6, 7, 8, 9, 10))
) + (
    (3, 4, 1),
    (10, 11, 1),
    (10, 12, 1),
    (11, 13, 2),
    (12, 14, 2),
)


R12_UD8A1X2_EDGES = tuple(
    (left, right, 1)
    for left, right in zip((1, 2, 3, 5, 6, 7, 8, 10), (2, 3, 5, 6, 7, 8, 10))
) + (
    (1, 12, 1),
    (1, 13, 1),
    (3, 4, 1),
    (8, 9, 1),
    (8, 11, 1),
    (16, 17, 2),
    (17, 18, 2),
    (19, 18, 2),
    (19, 14, 2),
    (14, 15, 2),
    (16, 15, 2),
    (12, 19, 2),
    (19, 10, 2),
    (16, 10, 2),
    (16, 13, 2),
    (12, 17, 2),
    (13, 14, 2),
    (12, 15, 2),
    (13, 18, 2),
    (11, 14, 2),
    (11, 17, 2),
    (9, 18, 2),
    (9, 15, 2),
)


STOCKS = (
    {
        "key": "R12_UE8A1x2",
        "lattice": "U+E8(-1)+A1(-1)^2",
        "walls": 14,
        "basis_rank": 12,
        "basis_determinant": -4,
        "source_lines": [8545, 8647],
        "edges": R12_UE8A1X2_EDGES,
    },
    {
        "key": "R12_UD8A1x2",
        "lattice": "U+D8(-1)+A1(-1)^2",
        "walls": 19,
        "basis_rank": 12,
        "basis_determinant": -16,
        "source_lines": [8648, 8803],
        "edges": R12_UD8A1X2_EDGES,
    },
)


EXPECTED = {
    "R12_UE8A1x2": {
        "group_order": 2,
        "group_orders": {"1": 1, "2": 1},
        "class_count": 2,
        "survivors": [2],
        "types": {
            2: {
                "facets": 12,
                "finite": 17,
                "cusps": 2,
                "edges": 99,
                "parabolic_audits": 3,
                "determinant": 4,
                "stable_orders": [1, 2],
                "maximal_count": 1,
                "index": 1,
                "euler": "-1/11147673600",
                "ordinary_matches": [
                    "graphs/Tumarkin pyramids n+2/10-tum04_12_01-k=2-l=4.coxiter",
                    "graphs/group units qf/10-1,1,1,1,1,1,1,1,1,1,1.coxiter",
                ],
                "doubling_match_count": 2,
            }
        },
    },
    "R12_UD8A1x2": {
        "group_order": 12,
        "group_orders": {"1": 1, "2": 7, "3": 2, "6": 2},
        "class_count": 10,
        "survivors": [2, 3],
        "types": {
            2: {
                "facets": 15,
                "finite": 42,
                "cusps": 10,
                "edges": 294,
                "parabolic_audits": 17,
                "determinant": 16,
                "stable_orders": [1, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4],
                "maximal_count": 3,
                "index": 1054,
                "euler": "-527/5573836800",
                "ordinary_matches": [],
                "doubling_match_count": 0,
            },
            3: {
                "facets": 16,
                "finite": 70,
                "cusps": 6,
                "edges": 426,
                "parabolic_audits": 14,
                "determinant": 16,
                "stable_orders": [1, 2, 2, 2, 4],
                "maximal_count": 1,
                "index": 992,
                "euler": "-31/348364800",
                "ordinary_matches": [],
                "doubling_match_count": 0,
            },
        },
    },
}


def _public_permutation(permutation: Permutation) -> list[int]:
    return [value + 1 for value in permutation]


def _relation_graph_isomorphic(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> bool:
    if len(left) != len(right):
        return False

    def graph(matrix: Sequence[Sequence[int]]) -> nx.Graph:
        result = nx.Graph()
        result.add_nodes_from(range(len(matrix)))
        for row in range(len(matrix)):
            for column in range(row):
                result.add_edge(row, column, relation=str(matrix[row][column]))
        return result

    return nx.is_isomorphic(
        graph(left),
        graph(right),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("relation", None),
    )


def _conjugate_subgroup(subgroup: Subgroup, action: Permutation) -> Subgroup:
    inverse = [0] * len(action)
    for source, target in enumerate(action):
        inverse[target] = source
    inverse_tuple = tuple(inverse)
    return frozenset(
        _compose(_compose(action, element), inverse_tuple) for element in subgroup
    )


def _subgroup_classes(group: tuple[Permutation, ...]) -> list[dict[str, Any]]:
    remaining = set(_all_subgroups(group))
    records = []
    while remaining:
        representative = min(remaining, key=lambda subgroup: (len(subgroup), sorted(subgroup)))
        conjugates = {_conjugate_subgroup(representative, action) for action in group}
        remaining.difference_update(conjugates)
        generators = []
        generated: Subgroup = frozenset({_identity(len(group[0]))})
        for element in sorted(representative):
            if element in generated:
                continue
            generators.append(element)
            generated = _generated_subgroup(generators)
        histogram = dict(sorted(Counter(_order(element) for element in representative).items()))
        records.append(
            {
                "order": len(representative),
                "conjugacy_class_size": len(conjugates),
                "element_order_histogram": {str(k): v for k, v in histogram.items()},
                "representative_generators": tuple(generators),
                "representative_subgroup": representative,
            }
        )
    records.sort(
        key=lambda row: (
            row["order"],
            tuple(sorted(row["element_order_histogram"].items())),
            sorted(row["representative_subgroup"]),
        )
    )
    for index, record in enumerate(records, start=1):
        record["class_index"] = index
    return records


def _lift_actions(
    group: Sequence[Permutation], wall_rows: sp.Matrix, basis_gram: sp.Matrix
) -> dict[Permutation, sp.Matrix]:
    """Lift every wall permutation through the first-rank wall basis."""

    rank = basis_gram.rows
    wall_count = wall_rows.rows
    actions = {}
    for permutation in group:
        action = wall_rows[[permutation[index] for index in range(rank)], :]
        expected = wall_rows[[permutation[index] for index in range(wall_count)], :]
        if not (
            all(value.is_Integer for value in action)
            and abs(int(action.det())) == 1
            and wall_rows * action == expected
            and action * basis_gram * action.T == basis_gram
        ):
            raise ArithmeticError("diagram symmetry failed its integral NS lift")
        actions[permutation] = action
    return actions


def _divisibility(gram: sp.Matrix, root: sp.Matrix) -> int:
    result = 0
    for value in root.T * gram:
        result = gcd(result, abs(int(value)))
    return result


def _fixed_lattice_data(
    generators: Sequence[Permutation],
    actions: dict[Permutation, sp.Matrix],
    basis_gram: sp.Matrix,
    wall_rows: sp.Matrix,
    orbits: tuple[tuple[int, ...], ...],
    orbit_gram: sp.Matrix,
) -> dict[str, Any]:
    rank = basis_gram.rows
    fixed_basis, kernel = _saturated_fixed_basis(
        [actions[generator] for generator in generators], rank
    )
    fixed_gram = fixed_basis * basis_gram * fixed_basis.T
    roots = []
    reflection_checks = []
    for orbit_index, orbit in enumerate(orbits, start=1):
        normal = sum((wall_rows.row(index) for index in orbit), sp.zeros(1, rank))
        solution = tuple(next(iter(sp.linsolve((fixed_basis.T, normal.T)))))
        root = sp.Matrix(solution)
        norm = int((root.T * fixed_gram * root)[0])
        reflection = _reflection(fixed_gram, root)
        check = {
            "orbit_index": orbit_index,
            "source_walls_one_based": [index + 1 for index in orbit],
            "fixed_coordinates": [int(value) for value in root],
            "norm": norm,
            "divisibility": _divisibility(fixed_gram, root),
            "integral_reflection": (
                root.T * fixed_basis == normal
                and all(value.is_Integer for value in root)
                and all(value.is_Integer for value in reflection)
                and reflection * reflection == sp.eye(fixed_gram.rows)
                and reflection.T * fixed_gram * reflection == fixed_gram
            ),
        }
        roots.append(root)
        reflection_checks.append(check)
    root_matrix = sp.Matrix.hstack(*roots)
    if root_matrix.T * fixed_gram * root_matrix != orbit_gram:
        raise ArithmeticError("fixed-lattice roots do not recover the orbit Gram")
    reflections = [_reflection(fixed_gram, root) for root in roots]
    stable = _stable_subgroups(fixed_gram, reflections)
    groups = stable.pop("groups_internal")
    left_inverse = stable.pop("left_smith_inverse")
    maximal_overlattices = []
    for position in stable["maximal_stable_subgroup_positions_zero_based"]:
        basis = _overlattice_basis(
            fixed_gram,
            left_inverse,
            stable["stable_subgroups"][position]["cosets_in_smith_coordinates"],
        )
        maximal_gram = sp.simplify(basis.T * fixed_gram * basis)
        reflection_actions = [
            sp.simplify(basis.inv() * reflection * basis) for reflection in reflections
        ]
        parity = (
            "even"
            if all(maximal_gram[index, index] % 2 == 0 for index in range(maximal_gram.rows))
            else "odd"
        )
        maximal_overlattices.append(
            {
                "stable_subgroup_position_zero_based": position,
                "index_over_fixed_lattice": len(groups[position]),
                "basis_columns_in_fixed_lattice": _payload(basis),
                "gram": _payload(maximal_gram),
                "determinant": int(maximal_gram.det()),
                "smith_invariants": _smith(maximal_gram),
                "parity": parity,
                "discriminant_q_histogram_mod_two": (
                    _q_histogram(maximal_gram) if parity == "even" else None
                ),
                "all_reflections_integral": all(
                    value.is_Integer for action in reflection_actions for value in action
                ),
                "all_reflections_preserve_gram": all(
                    action.T * maximal_gram * action == maximal_gram
                    for action in reflection_actions
                ),
            }
        )
    return {
        "basis_rows_in_source_basis": _matrix_payload(fixed_basis),
        "gram": _matrix_payload(fixed_gram),
        "rank": fixed_gram.rows,
        "determinant": int(fixed_gram.det()),
        "smith_invariants": _smith(fixed_gram),
        "even": all(fixed_gram[index, index] % 2 == 0 for index in range(fixed_gram.rows)),
        "discriminant_q_histogram_mod_two": _q_histogram(fixed_gram),
        "smith_kernel_audit": kernel,
        "orbit_root_coordinate_columns": _matrix_payload(root_matrix),
        "reflection_checks": reflection_checks,
        "all_reflections_integral": all(row["integral_reflection"] for row in reflection_checks),
        "stable_overlattice_audit": stable,
        "maximal_stable_overlattices": maximal_overlattices,
    }


def _source_formula_checks(key: str, wall_rows: sp.Matrix) -> dict[str, bool]:
    """Check the extra-wall coordinates against Roulleau's printed D2 relations."""

    if key == "R12_UE8A1x2":
        d2 = sp.Matrix([[2, 4, 6, 3, 5, 4, 3, 2, 2, 2, 1, 1]])
        decompositions = (
            wall_rows.row(8) + 2 * wall_rows.row(9) + 2 * wall_rows.row(10)
            + wall_rows.row(11) + wall_rows.row(12),
            wall_rows.row(8) + 2 * wall_rows.row(9) + wall_rows.row(10)
            + 2 * wall_rows.row(11) + wall_rows.row(13),
        )
    elif key == "R12_UD8A1x2":
        d2 = sp.Matrix([[0, 1, 2, 1, 2, 2, 2, 2, 1, 1, 1, 0]])
        decompositions = (
            wall_rows.row(10) + wall_rows.row(11) + wall_rows.row(16),
            wall_rows.row(10) + wall_rows.row(12) + wall_rows.row(13),
            wall_rows.row(9) + wall_rows.row(11) + wall_rows.row(18),
            wall_rows.row(9) + wall_rows.row(12) + wall_rows.row(15),
            wall_rows.row(8) + wall_rows.row(11) + wall_rows.row(14),
            wall_rows.row(8) + wall_rows.row(12) + wall_rows.row(17),
            2 * wall_rows.row(0) + 3 * wall_rows.row(1)
            + 4 * wall_rows.row(2) + 2 * wall_rows.row(3)
            + 3 * wall_rows.row(4) + 2 * wall_rows.row(5)
            + wall_rows.row(6) + wall_rows.row(11) + wall_rows.row(12),
        )
    else:
        raise ValueError(key)
    return {
        "printed_D2_coordinate_vector": all(value.is_Integer for value in d2),
        "every_printed_D2_decomposition": all(row == d2 for row in decompositions),
    }


def _stock_certificate(spec: dict[str, Any]) -> dict[str, Any]:
    wall_gram = _gram(spec["walls"], spec["edges"])
    rank = spec["basis_rank"]
    basis_gram = wall_gram[:rank, :rank]
    wall_rows = wall_gram[:, :rank] * basis_gram.inv()
    reconstruction_checks = {
        "basis_determinant": int(basis_gram.det()) == spec["basis_determinant"],
        "basis_inertia": exact_inertia(basis_gram) == (1, rank - 1, 0),
        "basis_rows_identity": wall_rows[:rank, :] == sp.eye(rank),
        "all_wall_coordinates_integral": all(value.is_Integer for value in wall_rows),
        "full_gram_reconstructed": wall_rows * basis_gram * wall_rows.T == wall_gram,
        "wall_gram_rank": wall_gram.rank() == rank,
        "wall_gram_inertia": exact_inertia(wall_gram) == (1, rank - 1, spec["walls"] - rank),
        **_source_formula_checks(spec["key"], wall_rows),
    }
    if not all(reconstruction_checks.values()):
        raise ArithmeticError(f"{spec['key']} source reconstruction failed: {reconstruction_checks}")

    group, search_nodes = _weighted_graph_automorphisms(wall_gram)
    if not all(_group_checks(group).values()):
        raise ArithmeticError(f"{spec['key']} graph automorphisms are not a group")
    actions = _lift_actions(group, wall_rows, basis_gram)
    classes = _subgroup_classes(group)
    prior_h10_15 = None
    if spec["key"] == "R12_UD8A1x2":
        prior_payload = json.loads(R14_27_ARTIFACT.read_text())
        prior_h10_15 = next(
            case
            for case in prior_payload["cases"]
            if case["geometry"]["hyperbolic_dimension"] == 10
            and case["geometry"]["facet_count"] == 15
        )
    public_classes = []
    survivors = []
    survivor_records = []
    for row in classes:
        subgroup: Subgroup = row.pop("representative_subgroup")
        generators = row.pop("representative_generators")
        orbits = _orbits(subgroup, spec["walls"])
        orbit_gram = _orbit_gram(wall_gram, orbits)
        fixed_rank = orbit_gram.rank()
        disjoint = all(_pairwise_disjoint(wall_gram, orbit) for orbit in orbits)
        negative_orbits = all(orbit_gram[index, index] < 0 for index in range(orbit_gram.rows))
        coxeter = _coxeter_audit(orbit_gram, range(orbit_gram.rows), [[i + 1 for i in orbit] for orbit in orbits])
        record = {
            **row,
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "wall_orbits_one_based": [[value + 1 for value in orbit] for orbit in orbits],
            "representative_generators_one_based": [_public_permutation(g) for g in generators],
            "all_orbits_pairwise_disjoint": disjoint,
            "all_orbit_normals_negative": negative_orbits,
            "exact_coxeter_pair_gate": coxeter,
        }
        public_classes.append(record)
        passes_gate = (
            len(subgroup) > 1
            and fixed_rank >= 11
            and disjoint
            and negative_orbits
            and coxeter["all_pairs_are_exact_coxeter_relations"]
        )
        if passes_gate:
            survivors.append(record["class_index"])
            cycles_one = [[value + 1 for value in orbit] for orbit in orbits]
            reflection_audit = _fold_reflection_audit(
                wall_gram, orbit_gram, list(range(orbit_gram.rows)), cycles_one
            )
            ns_checks = {
                "actions_integral": all(
                    value.is_Integer
                    for generator in generators
                    for value in actions[generator]
                ),
                "actions_unimodular": all(
                    abs(int(actions[generator].det())) == 1 for generator in generators
                ),
                "actions_preserve_NS": all(
                    actions[generator] * basis_gram * actions[generator].T == basis_gram
                    for generator in generators
                ),
            }
            ns_action = {
                "generator_permutations_one_based": [
                    _public_permutation(generator) for generator in generators
                ],
                "generator_action_matrices": [
                    _matrix_payload(actions[generator]) for generator in generators
                ],
                "checks": ns_checks,
                "all_integral_lattice_checks_pass": all(ns_checks.values()),
                "surface_realization_status": (
                    "not_realized_on_the_general_K3: Roulleau records Aut(X)=C2 "
                    "with tag t, hence the general surface involution acts trivially "
                    "on NS while this diagram action is nontrivial"
                ),
            }
            geometry = _geometry_certificate(
                rank, cycles_one, orbit_gram, coxeter, reflection_audit, ns_action
            )
            fixed_lattice = _fixed_lattice_data(
                generators,
                actions,
                basis_gram,
                wall_rows,
                orbits,
                orbit_gram,
            )
            positive = [
                [str(-orbit_gram[left, right]) for right in range(orbit_gram.rows)]
                for left in range(orbit_gram.rows)
            ]
            relations = _relation_matrix(positive)
            solomon = _euler_from_connected_components(
                relations, _connected_masks_direct(relations)
            )
            coxiter_lineage = _scan_coxiter(
                relations, fixed_rank - 1, orbit_gram.rows
            )
            expected_type = EXPECTED[spec["key"]]["types"][record["class_index"]]
            vinberg = geometry["vinberg_volume_criterion"]
            maximal = fixed_lattice["maximal_stable_overlattices"]
            ordinary_paths = sorted(
                match["path"] for match in coxiter_lineage["ordinary_exact_matches"]
            )
            full_euler = sp.Rational(-1, 11147673600)
            reflection_subgroup_index = int(
                sp.Rational(solomon["rational_euler_characteristic"]) / full_euler
            )
            lineage_checks = {
                "geometry": geometry["all_checks_pass"],
                "finite_volume_noncompact": (
                    vinberg["decision"] == "finite_volume_noncompact"
                    and not vinberg["finite_volume_bad_edges"]
                ),
                "volume_counts": (
                    geometry["hyperbolic_dimension"] == 10
                    and geometry["facet_count"] == expected_type["facets"]
                    and geometry["finite_vertex_count"] == expected_type["finite"]
                    and geometry["ideal_vertex_count"] == expected_type["cusps"]
                    and vinberg["elliptic_edge_count"] == expected_type["edges"]
                    and vinberg["parabolic_neighborhood_audit_count"]
                    == expected_type["parabolic_audits"]
                ),
                "fixed_lattice": (
                    fixed_lattice["determinant"] == expected_type["determinant"]
                    and fixed_lattice["rank"] == 11
                    and fixed_lattice["even"]
                    and fixed_lattice["all_reflections_integral"]
                ),
                "stable_overlattice_denominator": (
                    fixed_lattice["stable_overlattice_audit"][
                        "stable_integral_subgroup_orders"
                    ]
                    == expected_type["stable_orders"]
                    and len(maximal) == expected_type["maximal_count"]
                    and all(
                        row["index_over_fixed_lattice"]
                        == (2 if spec["key"] == "R12_UE8A1x2" else 4)
                        and row["determinant"] == 1
                        and row["parity"] == "odd"
                        and row["all_reflections_integral"]
                        and row["all_reflections_preserve_gram"]
                        for row in maximal
                    )
                ),
                "solomon": (
                    solomon["rational_euler_characteristic"] == expected_type["euler"]
                    and reflection_subgroup_index == expected_type["index"]
                ),
                "coxiter_bounded_comparison": (
                    ordinary_paths == sorted(expected_type["ordinary_matches"])
                    and len(coxiter_lineage["doubling_exact_matches"])
                    == expected_type["doubling_match_count"]
                ),
                "prior_project_same_size_comparison": (
                    spec["key"] != "R12_UD8A1x2"
                    or record["class_index"] != 2
                    or not _relation_graph_isomorphic(
                        relations, prior_h10_15["folded_walls"]["relation_matrix"]
                    )
                ),
            }
            if not all(lineage_checks.values()):
                raise ArithmeticError(
                    f"{spec['key']} class {record['class_index']} failed: {lineage_checks}"
                )
            survivor_records.append(
                {
                    "class_index": record["class_index"],
                    "subgroup_order": len(subgroup),
                    "orbit_gram": _matrix_payload(orbit_gram),
                    "reflection_audit": reflection_audit,
                    "neron_severi_action": ns_action,
                    "geometry": geometry,
                    "fixed_lattice": fixed_lattice,
                    "normalized_relation_matrix": relations,
                    "relation_histogram": {
                        str(label): count
                        for label, count in sorted(
                            Counter(
                                relations[left][right]
                                for left in range(len(relations))
                                for right in range(left + 1, len(relations))
                            ).items()
                        )
                    },
                    "solomon": solomon,
                    "coxiter_lineage": coxiter_lineage,
                    "maximal_reflection_lineage": {
                        "containing_lattice": "I_(1,10), odd unimodular",
                        "scharlau_kirschmer_list": "List10",
                        "entry_one_based": 1,
                        "full_chamber_facets": 12,
                        "full_chamber_finite_vertices": 17,
                        "full_chamber_cusps": 2,
                        "full_chamber_euler_characteristic": str(full_euler),
                        "reflection_subgroup_index": reflection_subgroup_index,
                        "status": (
                            "known_exact_full_chamber"
                            if reflection_subgroup_index == 1
                            else "finite_index_subgroup_of_known_full_chamber"
                        ),
                    },
                    "priority_status": (
                        "known_abstract_full_chamber"
                        if reflection_subgroup_index == 1
                        else "literal_diagram_unmatched_in_pinned_corpora_worldwide_priority_open"
                    ),
                    "prior_project_ledger_comparison": (
                        {
                            "candidate_id": "R14-E8A1x4-V4-H10-15",
                            "candidate_artifact": str(R14_27_ARTIFACT.relative_to(ROOT)),
                            "candidate_relation_matrix_sha256": hashlib.sha256(
                                json.dumps(
                                    prior_h10_15["folded_walls"]["relation_matrix"],
                                    separators=(",", ":"),
                                ).encode()
                            ).hexdigest(),
                            "exact_relation_graph_isomorphic": False,
                            "status": "distinct_literal_H10_15_diagram",
                        }
                        if spec["key"] == "R12_UD8A1x2"
                        and record["class_index"] == 2
                        else {
                            "same_dimension_and_facet_prior_ledger_candidate_count": 0
                        }
                    ),
                    "checks": lineage_checks,
                    "all_checks_pass": all(lineage_checks.values()),
                }
            )

    expected_stock = EXPECTED[spec["key"]]
    group_histogram = {
        str(order): count
        for order, count in sorted(Counter(_order(element) for element in group).items())
    }
    stock_checks = {
        "source_reconstruction": all(reconstruction_checks.values()),
        "weighted_group": (
            len(group) == expected_stock["group_order"]
            and group_histogram == expected_stock["group_orders"]
        ),
        "subgroup_denominator": len(public_classes) == expected_stock["class_count"],
        "survivor_partition": survivors == expected_stock["survivors"],
        "all_survivors_certified": all(row["all_checks_pass"] for row in survivor_records),
    }
    if not all(stock_checks.values()):
        raise ArithmeticError(f"{spec['key']} stock checks failed: {stock_checks}")
    return {
        "stock_key": spec["key"],
        "lattice": spec["lattice"],
        "source_lines": spec["source_lines"],
        "wall_count": spec["walls"],
        "source_edges_one_based": [list(edge) for edge in spec["edges"]],
        "wall_gram": _matrix_payload(wall_gram),
        "basis_gram": _matrix_payload(basis_gram),
        "wall_coordinate_rows": _matrix_payload(wall_rows),
        "reconstruction_checks": reconstruction_checks,
        "weighted_diagram_group": {
            "order": len(group),
            "structure": "C2" if len(group) == 2 else "D12 = S3 x C2",
            "element_order_histogram": group_histogram,
            "search_nodes": search_nodes,
            "checks": _group_checks(group),
        },
        "subgroup_class_count": len(public_classes),
        "proper_high_coxeter_survivor_class_indices": survivors,
        "survivors": survivor_records,
        "classes": public_classes,
        "checks": stock_checks,
        "all_checks_pass": all(stock_checks.values()),
    }


def build_certificate() -> dict[str, Any]:
    if hashlib.sha256(SOURCE_TEX.read_bytes()).hexdigest() != SOURCE_TEX_SHA256:
        raise ArithmeticError("pinned Roulleau TeX hash changed")
    stocks = [_stock_certificate(spec) for spec in STOCKS]
    checks = {
        "two_complete_stocks": len(stocks) == 2,
        "all_stock_checks_pass": all(stock["all_checks_pass"] for stock in stocks),
        "three_positive_types": sum(len(stock["survivors"]) for stock in stocks) == 3,
    }
    return {
        "schema": "coxeter.roulleau_rank12_graph_all_subgroup_folds.v1",
        "artifact_id": "roulleau_rank12_graph_all_subgroup_folds",
        "source_tex_sha256": SOURCE_TEX_SHA256,
        "stocks": stocks,
        "summary": {
            "stock_count": len(stocks),
            "subgroup_conjugacy_class_count": sum(
                stock["subgroup_class_count"] for stock in stocks
            ),
            "proper_high_coxeter_survivor_class_count": sum(
                len(stock["proper_high_coxeter_survivor_class_indices"])
                for stock in stocks
            ),
            "finite_volume_noncompact_type_count": 3,
            "known_abstract_full_chamber_count": 1,
            "bounded_literal_diagram_priority_candidate_count": 2,
            "new_arithmetic_commensurability_class_count": 0,
            "general_surface_geometric_action_count": 0,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        if certificate != json.loads(OUTPUT.read_text()):
            raise SystemExit("recomputed certificate differs from archive")
    elif args.output:
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    else:
        print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
