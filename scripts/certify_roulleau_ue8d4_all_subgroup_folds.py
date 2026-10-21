#!/usr/bin/env python3
"""Certify all subgroup folds of Roulleau's 15-wall U+E8+D4 K3 stock."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Any, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _canonical_sha256,
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
OUTPUT = ROOT / "data/enumerations/roulleau_ue8d4_all_subgroup_folds.json"
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
SOURCE_TEX_SHA256 = "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
SOURCE_INVENTORY = ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
SCHEMA = "coxeter.roulleau_ue8d4_all_subgroup_folds.v1"
RANK = 14
WALL_COUNT = 15


CHAIN = (1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 14)
SOURCE_EDGES = tuple((*edge, 1) for edge in zip(CHAIN, CHAIN[1:])) + (
    (3, 4, 1),
    (12, 13, 1),
    (12, 15, 1),
)


EXPECTED_TYPES = {
    "R14_UE8D4_C2": {
        "class_indices": [2],
        "dimension": 12,
        "facets": 14,
        "edges": 234,
        "finite_vertices": 35,
        "cusps": 2,
        "parabolic_audits": 3,
        "relations": {"2": 78, "3": 12, "4": 1},
        "fixed_determinant": 4,
        "fixed_smith": [1] * 12 + [4],
        "stable_orders": [1, 2],
        "stable_even_orders": [1],
        "maximal_index": 2,
        "maximal_determinant": 1,
        "maximal_parity": "odd",
        "classification": "I_(1,12), unique odd unimodular Lorentzian lattice",
        "sk_list": "List12",
        "sk_entry": 1,
        "sk_symbol": "I_{12,1}()",
        "coxiter_matches": [
            "graphs/Tumarkin pyramids n+2/12-tum04_12_01.coxiter",
            "graphs/group units qf/12-1,1,1,1,1,1,1,1,1,1,1,1,1.coxiter",
        ],
        "euler": "691/382588157952000",
    },
    "R14_UE8D4_C3_S3": {
        "class_indices": [3, 4],
        "dimension": 11,
        "facets": 13,
        "edges": 162,
        "finite_vertices": 27,
        "cusps": 1,
        "parabolic_audits": 2,
        "relations": {"2": 66, "3": 11, "6": 1},
        "fixed_determinant": -3,
        "fixed_smith": [1] * 11 + [3],
        "stable_orders": [1],
        "stable_even_orders": [1],
        "maximal_index": 1,
        "maximal_determinant": -3,
        "maximal_parity": "even",
        "classification": "U+E8(-1)+A2(-1), even determinant-three lattice",
        "sk_list": "List11",
        "sk_entry": 4,
        "sk_symbol": "II_{11,1}(3^{-1})",
        "coxiter_matches": [
            "graphs/Tumarkin pyramids n+2/11-tum04_12_05.coxiter",
        ],
        "euler": "0",
    },
}


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
    expected_last = [2, 4, 6, 3, 5, 4, 3, 2, 1, 0, -1, -2, -1, -1]
    if not (
        int(basis_gram.det()) == -4
        and exact_inertia(basis_gram) == (1, 13, 0)
        and wall_rows[:RANK, :] == sp.eye(RANK)
        and [int(value) for value in wall_rows.row(14)] == expected_last
        and all(value.is_Integer for value in wall_rows)
        and wall_rows * basis_gram * wall_rows.T == gram
    ):
        raise ArithmeticError("A1,...,A14 failed to reconstruct U+E8+D4")
    return basis_gram, wall_rows


def _lift_actions(
    group: Sequence[Permutation],
    wall_rows: sp.Matrix,
    basis_gram: sp.Matrix,
    *,
    rank: int = RANK,
    wall_count: int = WALL_COUNT,
) -> dict[Permutation, sp.Matrix]:
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


def _conjugate_subgroup(
    subgroup: Subgroup, action: Permutation
) -> Subgroup:
    inverse = [0] * len(action)
    for source, target in enumerate(action):
        inverse[target] = source
    inverse_tuple = tuple(inverse)
    return frozenset(
        _compose(_compose(action, element), inverse_tuple) for element in subgroup
    )


def _subgroup_classes(
    group: tuple[Permutation, ...], *, wall_count: int = WALL_COUNT
) -> list[dict[str, Any]]:
    remaining = set(_all_subgroups(group))
    records = []
    while remaining:
        representative = min(remaining, key=lambda subgroup: (len(subgroup), sorted(subgroup)))
        conjugates = {_conjugate_subgroup(representative, action) for action in group}
        remaining.difference_update(conjugates)
        generators = []
        generated: Subgroup = frozenset({_identity(wall_count)})
        for element in sorted(representative):
            if element in generated:
                continue
            generators.append(element)
            generated = _generated_subgroup(generators)
        records.append(
            {
                "order": len(representative),
                "conjugacy_class_size": len(conjugates),
                "element_order_histogram": {
                    str(order): count
                    for order, count in sorted(
                        Counter(_order(element) for element in representative).items()
                    )
                },
                "representative_generators": tuple(generators),
                "representative_subgroup": representative,
            }
        )
    records.sort(key=lambda row: (row["order"], row["element_order_histogram"]))
    for index, record in enumerate(records, start=1):
        record["class_index"] = index
    return records


def _divisibility(gram: sp.Matrix, root: sp.Matrix) -> int:
    result = 0
    for value in root.T * gram:
        result = gcd(result, abs(int(value)))
    return result


def _fixed_data(
    generators: Sequence[Permutation],
    subgroup: Subgroup,
    actions: dict[Permutation, sp.Matrix],
    basis_gram: sp.Matrix,
    wall_rows: sp.Matrix,
    wall_gram: sp.Matrix,
    *,
    rank: int = RANK,
    wall_count: int = WALL_COUNT,
    surface_realization_status: str | None = None,
) -> dict[str, Any]:
    fixed_basis, kernel = _saturated_fixed_basis(
        [actions[generator] for generator in generators], rank
    )
    fixed_gram = fixed_basis * basis_gram * fixed_basis.T
    cycles_zero = _orbits(subgroup, wall_count)
    roots = []
    reflection_checks = []
    for orbit_index, orbit in enumerate(cycles_zero, start=1):
        normal = sum((wall_rows.row(index) for index in orbit), sp.zeros(1, rank))
        root = sp.Matrix(tuple(next(iter(sp.linsolve((fixed_basis.T, normal.T))))))
        if not (root.T * fixed_basis == normal and all(value.is_Integer for value in root)):
            raise ArithmeticError("orbit root is not integral in the fixed lattice")
        norm = int((root.T * fixed_gram * root)[0])
        reflection = _reflection(fixed_gram, root)
        disjoint = _pairwise_disjoint(wall_gram, tuple(orbit))
        record = {
            "orbit_index": orbit_index,
            "source_walls_one_based": [index + 1 for index in orbit],
            "fixed_coordinates": [int(value) for value in root],
            "norm": norm,
            "divisibility": _divisibility(fixed_gram, root),
            "pairwise_disjoint_source_roots": disjoint,
            "integral_reflection_on_fixed_lattice": (
                disjoint
                and norm == -2 * len(orbit)
                and all(value.is_Integer for value in reflection)
                and reflection * reflection == sp.eye(fixed_gram.rows)
                and reflection.T * fixed_gram * reflection == fixed_gram
            ),
        }
        roots.append(root)
        reflection_checks.append(record)
    root_matrix = sp.Matrix.hstack(*roots)
    orbit_gram = root_matrix.T * fixed_gram * root_matrix
    if orbit_gram != _orbit_gram(wall_gram, cycles_zero):
        raise ArithmeticError("fixed-lattice and orbit Grams disagree")
    reflection_audit = {
        "all_folded_reflections_are_integral": all(
            row["integral_reflection_on_fixed_lattice"] for row in reflection_checks
        ),
        "checks": reflection_checks,
    }
    ns_checks = {
        "actions_integral": all(
            value.is_Integer for generator in generators for value in actions[generator]
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
            [value + 1 for value in generator] for generator in generators
        ],
        "generator_action_matrices": [
            _matrix_payload(actions[generator]) for generator in generators
        ],
        "fixed_rank": fixed_gram.rows,
        "checks": ns_checks,
        "all_integral_lattice_checks_pass": all(ns_checks.values()),
        "surface_realization_status": surface_realization_status
        or (
            "not_realized_on_the_general_K3: Roulleau records Aut(X)=C2 and "
            "tag t, so its nontrivial hyperelliptic involution acts trivially "
            "on NS; every nontrivial S3 diagram action here acts nontrivially"
        ),
    }
    return {
        "cycles_zero": cycles_zero,
        "cycles_one": [[value + 1 for value in orbit] for orbit in cycles_zero],
        "fixed_basis": fixed_basis,
        "fixed_gram": fixed_gram,
        "fixed_kernel": kernel,
        "roots": roots,
        "root_matrix": root_matrix,
        "orbit_gram": orbit_gram,
        "reflection_audit": reflection_audit,
        "ns_action": ns_action,
    }


def _stable_lineage(
    fixed_gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    expected: dict[str, Any],
) -> dict[str, Any]:
    reflections = [_reflection(fixed_gram, root) for root in roots]
    audit = _stable_subgroups(fixed_gram, reflections)
    groups = audit.pop("groups_internal")
    left_inverse = audit.pop("left_smith_inverse")
    if not (
        audit["stable_integral_subgroup_orders"] == expected["stable_orders"]
        and audit["stable_even_subgroup_orders"] == expected["stable_even_orders"]
        and len(audit["maximal_stable_subgroup_positions_zero_based"]) == 1
    ):
        raise ArithmeticError("stable-overlattice denominator changed")
    position = audit["maximal_stable_subgroup_positions_zero_based"][0]
    basis = _overlattice_basis(
        fixed_gram,
        left_inverse,
        audit["stable_subgroups"][position]["cosets_in_smith_coordinates"],
    )
    maximal_gram = sp.simplify(basis.T * fixed_gram * basis)
    actions = [sp.simplify(basis.inv() * reflection * basis) for reflection in reflections]
    parity = (
        "even"
        if all(maximal_gram[index, index] % 2 == 0 for index in range(maximal_gram.rows))
        else "odd"
    )
    checks = {
        "index": abs(sp.Rational(1, basis.det())) == expected["maximal_index"],
        "determinant": maximal_gram.det() == expected["maximal_determinant"],
        "parity": parity == expected["maximal_parity"],
        "actions_integral": all(value.is_Integer for action in actions for value in action),
        "actions_preserve_gram": all(
            action.T * maximal_gram * action == maximal_gram for action in actions
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"maximal stable lattice failed: {checks}")
    return {
        "audit": audit,
        "unique_maximal": {
            "index_over_fixed_lattice": len(groups[position]),
            "basis_columns_in_fixed_lattice": _payload(basis),
            "gram": _payload(maximal_gram),
            "determinant": int(maximal_gram.det()),
            "smith_invariants": _smith(maximal_gram),
            "parity": parity,
            "discriminant_q_histogram_mod_two": (
                _q_histogram(maximal_gram) if parity == "even" else None
            ),
            "facet_reflection_actions": [_payload(action) for action in actions],
            "facet_reflection_actions_sha256": _canonical_sha256(
                [_payload(action) for action in actions]
            ),
            "classification": expected["classification"],
            "checks": checks,
        },
    }


def build_certificate() -> dict[str, Any]:
    if _sha256(SOURCE_TEX) != SOURCE_TEX_SHA256:
        raise ArithmeticError("pinned Roulleau TeX hash changed")
    inventory = json.loads(SOURCE_INVENTORY.read_text())
    stock = next(row for row in inventory["stocks"] if row["key"] == "R14_UE8D4")
    excerpt = "\n".join(SOURCE_TEX.read_text().splitlines()[9585:9668])
    source_checks = {
        "inventory_certified": inventory["all_checks_pass"] is True,
        "stock_rank": stock["picard_rank"] == 14,
        "stock_walls": stock["wall_count"] == 15,
        "section_token": r"U\oplus\mathbf{E}_{8}\oplus\mathbf{D}_{4}" in excerpt,
        "wall_count_token": "contains $15$ $(-2)$-curves" in excerpt,
        "basis_token": "$A_{1},\\dots,A_{14}$ generate" in excerpt,
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
    actions = _lift_actions(group, wall_rows, basis_gram)
    classes_internal = _subgroup_classes(group)
    if [row["order"] for row in classes_internal] != [1, 2, 3, 6]:
        raise ArithmeticError("S3 subgroup conjugacy denominator changed")

    class_records = []
    high_internal = []
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
            "representative_generators_zero_based": [list(g) for g in generators],
        }
        class_records.append(public)
        if len(subgroup) > 1:
            data = _fixed_data(
                generators, subgroup, actions, basis_gram, wall_rows, wall_gram
            )
            high_internal.append({**public, "_data": data})

    type_specs = [
        ("R14_UE8D4_C2", [2]),
        ("R14_UE8D4_C3_S3", [3, 4]),
    ]
    types = []
    for identifier, class_indices in type_specs:
        expected = EXPECTED_TYPES[identifier]
        representative = next(
            row for row in high_internal if row["class_index"] == class_indices[0]
        )
        data = representative["_data"]
        orbit_gram = data["orbit_gram"]
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
        positive = [
            [str(-orbit_gram[left, right]) for right in range(orbit_gram.rows)]
            for left in range(orbit_gram.rows)
        ]
        relations = _relation_matrix(positive)
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
        stable = _stable_lineage(data["fixed_gram"], data["roots"], expected)
        sk_line = (Path("/private/tmp/Scharlau-Kirschmer-Lorentz") / expected["sk_list"]).read_text().splitlines()[expected["sk_entry"]]
        ordinary_paths = sorted(row["path"] for row in coxiter["ordinary_exact_matches"])
        checks = {
            "class_partition": class_indices == expected["class_indices"],
            "fixed_lattice": (
                int(data["fixed_gram"].det()) == expected["fixed_determinant"]
                and _smith(data["fixed_gram"]) == expected["fixed_smith"]
                and all(data["fixed_gram"][i, i] % 2 == 0 for i in range(data["fixed_gram"].rows))
            ),
            "geometry": geometry["all_checks_pass"],
            "volume_counts": (
                geometry["hyperbolic_dimension"] == expected["dimension"]
                and geometry["facet_count"] == expected["facets"]
                and geometry["vinberg_volume_criterion"]["elliptic_edge_count"] == expected["edges"]
                and geometry["finite_vertex_count"] == expected["finite_vertices"]
                and geometry["ideal_vertex_count"] == expected["cusps"]
                and geometry["vinberg_volume_criterion"]["parabolic_neighborhood_audit_count"] == expected["parabolic_audits"]
            ),
            "relations": relation_histogram == expected["relations"],
            "solomon": solomon["rational_euler_characteristic"] == expected["euler"],
            "coxiter_exact_matches": ordinary_paths == sorted(expected["coxiter_matches"])
            and not coxiter["doubling_exact_matches"],
            "sk_row": expected["sk_symbol"] in sk_line,
            "full_chamber_index_one": (
                expected["facets"] == len(relations)
                and bool(coxiter["ordinary_exact_matches"])
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(f"type checks failed for {identifier}: {checks}")
        types.append(
            {
                "type_id": identifier,
                "subgroup_class_indices_one_based": class_indices,
                "subgroup_structures": [
                    next(row for row in class_records if row["class_index"] == index)[
                        "element_order_histogram"
                    ]
                    for index in class_indices
                ],
                "fixed_lattice": {
                    "basis_rows_in_A1_to_A14_basis": _matrix_payload(data["fixed_basis"]),
                    "gram": _matrix_payload(data["fixed_gram"]),
                    "rank": data["fixed_gram"].rows,
                    "determinant": int(data["fixed_gram"].det()),
                    "smith_invariants": _smith(data["fixed_gram"]),
                    "even": True,
                    "smith_kernel_audit": data["fixed_kernel"],
                    "orbit_root_coordinate_columns": _matrix_payload(data["root_matrix"]),
                },
                "reflection_audit": data["reflection_audit"],
                "ns_action": data["ns_action"],
                "geometry": geometry,
                "normalized_relation_matrix": relations,
                "relation_histogram": relation_histogram,
                "solomon": solomon,
                "stable_overlattice_lineage": stable,
                "maximal_reflection_lineage": {
                    "scharlau_kirschmer_list": expected["sk_list"],
                    "entry_one_based": expected["sk_entry"],
                    "genus_symbol": expected["sk_symbol"],
                    "list_row": sk_line.strip(),
                    "coxiter": coxiter,
                    "reflection_subgroup_index": 1,
                    "status": "known_exact_full_reflection_chamber",
                },
                "checks": checks,
                "all_checks_pass": all(checks.values()),
            }
        )

    checks = {
        "source": all(source_checks.values()),
        "wall_gram": wall_gram.rank() == 14 and exact_inertia(wall_gram) == (1, 13, 1),
        "diagram_group_S3": len(group) == 6,
        "four_subgroup_classes": len(class_records) == 4,
        "three_proper_high_classes": len(high_internal) == 3,
        "two_index_one_types": len(types) == 2 and all(row["all_checks_pass"] for row in types),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global checks failed: {checks}")
    return {
        "schema": SCHEMA,
        "artifact_id": "roulleau_ue8d4_all_subgroup_folds",
        "source": {
            "stock_key": "R14_UE8D4",
            "lattice": "U+E8(-1)+D4(-1)",
            "tex_section_lines": [9586, 9668],
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
            "general_K3_surface_automorphism_status": (
                "Aut(X)=C2 with Roulleau tag t; the nontrivial surface "
                "automorphism acts trivially on NS, so none of the nontrivial "
                "S3 diagram actions is realized on the general surface"
            ),
        },
        "subgroup_census": {
            "conjugacy_class_count": len(class_records),
            "proper_high_class_count": len(high_internal),
            "classes": class_records,
        },
        "types": types,
        "summary": {
            "proper_high_subgroup_class_count": 3,
            "distinct_fixed_cone_count": 2,
            "finite_volume_noncompact_count": 2,
            "known_abstract_full_chamber_count": 2,
            "new_abstract_polytope_count": 0,
            "algebro_geometric_divisor_lattice_fold_priority": "open",
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
