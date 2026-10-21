"""Close the full weighted-diagram subgroup census in K3 ranks 15--17.

The published K3 wall stocks in ``k3_high_rank_nef_chambers.json`` are
complete nef-wall configurations.  Earlier work exhausted only their
involutions.  Here we enumerate every element conjugacy class and every
subgroup conjugacy class of each exact weighted wall-Gram automorphism group.

Only proper fixed spaces of rank at least 11 enter the geometric queue.  The
two involution classes already have frozen certificates and are referenced,
not recomputed.  The only non-involution class newly exposed by the full
subgroup denominator is an S3 in rank 15.  Its wall-orbit partition and orbit
Gram are proved identical to those of the already certified C3 triality fold,
so it gives the same H12, 17-wall chamber rather than a new polytope.

Diagram/lattice isometries and automorphisms of a chosen K3 surface are kept
logically separate throughout.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import sympy as sp

from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    _ambient_reconstruction,
    _gap_subgroup_classes,
    _greedy_generators,
    _lift_actions,
)
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _fixed_rank,
    _generated_subgroup,
    _orbits,
    _permutation_order,
)
from scripts.search_k3_high_rank_involution_folds import (
    _compose,
    _conjugacy_classes,
    _group_checks,
    _identity,
    _inverse,
    _weighted_graph_automorphisms,
)
from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
INVOLUTION_ARCHIVE = (
    ROOT / "data/enumerations/k3_high_rank_involution_folds_rho15_17.json"
)
ORDER3_INDEPENDENT = ROOT / "data/candidates/k3_rho15_h12_cyclic_fold_audit.json"
ORDER3_PRODUCTION = ROOT / "data/enumerations/k3_rank15_order3_h12_fold.json"
OUTPUT = ROOT / "data/enumerations/k3_high_rank_all_subgroup_completion.json"

EXPECTED_HASHES = {
    "data/known_examples/k3_high_rank_nef_chambers.json": (
        "0dcc7ec9c10f3437cd39713cd0d6036dc4fa8ad0c978ffdad4c2b34748269120"
    ),
    "data/enumerations/k3_high_rank_involution_folds_rho15_17.json": (
        "df6aa10a21faaca264fc120c2459ecc36d5e727b87c6e84411d6b0b47ef345c5"
    ),
    "data/candidates/k3_rho15_h12_cyclic_fold_audit.json": (
        "a93178362a063fc96c721e3dfd721655cf2093cde3ce7cdb4ef60ccc4e22a28d"
    ),
    "data/enumerations/k3_rank15_order3_h12_fold.json": (
        "54fea0698d48241cd50d6654ccc7ad687e227e9fa9d8908e011e7455782e097b"
    ),
}

TARGET_RHOS = (15, 16, 17)
MIN_FIXED_RANK = 11
EXPECTED = {
    15: {
        "wall_count": 21,
        "group_name": "S3 x C2",
        "group_order": 12,
        "element_order_histogram": {1: 1, 2: 7, 3: 2, 6: 2},
        "element_class_count": 6,
        "subgroup_class_count": 10,
    },
    16: {
        "wall_count": 19,
        "group_name": "C2 x C2",
        "group_order": 4,
        "element_order_histogram": {1: 1, 2: 3},
        "element_class_count": 4,
        "subgroup_class_count": 5,
    },
    17: {
        "wall_count": 19,
        "group_name": "C2",
        "group_order": 2,
        "element_order_histogram": {1: 1, 2: 1},
        "element_class_count": 2,
        "subgroup_class_count": 2,
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json_digest(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _one_based_orbits(orbits: Iterable[Iterable[int]]) -> list[list[int]]:
    return [[int(value) + 1 for value in orbit] for orbit in orbits]


def _orbit_gram(
    wall_gram: sp.Matrix, orbits: tuple[tuple[int, ...], ...]
) -> sp.Matrix:
    incidence = sp.zeros(wall_gram.rows, len(orbits))
    for column, orbit in enumerate(orbits):
        for wall in orbit:
            incidence[wall, column] = 1
    return incidence.T * wall_gram * incidence


def _conjugate_subgroup(
    action: Permutation, subgroup: Subgroup
) -> Subgroup:
    inverse = _inverse(action)
    return frozenset(
        _compose(_compose(action, element), inverse) for element in subgroup
    )


def _subgroup_conjugates(
    group: tuple[Permutation, ...], subgroup: Subgroup
) -> tuple[Subgroup, ...]:
    return tuple(
        sorted({_conjugate_subgroup(action, subgroup) for action in group})
    )


def _group_type(subgroup: Subgroup) -> str:
    histogram = Counter(_permutation_order(element) for element in subgroup)
    signature = (len(subgroup), tuple(sorted(histogram.items())))
    names = {
        (1, ((1, 1),)): "1",
        (2, ((1, 1), (2, 1))): "C2",
        (3, ((1, 1), (3, 2))): "C3",
        (4, ((1, 1), (2, 3))): "C2 x C2",
        (4, ((1, 1), (2, 1), (4, 2))): "C4",
        (6, ((1, 1), (2, 1), (3, 2), (6, 2))): "C6",
        (6, ((1, 1), (2, 3), (3, 2))): "S3",
        (12, ((1, 1), (2, 7), (3, 2), (6, 2))): "S3 x C2",
    }
    if signature not in names:
        raise ArithmeticError(f"unrecognized small subgroup signature {signature}")
    return names[signature]


def _surface_status(rho: int, group_type: str, fixed_rank: int) -> dict[str, str]:
    if rho == 15 and group_type in {"C3", "S3"} and fixed_rank == 13:
        return {
            "status": "diagram_and_integral_NS_lattice_action_only",
            "general_surface": "not_realized_by_the_published_general_K3_action",
            "special_period": "unassessed",
            "reason": (
                "Kondo's general rank-15 surface realizes only the central C2 "
                "of the weighted graph on NS; the S3 factor is not its NS action."
            ),
        }
    if rho == 15 and group_type == "C2" and fixed_rank == 14:
        return {
            "status": "diagram_and_integral_NS_lattice_action_only",
            "general_surface": "not_realized_by_the_published_general_K3_action",
            "special_period": "unassessed",
            "reason": "This is the already archived noncentral S3-transposition class.",
        }
    if rho == 16 and group_type == "C2" and fixed_rank == 15:
        return {
            "status": "diagram_and_integral_NS_lattice_action_only",
            "general_surface": "excluded_as_an_involution_for_this_NS_lattice",
            "special_period": "higher_order_realization_not_classified",
            "reason": "This is the already archived tau_2 graph fold.",
        }
    return {
        "status": "not_reclassified_here",
        "general_surface": "see_the_frozen_involution_inventory_or_Kondo",
        "special_period": "not_applicable_to_this_rank_gate_record",
        "reason": "Surface realization is not used in the subgroup denominator.",
    }


def _element_classes(
    group: tuple[Permutation, ...],
    actions: dict[Permutation, sp.Matrix],
    ambient_rank: int,
    wall_count: int,
) -> list[dict[str, Any]]:
    records = []
    for class_index, conjugacy_class in enumerate(_conjugacy_classes(group, group)):
        representative = conjugacy_class[0]
        cyclic_subgroup = _generated_subgroup((representative,))
        fixed_rank = _fixed_rank((actions[representative],), ambient_rank)
        orbit_partition = _orbits(cyclic_subgroup, wall_count)
        records.append(
            {
                "class_index": class_index,
                "element_order": _permutation_order(representative),
                "conjugacy_class_size": len(conjugacy_class),
                "representative_images_zero_based": list(representative),
                "cyclic_subgroup_order": len(cyclic_subgroup),
                "fixed_rank": fixed_rank,
                "hyperbolic_dimension": fixed_rank - 1,
                "proper_fixed_rank_at_least_11": (
                    representative != _identity(wall_count)
                    and fixed_rank >= MIN_FIXED_RANK
                ),
                "cyclic_wall_orbits_one_based": _one_based_orbits(orbit_partition),
            }
        )
    return records


def _subgroup_classes(
    group: tuple[Permutation, ...],
    actions: dict[Permutation, sp.Matrix],
    wall_gram: sp.Matrix,
    ambient_rank: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    generators = _greedy_generators(group)
    gap_rows = _gap_subgroup_classes(group, generators)
    records = []
    high_records = []
    for class_index, gap_row in enumerate(gap_rows):
        subgroup: Subgroup = gap_row["representative_subgroup"]
        representative_generators: tuple[Permutation, ...] = gap_row[
            "representative_generators"
        ]
        matrices = tuple(actions[element] for element in representative_generators)
        fixed_rank = _fixed_rank(matrices, ambient_rank) if matrices else ambient_rank
        orbits = _orbits(subgroup, wall_gram.rows)
        orbit_gram = _orbit_gram(wall_gram, orbits)
        if orbit_gram.rank() != fixed_rank:
            raise ArithmeticError("wall-orbit span does not equal the fixed subspace")
        group_type = _group_type(subgroup)
        cyclic = any(
            _permutation_order(element) == len(subgroup) for element in subgroup
        )
        proper_high = len(subgroup) > 1 and fixed_rank >= MIN_FIXED_RANK
        conjugates = _subgroup_conjugates(group, subgroup)
        if len(conjugates) != int(gap_row["conjugacy_class_size"]):
            raise ArithmeticError("GAP subgroup class size disagrees with direct conjugation")
        conjugate_partitions = sorted(
            {
                tuple(_orbits(conjugate, wall_gram.rows))
                for conjugate in conjugates
            }
        )
        record = {
            "class_index": class_index,
            "group_type": group_type,
            "order": len(subgroup),
            "conjugacy_class_size": int(gap_row["conjugacy_class_size"]),
            "cyclic": cyclic,
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "proper_fixed_rank_at_least_11": proper_high,
            "element_order_histogram": {
                str(order): count
                for order, count in sorted(
                    Counter(
                        _permutation_order(element) for element in subgroup
                    ).items()
                )
            },
            "representative_generators_zero_based": [
                list(element) for element in representative_generators
            ],
            "representative_elements_zero_based": [
                list(element) for element in sorted(subgroup)
            ],
            "representative_wall_orbits_one_based": _one_based_orbits(orbits),
            "conjugate_wall_orbit_partitions_one_based": [
                _one_based_orbits(partition) for partition in conjugate_partitions
            ],
            "orbit_gram_rank": orbit_gram.rank(),
            "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
            "orbit_gram_sha256": _json_digest(_matrix_payload(orbit_gram)),
        }
        records.append(record)
        if proper_high:
            high_records.append(
                {
                    **record,
                    "surface_realization": _surface_status(
                        ambient_rank, group_type, fixed_rank
                    ),
                    "_subgroup": subgroup,
                    "_orbit_gram": orbit_gram,
                }
            )
    return records, high_records


def _load_and_verify_sources() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    paths = (SOURCE, INVOLUTION_ARCHIVE, ORDER3_INDEPENDENT, ORDER3_PRODUCTION)
    for path in paths:
        relative = str(path.relative_to(ROOT))
        digest = _sha256(path)
        if digest != EXPECTED_HASHES[relative]:
            raise ArithmeticError(f"unexpected SHA-256 for {relative}: {digest}")
    return tuple(json.loads(path.read_text()) for path in paths)  # type: ignore[return-value]


def build() -> dict[str, Any]:
    source, involution_archive, order3_independent, order3_production = (
        _load_and_verify_sources()
    )
    source_cases = {int(case["rho"]): case for case in source["cases"]}
    old_geometries = {
        int(case["rho"]): case["geometry_certificates"]
        for case in involution_archive["cases"]
    }
    known_order3_orbits = order3_independent["order_three_H12_fold"][
        "wall_orbits"
    ]
    known_order3_gram = sp.Matrix(order3_independent["matrices"]["orbit_gram_17"])

    cases = []
    high_queue: list[dict[str, Any]] = []
    internal_high: list[dict[str, Any]] = []
    for rho in TARGET_RHOS:
        case = source_cases[rho]
        expected = EXPECTED[rho]
        wall_gram = build_gram(case)
        if not (
            wall_gram.rows == expected["wall_count"]
            and wall_gram.rank() == rho
            and exact_inertia(wall_gram) == (1, rho - 1, wall_gram.rows - rho)
        ):
            raise ArithmeticError(f"rank-{rho} source wall Gram failed exact checks")
        ambient_gram, wall_rows = _ambient_reconstruction(wall_gram, rho)
        if not (
            int(ambient_gram.det()) == int(case["expected_basis_determinant"])
            and exact_inertia(ambient_gram) == (1, rho - 1, 0)
        ):
            raise ArithmeticError(f"rank-{rho} leading wall basis has wrong lattice")

        group, search_node_count = _weighted_graph_automorphisms(wall_gram)
        order_histogram = Counter(_permutation_order(element) for element in group)
        if not (
            len(group) == expected["group_order"]
            and dict(order_histogram) == expected["element_order_histogram"]
            and all(_group_checks(group).values())
        ):
            raise ArithmeticError(f"rank-{rho} weighted diagram group mismatch")
        actions = _lift_actions(group, wall_rows, ambient_gram)
        element_classes = _element_classes(group, actions, rho, wall_gram.rows)
        subgroup_classes, high = _subgroup_classes(
            group, actions, wall_gram, rho
        )
        if not (
            len(element_classes) == expected["element_class_count"]
            and len(subgroup_classes) == expected["subgroup_class_count"]
        ):
            raise ArithmeticError(f"rank-{rho} class-count mismatch")

        for record in high:
            group_type = record["group_type"]
            partitions = record["conjugate_wall_orbit_partitions_one_based"]
            if rho == 15 and group_type == "C2" and record["fixed_rank"] == 14:
                old_orbits = old_geometries[15][0]["facet_orbits"]
                if old_orbits not in partitions:
                    raise ArithmeticError("rank-15 high involution did not match its archive")
                status = "already_archived_involution_H13_19_facets"
                reference = "results/k3_high_rank_involution_folds.md"
            elif rho == 16 and group_type == "C2" and record["fixed_rank"] == 15:
                old_orbits = old_geometries[16][0]["facet_orbits"]
                if old_orbits not in partitions:
                    raise ArithmeticError("rank-16 high involution did not match its archive")
                status = "already_archived_involution_H14_17_facets"
                reference = "results/k3_high_rank_involution_folds.md"
            elif rho == 15 and group_type == "C3" and record["fixed_rank"] == 13:
                if record["representative_wall_orbits_one_based"] != known_order3_orbits:
                    raise ArithmeticError("rank-15 C3 did not match the triality archive")
                if record["_orbit_gram"] != known_order3_gram:
                    raise ArithmeticError("rank-15 C3 orbit Gram changed")
                status = "already_archived_C3_triality_H12_17_facets"
                reference = "results/k3_rank15_order3_h12_fold.md"
            elif rho == 15 and group_type == "S3" and record["fixed_rank"] == 13:
                if record["representative_wall_orbits_one_based"] != known_order3_orbits:
                    raise ArithmeticError("rank-15 S3 fixed cone differs from C3")
                if record["_orbit_gram"] != known_order3_gram:
                    raise ArithmeticError("rank-15 S3 orbit Gram differs from C3")
                status = "new_subgroup_class_but_same_fixed_cone_as_C3_triality"
                reference = "results/k3_rank15_order3_h12_fold.md"
            else:
                raise ArithmeticError(f"unexpected unscreened high class: rho={rho}, {group_type}")
            public_record = {
                key: value
                for key, value in record.items()
                if not key.startswith("_")
            }
            public_record.update(
                {
                    "rho": rho,
                    "geometric_status": status,
                    "certificate_reference": reference,
                }
            )
            high_queue.append(public_record)
            internal_high.append(record)

        cases.append(
            {
                "rho": rho,
                "ambient_lattice": case["lattice"],
                "wall_count": wall_gram.rows,
                "wall_gram_rank": wall_gram.rank(),
                "wall_gram_inertia": list(exact_inertia(wall_gram)),
                "ambient_basis_gram": _matrix_payload(ambient_gram),
                "ambient_basis_determinant": int(ambient_gram.det()),
                "weighted_diagram_group": {
                    "name": expected["group_name"],
                    "order": len(group),
                    "search_node_count": search_node_count,
                    "checks": _group_checks(group),
                    "all_actions_integral_unimodular_Lorentz": True,
                    "element_order_histogram": {
                        str(order): count
                        for order, count in sorted(order_histogram.items())
                    },
                    "greedy_generators_zero_based": [
                        list(element) for element in _greedy_generators(group)
                    ],
                    "elements_zero_based": [list(element) for element in group],
                },
                "element_conjugacy_classes": element_classes,
                "subgroup_conjugacy_classes": subgroup_classes,
                "summary": {
                    "element_conjugacy_class_count": len(element_classes),
                    "subgroup_conjugacy_class_count": len(subgroup_classes),
                    "proper_cyclic_fixed_rank_at_least_11_count": sum(
                        row["proper_fixed_rank_at_least_11"]
                        for row in element_classes
                    ),
                    "proper_subgroup_fixed_rank_at_least_11_count": len(high),
                },
            }
        )

    if len(high_queue) != 4:
        raise ArithmeticError(f"expected four high subgroup classes, got {len(high_queue)}")
    c3 = next(row for row in internal_high if row["group_type"] == "C3")
    s3 = next(row for row in internal_high if row["group_type"] == "S3")
    if not (
        c3["_subgroup"].issubset(s3["_subgroup"])
        and c3["representative_wall_orbits_one_based"]
        == s3["representative_wall_orbits_one_based"]
        and c3["_orbit_gram"] == s3["_orbit_gram"]
    ):
        raise ArithmeticError("S3/C3 fixed-cone equality failed")

    independent_fold = order3_independent["order_three_H12_fold"]
    independent_literal = order3_independent["literal_vinberg_coxiter_prop_1_10"]
    production_literal = order3_production["literal_coxiter_prop_1_10"]
    h12_replay_checks = {
        "C3_is_a_subgroup_of_the_high_S3": c3["_subgroup"].issubset(s3["_subgroup"]),
        "C3_and_S3_wall_partitions_are_identical": (
            c3["representative_wall_orbits_one_based"]
            == s3["representative_wall_orbits_one_based"]
            == known_order3_orbits
        ),
        "C3_and_S3_orbit_Grams_are_identical_to_archive": (
            c3["_orbit_gram"] == s3["_orbit_gram"] == known_order3_gram
        ),
        "fixed_rank_13_and_H12": (
            c3["fixed_rank"] == s3["fixed_rank"] == 13
            and independent_fold["hyperbolic_dimension"] == 12
        ),
        "all_17_walls_have_strict_irredundancy_witnesses": all(
            audit["irredundant_with_strict_relative_interior"]
            for audit in independent_fold["facet_audits"]
        )
        and len(independent_fold["facet_audits"]) == 17,
        "exact_Coxeter_relation_counts_replayed": (
            independent_fold["coxeter_relation_counts"]
            == {
                "2": 116,
                "3": 12,
                "6": 2,
                "infinity": 3,
                "ultraparallel": 3,
            }
        ),
        "exact_extreme_rays_93_finite_3_ideal_0_hyperideal": (
            independent_fold["finite_vertex_count"],
            independent_fold["ideal_vertex_count"],
            independent_fold["hyperideal_vertex_count"],
        )
        == (93, 3, 0),
        "literal_623_edge_Vinberg_test_has_no_bad_edges": (
            independent_literal["elliptic_edge_count"] == 623
            and independent_literal["extension_histogram"] == {"2": 623}
            and independent_literal["bad_edges"] == []
            and independent_literal["finite_volume_noncompact"] is True
        ),
        "production_literal_certificate_agrees": (
            production_literal["spherical_rank_11_edge_count"] == 623
            and production_literal["spherical_rank_12_vertex_count"] == 93
            and production_literal["parabolic_rank_11_vertex_count"] == 3
            and production_literal["bad_edges"] == []
            and production_literal["finite_covolume"] is True
            and order3_production["all_checks_pass"] is True
        ),
        "independent_certificate_all_checks_pass": all(
            order3_independent["checks"].values()
        ),
    }
    if not all(h12_replay_checks.values()):
        raise ArithmeticError("the archived H12 triality certificate did not replay")

    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "status_date": "2026-08-06",
        "artifact_id": "k3_high_rank_all_subgroup_completion",
        "certificate_kind": (
            "exhaustive_exact_element_and_subgroup_conjugacy_class_completion"
        ),
        "source_hashes": EXPECTED_HASHES,
        "scope": {
            "ambient_picard_ranks": list(TARGET_RHOS),
            "wall_stock": "published complete finite-automorphism K3 nef chambers",
            "group_denominator": (
                "full exact automorphism group of each complete weighted wall Gram"
            ),
            "element_classes": "all conjugacy classes, including identity",
            "subgroup_classes": "all conjugacy classes from GAP",
            "geometric_rank_gate": f"proper fixed rank >= {MIN_FIXED_RANK}",
            "involution_policy": (
                "match and reference the frozen exhaustive involution archive; do not duplicate geometry"
            ),
        },
        "method": {
            "weighted_group": (
                "exact complete-Gram row-colour backtracking, preserving zero and nonzero entries"
            ),
            "lattice_lift": (
                "reconstruct every wall in the leading integral curve basis and verify every diagram action is an integral unimodular Lorentz isometry"
            ),
            "subgroup_exhaustion": "GAP ConjugacyClassesSubgroups plus direct conjugacy-size check",
            "fixed_cone": "wall-orbit sums; orbit-span rank checked against the common fixed rank",
            "geometry_reuse_rule": (
                "an identical complete wall-orbit partition and identical orbit Gram define the same fixed cone, so the existing strict-facet and Vinberg certificates transfer exactly"
            ),
        },
        "cases": cases,
        "proper_fixed_rank_at_least_11_queue": high_queue,
        "rank15_C3_S3_same_fixed_cone_certificate": {
            "wall_orbits_one_based": known_order3_orbits,
            "orbit_gram": _matrix_payload(known_order3_gram),
            "orbit_gram_sha256": _json_digest(_matrix_payload(known_order3_gram)),
            "checks": h12_replay_checks,
            "transferred_geometry": {
                "hyperbolic_dimension": 12,
                "facet_count": 17,
                "finite_vertex_count": 93,
                "ideal_vertex_count": 3,
                "hyperideal_vertex_count": 0,
                "elliptic_edge_count": 623,
                "endpoint_type_histogram": {
                    "0,2": 3,
                    "1,1": 124,
                    "2,0": 496,
                },
                "volume": "finite_volume_noncompact",
            },
            "lineage": (
                "Scharlau--Kirschmer List12/Gram12 entry 6; full 17-wall chamber, index 1"
            ),
        },
        "surface_realization_boundary": {
            "search_output": "exact diagram and integral Neron-Severi lattice actions",
            "not_automatic": "automorphism of a chosen algebraic K3 surface",
            "rank15_general_surface": (
                "Kondo realizes only the central C2 wall-graph action on NS; neither the C3 triality nor the high S3 subgroup is asserted as its automorphism action"
            ),
            "special_periods": "unassessed and unnecessary for the abstract Coxeter-polytope statement",
        },
        "summary": {
            "weighted_group_orders": {"15": 12, "16": 4, "17": 2},
            "element_conjugacy_class_counts": {"15": 6, "16": 4, "17": 2},
            "subgroup_conjugacy_class_counts": {"15": 10, "16": 5, "17": 2},
            "total_subgroup_conjugacy_classes": 17,
            "proper_fixed_rank_at_least_11_subgroup_classes": 4,
            "already_archived_involution_classes": 2,
            "already_archived_C3_class": 1,
            "new_high_subgroup_classes": 1,
            "new_high_subgroup_classes_with_new_fixed_cone": 0,
            "new_hyperbolic_Coxeter_polytopes": 0,
            "closure_statement": (
                "The full rank-15--17 weighted-diagram subgroup denominator is closed. The sole new high subgroup class is S3, and it fixes exactly the known C3 triality H12 chamber."
            ),
        },
    }
    payload["all_checks_pass"] = all(h12_replay_checks.values())
    payload["certificate_payload_sha256"] = _json_digest(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text() != encoded:
            raise SystemExit(f"archive mismatch: {args.output}")
        print(f"archive verified: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
        print(args.output)
    print(json.dumps(payload["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
