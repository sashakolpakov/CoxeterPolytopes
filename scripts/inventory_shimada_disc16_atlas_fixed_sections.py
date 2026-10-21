"""Inventory archived-stabilizer fixed sections in Shimada's disc16 atlas.

This is the cheap, exact first gate for the next algebraic-stock fishing cast.
Shimada's atlas contains 4,539 published induced chambers.  Most have trivial
archived ``Aut(D)`` and hence no proper fixed section in this denominator.
For every chamber with nontrivial archived lattice/chamber stabilizer we
enumerate *all* subgroups, quotient by conjugacy in that archived group, and
compute the fixed rank on the rank-20 Neron--Severi lattice.  No wall-angle or
volume claim is made here.

This is deliberately not the full weighted-wall graph automorphism inventory.
That larger symmetry computation is a separate fishing lane and may expose
additional casts not represented by Shimada's archived ``Aut(D)`` matrices.

The finite matrix group is converted to its faithful regular permutation
action only to reuse the audited subgroup enumerator.  Fixed ranks are always
computed from the original integral 20 by 20 matrices.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    _subgroup_conjugacy_classes,
)
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _all_subgroups,
    _fixed_rank,
    _identity,
    _permutation_order,
)
from scripts.shimada_disc11_atlas import (
    ShimadaDisc11Atlas,
    _matrix_product,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path("/private/tmp/threediscs/disc16.txt")
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_atlas_fixed_section_inventory.json"
)
EXPECTED_SOURCE_SHA256 = (
    "a19b797d6c314bc85e1ac4991fe3b190fe29a18a4746103bc6fa9b05e26a7db6"
)
AMBIENT_RANK = 20
MIN_FIXED_RANK = 11  # hyperbolic dimension fixed_rank - 1 is at least 10


MatrixTuple = tuple[tuple[int, ...], ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _regular_action(
    matrices: tuple[MatrixTuple, ...],
) -> tuple[tuple[Permutation, ...], dict[Permutation, MatrixTuple]]:
    """Return the left regular action and its exact matrix lookup."""

    index = {matrix: position for position, matrix in enumerate(matrices)}
    if len(index) != len(matrices):
        raise ArithmeticError("the archived stabilizer contains duplicates")
    permutations: list[Permutation] = []
    lookup: dict[Permutation, MatrixTuple] = {}
    for left in matrices:
        permutation = tuple(
            index[_matrix_product(left, right)] for right in matrices
        )
        permutations.append(permutation)
        lookup[permutation] = left
    group = tuple(permutations)
    if len(lookup) != len(matrices):
        raise ArithmeticError("the regular action is not faithful")
    identity = _identity(len(matrices))
    if identity not in lookup:
        raise ArithmeticError("the regular action has no identity")
    return group, lookup


def _subgroup_payload(
    subgroup: Subgroup,
    lookup: dict[Permutation, MatrixTuple],
) -> dict[str, Any]:
    matrices = tuple(sp.Matrix(lookup[element]) for element in subgroup)
    fixed_rank = _fixed_rank(matrices, AMBIENT_RANK)
    element_orders = Counter(_permutation_order(element) for element in subgroup)
    return {
        "order": len(subgroup),
        "fixed_rank": fixed_rank,
        "hyperbolic_dimension": fixed_rank - 1,
        "element_order_histogram": {
            str(order): count for order, count in sorted(element_orders.items())
        },
        "regular_action_elements": [list(element) for element in sorted(subgroup)],
    }


def build(source: Path = SOURCE) -> dict[str, Any]:
    source_hash = _sha256(source)
    if source_hash != EXPECTED_SOURCE_SHA256:
        raise ValueError(
            f"unexpected disc16.txt SHA-256 {source_hash}; expected "
            f"{EXPECTED_SOURCE_SHA256}"
        )

    atlas = ShimadaDisc11Atlas(source)
    order_by_chamber = tuple(
        int(atlas.assignment(atlas._chamber_name(index, "orderaut")))
        for index in range(atlas.chamber_count)
    )
    nontrivial = tuple(
        index for index, order in enumerate(order_by_chamber) if order > 1
    )

    chamber_rows: list[dict[str, Any]] = []
    total_subgroups = 0
    total_classes = 0
    high_classes = 0
    proper_classes = 0
    proper_high_classes = 0
    proper_high_chambers: set[int] = set()
    proper_high_dimension_histogram: Counter[int] = Counter()
    fixed_rank_histogram: Counter[int] = Counter()

    for chamber_index in nontrivial:
        matrices = atlas.stabilizer_group(chamber_index)
        expected_order = order_by_chamber[chamber_index]
        if len(matrices) != expected_order:
            raise ArithmeticError("stabilizer order does not replay")
        group, lookup = _regular_action(matrices)
        subgroups = _all_subgroups(group)
        classes = _subgroup_conjugacy_classes(group, subgroups)
        class_rows: list[dict[str, Any]] = []
        for class_index, conjugacy_class in enumerate(classes):
            representative_index = conjugacy_class[0]
            payload = _subgroup_payload(subgroups[representative_index], lookup)
            fixed_rank_histogram[payload["fixed_rank"]] += 1
            above_h9 = payload["fixed_rank"] >= MIN_FIXED_RANK
            proper = payload["order"] > 1
            proper_above_h9 = proper and above_h9
            if above_h9:
                high_classes += 1
            if proper:
                proper_classes += 1
            if proper_above_h9:
                proper_high_classes += 1
                proper_high_chambers.add(chamber_index)
                proper_high_dimension_histogram[
                    payload["hyperbolic_dimension"]
                ] += 1
            class_rows.append(
                {
                    "class_index": class_index,
                    "representative_subgroup_index": representative_index,
                    "subgroup_indices": list(conjugacy_class),
                    "conjugacy_class_size": len(conjugacy_class),
                    **payload,
                    "above_h9": above_h9,
                    "proper_subgroup": proper,
                    "proper_fixed_section_above_h9": proper_above_h9,
                }
            )
        total_subgroups += len(subgroups)
        total_classes += len(classes)
        chamber = atlas.chamber(chamber_index)
        chamber_rows.append(
            {
                "chamber_index": chamber_index,
                "wall_count": chamber.wall_count,
                "stabilizer_order": expected_order,
                "subgroup_count": len(subgroups),
                "subgroup_conjugacy_class_count": len(classes),
                "high_dimensional_class_count": sum(
                    row["above_h9"] for row in class_rows
                ),
                "proper_high_dimensional_class_count": sum(
                    row["proper_fixed_section_above_h9"] for row in class_rows
                ),
                "classes": class_rows,
            }
        )

    stabilizer_histogram = Counter(order_by_chamber)
    checks = {
        "pinned_atlas_has_4539_chambers": atlas.chamber_count == 4_539,
        "exactly_98_chambers_have_archived_Aut_D_order_gt_1": len(nontrivial)
        == 98,
        "one_identity_subgroup_class_per_enumerated_chamber": (
            high_classes - proper_high_classes == len(nontrivial) == 98
        ),
        "all_proper_archived_subgroup_classes_act_nontrivially": all(
            row["fixed_rank"] < AMBIENT_RANK
            for chamber in chamber_rows
            for row in chamber["classes"]
            if row["proper_subgroup"]
        ),
        "proper_above_H9_count_is_107": proper_high_classes == 107,
        "proper_above_H9_dimension_histogram_is_expected": (
            proper_high_dimension_histogram == Counter({11: 101, 10: 4, 14: 2})
        ),
        "all_346_subgroups_enumerated": total_subgroups == 346,
        "all_257_inner_conjugacy_classes_enumerated": total_classes == 257,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"disc16 atlas inventory checks failed: {checks}")

    return {
        "schema": "coxeter.shimada_disc16_atlas_fixed_section_inventory.v1",
        "source": {
            "path": str(source),
            "sha256": source_hash,
            "chamber_count": atlas.chamber_count,
            "ambient_rank": AMBIENT_RANK,
        },
        "scope": {
            "all_published_chambers": True,
            "stabilizer_source": (
                "lattice/chamber stabilizer matrices Aut(D) archived in "
                "Shimada disc16.txt"
            ),
            "all_nontrivial_archived_Aut_D_stabilizers": True,
            "all_subgroups_of_each_archived_Aut_D": True,
            "subgroups_quotiented_by_inner_conjugacy_in_archived_Aut_D": True,
            "full_weighted_wall_graph_automorphism_groups_computed": False,
            "full_weighted_wall_graph_symmetry_inventory_is_separate": True,
            "does_not_claim_archived_Aut_D_is_full_weighted_wall_symmetry": True,
            "does_not_deduplicate_isomorphic_casts_across_chambers": True,
            "minimum_fixed_rank": MIN_FIXED_RANK,
            "minimum_hyperbolic_dimension": MIN_FIXED_RANK - 1,
            "wall_angle_gate_performed": False,
            "volume_gate_performed": False,
        },
        "summary": {
            "chamber_count": atlas.chamber_count,
            "nontrivial_stabilizer_chamber_count": len(nontrivial),
            "stabilizer_order_histogram": {
                str(order): count
                for order, count in sorted(stabilizer_histogram.items())
            },
            "total_subgroup_count_in_nontrivial_chambers": total_subgroups,
            "total_subgroup_conjugacy_class_count": total_classes,
            "high_dimensional_subgroup_class_count": high_classes,
            "identity_ambient_class_count": len(nontrivial),
            "proper_subgroup_conjugacy_class_count": proper_classes,
            "proper_high_dimensional_subgroup_class_count": proper_high_classes,
            "chambers_with_a_proper_class_above_h9": len(proper_high_chambers),
            "proper_high_dimensional_class_histogram": {
                f"H{dimension}": count
                for dimension, count in sorted(
                    proper_high_dimension_histogram.items()
                )
            },
            "class_fixed_rank_histogram": {
                str(rank): count
                for rank, count in sorted(fixed_rank_histogram.items())
            },
        },
        "chambers": chamber_rows,
        "decision": {
            "archived_stabilizer_denominator_complete": True,
            "proper_fixed_section_classes_above_h9": proper_high_classes,
            "requires_wall_angle_reduction_and_volume_audit": True,
            "full_weighted_wall_graph_denominator_complete": False,
            "reason_full_graph_denominator_remains_open": (
                "the full weighted-wall graph automorphism groups are being "
                "inventoried separately and can be larger than archived Aut(D)"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build(arguments.source)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(json.dumps(payload["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
