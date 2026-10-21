"""Compare the old and orientation-first full disc16 wall inventories exactly.

The original complete weighted-wall inventory was built from Shimada's raw
mirror rows.  The determinant-parameterized successor first orients every
wall toward the published Weyl chamber, resolving zero Weyl pairings from the
archived adjacent chamber, and only then computes the exact weighted Gram
automorphism group.  This audit proves that the stricter normalization does
not alter the fishing denominator or reveal an unscreened high-rank class.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OLD_INVENTORY = (
    ROOT / "data" / "enumerations" / "shimada_disc16_weighted_wall_symmetry_inventory.json"
)
ORIENTED_INVENTORY = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_oriented_weighted_wall_symmetry_inventory.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_oriented_inventory_normalization_audit.json"
)
SCHEMA = "coxeter.shimada_disc16_oriented_inventory_normalization_audit.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _class_signature(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        int(record["order"]),
        int(record["fixed_rank"]),
        int(record["conjugacy_class_size"]),
        tuple(
            sorted(
                (int(key), int(value))
                for key, value in record["element_order_histogram"].items()
            )
        ),
    )


def _action_signature(record: dict[str, Any], *, old: bool) -> tuple[Any, ...]:
    matrix_key = "dual_action" if old else "dual_row_action"
    return (
        tuple(map(int, record["permutation"])),
        tuple(tuple(row) for row in record[matrix_key]),
        int(record["order"]),
        bool(record["primal_action_integral"]),
    )


def main() -> None:
    old_payload = json.loads(OLD_INVENTORY.read_text())
    oriented_payload = json.loads(ORIENTED_INVENTORY.read_text())
    old_all = {
        int(record["chamber_index"]): record for record in old_payload["all_chambers"]
    }
    old_nontrivial = {
        int(record["chamber_index"]): record
        for record in old_payload["nontrivial_weighted_group_chambers"]
    }
    oriented = {
        int(record["chamber_index"]): record for record in oriented_payload["chambers"]
    }
    if set(old_all) != set(oriented) or len(oriented) != 4_539:
        raise ArithmeticError("the two inventories do not cover the same 4,539 chambers")

    summary_mismatches = []
    class_mismatches = []
    strict_action_mismatches = []
    strict_indices = []
    for chamber_index in sorted(oriented):
        old_summary = old_all[chamber_index]
        new_record = oriented[chamber_index]
        new_archived = new_record["archived_AutD"]
        new_group = new_record["oriented_weighted_wall_group"]
        comparable = {
            "wall_count": (
                int(old_summary["wall_count"]),
                int(new_record["wall_count"]),
            ),
            "archived_order": (
                int(old_summary["archived_stabilizer_order"]),
                int(new_archived["order"]),
            ),
            "weighted_order": (
                int(old_summary["weighted_wall_group_order"]),
                int(new_group["order"]),
            ),
            "groups_equal": (
                bool(old_summary["groups_equal"]),
                not bool(new_group["strictly_extends_archived_AutD"]),
            ),
        }
        failed = {
            key: values for key, values in comparable.items() if values[0] != values[1]
        }
        if failed:
            summary_mismatches.append({"chamber_index": chamber_index, "failed": failed})

        if int(new_group["order"]) > 1:
            old_record = old_nontrivial[chamber_index]
            old_classes = Counter(_class_signature(row) for row in old_record["classes"])
            new_classes = Counter(_class_signature(row) for row in new_group["classes"])
            if old_classes != new_classes:
                class_mismatches.append(
                    {
                        "chamber_index": chamber_index,
                        "old": [[list(key), value] for key, value in old_classes.items()],
                        "oriented": [[list(key), value] for key, value in new_classes.items()],
                    }
                )
            if not old_summary["groups_equal"]:
                strict_indices.append(chamber_index)
                old_actions = {
                    _action_signature(row, old=True)
                    for row in old_record["weighted_group_elements"]
                }
                new_actions = {
                    _action_signature(row, old=False)
                    for row in new_group["all_actions"]
                }
                if old_actions != new_actions:
                    strict_action_mismatches.append(chamber_index)

    old_summary = old_payload["summary"]
    new_summary = oriented_payload["summary"]
    high_dimensions = {
        str(key): value
        for key, value in sorted(
            Counter(
                int(row["hyperbolic_dimension"])
                for chamber in oriented.values()
                for row in chamber["oriented_weighted_wall_group"]["classes"]
                if int(row["order"]) > 1 and int(row["fixed_rank"]) >= 11
            ).items()
        )
    }
    checks = {
        "both_inventories_complete": bool(oriented_payload["complete"])
        and len(old_all) == 4_539,
        "same_source_hash": old_payload["source"]["sha256"]
        == oriented_payload["source"]["sha256"],
        "all_chamber_summaries_equal": not summary_mismatches,
        "all_nontrivial_class_invariant_multisets_equal": not class_mismatches,
        "all_strict_extension_exact_action_sets_equal": not strict_action_mismatches,
        "same_group_order_histogram": old_summary["weighted_group_order_histogram"]
        == new_summary["oriented_weighted_group_order_histogram"],
        "same_strict_extension_count": old_summary["weighted_group_strict_extension_count"]
        == new_summary["weighted_group_strict_extension_count"]
        == 575,
        "same_subgroup_count": old_summary[
            "total_subgroup_count_in_nontrivial_weighted_groups"
        ]
        == new_summary["weighted_total_subgroup_count_nontrivial_groups"],
        "same_subgroup_class_count": old_summary[
            "total_subgroup_conjugacy_class_count"
        ]
        == new_summary["weighted_total_subgroup_conjugacy_class_count"],
        "same_high_dimension_histogram": high_dimensions
        == old_summary["proper_high_dimension_histogram"]
        == new_summary[
            "proper_weighted_fixed_rank_at_least_11_dimension_histogram"
        ],
        "high_dimension_class_count_is_682": sum(high_dimensions.values()) == 682,
        "strict_extension_count_is_575": len(strict_indices) == 575,
        "zero_Weyl_orbits_resolved": new_summary[
            "zero_published_Weyl_wall_orbit_count"
        ]
        == 105,
    }
    failed_checks = [key for key, value in checks.items() if not value]
    if failed_checks:
        raise ArithmeticError(f"normalization audit failed: {failed_checks}")

    payload = {
        "schema": SCHEMA,
        "scope": {
            "statement": "orientation-before-symmetry changes neither the complete disc16 group denominator nor any high-rank subgroup class",
            "facet_angle_or_volume_recomputed": False,
            "reason_no_rescreen_is_needed": "all 107 archived classes are the already screened archived actions, while all 575 strict extensions have exact action sets identical to the already screened fresh inventory",
        },
        "sources": {
            "raw_weighted_inventory": {
                "path": str(OLD_INVENTORY.relative_to(ROOT)),
                "sha256": _sha256(OLD_INVENTORY),
            },
            "oriented_weighted_inventory": {
                "path": str(ORIENTED_INVENTORY.relative_to(ROOT)),
                "sha256": _sha256(ORIENTED_INVENTORY),
            },
            "primary_atlas_sha256": oriented_payload["source"]["sha256"],
        },
        "summary": {
            "chamber_count": len(oriented),
            "nontrivial_group_chamber_count": sum(
                chamber["oriented_weighted_wall_group"]["order"] > 1
                for chamber in oriented.values()
            ),
            "strict_extension_chamber_count": len(strict_indices),
            "proper_high_dimension_class_count": sum(high_dimensions.values()),
            "proper_high_dimension_histogram": high_dimensions,
            "zero_published_Weyl_wall_orbit_count": 105,
            "primal_action_denominator_histogram": new_summary[
                "weighted_primal_action_denominator_histogram"
            ],
            "all_checks_pass": all(checks.values()),
        },
        "checks": checks,
        "mismatches": {
            "chamber_summaries": summary_mismatches,
            "class_invariant_multisets": class_mismatches,
            "strict_extension_action_sets": strict_action_mismatches,
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
