from __future__ import annotations

import json

from scripts.certify_shimada_disc16_oriented_inventory_normalization import (
    OLD_INVENTORY,
    ORIENTED_INVENTORY,
    OUTPUT,
    SCHEMA,
    _sha256,
)


def _load(path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_orientation_first_inventory_has_the_complete_expected_denominator() -> None:
    inventory = _load(ORIENTED_INVENTORY)
    assert inventory["complete"] is True
    summary = inventory["summary"]
    assert summary["completed_chamber_count"] == 4539
    assert summary["pending_chamber_count"] == 0
    assert summary["oriented_weighted_group_order_histogram"] == {
        "1": 3866,
        "2": 659,
        "3": 4,
        "4": 2,
        "6": 5,
        "8": 1,
        "20": 1,
        "48": 1,
    }
    assert summary["weighted_group_strict_extension_count"] == 575
    assert summary[
        "proper_weighted_fixed_rank_at_least_11_dimension_histogram"
    ] == {"10": 4, "11": 101, "13": 1, "14": 2, "18": 574}
    assert summary["weighted_primal_action_denominator_histogram"] == {
        "1": 4735,
        "2": 574,
        "7": 1,
    }
    assert summary["zero_published_Weyl_wall_orbit_count"] == 105


def test_exact_normalization_audit_has_no_mismatch() -> None:
    audit = _load(OUTPUT)
    assert audit["schema"] == SCHEMA
    assert audit["summary"] == {
        "all_checks_pass": True,
        "chamber_count": 4539,
        "nontrivial_group_chamber_count": 673,
        "primal_action_denominator_histogram": {"1": 4735, "2": 574, "7": 1},
        "proper_high_dimension_class_count": 682,
        "proper_high_dimension_histogram": {
            "10": 4,
            "11": 101,
            "13": 1,
            "14": 2,
            "18": 574,
        },
        "strict_extension_chamber_count": 575,
        "zero_published_Weyl_wall_orbit_count": 105,
    }
    assert all(audit["checks"].values())
    assert audit["mismatches"] == {
        "chamber_summaries": [],
        "class_invariant_multisets": [],
        "strict_extension_action_sets": [],
    }


def test_normalization_audit_is_tied_to_the_frozen_inputs() -> None:
    audit = _load(OUTPUT)
    assert _sha256(OLD_INVENTORY) == audit["sources"]["raw_weighted_inventory"][
        "sha256"
    ]
    assert _sha256(ORIENTED_INVENTORY) == audit["sources"][
        "oriented_weighted_inventory"
    ]["sha256"]

