#!/usr/bin/env python3
"""Verify the partitioned terminal Record-31 a=1 source-orbit-4 no-go."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from scripts.verify_h9_14_record31_a1_root_class_split import (
    verify as verify_split_engine,
)


ROOT = Path(__file__).resolve().parents[1]
BOUNDARY = (
    ROOT
    / "data/enumerations"
    / "h9_14_record31_a1_orbit4_inline_root_boundaries.json"
)
FRAGMENT = (
    ROOT
    / "data/enumerations"
    / "h9_14_record31_a1_orbit4_root_classes4_6.json"
)
EXPECTED_BOUNDARY_SHA256 = (
    "5644ed22f418c767849a2e9a9c4f57f0b2088e2d2d6aa4eb9a4c66a044c3c209"
)
EXPECTED_FRAGMENT_SHA256 = (
    "dc463688433307ad8fe04c0a03af73b2080ef0fa49e71732743eaaca80675acc"
)
EXPECTED_ORDER_FIVE = {
    "120": 4749,
    "121": 2643,
    "122": 1059,
    "123": 816,
    "124": 457,
    "125": 188,
    "126": 40,
    "127": 62,
    "128": 3,
    "130": 11,
}
EXPECTED_SCHEDULED = {
    "120": 3184,
    "121": 667,
    "122": 108,
    "123": 258,
    "124": 335,
    "125": 129,
    "126": 19,
    "127": 9,
    "128": 15,
    "129": 22,
    "130": 1,
    "132": 2,
}


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AssertionError(f"{label} is not an object")
    return value


def _load(path: Path, expected_digest: str) -> tuple[Mapping[str, object], str]:
    raw = path.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != expected_digest:
        raise AssertionError(f"partitioned source-4 bytes changed: {path.name}")
    return _mapping(json.loads(raw), path.name), digest


def _histogram(value: object, label: str) -> dict[str, int]:
    return {str(key): int(count) for key, count in _mapping(value, label).items()}


def verify() -> dict[str, object]:
    split_result = verify_split_engine()
    if (
        split_result.get("unrestricted_transcript_core_identical") is not True
        or split_result.get("bounded_target_orbit_union_identical") is not True
        or split_result.get("split_internal_counts_nonadditive") is not True
        or split_result.get("prepared_source4_root_class_range") != [4, 6]
    ):
        raise AssertionError("root-class split-engine certificate changed")

    boundary, boundary_digest = _load(BOUNDARY, EXPECTED_BOUNDARY_SHA256)
    fragment, fragment_digest = _load(FRAGMENT, EXPECTED_FRAGMENT_SHA256)
    if (
        fragment.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
        or fragment.get("compact_only") is not True
        or fragment.get("enumeration_mode") != "nauty_partial_quotient"
        or fragment.get("orbit_count") != 1
        or fragment.get("pure_120_facet_target_count_across_representatives")
        != 3184
    ):
        raise AssertionError("terminal source-4 fragment header changed")
    rows = fragment.get("orbits")
    if not isinstance(rows, list) or len(rows) != 1:
        raise AssertionError("terminal source-4 fragment row changed")
    row = _mapping(rows[0], "terminal source-4 row")
    expected_scalars = {
        "orbit_index": 4,
        "forced_A5_shadow_count": 8,
        "additional_A5_count": 7,
        "candidate_count": 952,
        "base_small_transversal_count": 74,
        "search_node_count": 105_066_986,
        "quotient_isomorphic_prune_count": 33_570_325,
        "dead_uncoverable_count": 10_104_203,
        "dead_five_count": 0,
        "dead_slot_bound_count": 59_895_867,
        "no_small_transversal_count": 10_028,
        "pure_120_facet_target_count": 3_184,
        "ridge_valid_pure_target_count": 0,
        "ridge_valid_target_orbit_count": 0,
        "root_class_range_restricted": True,
        "requested_root_class_begin": 4,
        "requested_root_class_end": 6,
        "root_class_discovered_count": 6,
        "root_class_selected_count": 3,
        "root_class_skipped_count": 3,
    }
    for key, expected in expected_scalars.items():
        if row.get(key) != expected:
            raise AssertionError(f"terminal source-4 field {key} changed")
    if (
        row.get("positive_A4_masks") != [[2, 3, 11, 13]]
        or row.get("fixed_local_failure") is not False
        or row.get("fixed_local_failure_witnesses") != []
        or row.get("ridge_valid_target_minimal_nonfaces") != []
        or row.get("ridge_valid_target_orbit_literal_multiplicities") != []
    ):
        raise AssertionError("terminal source-4 identity or target arrays changed")

    depth_nodes = [int(value) for value in row["partition_depth_node_counts"]]
    seen = [int(value) for value in row["quotient_seen_counts"]]
    if depth_nodes != [1, 196, 588, 11704, 337227, 11025632, 93591394, 100244]:
        raise AssertionError("terminal source-4 depth-node partition changed")
    if seen != [1, 3, 69, 2145, 106155, 5683652, 65673348, 31285]:
        raise AssertionError("terminal source-4 quotient-seen partition changed")
    if sum(depth_nodes) != row["search_node_count"]:
        raise AssertionError("terminal source-4 node sum failed")
    if row["search_node_count"] - sum(seen) != (
        row["quotient_isomorphic_prune_count"] + row["root_class_skipped_count"]
    ):
        raise AssertionError("terminal source-4 prune/skip identity failed")

    order_five = _histogram(
        row.get("order_five_transversal_count_histogram"), "order-five histogram"
    )
    scheduled = _histogram(
        row.get("scheduled_blocker_size_histogram"), "scheduled histogram"
    )
    pure_ridges = _histogram(
        row.get("pure_target_ridge_histogram"), "pure-ridge histogram"
    )
    if order_five != EXPECTED_ORDER_FIVE or scheduled != EXPECTED_SCHEDULED:
        raise AssertionError("terminal source-4 blocker histograms changed")
    if sum(order_five.values()) != 10_028:
        raise AssertionError("terminal source-4 order-five total changed")
    if sum(scheduled.values()) != order_five["120"]:
        raise AssertionError("terminal source-4 scheduled total changed")
    if scheduled["120"] != 3_184:
        raise AssertionError("terminal source-4 pure-120 identity changed")
    if len(pure_ridges) != 818 or sum(pure_ridges.values()) != 3_184:
        raise AssertionError("terminal source-4 pure-ridge census changed")
    if pure_ridges.get("2:540", 0) != 0:
        raise AssertionError("terminal source-4 fragment gained a ridge target")

    completed = boundary.get("completed_class_boundaries")
    if not isinstance(completed, list) or [item.get("completed_class_id") for item in completed] != [1, 2, 3]:
        raise AssertionError("source-4 first-fragment class partition changed")
    if any(item.get("target_orbit_lower_bound") != 0 for item in completed):
        raise AssertionError("source-4 first fragment gained a target")
    clean_interrupt = _mapping(boundary.get("clean_interrupt"), "clean interrupt")
    if (
        clean_interrupt.get("terminal_stdout_bytes") != 0
        or clean_interrupt.get("completed_class_ids_retained_as_exact_boundary")
        != [1, 2, 3]
        or clean_interrupt.get("partial_class_data_is_diagnostic_only") is not True
    ):
        raise AssertionError("source-4 clean boundary changed")
    restart = _mapping(
        boundary.get("planned_selective_restart"), "selective restart"
    )
    terminal = _mapping(restart.get("terminal_completion"), "terminal fragment")
    if (
        boundary.get("status")
        != "terminal_complete_source_orbit_certificate_by_root_class_partition"
        or restart.get("status") != "terminal_complete"
        or terminal.get("selected_root_class_ids") != [4, 5, 6]
        or terminal.get("ridge_valid_target_orbit_count") != 0
        or terminal.get("terminal_transcript_sha256") != fragment_digest
    ):
        raise AssertionError("source-4 terminal partition metadata changed")
    claims = _mapping(boundary.get("claims"), "source-4 claims")
    if (
        claims.get("source_orbit_census_complete") is not True
        or claims.get("target_orbit_count_is_final") is not True
        or claims.get("source_orbit_ridge_valid_target_orbit_count") != 0
        or claims.get("record31_a1_schedule_theorem") is not False
    ):
        raise AssertionError("source-4 scope was changed or over-promoted")

    return {
        "kind": "verified_h9_14_record31_a1_orbit4_partitioned_incidence_no_go",
        "boundary_archive_sha256": boundary_digest,
        "terminal_fragment_sha256": fragment_digest,
        "split_engine_generated_source_sha256": split_result[
            "generated_source_sha256"
        ],
        "completed_original_root_class_ids": [1, 2, 3, 4, 5, 6],
        "first_fragment_ridge_valid_target_orbit_count": 0,
        "complementary_fragment_ridge_valid_target_orbit_count": 0,
        "canonical_target_orbit_union_count": 0,
        "source_orbit_incidence_impossible": True,
        "fragment_internal_counts_combined": False,
        "scope": (
            "Complete compact-incidence no-go for Record-31 a=1 source orbit "
            "4 only; source 2 and independent-incidence gaps remain, and no "
            "finite-volume noncompact statement is made."
        ),
    }


def main() -> None:
    print(json.dumps(verify(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
