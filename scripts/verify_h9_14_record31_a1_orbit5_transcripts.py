#!/usr/bin/env python3
"""Verify the two terminal corrected Record-31 a=1 source-orbit-5 transcripts."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from scripts.generate_h9_14_record31_a1_inline_key_engine import generate_inline


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record31_a1_orbit5_corrected_quotient.json"
)
INDEPENDENT = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record31_a1_orbit5_independent_matrix_quotient.json"
)

EXPECTED_PRODUCTION_SHA256 = (
    "90fb08a5889d6c3fa0dd30aae35cb33feee334695482a31a29513a64ad92d294"
)
EXPECTED_INDEPENDENT_SHA256 = (
    "7dee139c3b5740635df2f7814e8e34f58983a9a6469f40fa7b183d586d9c395b"
)
EXPECTED_INLINE_GENERATED_SOURCE_SHA256 = (
    "ec577df47018f1d03d48e8304c56a2541ef02180144937600569ec0a9a7e68c7"
)
EXPECTED_PRODUCTION_BINARY_SHA256 = (
    "0cb72b4b62c4ec8f66f3440f87d2eb4b1771a67dab6203ab67b0fe611e5c3bf4"
)

EXPECTED_SCHEDULED_HISTOGRAM = {
    "120": 735,
    "121": 305,
    "122": 29,
    "123": 144,
    "124": 138,
    "125": 34,
    "126": 5,
    "127": 12,
    "128": 2,
    "129": 7,
}


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AssertionError(f"{label} is not an object")
    return value


def _integer_histogram(value: object, label: str) -> dict[str, int]:
    mapping = _mapping(value, label)
    return {str(key): int(count) for key, count in mapping.items()}


def _load(path: Path, expected_digest: str) -> tuple[Mapping[str, object], str]:
    raw = path.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != expected_digest:
        raise AssertionError(f"terminal archive bytes changed: {path.name}")
    return _mapping(json.loads(raw), path.name), digest


def _row(
    document: Mapping[str, object],
    mode: str,
    expected: Mapping[str, object],
) -> Mapping[str, object]:
    if (
        document.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
        or document.get("compact_only") is not True
        or document.get("enumeration_mode") != mode
        or document.get("orbit_count") != 1
        or document.get("pure_120_facet_target_count_across_representatives")
        != 735
    ):
        raise AssertionError(f"{mode} transcript header changed")
    raw_rows = document.get("orbits")
    if not isinstance(raw_rows, list) or len(raw_rows) != 1:
        raise AssertionError(f"{mode} transcript row count changed")
    row = _mapping(raw_rows[0], f"{mode} source row")

    common = {
        "orbit_index": 5,
        "forced_A5_shadow_count": 8,
        "additional_A5_count": 7,
        "candidate_count": 924,
        "base_small_transversal_count": 74,
        "dead_five_count": 0,
        "pure_120_facet_target_count": 735,
        "ridge_valid_pure_target_count": 0,
        "ridge_valid_target_orbit_count": 0,
    }
    for key, value in (*common.items(), *expected.items()):
        if row.get(key) != value:
            raise AssertionError(f"{mode} scalar {key} changed")
    if (
        row.get("positive_A4_masks") != [[2, 11, 12, 13]]
        or row.get("fixed_local_failure") is not False
        or row.get("fixed_local_failure_witnesses") != []
        or row.get("ridge_valid_target_minimal_nonfaces") != []
        or row.get("ridge_valid_target_orbit_literal_multiplicities") != []
    ):
        raise AssertionError(f"{mode} source identity or terminal target arrays changed")

    depth_nodes = [int(value) for value in row["partition_depth_node_counts"]]
    seen = [int(value) for value in row["quotient_seen_counts"]]
    if depth_nodes != expected["partition_depth_node_counts"]:
        raise AssertionError(f"{mode} depth-node partition changed")
    if seen != expected["quotient_seen_counts"]:
        raise AssertionError(f"{mode} quotient-seen partition changed")
    if sum(depth_nodes) != row["search_node_count"]:
        raise AssertionError(f"{mode} search-node sum failed")
    if row["search_node_count"] - sum(seen) != row[
        "quotient_isomorphic_prune_count"
    ]:
        raise AssertionError(f"{mode} quotient-prune identity failed")

    order_five = _integer_histogram(
        row.get("order_five_transversal_count_histogram"),
        f"{mode} order-five histogram",
    )
    scheduled = _integer_histogram(
        row.get("scheduled_blocker_size_histogram"),
        f"{mode} blocker-size histogram",
    )
    pure_ridges = _integer_histogram(
        row.get("pure_target_ridge_histogram"),
        f"{mode} pure-ridge histogram",
    )
    if order_five != expected["order_five_transversal_count_histogram"]:
        raise AssertionError(f"{mode} order-five histogram changed")
    if scheduled != EXPECTED_SCHEDULED_HISTOGRAM:
        raise AssertionError(f"{mode} scheduled-blocker histogram changed")
    if sum(order_five.values()) != row["no_small_transversal_count"]:
        raise AssertionError(f"{mode} order-five total failed")
    if sum(scheduled.values()) != order_five["120"]:
        raise AssertionError(f"{mode} scheduled-blocker total failed")
    if scheduled["120"] != row["pure_120_facet_target_count"]:
        raise AssertionError(f"{mode} pure-120 identity failed")
    if len(pure_ridges) != 231 or sum(pure_ridges.values()) != 735:
        raise AssertionError(f"{mode} pure-ridge census changed")
    if pure_ridges.get("2:540", 0) != 0:
        raise AssertionError(f"{mode} unexpectedly has a corrected ridge target")
    return row


def verify() -> dict[str, object]:
    production, production_digest = _load(
        PRODUCTION, EXPECTED_PRODUCTION_SHA256
    )
    independent, independent_digest = _load(
        INDEPENDENT, EXPECTED_INDEPENDENT_SHA256
    )
    source_digest = sha256(generate_inline().encode()).hexdigest()
    if source_digest != EXPECTED_INLINE_GENERATED_SOURCE_SHA256:
        raise AssertionError("equivalence-certified inline source changed")

    production_expected = {
        "search_node_count": 41_488_368,
        "partition_depth_node_counts": [
            1, 196, 518, 8456, 165396, 3923619, 37353399, 36783
        ],
        "quotient_seen_counts": [
            1, 3, 51, 1099, 41169, 1713317, 25378215, 8231
        ],
        "quotient_isomorphic_prune_count": 14_346_282,
        "dead_uncoverable_count": 5_081_805,
        "dead_slot_bound_count": 21_487_754,
        "no_small_transversal_count": 3_282,
        "order_five_transversal_count_histogram": {
            "120": 1411,
            "121": 962,
            "122": 315,
            "123": 190,
            "124": 244,
            "125": 98,
            "126": 9,
            "127": 15,
            "128": 24,
            "129": 14,
        },
    }
    independent_expected = {
        "search_node_count": 6_042_612,
        "partition_depth_node_counts": [
            1, 196, 644, 5068, 60477, 1032822, 4905830, 37574
        ],
        "quotient_seen_counts": [
            1, 4, 38, 477, 11759, 412438, 2934227, 7995
        ],
        "quotient_isomorphic_prune_count": 2_675_673,
        "dead_uncoverable_count": 878_415,
        "dead_slot_bound_count": 2_381_155,
        "no_small_transversal_count": 3_283,
        "order_five_transversal_count_histogram": {
            "120": 1411,
            "121": 960,
            "122": 317,
            "123": 190,
            "124": 244,
            "125": 99,
            "126": 9,
            "127": 15,
            "128": 24,
            "129": 14,
        },
    }
    production_row = _row(
        production, "nauty_partial_quotient", production_expected
    )
    independent_row = _row(
        independent,
        "direct_binary_matrix_partial_quotient",
        independent_expected,
    )

    # Quotient-representative counts differ before the pure-120 gate, as
    # expected for the two exact state keys.  The theorem-facing terminal
    # census agrees exactly.
    if production_row["search_node_count"] == independent_row["search_node_count"]:
        raise AssertionError("independent traversal collapsed to production")
    if (
        production_row["scheduled_blocker_size_histogram"]
        != independent_row["scheduled_blocker_size_histogram"]
        or production_row["pure_target_ridge_histogram"]
        != independent_row["pure_target_ridge_histogram"]
    ):
        raise AssertionError("the two terminal pure-target censuses differ")

    return {
        "kind": "verified_h9_14_record31_a1_orbit5_corrected_quotients",
        "production_archive_sha256": production_digest,
        "independent_archive_sha256": independent_digest,
        "inline_generated_source_sha256": source_digest,
        "production_binary_sha256": EXPECTED_PRODUCTION_BINARY_SHA256,
        "production_search_node_count": production_row["search_node_count"],
        "independent_search_node_count": independent_row["search_node_count"],
        "pure_120_facet_target_count": 735,
        "pure_ridge_signature_count": 231,
        "ridge_valid_target_orbit_count": 0,
        "terminal_pure_target_censuses_identical": True,
        "scope": (
            "Complete corrected source-orbit-5 compact incidence no-go in two "
            "exact quotient traversals; not a full Record-31 a=1 theorem and "
            "not a finite-volume noncompact inference."
        ),
    }


def main() -> None:
    print(json.dumps(verify(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
