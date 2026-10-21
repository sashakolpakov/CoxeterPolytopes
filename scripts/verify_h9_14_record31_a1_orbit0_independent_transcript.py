#!/usr/bin/env python3
"""Verify the independent matrix-quotient replay of Record-31 a=1 source 0."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from scripts.explore_h9_14_record31_a1_downstream import canonical_family_key
from scripts.explore_h9_14_record31_a3_downstream import (
    EDGE_MASKS,
    mask,
    minimal_transversals,
    ridge_degree_histogram,
)


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION = (
    ROOT / "data/enumerations/h9_14_record31_a1_orbit0_corrected.json"
)
INDEPENDENT = (
    ROOT
    / "data/enumerations"
    / "h9_14_record31_a1_orbit0_independent_matrix_quotient.json"
)
EXPECTED_PRODUCTION_SHA256 = (
    "6f771cbaf9b1b2249217c0e25aa3c4a5217b63a6f7abd1f1328155e700bfaa39"
)
EXPECTED_INDEPENDENT_SHA256 = (
    "f39712bdadeee4ca1a7852f52efa10ccaa75ac92b9409949ec2285fd84384325"
)
EXPECTED_CANONICAL_FAMILY_KEY = (
    0, 0, 12, 12, 16, 16, 5, 5, 14, 17, 17, 42, 42, 50, 50
)


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AssertionError(f"{label} is not an object")
    return value


def _load(path: Path, expected_digest: str) -> tuple[Mapping[str, object], str]:
    raw = path.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != expected_digest:
        raise AssertionError(f"terminal archive bytes changed: {path.name}")
    return _mapping(json.loads(raw), path.name), digest


def _row(
    document: Mapping[str, object],
    mode: str,
    expected_nodes: int,
    expected_depth: list[int],
    expected_seen: list[int],
) -> Mapping[str, object]:
    if (
        document.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
        or document.get("compact_only") is not True
        or document.get("enumeration_mode") != mode
        or document.get("orbit_count") != 1
        or document.get("pure_120_facet_target_count_across_representatives")
        != 491
    ):
        raise AssertionError(f"{mode} header changed")
    rows = document.get("orbits")
    if not isinstance(rows, list) or len(rows) != 1:
        raise AssertionError(f"{mode} row count changed")
    row = _mapping(rows[0], f"{mode} row")
    expected_common = {
        "orbit_index": 0,
        "positive_A4_masks": [[2, 3, 4, 5]],
        "additional_A5_count": 5,
        "search_node_count": expected_nodes,
        "partition_depth_node_counts": expected_depth,
        "quotient_seen_counts": expected_seen,
        "no_small_transversal_count": 711,
        "pure_120_facet_target_count": 491,
        "ridge_valid_pure_target_count": 1,
        "ridge_valid_target_orbit_count": 1,
        "ridge_valid_target_orbit_literal_multiplicities": [1],
    }
    for key, expected in expected_common.items():
        if row.get(key) != expected:
            raise AssertionError(f"{mode} field {key} changed")
    if sum(expected_depth) != expected_nodes:
        raise AssertionError(f"{mode} depth-node sum failed")
    if expected_nodes - sum(expected_seen) != row[
        "quotient_isomorphic_prune_count"
    ]:
        raise AssertionError(f"{mode} quotient-prune identity failed")
    if sum(int(value) for value in row["order_five_transversal_count_histogram"].values()) != 711:
        raise AssertionError(f"{mode} order-five total failed")
    if sum(int(value) for value in row["scheduled_blocker_size_histogram"].values()) != 525:
        raise AssertionError(f"{mode} scheduled-blocker total failed")
    pure_ridges = _mapping(
        row.get("pure_target_ridge_histogram"), f"{mode} ridge histogram"
    )
    if len(pure_ridges) != 42 or sum(map(int, pure_ridges.values())) != 491:
        raise AssertionError(f"{mode} pure-ridge census changed")
    if pure_ridges.get("2:540") != 1:
        raise AssertionError(f"{mode} corrected ridge count changed")
    return row


def _family(row: Mapping[str, object]) -> tuple[int, ...]:
    targets = row.get("ridge_valid_target_minimal_nonfaces")
    if not isinstance(targets, list) or len(targets) != 1:
        raise AssertionError("source-0 target list changed")
    raw_family = targets[0]
    if not isinstance(raw_family, list):
        raise AssertionError("source-0 target family changed")
    return tuple(mask(map(int, face)) for face in raw_family)


def verify() -> dict[str, object]:
    production, production_digest = _load(
        PRODUCTION, EXPECTED_PRODUCTION_SHA256
    )
    independent, independent_digest = _load(
        INDEPENDENT, EXPECTED_INDEPENDENT_SHA256
    )
    production_row = _row(
        production,
        "nauty_partial_quotient",
        231_773,
        [1, 56, 420, 8495, 195250, 27551],
        [1, 4, 47, 1393, 51035, 8396],
    )
    independent_row = _row(
        independent,
        "direct_binary_matrix_partial_quotient",
        242_124,
        [1, 56, 420, 8495, 205835, 27317],
        [1, 4, 47, 1523, 55460, 8447],
    )
    if (
        production_row["scheduled_blocker_size_histogram"]
        != independent_row["scheduled_blocker_size_histogram"]
        or production_row["pure_target_ridge_histogram"]
        != independent_row["pure_target_ridge_histogram"]
    ):
        raise AssertionError("source-0 terminal pure-target censuses differ")

    production_family = _family(production_row)
    independent_family = _family(independent_row)
    production_key = canonical_family_key(production_family)
    independent_key = canonical_family_key(independent_family)
    if (
        production_key != EXPECTED_CANONICAL_FAMILY_KEY
        or independent_key != EXPECTED_CANONICAL_FAMILY_KEY
    ):
        raise AssertionError("source-0 target orbit changed")
    for family in (production_family, independent_family):
        if Counter(face.bit_count() for face in family) != {4: 1, 5: 5}:
            raise AssertionError("source-0 target order profile changed")
        blocker = minimal_transversals((*EDGE_MASKS, *family))
        if Counter(face.bit_count() for face in blocker) != {5: 120}:
            raise AssertionError("source-0 literal blocker replay failed")
        if ridge_degree_histogram(blocker) != {"2": 540}:
            raise AssertionError("source-0 literal ridge replay failed")

    return {
        "kind": "verified_h9_14_record31_a1_orbit0_independent_incidence",
        "production_archive_sha256": production_digest,
        "independent_archive_sha256": independent_digest,
        "production_search_node_count": production_row["search_node_count"],
        "independent_search_node_count": independent_row["search_node_count"],
        "ridge_valid_target_orbit_count": 1,
        "canonical_family_key": list(EXPECTED_CANONICAL_FAMILY_KEY),
        "terminal_pure_target_censuses_identical": True,
        "independent_literal_blocker_and_ridge_replay": True,
        "scope": (
            "Complete independent incidence replay for source orbit 0; its "
            "separate frozen downstream compact-label no-go remains the reason "
            "the incidence target yields no Gram candidate."
        ),
    }


def main() -> None:
    print(json.dumps(verify(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
