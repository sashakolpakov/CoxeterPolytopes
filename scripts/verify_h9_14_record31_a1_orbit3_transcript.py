#!/usr/bin/env python3
"""Verify the terminal corrected Record-31 a=1 source-orbit-3 transcript."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from scripts.explore_h9_14_record31_a3_downstream import (
    EDGE_MASKS,
    minimal_transversals,
    ridge_degree_histogram,
)
from scripts.generate_h9_14_record31_a1_compact_key_engine import generate_compact


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record31_a1_orbit3_corrected_quotient.json"
)
EXPECTED_ARCHIVE_SHA256 = (
    "b4ea149c1259d713dd8c47910c5052d2f00c51a8b4b55c988bf22c22ec035345"
)
EXPECTED_GENERATED_SOURCE_SHA256 = (
    "186418cf82872f3145c8ebcf4e1fe953a61cb9836eb04c4734ee0e1c35fd8443"
)
EXPECTED_BINARY_SHA256 = (
    "b740f0d8ba7fc357df6db1e50b70629533b9375a585ea2fc62f09846907a86cc"
)

EXPECTED_FAMILIES = (
    (8220, 121, 122, 1144, 1924, 7040, 7264, 9092),
    (8220, 481, 482, 1504, 1564, 6680, 7744, 8612),
    (8220, 481, 482, 1504, 1564, 6680, 6848, 9508),
    (8220, 233, 234, 248, 1812, 6368, 7936, 9988),
    (8220, 233, 248, 1812, 3844, 4328, 6368, 7936),
    (8220, 481, 992, 1564, 3100, 7192, 7744, 8612),
    (8220, 481, 992, 1564, 3100, 7192, 7360, 8996),
    (8220, 233, 248, 1812, 3844, 6368, 7936, 9988),
    (8220, 121, 248, 1924, 3844, 6368, 7936, 9988),
    (8220, 481, 488, 1564, 3604, 4576, 7744, 8616),
    (8220, 121, 124, 1924, 3968, 4216, 6496, 9860),
)


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AssertionError(f"{label} is not an object")
    return value


def _integer_histogram(value: object, label: str) -> dict[str, int]:
    mapping = _mapping(value, label)
    return {str(key): int(count) for key, count in mapping.items()}


def _mask(vertices: Sequence[object]) -> int:
    return sum(1 << int(vertex) for vertex in vertices)


def verify() -> dict[str, object]:
    raw = ARCHIVE.read_bytes()
    archive_digest = sha256(raw).hexdigest()
    if archive_digest != EXPECTED_ARCHIVE_SHA256:
        raise AssertionError("terminal source-orbit-3 archive bytes changed")
    source_digest = sha256(generate_compact().encode()).hexdigest()
    if source_digest != EXPECTED_GENERATED_SOURCE_SHA256:
        raise AssertionError("equivalence-certified compact source changed")

    document = _mapping(json.loads(raw), "terminal transcript")
    if (
        document.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
        or document.get("compact_only") is not True
        or document.get("enumeration_mode") != "nauty_partial_quotient"
        or document.get("orbit_count") != 1
        or document.get("pure_120_facet_target_count_across_representatives")
        != 4445
    ):
        raise AssertionError("terminal transcript header changed")
    raw_rows = document.get("orbits")
    if not isinstance(raw_rows, list) or len(raw_rows) != 1:
        raise AssertionError("terminal transcript row count changed")
    row = _mapping(raw_rows[0], "terminal source row")
    expected_scalars = {
        "orbit_index": 3,
        "forced_A5_shadow_count": 8,
        "additional_A5_count": 7,
        "candidate_count": 952,
        "base_small_transversal_count": 74,
        "search_node_count": 469_990_577,
        "dead_uncoverable_count": 39_299_444,
        "dead_five_count": 0,
        "dead_slot_bound_count": 265_737_338,
        "quotient_isomorphic_prune_count": 159_523_605,
        "no_small_transversal_count": 13_273,
        "pure_120_facet_target_count": 4_445,
        "ridge_valid_pure_target_count": 11,
        "ridge_valid_target_orbit_count": 11,
    }
    for key, expected in expected_scalars.items():
        if row.get(key) != expected:
            raise AssertionError(f"terminal scalar {key} changed")
    if (
        row.get("positive_A4_masks") != [[2, 3, 4, 13]]
        or row.get("fixed_local_failure") is not False
        or row.get("fixed_local_failure_witnesses") != []
    ):
        raise AssertionError("terminal source identity changed")

    depth_nodes = [int(value) for value in row["partition_depth_node_counts"]]
    seen = [int(value) for value in row["quotient_seen_counts"]]
    if depth_nodes != [1, 196, 1554, 25284, 691780, 29702118, 439208906, 360738]:
        raise AssertionError("terminal depth-node partition changed")
    if seen != [1, 9, 154, 4388, 235670, 14715798, 295348921, 162031]:
        raise AssertionError("terminal quotient-seen partition changed")
    if sum(depth_nodes) != row["search_node_count"]:
        raise AssertionError("search nodes do not equal depth-node sum")
    if row["search_node_count"] - sum(seen) != row[
        "quotient_isomorphic_prune_count"
    ]:
        raise AssertionError("quotient prune identity failed")

    order_five = _integer_histogram(
        row.get("order_five_transversal_count_histogram"), "order-five histogram"
    )
    scheduled = _integer_histogram(
        row.get("scheduled_blocker_size_histogram"), "blocker-size histogram"
    )
    pure_ridges = _integer_histogram(
        row.get("pure_target_ridge_histogram"), "pure ridge histogram"
    )
    if order_five != {
        "120": 6193,
        "121": 3552,
        "122": 1896,
        "123": 612,
        "124": 547,
        "125": 290,
        "126": 60,
        "127": 60,
        "128": 35,
        "129": 28,
    }:
        raise AssertionError("order-five histogram changed")
    if scheduled != {
        "120": 4445,
        "121": 554,
        "122": 79,
        "123": 480,
        "124": 400,
        "125": 87,
        "126": 79,
        "127": 4,
        "128": 28,
        "129": 33,
        "130": 3,
        "132": 1,
    }:
        raise AssertionError("scheduled blocker histogram changed")
    if sum(order_five.values()) != row["no_small_transversal_count"]:
        raise AssertionError("order-five histogram total changed")
    if sum(scheduled.values()) != order_five["120"]:
        raise AssertionError("scheduled blocker histogram total changed")
    if scheduled["120"] != row["pure_120_facet_target_count"]:
        raise AssertionError("pure-120 target identity failed")
    if sum(pure_ridges.values()) != row["pure_120_facet_target_count"]:
        raise AssertionError("pure ridge histogram total changed")
    if pure_ridges.get("2:540") != row["ridge_valid_pure_target_count"]:
        raise AssertionError("corrected ridge gate count changed")

    raw_families = row.get("ridge_valid_target_minimal_nonfaces")
    raw_multiplicities = row.get("ridge_valid_target_orbit_literal_multiplicities")
    if not isinstance(raw_families, list) or not isinstance(raw_multiplicities, list):
        raise AssertionError("terminal target arrays changed")
    families = tuple(
        tuple(_mask(face) for face in raw_family) for raw_family in raw_families
    )
    multiplicities = tuple(int(value) for value in raw_multiplicities)
    if families != EXPECTED_FAMILIES:
        raise AssertionError("terminal minimal-nonface families changed")
    if multiplicities != (1,) * len(EXPECTED_FAMILIES):
        raise AssertionError("terminal literal multiplicities changed")
    if sum(multiplicities) != row["ridge_valid_pure_target_count"]:
        raise AssertionError("literal target multiplicity total changed")

    for family in families:
        if Counter(face.bit_count() for face in family) != {4: 1, 5: 7}:
            raise AssertionError("terminal family order profile changed")
        blocker = minimal_transversals((*EDGE_MASKS, *family))
        if Counter(face.bit_count() for face in blocker) != {5: 120}:
            raise AssertionError("independent literal blocker replay failed")
        if ridge_degree_histogram(blocker) != {"2": 540}:
            raise AssertionError("independent corrected ridge replay failed")

    return {
        "kind": "verified_h9_14_record31_a1_orbit3_corrected_quotient",
        "archive_sha256": archive_digest,
        "generated_source_sha256": source_digest,
        "production_binary_sha256": EXPECTED_BINARY_SHA256,
        "search_node_count": row["search_node_count"],
        "pure_120_facet_target_count": row["pure_120_facet_target_count"],
        "ridge_valid_literal_target_count": row["ridge_valid_pure_target_count"],
        "ridge_valid_target_orbit_count": row["ridge_valid_target_orbit_count"],
        "target_literal_multiplicities": list(multiplicities),
        "independent_blocker_and_ridge_replay": True,
        "scope": (
            "Complete corrected production source-orbit-3 incidence transcript; "
            "not a topology, Coxeter-label, Gram, hyperbolic-example, or full "
            "Record-31 a=1 schedule theorem."
        ),
    }


def main() -> None:
    print(json.dumps(verify(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
