#!/usr/bin/env python3
"""Verify the root-class fragment engine and its complete orbit-0 split replay."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from scripts.explore_h9_14_record31_a1_downstream import canonical_family_key
from scripts.explore_h9_14_record31_a3_downstream import mask
from scripts.generate_h9_14_record31_a1_root_class_engine import (
    generate_root_class_engine,
)


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data/enumerations/h9_14_record31_a1_orbit0_corrected.json"
FULL = ROOT / "data/enumerations/h9_14_record31_a1_root_class_orbit0_full.json"
LEFT = ROOT / "data/enumerations/h9_14_record31_a1_root_class_orbit0_1_2.json"
RIGHT = ROOT / "data/enumerations/h9_14_record31_a1_root_class_orbit0_3_4.json"
BOUNDARY = (
    ROOT
    / "data/enumerations"
    / "h9_14_record31_a1_orbit4_inline_root_boundaries.json"
)

EXPECTED_SHA256 = {
    REFERENCE: "6f771cbaf9b1b2249217c0e25aa3c4a5217b63a6f7abd1f1328155e700bfaa39",
    FULL: "077476c0b8390d6b2b9bae0dcee1957f7078f0114a17a646b1c83d8cde31ade1",
    LEFT: "2632f473a7f2f318d8ef13fc44c2585fab6d7b2fbd0b6b6cbdcd58e9f8d56fbe",
    RIGHT: "9160190fc43504ba7f789bdd46a00c4d4ed3b9a668d17df0b8b977897d205ec4",
}
EXPECTED_GENERATED_SOURCE_SHA256 = (
    "b2e34b810432db09fe685b52b39ff97ade0c7c1b3250deea0fbe8904713ba945"
)
ROOT_FIELDS = (
    "root_class_range_restricted",
    "requested_root_class_begin",
    "requested_root_class_end",
    "root_class_discovered_count",
    "root_class_selected_count",
    "root_class_skipped_count",
)


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise AssertionError(f"{label} is not an object")
    return value


def _load(path: Path) -> tuple[dict[str, object], str]:
    raw = path.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != EXPECTED_SHA256[path]:
        raise AssertionError(f"archived bytes changed: {path.name}")
    document = json.loads(raw)
    if not isinstance(document, dict):
        raise AssertionError(f"archive is not an object: {path.name}")
    return document, digest


def _row(document: Mapping[str, object], label: str) -> Mapping[str, object]:
    if (
        document.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
        or document.get("compact_only") is not True
        or document.get("enumeration_mode") != "nauty_partial_quotient"
        or document.get("orbit_count") != 1
    ):
        raise AssertionError(f"{label} header changed")
    rows = document.get("orbits")
    if not isinstance(rows, list) or len(rows) != 1:
        raise AssertionError(f"{label} row count changed")
    row = _mapping(rows[0], f"{label} row")
    if row.get("orbit_index") != 0 or row.get("positive_A4_masks") != [[2, 3, 4, 5]]:
        raise AssertionError(f"{label} source changed")
    return row


def _target_keys(row: Mapping[str, object]) -> set[tuple[int, ...]]:
    raw_targets = row.get("ridge_valid_target_minimal_nonfaces")
    if not isinstance(raw_targets, list):
        raise AssertionError("root-fragment target array changed")
    result = set()
    for raw_family in raw_targets:
        if not isinstance(raw_family, list):
            raise AssertionError("root-fragment target changed")
        family = tuple(mask(map(int, face)) for face in raw_family)
        result.add(canonical_family_key(family))
    return result


def verify() -> dict[str, object]:
    documents: dict[Path, dict[str, object]] = {}
    digests: dict[str, str] = {}
    for path in EXPECTED_SHA256:
        document, digest = _load(path)
        documents[path] = document
        digests[path.name] = digest

    source = generate_root_class_engine()
    source_digest = sha256(source.encode()).hexdigest()
    if source_digest != EXPECTED_GENERATED_SOURCE_SHA256:
        raise AssertionError("root-class generated C++ changed")
    required_source_fragments = (
        "std::unordered_set<std::string> root_class_codes",
        "if (!audit->root_class_codes.insert(canonical_code).second)",
        "const int root_class_id = ++audit->root_class_discovered",
        "root_class_id < gRootClassBegin || root_class_id > gRootClassEnd",
        'std::string(argv[2]) == "--root-classes"',
    )
    for fragment in required_source_fragments:
        if fragment not in source:
            raise AssertionError("root-class range guard changed")

    reference = documents[REFERENCE]
    full = documents[FULL]
    full_core = json.loads(json.dumps(full))
    full_rows = full_core.get("orbits")
    if not isinstance(full_rows, list) or len(full_rows) != 1:
        raise AssertionError("full root-class replay row changed")
    full_core_row = _mapping(full_rows[0], "full core row")
    for key in ROOT_FIELDS:
        full_core_row.pop(key, None)  # type: ignore[attr-defined]
    full_core["scope"] = reference["scope"]
    if full_core != reference:
        raise AssertionError("unrestricted root-class engine changed the exact transcript")

    full_row = _row(full, "full")
    left_row = _row(documents[LEFT], "left fragment")
    right_row = _row(documents[RIGHT], "right fragment")
    expected_metadata = (
        (full_row, False, 1, 1_000_000, 4, 4, 0),
        (left_row, True, 1, 2, 4, 2, 2),
        (right_row, True, 3, 4, 4, 2, 2),
    )
    for row, restricted, begin, end, discovered, selected, skipped in expected_metadata:
        if (
            row.get("root_class_range_restricted") is not restricted
            or row.get("requested_root_class_begin") != begin
            or row.get("requested_root_class_end") != end
            or row.get("root_class_discovered_count") != discovered
            or row.get("root_class_selected_count") != selected
            or row.get("root_class_skipped_count") != skipped
        ):
            raise AssertionError("root-class fragment metadata changed")

    full_keys = _target_keys(full_row)
    left_keys = _target_keys(left_row)
    right_keys = _target_keys(right_row)
    if len(full_keys) != 1 or left_keys | right_keys != full_keys:
        raise AssertionError("root-fragment target-orbit union differs from full run")
    if left_keys or right_keys != full_keys:
        raise AssertionError("bounded split fixture partition changed")
    full_nodes = int(full_row["search_node_count"])
    fragment_nodes = int(left_row["search_node_count"]) + int(
        right_row["search_node_count"]
    )
    if full_nodes != 231_773 or fragment_nodes != 266_434 or fragment_nodes <= full_nodes:
        raise AssertionError("root-fragment nonadditivity witness changed")

    boundary = _mapping(json.loads(BOUNDARY.read_bytes()), "source-4 boundary")
    if (
        boundary.get("status")
        != "terminal_complete_source_orbit_certificate_by_root_class_partition"
        or boundary.get("source_orbit_index") != 4
        or boundary.get("positive_A4_mask") != 10252
        or boundary.get("run", {}).get("full_prefix_depth_one_class_count") != 6
        or boundary.get("latest_observed_checkpoint", {}).get("completed_class_ids")
        != [1, 2, 3]
        or boundary.get("planned_selective_restart", {}).get("root_class_begin") != 4
        or boundary.get("planned_selective_restart", {}).get("root_class_end") != 6
        or boundary.get("claims", {}).get("source_orbit_census_complete") is not True
    ):
        raise AssertionError("source-4 selective-restart boundary changed")

    return {
        "kind": "verified_h9_14_record31_a1_root_class_split",
        "archive_sha256": digests,
        "generated_source_sha256": source_digest,
        "unrestricted_transcript_core_identical": True,
        "bounded_full_search_node_count": full_nodes,
        "bounded_split_search_node_count": fragment_nodes,
        "bounded_target_orbit_union_identical": True,
        "split_internal_counts_nonadditive": True,
        "prepared_source4_root_class_range": [4, 6],
        "scope": (
            "The range engine exhausts selected canonical root classes but "
            "restarts deeper memoization, so target-orbit unions are exact "
            "while internal fragment counts are not additive."
        ),
    }


def main() -> None:
    print(json.dumps(verify(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
