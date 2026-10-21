#!/usr/bin/env python3
"""Verify the terminal Record-31 ``a=1`` source-2 matrix root-1 tranche."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from scripts.explore_h9_14_record31_a1_downstream import canonical_family_key
from scripts.generate_h9_14_record31_a1_depth_two_class_engine import (
    generate_depth_two_class_engine,
)


ROOT = Path(__file__).resolve().parents[1]
INCIDENCE = (
    ROOT / "data/enumerations/h9_14_record31_a1_orbit2_matrix_root1.json"
)
DOWNSTREAM = (
    ROOT
    / "data/enumerations/h9_14_record31_a1_orbit2_matrix_root1_downstream.json"
)
HOSTILE = (
    ROOT
    / "data/enumerations/"
    "h9_14_record31_a1_orbit2_matrix_root1_independent_downstream.json"
)

EXPECTED_FILE_HASHES = {
    INCIDENCE: "5da7b0f18b7c09236d787da89eaa81dff6d05f7df883b856f481c012b4b83eea",
    DOWNSTREAM: "da269ae7215ca03cea8f5fc07b19f477a42ceb1e8d5352788e0bfd1d115a23c9",
    HOSTILE: "19b0ddc3db5d7265f7ce06a98039bc6dfb154410c0ba0dad1fd884bde3860462",
}
EXPECTED_GENERATED_CPP_SHA256 = (
    "0bf66b0bd22827f90e4034ba80539eb1cc57b712cdd59e40070b4150e074ddfc"
)


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"malformed {name}")
    return value


def _mask(vertices: object) -> int:
    if not isinstance(vertices, list):
        raise ValueError("malformed target face")
    return sum(1 << int(vertex) for vertex in vertices)


def _raw_target_keys(row: Mapping[str, object]) -> set[tuple[int, ...]]:
    families = row.get("ridge_valid_target_minimal_nonfaces")
    if not isinstance(families, list):
        raise ValueError("source-2 target list changed")
    return {
        canonical_family_key(tuple(_mask(face) for face in family))
        for family in families
    }


def verify() -> dict[str, object]:
    for path, expected in EXPECTED_FILE_HASHES.items():
        if sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"frozen source-2 root-1 artifact changed: {path.name}")
    generated_hash = sha256(generate_depth_two_class_engine().encode()).hexdigest()
    if generated_hash != EXPECTED_GENERATED_CPP_SHA256:
        raise ValueError("source-2 root-1 generated engine changed")

    incidence = _mapping(json.loads(INCIDENCE.read_text()), "incidence transcript")
    if (
        incidence.get("kind") != "h9_14_record31_a1_pivoted_addition_scout"
        or incidence.get("compact_only") is not True
        or incidence.get("enumeration_mode")
        != "direct_binary_matrix_partial_quotient"
        or incidence.get("orbit_count") != 1
    ):
        raise ValueError("source-2 root-1 incidence header changed")
    rows = incidence.get("orbits")
    if not isinstance(rows, list) or len(rows) != 1:
        raise ValueError("source-2 root-1 incidence row changed")
    row = _mapping(rows[0], "source-2 incidence row")
    expected_fields = {
        "orbit_index": 2,
        "positive_A4_masks": [[2, 3, 11, 12]],
        "additional_A5_count": 7,
        "search_node_count": 11_053_340,
        "partition_depth_node_counts": [
            1,
            196,
            196,
            2_772,
            53_831,
            1_301_191,
            9_658_556,
            36_597,
        ],
        "quotient_seen_counts": [
            1,
            1,
            17,
            382,
            13_050,
            614_326,
            6_307_913,
            12_201,
        ],
        "root_class_discovered_count": 6,
        "root_class_selected_count": 1,
        "root_class_skipped_count": 5,
        "depth_two_class_discovered_count": 17,
        "depth_two_class_selected_count": 17,
        "depth_two_class_skipped_count": 0,
        "pure_120_facet_target_count": 1_407,
        "ridge_valid_pure_target_count": 8,
        "ridge_valid_target_orbit_count": 8,
        "ridge_valid_target_orbit_literal_multiplicities": [1] * 8,
    }
    for key, expected in expected_fields.items():
        if row.get(key) != expected:
            raise ValueError(f"source-2 root-1 field changed: {key}")
    raw_keys = _raw_target_keys(row)
    if len(raw_keys) != 8:
        raise ValueError("source-2 root-1 canonical target count changed")

    downstream = _mapping(json.loads(DOWNSTREAM.read_text()), "downstream audit")
    expected_downstream = {
        "input_target_orbit_count": 8,
        "connected_mod2_homology_sphere_orbit_count": 8,
        "unrestricted_compact_label_no_go_orbit_count": 8,
        "compact_label_survivor_orbit_count": 0,
    }
    for key, expected in expected_downstream.items():
        if downstream.get(key) != expected:
            raise ValueError(f"source-2 root-1 downstream field changed: {key}")
    target_rows = downstream.get("target_orbits")
    if not isinstance(target_rows, list) or len(target_rows) != 8:
        raise ValueError("source-2 root-1 downstream targets changed")
    downstream_keys = {
        tuple(_mapping(target, "downstream target")["canonical_family_key"])
        for target in target_rows
    }
    if downstream_keys != raw_keys:
        raise ValueError("source-2 root-1 downstream target union changed")
    for target in target_rows:
        target = _mapping(target, "downstream target")
        topology = _mapping(target.get("global_topology_audit"), "topology")
        label = _mapping(target.get("sparse_compact_label_gate"), "label gate")
        if topology.get("passes_connected_homology_audit") is not True:
            raise ValueError("source-2 root-1 target failed topology replay")
        if (
            label.get("state_cap_used") is not False
            or label.get("terminal_sparse_state_count") != 0
            or label.get("unrestricted_compact_label_no_go") is not True
        ):
            raise ValueError("source-2 root-1 target label disposition changed")

    hostile = _mapping(json.loads(HOSTILE.read_text()), "hostile downstream audit")
    if (
        hostile.get("kind")
        != "independent_h9_14_record31_a1_downstream_audit"
        or hostile.get("literal_target_orbit_count") != 8
        or hostile.get("compact_label_survivor_count") != 0
        or hostile.get("all_targets_have_unrestricted_compact_label_no_go")
        is not True
        or hostile.get("classical_lanner_placement_counts") != {"4": 108, "5": 420}
    ):
        raise ValueError("source-2 root-1 hostile audit changed")
    hostile_rows = hostile.get("ridge_survivor_label_audits")
    if not isinstance(hostile_rows, list) or len(hostile_rows) != 8:
        raise ValueError("source-2 root-1 hostile target rows changed")
    hostile_families = {
        canonical_family_key(tuple(_mapping(target, "hostile target")["minimal_nonface_masks"]))
        for target in hostile_rows
    }
    if hostile_families != raw_keys:
        raise ValueError("source-2 root-1 hostile target union changed")
    return {
        "root_class_id": 1,
        "root_class_count": 6,
        "incidence_target_orbit_count": 8,
        "topology_survivor_count": 8,
        "compact_label_survivor_count": 0,
        "generated_cpp_sha256": generated_hash,
        "scope": "one of six source-2 root classes only",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
