#!/usr/bin/env python3
"""Verify the complete partitioned source-2 matrix root-3 census."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from scripts.assemble_h9_14_record31_a1_orbit2_matrix_partition import assemble
from scripts.explore_h9_14_record31_a1_downstream import audit as production_audit
from scripts.generate_h9_14_record31_a1_depth_two_class_engine import (
    generate_depth_two_class_engine,
)
from scripts.independent_h9_14_record31_a1_downstream_audit import (
    audit as hostile_audit,
)


ROOT = Path(__file__).resolve().parents[1]
FRAGMENTS = sorted(
    (ROOT / "data/enumerations").glob(
        "h9_14_record31_a1_orbit2_matrix_root3_d2_*.json"
    )
)
TARGET_UNION = (
    ROOT
    / "data/enumerations/h9_14_record31_a1_orbit2_matrix_root3_target_union.json"
)
DOWNSTREAM = (
    ROOT
    / "data/enumerations/h9_14_record31_a1_orbit2_matrix_root3_downstream.json"
)
HOSTILE = (
    ROOT
    / "data/enumerations/"
    "h9_14_record31_a1_orbit2_matrix_root3_independent_downstream.json"
)
EXPECTED_GENERATED_CPP_SHA256 = (
    "0bf66b0bd22827f90e4034ba80539eb1cc57b712cdd59e40070b4150e074ddfc"
)
EXPECTED_BINARY_SHA256 = (
    "9d4d9cc71e1f7a70bc94169e03d31f4af507dd3fa9ee8033900604c6a07e261d"
)
EXPECTED_WATCHDOG_SHA256 = (
    "f118b000a60db7f04a9e24cabc5c191344ac61d2d8ab743841a6ee7d425f0b90"
)


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"malformed {name}")
    return value


def _verify_terminal_fragments() -> dict[int, int]:
    if len(FRAGMENTS) != 19:
        raise ValueError("source-2 root-3 fragment count changed")
    literal_by_class: dict[int, int] = {}
    for path in FRAGMENTS:
        document = _mapping(json.loads(path.read_text()), path.name)
        rows = document.get("orbits")
        if not isinstance(rows, list) or len(rows) != 1:
            raise ValueError(f"source-2 root-3 fragment row changed: {path.name}")
        row = _mapping(rows[0], path.name)
        begin = int(row["requested_depth_two_class_begin"])
        end = int(row["requested_depth_two_class_end"])
        if (
            int(row["requested_root_class_begin"]) != 3
            or int(row["requested_root_class_end"]) != 3
            or int(row["depth_two_class_discovered_count"]) != 19
            or begin != end
            or begin in literal_by_class
        ):
            raise ValueError(f"source-2 root-3 partition changed: {path.name}")
        literal_by_class[begin] = int(row["ridge_valid_target_orbit_count"])

        provenance = _mapping(
            document.get("terminal_fragment_provenance"),
            f"{path.name} terminal provenance",
        )
        if (
            provenance.get("watchdog_disposition") != "child_exit"
            or provenance.get("watchdog_scope") != "direct"
            or float(provenance["watchdog_memory_ceiling_GiB"]) != 0.6
            or not 0.0 < float(provenance["watchdog_peak_RSS_GiB"]) <= 0.6
            or not 0.0 < float(provenance["watchdog_elapsed_seconds"]) <= 600.0
            or provenance.get("watchdog_timeout_seconds") != 600.0
            or provenance.get("generated_cpp_sha256")
            != EXPECTED_GENERATED_CPP_SHA256
            or provenance.get("binary_sha256") != EXPECTED_BINARY_SHA256
            or provenance.get("watchdog_sha256") != EXPECTED_WATCHDOG_SHA256
        ):
            raise ValueError(f"source-2 root-3 fragment was not terminal: {path.name}")

    if set(literal_by_class) != set(range(1, 20)):
        raise ValueError("source-2 root-3 singleton partition is incomplete")
    if literal_by_class != {**{index: 0 for index in range(1, 20)}, 14: 5}:
        raise ValueError("source-2 root-3 literal target distribution changed")
    return literal_by_class


def verify() -> dict[str, object]:
    literal_by_class = _verify_terminal_fragments()
    generated_hash = sha256(generate_depth_two_class_engine().encode()).hexdigest()
    if generated_hash != EXPECTED_GENERATED_CPP_SHA256:
        raise ValueError("source-2 root-3 generated engine changed")

    archived_union = json.loads(TARGET_UNION.read_text())
    recomputed_union = assemble(FRAGMENTS, {3})
    if archived_union != recomputed_union:
        raise ValueError("source-2 root-3 canonical target union changed")
    rows = archived_union.get("orbits")
    if not isinstance(rows, list) or len(rows) != 1:
        raise ValueError("source-2 root-3 union row changed")
    row = rows[0]
    if (
        row.get("root_class_ids") != [3]
        or row.get("fragment_count") != 19
        or row.get("fragment_internal_counts_combined") is not False
        or row.get("raw_fragment_target_orbit_count") != 5
        or row.get("ridge_valid_target_orbit_count") != 5
        or len(row.get("canonical_target_family_keys", [])) != 5
    ):
        raise ValueError("source-2 root-3 union counts changed")

    archived_downstream = json.loads(DOWNSTREAM.read_text())
    if archived_downstream != production_audit(archived_union):
        raise ValueError("source-2 root-3 production downstream replay changed")
    if (
        archived_downstream.get("input_target_orbit_count") != 5
        or archived_downstream.get("connected_mod2_homology_sphere_orbit_count")
        != 5
        or archived_downstream.get("unrestricted_compact_label_no_go_orbit_count")
        != 5
        or archived_downstream.get("compact_label_survivor_orbit_count") != 0
    ):
        raise ValueError("source-2 root-3 production downstream counts changed")

    archived_hostile = json.loads(HOSTILE.read_text())
    if archived_hostile != hostile_audit(archived_union):
        raise ValueError("source-2 root-3 hostile downstream replay changed")
    if (
        archived_hostile.get("literal_target_orbit_count") != 5
        or archived_hostile.get("compact_label_survivor_count") != 0
        or archived_hostile.get("all_targets_have_unrestricted_compact_label_no_go")
        is not True
        or archived_hostile.get("classical_lanner_placement_counts")
        != {"4": 108, "5": 420}
    ):
        raise ValueError("source-2 root-3 hostile downstream counts changed")

    return {
        "root_class_id": 3,
        "depth_two_class_count": 19,
        "terminal_fragment_count": 19,
        "literal_target_distribution": {
            str(index): count for index, count in sorted(literal_by_class.items())
        },
        "raw_fragment_target_orbit_count": 5,
        "canonical_target_union_count": 5,
        "topology_survivor_count": 5,
        "compact_label_survivor_count": 0,
        "generated_cpp_sha256": generated_hash,
        "target_union_sha256": sha256(TARGET_UNION.read_bytes()).hexdigest(),
        "production_downstream_sha256": sha256(DOWNSTREAM.read_bytes()).hexdigest(),
        "hostile_downstream_sha256": sha256(HOSTILE.read_bytes()).hexdigest(),
        "scope": "root class 3 only; source orbit 2 remains open",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
