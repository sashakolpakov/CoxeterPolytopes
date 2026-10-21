#!/usr/bin/env python3
"""Verify the complete partitioned source-2 matrix root-2 census."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

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
        "h9_14_record31_a1_orbit2_matrix_root2_d2_*.json"
    )
)
TARGET_UNION = (
    ROOT
    / "data/enumerations/h9_14_record31_a1_orbit2_matrix_root2_target_union.json"
)
DOWNSTREAM = (
    ROOT
    / "data/enumerations/h9_14_record31_a1_orbit2_matrix_root2_downstream.json"
)
HOSTILE = (
    ROOT
    / "data/enumerations/"
    "h9_14_record31_a1_orbit2_matrix_root2_independent_downstream.json"
)
EXPECTED_GENERATED_CPP_SHA256 = (
    "0bf66b0bd22827f90e4034ba80539eb1cc57b712cdd59e40070b4150e074ddfc"
)


def verify() -> dict[str, object]:
    if len(FRAGMENTS) != 13:
        raise ValueError("source-2 root-2 fragment count changed")
    generated_hash = sha256(generate_depth_two_class_engine().encode()).hexdigest()
    if generated_hash != EXPECTED_GENERATED_CPP_SHA256:
        raise ValueError("source-2 root-2 generated engine changed")

    archived_union = json.loads(TARGET_UNION.read_text())
    recomputed_union = assemble(FRAGMENTS, {2})
    if archived_union != recomputed_union:
        raise ValueError("source-2 root-2 canonical target union changed")
    rows = archived_union.get("orbits")
    if not isinstance(rows, list) or len(rows) != 1:
        raise ValueError("source-2 root-2 union row changed")
    row = rows[0]
    if (
        row.get("root_class_ids") != [2]
        or row.get("fragment_count") != 13
        or row.get("fragment_internal_counts_combined") is not False
        or row.get("raw_fragment_target_orbit_count") != 7
        or row.get("ridge_valid_target_orbit_count") != 3
        or len(row.get("canonical_target_family_keys", [])) != 3
    ):
        raise ValueError("source-2 root-2 union counts changed")

    archived_downstream = json.loads(DOWNSTREAM.read_text())
    if archived_downstream != production_audit(archived_union):
        raise ValueError("source-2 root-2 production downstream replay changed")
    if (
        archived_downstream.get("input_target_orbit_count") != 3
        or archived_downstream.get("connected_mod2_homology_sphere_orbit_count")
        != 3
        or archived_downstream.get("unrestricted_compact_label_no_go_orbit_count")
        != 3
        or archived_downstream.get("compact_label_survivor_orbit_count") != 0
    ):
        raise ValueError("source-2 root-2 production downstream counts changed")

    archived_hostile = json.loads(HOSTILE.read_text())
    if archived_hostile != hostile_audit(archived_union):
        raise ValueError("source-2 root-2 hostile downstream replay changed")
    if (
        archived_hostile.get("literal_target_orbit_count") != 3
        or archived_hostile.get("compact_label_survivor_count") != 0
        or archived_hostile.get("all_targets_have_unrestricted_compact_label_no_go")
        is not True
        or archived_hostile.get("classical_lanner_placement_counts")
        != {"4": 108, "5": 420}
    ):
        raise ValueError("source-2 root-2 hostile downstream counts changed")

    return {
        "root_class_id": 2,
        "depth_two_class_count": 26,
        "terminal_fragment_count": 13,
        "raw_fragment_target_orbit_count": 7,
        "canonical_target_union_count": 3,
        "topology_survivor_count": 3,
        "compact_label_survivor_count": 0,
        "generated_cpp_sha256": generated_hash,
        "target_union_sha256": sha256(TARGET_UNION.read_bytes()).hexdigest(),
        "scope": "root class 2 only; source orbit 2 remains open",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
