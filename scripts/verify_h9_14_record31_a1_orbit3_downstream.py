#!/usr/bin/env python3
"""Replay and verify both downstream audits for Record-31 a=1 source 3."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from scripts.explore_h9_14_record31_a1_downstream import audit as production_audit
from scripts.independent_h9_14_record31_a1_downstream_audit import (
    audit as independent_audit,
)
from scripts.verify_h9_14_record31_a1_orbit3_transcript import (
    ARCHIVE as INCIDENCE,
    EXPECTED_FAMILIES,
)


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record31_a1_orbit3_downstream.json"
)
INDEPENDENT = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record31_a1_orbit3_independent_downstream.json"
)
EXPECTED_PRODUCTION_SHA256 = (
    "1562231aedd9f3fc138aeb31c6af48b01cd70eee30fa32e12bab1987c9c5184b"
)
EXPECTED_INDEPENDENT_SHA256 = (
    "efd2f975f79bc3d28215e952269165b934d3744fd6e9e7a084810b93f07542af"
)


def _load(path: Path) -> Mapping[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise AssertionError(f"{path} is not a JSON object")
    return value


def verify(*, replay: bool = True) -> dict[str, object]:
    incidence = _load(INCIDENCE)
    production = _load(PRODUCTION)
    independent = _load(INDEPENDENT)
    production_digest = sha256(PRODUCTION.read_bytes()).hexdigest()
    independent_digest = sha256(INDEPENDENT.read_bytes()).hexdigest()
    if production_digest != EXPECTED_PRODUCTION_SHA256:
        raise AssertionError("production downstream archive bytes changed")
    if independent_digest != EXPECTED_INDEPENDENT_SHA256:
        raise AssertionError("independent downstream archive bytes changed")

    if replay:
        if production_audit(incidence) != production:
            raise AssertionError("production downstream replay differs")
        if independent_audit(incidence) != independent:
            raise AssertionError("independent downstream replay differs")

    expected_summary = {
        "input_target_orbit_count": 11,
        "connected_mod2_homology_sphere_orbit_count": 11,
        "unrestricted_compact_label_no_go_orbit_count": 11,
        "compact_label_survivor_orbit_count": 0,
    }
    for key, expected in expected_summary.items():
        if production.get(key) != expected:
            raise AssertionError(f"production downstream summary {key} changed")
    if (
        independent.get("literal_target_orbit_count") != 11
        or independent.get("all_targets_have_unrestricted_compact_label_no_go")
        is not True
        or independent.get("compact_label_survivor_count") != 0
        or independent.get("classical_lanner_placement_counts")
        != {"4": 108, "5": 420}
    ):
        raise AssertionError("independent downstream summary changed")

    production_rows = production.get("target_orbits")
    independent_rows = independent.get("ridge_survivor_label_audits")
    if not isinstance(production_rows, list) or not isinstance(independent_rows, list):
        raise AssertionError("downstream target rows changed")
    if not (
        len(production_rows) == len(independent_rows) == len(EXPECTED_FAMILIES)
    ):
        raise AssertionError("downstream row count changed")
    for index, expected_family in enumerate(EXPECTED_FAMILIES):
        production_row = production_rows[index]
        independent_row = independent_rows[index]
        if not isinstance(production_row, Mapping) or not isinstance(
            independent_row, Mapping
        ):
            raise AssertionError("downstream row is malformed")
        if production_row.get("minimal_nonface_masks") != list(expected_family):
            raise AssertionError("production downstream family changed")
        if independent_row.get("minimal_nonface_masks") != list(expected_family):
            raise AssertionError("independent downstream family changed")
        topology = production_row.get("global_topology_audit")
        labels = production_row.get("sparse_compact_label_gate")
        if not isinstance(topology, Mapping) or not isinstance(labels, Mapping):
            raise AssertionError("production downstream gates are missing")
        if topology.get("passes_connected_homology_audit") is not True:
            raise AssertionError("a topology-valid source-3 target changed")
        if (
            labels.get("state_cap_used") is not False
            or labels.get("terminal_sparse_state_count") != 0
            or labels.get("unrestricted_compact_label_no_go") is not True
        ):
            raise AssertionError("production compact-label no-go changed")
        if independent_row.get("terminal_sparse_state_count") != 0:
            raise AssertionError("independent compact-label no-go changed")

    return {
        "kind": "verified_h9_14_record31_a1_orbit3_downstream",
        "production_downstream_sha256": production_digest,
        "independent_downstream_sha256": independent_digest,
        "target_orbit_count": 11,
        "topology_survivor_count": 11,
        "compact_label_survivor_count": 0,
        "gram_candidate_count": 0,
        "production_and_independent_replayed": replay,
        "scope": (
            "Complete compact downstream no-go for source orbit 3 only; not "
            "a full Record-31 a=1 theorem or a noncompact finite-volume claim."
        ),
    }


def main() -> None:
    print(json.dumps(verify(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
