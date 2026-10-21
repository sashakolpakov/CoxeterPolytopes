#!/usr/bin/env python3
"""Freeze and verify the ``a=0`` completion of H9/14 record 32."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_record32_a0 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_RECORD_INDEX,
    TARGET_SCHEDULE_INDEX,
    audit_record32_a0_common_four_core,
    audit_record32_a0_incidence,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record32_a0_selection.json"
DEPENDENCIES = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record32_a0.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record35_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
)
EXPECTED_CERTIFICATE_SHA256 = (
    "14c3f9f5a425eb6934d7f53a7a7f071467eba7cf2e4a64b7cf9a194d3490d0bb"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _critical_digest(payload: Mapping[str, object]) -> str:
    critical = {
        key: payload[key]
        for key in (
            "kind",
            "claim",
            "nonclaims",
            "source",
            "incidence",
            "global_selection_audit",
            "common_four_core_coxeter_gate",
        )
    }
    encoded = json.dumps(critical, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def _support_record() -> tuple[dict[str, object], str]:
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    verify_support_manifest(support)
    survivors = [
        record for record in support["records"] if record["survives_second_layer"]
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("second-layer survivor record 32 changed")
    return record, _file_sha256(SUPPORT)


def build() -> dict[str, object]:
    source, support_sha = _support_record()
    incidence = audit_record32_a0_incidence(TARGET_GRAPH6)
    facets = incidence.facet_complement_masks
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        facets,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("record-32 a=0 completion failed the full global audit")
    if global_audit.face_vector != TARGET_F_VECTOR:
        raise AssertionError("record-32 a=0 face vector changed")
    if global_audit.matching_schedule_index != TARGET_SCHEDULE_INDEX:
        raise AssertionError("record-32 a=0 Dehn--Sommerville branch changed")
    expected_nonfaces = {
        2: tuple(
            sorted(sum(1 << vertex for vertex in face) for face in TARGET_MISSING_EDGES)
        ),
        5: tuple(
            sorted(sum(1 << vertex for vertex in face) for face in TARGET_MISSING_FIVES)
        ),
    }
    if global_audit.minimal_nonfaces_by_size != expected_nonfaces:
        raise AssertionError("record-32 a=0 generated face table has wrong clutter")

    coxeter_gate = audit_record32_a0_common_four_core()
    if not coxeter_gate.unrestricted_external_label_no_go:
        raise AssertionError("record-32 a=0 compact-Coxeter gate failed")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record32_a0_abstract_completion_and_coxeter_no_go",
        "claim": (
            "The displayed a=0 clutter over H9/14 record 32 has a pure "
            "114-member blocker and passes the complete incidence, ridge, "
            "connectedness, Euler-link, and mod-2 homology-link audit.  This "
            "selected nerve cannot be the nerve of a compact hyperbolic "
            "Coxeter polytope: its common-four-core Lanner cluster has no "
            "locally elliptic labeling, even with arbitrary remaining labels."
        ),
        "nonclaims": [
            "No ordinary convex realization of this abstract completion is asserted.",
            "No Lorentzian Gram realization or finite-volume Coxeter realization is supplied.",
            "The no-go concerns this selected nerve, not every completion of the same missing-edge graph shadow.",
            "No uniqueness claim is made for the a=0 or record-32 completion problem.",
        ],
        "source": {
            "support_path": "data/enumerations/h9_14_nerve_support.json",
            "support_sha256": support_sha,
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "dehn_sommerville_schedule_index": TARGET_SCHEDULE_INDEX,
            "graph6": TARGET_GRAPH6,
            "dependency_sha256": {
                str(path.relative_to(ROOT)): _file_sha256(path)
                for path in DEPENDENCIES
            },
            "search_provenance": (
                "A finite blocker search found the displayed a=0 clutter; "
                "every archived assertion replays from that explicit clutter."
            ),
        },
        "incidence": incidence.to_record(),
        "global_selection_audit": global_audit.to_record(),
        "common_four_core_coxeter_gate": coxeter_gate.to_record(),
    }
    manifest["certificate_sha256"] = _critical_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != (
        "h9_14_record32_a0_abstract_completion_and_coxeter_no_go"
    ):
        raise ValueError("wrong record-32 a=0 certificate kind")
    rebuilt = build()
    for key in (
        "claim",
        "nonclaims",
        "source",
        "incidence",
        "global_selection_audit",
        "common_four_core_coxeter_gate",
        "certificate_sha256",
    ):
        if manifest.get(key) != rebuilt.get(key):
            raise ValueError(f"record-32 a=0 field does not verify: {key}")
    digest = _critical_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-32 a=0 critical digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-32 a=0 certificate digest is not frozen")
    gate = rebuilt["common_four_core_coxeter_gate"]
    return {
        "facet_count": len(rebuilt["incidence"]["facet_complement_masks"]),
        "face_vector": rebuilt["global_selection_audit"]["face_vector"],
        "rank_six_common_core_survivors": gate["sparse_shadow_passing_state_count"],
        "unrestricted_compact_coxeter_no_go": gate[
            "unrestricted_external_label_no_go"
        ],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-certificate",
        action="store_true",
        help="rebuild and write the deterministic JSON certificate",
    )
    args = parser.parse_args(argv)
    if args.write_certificate:
        manifest = build()
        OUTPUT.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps({"output": str(OUTPUT), **verify(manifest)}, indent=2))


if __name__ == "__main__":
    main()
