#!/usr/bin/env python3
"""Build and verify the solver-independent H9/14 record-37 exclusion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest
from scripts.audit_h9_14_record37_targets import audit_target_families
from scripts.enumerate_h9_14_record37_a0 import (
    TARGET_GRAPH6,
    TARGET_INDEX,
    enumerate_a0,
)
from scripts.explore_h9_14_record37 import enumerate_positive_branches


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
BACKEND_SOURCE = ROOT / "data/enumerations/h9_14_facet_completion.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record37_certificate.json"
PROOF_SOURCES = (
    ROOT / "scripts/explore_h9_14_record37.py",
    ROOT / "scripts/enumerate_h9_14_record37_a0.py",
    ROOT / "scripts/audit_h9_14_record37_targets.py",
    ROOT / "scripts/audit_h9_14_four_a5_candidate.py",
)
EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_BACKEND_SOURCE_SHA256 = (
    "2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "4da8617bd4dfc77889c3c3f987e7988efd87d6d366d7efd29c7403370536cc5c"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    # Canonicalize through JSON first: JSON object keys are strings, while
    # census maps naturally use integer keys before the archive is written.
    # Without this round trip the embedded digest changes after reloading the
    # very JSON file it is meant to authenticate.
    json_normal = json.loads(json.dumps(critical))
    return _canonical_sha256(json_normal)


def _source_record(source: Mapping[str, object]) -> Mapping[str, object]:
    records = source.get("records")
    if not isinstance(records, list):
        raise ValueError("nerve-support source lacks records")
    survivors = [
        record
        for record in records
        if isinstance(record, Mapping) and record.get("survives_second_layer")
    ]
    record = survivors[TARGET_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("survivor index 37 has changed")
    return record


def _backend_corroboration() -> dict[str, object]:
    backend = json.loads(BACKEND_SOURCE.read_text(encoding="utf-8"))
    matching = [
        record
        for record in backend["records"]
        if record.get("index") == TARGET_INDEX
    ]
    if len(matching) != 1 or matching[0]["solver"]["status"] != "unknown":
        raise ValueError("record-37 bounded backend disposition changed")
    solver = matching[0]["solver"]
    return {
        "path": "data/enumerations/h9_14_facet_completion.json",
        "sha256": _file_sha256(BACKEND_SOURCE),
        "status": solver["status"],
        "reason_unknown": solver["reason_unknown"],
        "backend": solver["backend"],
        "backend_version": solver["backend_version"],
        "assertion_count": solver["assertion_count"],
        "elapsed_seconds": solver["elapsed_seconds"],
        "trusted_for_proof": False,
    }


def _proof_source_records() -> list[dict[str, str]]:
    return [
        {
            "path": str(path.relative_to(ROOT)),
            "sha256": _file_sha256(path),
        }
        for path in PROOF_SOURCES
    ]


def _compact_a0_census(full: Mapping[str, object]) -> dict[str, object]:
    compact = dict(full)
    compact.pop("target-states", None)
    return compact


def _audit_summary(manifest: Mapping[str, object]) -> None:
    positive = manifest.get("positive_a_branches")
    census = manifest.get("a0_blocker_census")
    target = manifest.get("target_count_state_audits")
    if not isinstance(positive, Mapping) or not isinstance(census, Mapping):
        raise ValueError("record-37 certificate lacks branch censuses")
    if not isinstance(target, Mapping):
        raise ValueError("record-37 certificate lacks target-state audits")
    branches = positive.get("branches")
    if not isinstance(branches, Mapping):
        raise ValueError("record-37 positive census lacks branches")
    a1 = branches.get("1")
    a2 = branches.get("2")
    if not isinstance(a1, Mapping) or not isinstance(a2, Mapping):
        raise ValueError("record-37 positive census has wrong branches")
    if (
        positive.get("individually_interval_feasible_a4") != 462
        or a1.get("total_a4_families") != 462
        or a1.get("locally_feasible_a4_families") != 462
        or a1.get("shadow_strict_families") != 0
        or a1.get("forced_shadow_count_overrun_families") != 0
        or a2.get("total_a4_families") != 106_491
        or a2.get("locally_feasible_a4_families") != 96_033
        or a2.get("forced_shadow_count_overrun_families") != 51_429
        or a2.get("shadow_strict_families") != 0
    ):
        raise AssertionError("record-37 positive branch census changed")
    histogram = {
        int(key): int(value)
        for key, value in census["facet-count-histogram"].items()
    }
    observed_a0 = {
        "quotient": census.get("row-normalized-isolate-quotient-states"),
        "distinct": census.get("distinct-row-states"),
        "blocker": census.get("blocker-valid-states"),
        "required": census.get("required-facet-count"),
        "targets": census.get("target-state-count"),
        "histogram_target": histogram.get(102),
        "histogram_sum": sum(histogram.values()),
        "histogram_min": min(histogram),
        "histogram_max": max(histogram),
    }
    expected_a0 = {
        "quotient": 2_865_450,
        "distinct": 2_774_810,
        "blocker": 522_427,
        "required": 102,
        "targets": 12,
        "histogram_target": 12,
        "histogram_sum": 522_427,
        "histogram_min": 40,
        "histogram_max": 105,
    }
    if observed_a0 != expected_a0:
        raise AssertionError(
            f"record-37 a=0 census changed: {observed_a0!r}"
        )
    families = target.get("families")
    raw_map = target.get("raw_to_unique_family_map")
    if (
        target.get("raw_target_state_count") != 12
        or target.get("unique_target_family_count") != 5
        or target.get("passing_family_count") != 0
        or not isinstance(families, list)
        or len(families) != 5
        or not isinstance(raw_map, list)
        or len(raw_map) != 12
    ):
        raise AssertionError("record-37 target-state deduplication changed")
    expected_ridges = (
        {1: 28, 2: 403, 3: 28},
        {1: 24, 2: 411, 3: 24},
    )
    for family in families:
        blocker_histogram = {
            int(key): int(value)
            for key, value in family["minimal_transversal_size_histogram"].items()
        }
        ridge_histogram = {
            int(key): int(value)
            for key, value in family["ridge_degree_histogram"].items()
        }
        if (
            blocker_histogram != {5: 102}
            or len(family["facet_complement_masks"]) != 102
            or ridge_histogram not in expected_ridges
            or family["all_ridge_degrees_two"] is not False
            or family["passes_selection_model"] is not False
        ):
            raise AssertionError("record-37 explicit blocker audit changed")


def build() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    schedules = {
        schedule["missing_four_count"]: {
            "missing_five_count_offset": schedule["missing_five_count_offset"],
            "f_vector": list(schedule["f_vector"]),
        }
        for schedule in record["dehn_sommerville_schedule"]
    }
    positive = enumerate_positive_branches()
    a0_full = enumerate_a0()
    target_audits = audit_target_families(a0_full)
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record37_solver_independent_exclusion",
        "claim": (
            "H9/14 survivor record 37 cannot be the nerve of a compact "
            "hyperbolic Coxeter 9-polytope."
        ),
        "record_index": TARGET_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
        "proof_backend": (
            "exact local g-vectors, five-row isolate-column compression, "
            "blocker duality, and complete ridge replay; no SAT/SMT"
        ),
        "source": {
            "path": "data/enumerations/h9_14_nerve_support.json",
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_INDEX,
        },
        "proof_sources": _proof_source_records(),
        "proof_dependencies": [
            "Burcroff's no-missing-triples theorem for compact H9/14 nerves",
            "every compact Coxeter missing face is Lanner of order at most five",
            "strict generalized lower-bound and local Macaulay inequalities",
            "Dehn--Sommerville face-vector schedules",
            (
                "maximal faces are complementary to inclusion-minimal "
                "transversals of the complete minimal-nonface clutter"
            ),
            "every ridge of a simplicial polytope boundary has degree two",
        ],
        "schedules": schedules,
        "high_a_exclusion": {
            "excluded_missing_four_counts": [3, 4],
            "summed_vertex_link_g4": {"3": 12, "4": 6},
            "required_minimum": 14,
        },
        "positive_a_branches": positive,
        "a0_blocker_census": _compact_a0_census(a0_full),
        "a0_full_census_sha256": _canonical_sha256(a0_full),
        "target_count_state_audits": target_audits,
        "target_count_audit_sha256": _canonical_sha256(target_audits),
        "backend_corroboration": _backend_corroboration(),
        "updated_theorem_partition": {
            "independently_certified_impossible": 120,
            "certified_explicit_ordinary_polytopal": 5,
            "theorem_undecided": 17,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "This excludes one graph shadow, not all compact H9/14 nerves.",
            "The archived bounded Z3 timeout contributes no logical step.",
            "No Coxeter labels or Lorentzian Gram matrices are classified.",
            "This is not a general non-existence theorem in dimension 9.",
        ],
    }
    _audit_summary(manifest)
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(
    manifest: Mapping[str, object], *, recompute: bool = False
) -> dict[str, object]:
    if manifest.get("schema_version") != "1.0.0" or manifest.get("kind") != (
        "h9_14_record37_solver_independent_exclusion"
    ):
        raise ValueError("unexpected record-37 certificate metadata")
    if manifest.get("record_index") != TARGET_INDEX or manifest.get(
        "graph6"
    ) != TARGET_GRAPH6:
        raise ValueError("record-37 target mismatch")
    if manifest.get("theorem_status") != "certified_unsat":
        raise ValueError("record-37 theorem status was not promoted")
    if _file_sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ValueError("nerve-support source is not frozen")
    if _file_sha256(BACKEND_SOURCE) != EXPECTED_BACKEND_SOURCE_SHA256:
        raise ValueError("backend source is not frozen")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    _source_record(source)
    if manifest.get("proof_sources") != _proof_source_records():
        raise ValueError("record-37 proof-source provenance mismatch")
    if manifest.get("backend_corroboration") != _backend_corroboration():
        raise ValueError("record-37 backend corroboration mismatch")
    backend = manifest["backend_corroboration"]
    if not isinstance(backend, Mapping) or backend.get("trusted_for_proof") is not False:
        raise ValueError("bounded backend timeout was incorrectly trusted")
    _audit_summary(manifest)
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-37 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-37 certificate is not frozen")
    if recompute:
        rebuilt = build()
        # JSON object keys are strings after archive reload, whereas several
        # exact census maps naturally use integer keys in memory.  Compare
        # their JSON normal forms so a faithful round trip is not mistaken
        # for census drift.
        rebuilt_json = json.loads(json.dumps(rebuilt, sort_keys=True))
        manifest_json = json.loads(json.dumps(manifest, sort_keys=True))
        if rebuilt_json != manifest_json:
            raise ValueError("record-37 full census does not reproduce")
    return {
        "record_index": TARGET_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": manifest["theorem_status"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    else:
        manifest = build()
        OUTPUT.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    proof = verify(manifest, recompute=args.recompute)
    print(json.dumps({"output": str(OUTPUT), **proof}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
