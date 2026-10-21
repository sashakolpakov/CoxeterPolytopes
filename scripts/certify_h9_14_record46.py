"""Build and verify the solver-independent H9/14 record-46 exclusion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record46 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record46_structure,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest
from scripts.enumerate_h9_14_four_a5_blocker import enumerate_record


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
BACKEND_SOURCE = ROOT / "data/enumerations/h9_14_facet_completion.json"
STRUCTURE_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record46.py"
CENSUS_SOURCE = ROOT / "scripts/enumerate_h9_14_four_a5_blocker.py"
OUTPUT = ROOT / "data/enumerations/h9_14_record46_certificate.json"
EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_BACKEND_SOURCE_SHA256 = (
    "2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "5f28acdd31ae4c043e799c4dfb62c63dd4b22f95af7032c6b8af83a5fa04b3f3"
)
EXPECTED_HISTOGRAM = {
    48: 1,
    51: 1,
    52: 4,
    54: 2,
    55: 5,
    56: 7,
    57: 2,
    58: 16,
    59: 1,
    60: 24,
    61: 8,
    62: 14,
    63: 10,
    64: 54,
    65: 8,
    66: 22,
    67: 23,
    68: 38,
    69: 6,
    70: 23,
    71: 5,
    72: 22,
    73: 7,
    74: 5,
    75: 1,
    76: 8,
    77: 2,
    79: 2,
    80: 3,
    81: 1,
    82: 2,
    83: 3,
}


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
    return _canonical_sha256(critical)


def _source_record(source: Mapping[str, object]) -> Mapping[str, object]:
    records = source.get("records")
    if not isinstance(records, list):
        raise ValueError("nerve-support source lacks records")
    survivors = [
        record
        for record in records
        if isinstance(record, Mapping) and record.get("survives_second_layer")
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("survivor index 46 has changed")
    return record


def _backend_corroboration() -> dict[str, object]:
    backend = json.loads(BACKEND_SOURCE.read_text(encoding="utf-8"))
    matching = [
        record
        for record in backend["records"]
        if record.get("index") == TARGET_RECORD_INDEX
    ]
    if len(matching) != 1 or matching[0]["solver"]["status"] != "unknown":
        raise ValueError("record-46 bounded backend disposition changed")
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


def _compact_census(full: Mapping[str, object]) -> dict[str, object]:
    fields = (
        "record_index",
        "graph6",
        "missing_edges",
        "distinguished_allowed_vertices",
        "forbidden_nonisolated_vertices",
        "isolated_vertices",
        "active_assignment_row_count_cells",
        "raw_active_isolate_joins",
        "full_local_states_before_row_quotient",
        "distinct_row_states_before_row_quotient",
        "canonical_local_states",
        "blocker_valid_states",
        "required_facet_count",
        "facet_count_histogram",
        "target_state_count",
        "target_states",
    )
    return {field: full[field] for field in fields}


def _audit_census(census: Mapping[str, object]) -> None:
    histogram = {
        int(key): int(value)
        for key, value in census["facet_count_histogram"].items()
    }
    if (
        census["record_index"] != TARGET_RECORD_INDEX
        or census["graph6"] != TARGET_GRAPH6
        or census["distinguished_allowed_vertices"] != [1, 11]
        or census["forbidden_nonisolated_vertices"] != [0, 12, 13]
        or census["isolated_vertices"] != list(range(2, 11))
        or census["active_assignment_row_count_cells"] != 15
        or census["raw_active_isolate_joins"] != 55_458
        or census["full_local_states_before_row_quotient"] != 55_458
        or census["distinct_row_states_before_row_quotient"] != 53_902
        or census["canonical_local_states"] != 2_523
        or census["blocker_valid_states"] != 330
        or census["required_facet_count"] != 90
        or histogram != EXPECTED_HISTOGRAM
        or census["target_state_count"] != 0
        or census["target_states"] != []
    ):
        raise AssertionError("record-46 blocker census did not reproduce")


def build() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    structure = audit_record46_structure(
        TARGET_GRAPH6, record["dehn_sommerville_schedule"]
    ).to_record()
    full_census = enumerate_record(TARGET_RECORD_INDEX)
    census = _compact_census(full_census)
    _audit_census(census)
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record46_solver_independent_exclusion",
        "claim": (
            "H9/14 survivor record 46 cannot be the nerve of a compact "
            "hyperbolic Coxeter 9-polytope."
        ),
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
        "proof_backend": (
            "exact local g-vectors, four-row/16-column-pattern symmetry "
            "compression, "
            "and blocker enumeration; no SAT/SMT"
        ),
        "source": {
            "path": "data/enumerations/h9_14_nerve_support.json",
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
        },
        "proof_sources": [
            {
                "path": "src/coxeter_polytopes/combinatorics/h9_14_record46.py",
                "sha256": _file_sha256(STRUCTURE_SOURCE),
            },
            {
                "path": "scripts/enumerate_h9_14_four_a5_blocker.py",
                "sha256": _file_sha256(CENSUS_SOURCE),
            },
        ],
        "proof_dependencies": [
            "Burcroff's no-missing-triples theorem for compact H9/14 nerves",
            "every compact Coxeter missing face is Lanner of order at most five",
            "strict generalized lower-bound and local Macaulay inequalities",
            "Dehn--Sommerville face-vector schedules",
            (
                "maximal faces are complementary to inclusion-minimal "
                "transversals of the complete minimal-nonface clutter"
            ),
        ],
        "structure": structure,
        "a0_blocker_census": census,
        "full_census_sha256": _canonical_sha256(full_census),
        "backend_corroboration": _backend_corroboration(),
        "updated_theorem_partition": {
            "independently_certified_impossible": 119,
            "certified_explicit_ordinary_polytopal": 2,
            "theorem_undecided": 21,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "This excludes one graph shadow, not all compact H9/14 nerves.",
            "The archived bounded Z3 timeout contributes no logical step.",
            "No Coxeter labels or Lorentzian Gram matrices are classified.",
            "This is not a general non-existence theorem in dimension 9.",
        ],
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(
    manifest: Mapping[str, object], *, recompute: bool = False
) -> dict[str, object]:
    if manifest.get("schema_version") != "1.0.0" or manifest.get("kind") != (
        "h9_14_record46_solver_independent_exclusion"
    ):
        raise ValueError("unexpected record-46 certificate metadata")
    if manifest.get("record_index") != TARGET_RECORD_INDEX or manifest.get(
        "graph6"
    ) != TARGET_GRAPH6:
        raise ValueError("record-46 target mismatch")
    if manifest.get("theorem_status") != "certified_unsat":
        raise ValueError("record-46 theorem status was not promoted")
    if _file_sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ValueError("nerve-support source is not frozen")
    if _file_sha256(BACKEND_SOURCE) != EXPECTED_BACKEND_SOURCE_SHA256:
        raise ValueError("backend source is not frozen")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    expected_structure = audit_record46_structure(
        TARGET_GRAPH6, record["dehn_sommerville_schedule"]
    ).to_record()
    if manifest.get("structure") != expected_structure:
        raise ValueError("record-46 structural proof does not verify")
    expected_sources = [
        {
            "path": "src/coxeter_polytopes/combinatorics/h9_14_record46.py",
            "sha256": _file_sha256(STRUCTURE_SOURCE),
        },
        {
            "path": "scripts/enumerate_h9_14_four_a5_blocker.py",
            "sha256": _file_sha256(CENSUS_SOURCE),
        },
    ]
    if manifest.get("proof_sources") != expected_sources:
        raise ValueError("record-46 proof-source provenance mismatch")
    census = manifest.get("a0_blocker_census")
    if not isinstance(census, Mapping):
        raise ValueError("record-46 certificate lacks blocker census")
    _audit_census(census)
    if manifest.get("backend_corroboration") != _backend_corroboration():
        raise ValueError("record-46 backend corroboration mismatch")
    backend = manifest["backend_corroboration"]
    if not isinstance(backend, Mapping) or backend.get("trusted_for_proof") is not False:
        raise ValueError("bounded backend timeout was incorrectly trusted")
    if recompute:
        full_census = enumerate_record(TARGET_RECORD_INDEX)
        if _canonical_sha256(_compact_census(full_census)) != (
            _canonical_sha256(census)
        ):
            raise ValueError("record-46 blocker census does not reproduce")
        if _canonical_sha256(full_census) != manifest.get("full_census_sha256"):
            raise ValueError("record-46 full-census digest mismatch")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-46 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-46 certificate is not frozen")
    return expected_structure


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
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "certificate_sha256": manifest["certificate_sha256"],
                "verified": proof["verified"],
                "recomputed": bool(args.recompute or not args.verify_only),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
