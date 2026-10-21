#!/usr/bin/env python3
"""Freeze and verify the positive abstract completion of H9/14 record 43.

Certificate verification is solver-independent.  The optional build path uses
Z3 only to find 1001 Gale chirotope signs; the archived signs are then checked
directly against all 45,045 Grassmann--Pluecker relations and all 2,002
five-set circuit predicates.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.gale_chirotope import (
    EXPECTED_CHIROTOPE_SIGN_COUNT,
    EXPECTED_GRASSMANN_PLUECKER_COUNT,
    audit_uniform_rank_four_gale,
    solve_uniform_rank_four_gale,
)
from coxeter_polytopes.combinatorics.h9_14_record43 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_RECORD_INDEX,
    audit_record43_incidence,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record43_selection.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "af521ed4e7476c95945df9b770f36a343dcf3fc79dabed557aaab0fd97a671a3"
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def critical_digest(payload: Mapping[str, object]) -> str:
    critical = {
        "kind": payload["kind"],
        "claim": payload["claim"],
        "source": payload["source"],
        "incidence": payload["incidence"],
        "global_selection_audit": payload["global_selection_audit"],
        "gale_chirotope_signs": payload["gale_chirotope_signs"],
        "gale_audit": payload["gale_audit"],
    }
    encoded = json.dumps(critical, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def support_record() -> tuple[dict[str, object], str]:
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    verify_support_manifest(support)
    survivors = [
        record for record in support["records"] if record["survives_second_layer"]
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("second-layer survivor record 43 changed")
    return record, file_sha256(SUPPORT)


def gale_audit_record(facets: tuple[int, ...], signs: Sequence[int]) -> dict[str, object]:
    audit = audit_uniform_rank_four_gale(facets, signs)
    if not audit.passes:
        raise AssertionError("archived record-43 Gale signs failed independent audit")
    return {
        "rank": 4,
        "uniform": True,
        "sign_count": audit.sign_count,
        "grassmann_pluecker_relation_count": EXPECTED_GRASSMANN_PLUECKER_COUNT,
        "first_grassmann_pluecker_failure": audit.first_grassmann_pluecker_failure,
        "positive_circuit_count": len(audit.positive_circuit_masks),
        "positive_circuit_masks": list(audit.positive_circuit_masks),
        "positive_circuits_match_facets": audit.positive_circuits_match_facets,
        "every_element_in_positive_circuit": audit.every_element_in_positive_circuit,
        "every_element_avoided_by_positive_circuit": (
            audit.every_element_avoided_by_positive_circuit
        ),
        "independent_sign_audit_passes": audit.passes,
        "scope": (
            "Necessary uniform oriented-matroid feasibility only; neither real "
            "realizability nor polytopality follows."
        ),
    }


def build(signs: Sequence[int] | None = None) -> dict[str, object]:
    source, support_sha = support_record()
    incidence = audit_record43_incidence(TARGET_GRAPH6)
    facets = incidence.facet_complement_masks
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        facets,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("record-43 completion failed the full independent audit")
    if global_audit.face_vector != TARGET_F_VECTOR:
        raise AssertionError("record-43 face vector changed")
    expected_nonfaces = {
        2: tuple(sorted(sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES)),
        5: tuple(sorted(sum(1 << vertex for vertex in face) for face in TARGET_MISSING_FIVES)),
    }
    if global_audit.minimal_nonfaces_by_size != expected_nonfaces:
        raise AssertionError("generated face table has unexpected minimal nonfaces")

    search: dict[str, object]
    if signs is None:
        result = solve_uniform_rank_four_gale(facets, timeout_ms=300_000)
        if result.status != "sat" or len(result.signs) != EXPECTED_CHIROTOPE_SIGN_COUNT:
            raise RuntimeError(
                f"record-43 Gale search did not finish SAT: {result.status} "
                f"({result.reason_unknown})"
            )
        normalized_signs = result.signs
        search = {
            "backend": result.backend,
            "backend_version": result.backend_version,
            "status": result.status,
            "assertion_count": result.assertion_count,
            "timeout_ms": 300_000,
            "role": "search provenance only; the signs are audited independently",
        }
    else:
        normalized_signs = tuple(int(sign) for sign in signs)
        search = {
            "backend": "archived-sign replay",
            "status": "not rerun",
            "role": "verification imports no solver",
        }
    gale_audit = gale_audit_record(facets, normalized_signs)

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record43_positive_abstract_completion",
        "claim": (
            "H9/14 survivor record 43 has an explicit 90-facet abstract "
            "completion passing the full incidence and mod-two link-homology "
            "audits and the uniform rank-four Gale oriented-matroid gate."
        ),
        "nonclaims": [
            "The abstract completion is not proved to be a PL sphere.",
            "The archived Gale oriented matroid is not proved realizable over R.",
            "No convex 9-polytope realization is supplied.",
            "No Coxeter labeling or Lorentzian Gram realization is supplied.",
        ],
        "source": {
            "support_path": "data/enumerations/h9_14_nerve_support.json",
            "support_sha256": support_sha,
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
        },
        "incidence": incidence.to_record(),
        "global_selection_audit": global_audit.to_record(),
        "gale_search": search,
        "gale_chirotope_signs": list(normalized_signs),
        "gale_sign_sha256": sha256(
            bytes(1 if sign > 0 else 0 for sign in normalized_signs)
        ).hexdigest(),
        "gale_audit": gale_audit,
    }
    manifest["certificate_sha256"] = critical_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_record43_positive_abstract_completion":
        raise ValueError("wrong record-43 certificate kind")
    signs = manifest.get("gale_chirotope_signs")
    if not isinstance(signs, list):
        raise ValueError("record-43 certificate lacks Gale signs")
    rebuilt = build(signs)
    for key in (
        "claim",
        "nonclaims",
        "source",
        "incidence",
        "global_selection_audit",
        "gale_chirotope_signs",
        "gale_sign_sha256",
        "gale_audit",
        "certificate_sha256",
    ):
        if manifest.get(key) != rebuilt.get(key):
            raise ValueError(f"record-43 certificate field does not verify: {key}")
    digest = critical_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-43 critical digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-43 certificate digest is not frozen")
    return {
        "facet_count": len(rebuilt["incidence"]["facet_complement_masks"]),
        "face_vector": rebuilt["global_selection_audit"]["face_vector"],
        "gale_sign_count": len(signs),
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-certificate",
        action="store_true",
        help="run the Gale search and write a new certificate",
    )
    args = parser.parse_args(argv)
    if args.write_certificate:
        manifest = build()
        OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    else:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    print(json.dumps({"output": str(OUTPUT), **summary}, indent=2))


if __name__ == "__main__":
    main()
