#!/usr/bin/env python3
"""Freeze and verify the audited abstract H9/14 record-18 candidate nerve.

The optional build path asks Z3 only for a rank-four Gale sign word.  Normal
verification is solver-independent: it rebuilds the complete blocker and
face table, checks all links, and audits every archived chirotope relation.
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
from coxeter_polytopes.combinatorics.h9_14_record18 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_RECORD_INDEX,
    audit_record18_incidence,
    audit_record18_structure,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record18_selection.json"

PROOF_SOURCES = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record18.py",
    ROOT / "src/coxeter_polytopes/combinatorics/facet_completion.py",
    ROOT / "src/coxeter_polytopes/combinatorics/gale_chirotope.py",
    ROOT / "src/coxeter_polytopes/combinatorics/local_link.py",
    ROOT / "src/coxeter_polytopes/combinatorics/missing_edges.py",
    ROOT / "src/coxeter_polytopes/combinatorics/nerve_support.py",
    ROOT / "src/coxeter_polytopes/combinatorics/simplicial_homology.py",
    ROOT / "scripts/independent_h9_14_record18_audit.py",
)

EXPECTED_SUPPORT_SHA256 = "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
EXPECTED_PROOF_SOURCE_SHA256 = {
    "src/coxeter_polytopes/combinatorics/h9_14_record18.py": "c96d21616004da447a43c6c803a1b27583a2c997454eba64ccafc2eeee512e15",
    "src/coxeter_polytopes/combinatorics/facet_completion.py": "5f2d078cfe400b39f3e194352f7a0c6d3bd8cb6bb8e7cee37be159e8766d9400",
    "src/coxeter_polytopes/combinatorics/gale_chirotope.py": "c051432fc182d00079633998bf554b798b6625e2aae6c0f4d210d6a6034017db",
    "src/coxeter_polytopes/combinatorics/local_link.py": "8f8169edab8b87aff682785667c2a8eab5c52b92a7c68cef130aeec4d3ae807d",
    "src/coxeter_polytopes/combinatorics/missing_edges.py": "231e85fa98e27e864f997cd4ebe59136707f6aa8dd2d8136b80cff543dfbfbb8",
    "src/coxeter_polytopes/combinatorics/nerve_support.py": "aad8fb25c2f9a56e45121fa1d1cbbdedfdff8ba1f67e379be74eed73940ad2a2",
    "src/coxeter_polytopes/combinatorics/simplicial_homology.py": "547f51b58c1146aac80c6e08910754bf93ace5fa808dae5a342da7e042958110",
    "scripts/independent_h9_14_record18_audit.py": "a6be2665517ed97ae8c76b19d0063ff7f26a6b1bbcdfa7f5e2c5398a40c00d7a",
}
EXPECTED_CERTIFICATE_SHA256 = "27193c530804237c8505b7ce7400e7551a6b152d5a4391f2dc624371d8f54406"


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def certificate_digest(payload: Mapping[str, object]) -> str:
    normalized = dict(payload)
    normalized.pop("certificate_sha256", None)
    encoded = json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def proof_source_records() -> list[dict[str, str]]:
    return [
        {"path": str(path.relative_to(ROOT)), "sha256": file_sha256(path)}
        for path in PROOF_SOURCES
    ]


def support_record() -> tuple[dict[str, object], str]:
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    verify_support_manifest(support)
    survivors = [record for record in support["records"] if record["survives_second_layer"]]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("second-layer survivor record 18 changed")
    return record, file_sha256(SUPPORT)


def gale_audit_record(facets: tuple[int, ...], signs: Sequence[int]) -> dict[str, object]:
    audit = audit_uniform_rank_four_gale(facets, signs)
    if not audit.passes:
        raise AssertionError("archived record-18 Gale signs failed exact audit")
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


def build(
    signs: Sequence[int] | None = None,
    *,
    search_provenance: Mapping[str, object] | None = None,
) -> dict[str, object]:
    source, support_sha = support_record()
    incidence = audit_record18_incidence(TARGET_GRAPH6)
    structure = audit_record18_structure()
    facets = incidence.facet_complement_masks
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        facets,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("record-18 candidate failed the complete face/link audit")
    if global_audit.face_vector != TARGET_F_VECTOR:
        raise AssertionError("record-18 face vector changed")
    expected_nonfaces = {
        2: tuple(sorted(sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES)),
        5: tuple(sorted(sum(1 << vertex for vertex in face) for face in TARGET_MISSING_FIVES)),
    }
    if global_audit.minimal_nonfaces_by_size != expected_nonfaces:
        raise AssertionError("generated face table has unexpected minimal nonfaces")

    if signs is None:
        result = solve_uniform_rank_four_gale(facets, timeout_ms=300_000)
        if result.status != "sat" or len(result.signs) != EXPECTED_CHIROTOPE_SIGN_COUNT:
            raise RuntimeError(
                f"record-18 Gale search did not finish SAT: {result.status} "
                f"({result.reason_unknown})"
            )
        normalized_signs = result.signs
        search: dict[str, object] = {
            "backend": result.backend,
            "backend_version": result.backend_version,
            "status": result.status,
            "assertion_count": result.assertion_count,
            "elapsed_seconds": result.elapsed_seconds,
            "timeout_ms": 300_000,
            "trusted_for_proof": False,
            "role": "search provenance only; archived signs are audited directly",
        }
    else:
        normalized_signs = tuple(int(sign) for sign in signs)
        if search_provenance is None:
            search = {
                "backend": "archived-sign replay",
                "status": "not rerun",
                "trusted_for_proof": False,
                "role": "verification imports no solver",
            }
        else:
            search = dict(search_provenance)
    gale_audit = gale_audit_record(facets, normalized_signs)

    schedules = [
        {
            "missing_four_count": schedule["missing_four_count"],
            "missing_five_count_offset": schedule["missing_five_count_offset"],
            "facet_count": schedule["f_vector"][-1],
        }
        for schedule in source["dehn_sommerville_schedule"]
    ]
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record18_abstract_audited_nerve_candidate",
        "claim": (
            "H9/14 survivor record 18 has an explicit 110-facet abstract nerve "
            "candidate passing the complete blocker, pseudomanifold, local-g, "
            "Macaulay, Euler-link, mod-two link-homology, and uniform rank-four "
            "oriented-matroid audits."
        ),
        "status": "abstract_audited_nerve_candidate",
        "nonclaims": [
            "The generated complex is not proved to be a PL sphere.",
            "The archived Gale oriented matroid is not proved realizable over R.",
            "No convex 9-polytope realization is supplied.",
            "No Coxeter labels or Lorentzian Gram realization is supplied.",
            "No compactness or finite-volume conclusion is made.",
        ],
        "source": {
            "support_path": "data/enumerations/h9_14_nerve_support.json",
            "support_sha256": support_sha,
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
        },
        "proof_sources": proof_source_records(),
        "schedule_table": schedules,
        "incidence": incidence.to_record(),
        "structural_reduction": structure.to_record(),
        "global_selection_audit": global_audit.to_record(),
        "gale_search": search,
        "gale_chirotope_signs": list(normalized_signs),
        "gale_sign_sha256": sha256(
            bytes(1 if sign > 0 else 0 for sign in normalized_signs)
        ).hexdigest(),
        "gale_audit": gale_audit,
        "vinberg_scope": {
            "compact_combinatorial_analogue_checked": (
                "every selected abstract ridge lies in exactly two facets"
            ),
            "finite_volume_criterion_checked": False,
            "reason": (
                "Vinberg's finite-volume criterion requires an actual Coxeter "
                "Gram diagram and its elliptic/parabolic subdiagrams; neither "
                "labels nor a Gram matrix are present here."
            ),
        },
    }
    manifest = json.loads(json.dumps(manifest, sort_keys=True))
    manifest["certificate_sha256"] = certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_record18_abstract_audited_nerve_candidate":
        raise ValueError("wrong record-18 certificate kind")
    if manifest.get("status") != "abstract_audited_nerve_candidate":
        raise ValueError("record-18 certificate overclaims its status")
    support_hash = file_sha256(SUPPORT)
    if EXPECTED_SUPPORT_SHA256 != "TO_BE_FROZEN" and support_hash != EXPECTED_SUPPORT_SHA256:
        raise ValueError("record-18 support source changed")
    for path in PROOF_SOURCES:
        relative = str(path.relative_to(ROOT))
        observed = file_sha256(path)
        expected = EXPECTED_PROOF_SOURCE_SHA256[relative]
        if expected != "TO_BE_FROZEN" and observed != expected:
            raise ValueError(f"record-18 frozen proof source changed: {relative}")
    signs = manifest.get("gale_chirotope_signs")
    if not isinstance(signs, list):
        raise ValueError("record-18 certificate lacks Gale signs")
    provenance = manifest.get("gale_search")
    if not isinstance(provenance, Mapping):
        raise ValueError("record-18 certificate lacks search provenance")
    rebuilt = build(signs, search_provenance=provenance)
    if manifest != rebuilt:
        raise ValueError("record-18 certificate did not reproduce")
    digest = certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-18 critical digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-18 certificate digest is not frozen")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "status": "abstract_audited_nerve_candidate",
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
    parser.add_argument(
        "--refresh-archived-signs",
        action="store_true",
        help="rebuild all non-solver fields while retaining the archived signs",
    )
    args = parser.parse_args(argv)
    if args.write_certificate and args.refresh_archived_signs:
        parser.error("choose at most one certificate-writing mode")
    if args.write_certificate:
        manifest = build()
        OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    elif args.refresh_archived_signs:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        manifest = build(
            archived["gale_chirotope_signs"],
            search_provenance=archived["gale_search"],
        )
        OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    else:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps({"output": str(OUTPUT), **verify(manifest)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
