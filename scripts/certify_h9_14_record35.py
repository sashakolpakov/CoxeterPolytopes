#!/usr/bin/env python3
"""Freeze and verify the ordinary-polytopal H9/14 record-35 completion."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_record35 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_RECORD_INDEX,
    audit_record35_base_realization,
    audit_record35_incidence,
    audit_record35_multiwedge,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
BASELINE = ROOT / "data/enumerations/h9_14_facet_completion.json"
PROOF_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record35.py"
)
OUTPUT = ROOT / "data/enumerations/h9_14_record35_polytopality.json"
EXPECTED_SUPPORT_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "bc535f94895c6e151c5bd05db142a7735e1260bc736012daa6ccbe5787d00db7"
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
    return _canonical_sha256(critical)


def _support_record() -> tuple[dict[str, object], str]:
    support_sha = _file_sha256(SUPPORT)
    if support_sha != EXPECTED_SUPPORT_SHA256:
        raise ValueError("record-35 support source is not frozen")
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    verify_support_manifest(support)
    survivors = [
        record for record in support["records"] if record["survives_second_layer"]
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("second-layer survivor record 35 changed")
    return record, support_sha


def _backend_corroboration() -> dict[str, object]:
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    records = [
        record
        for record in baseline["records"]
        if record.get("index") == TARGET_RECORD_INDEX
    ]
    if len(records) != 1 or records[0]["solver"]["status"] != "unknown":
        raise ValueError("record-35 baseline timeout changed")
    solver = records[0]["solver"]
    return {
        "path": "data/enumerations/h9_14_facet_completion.json",
        "sha256": _file_sha256(BASELINE),
        "status": solver["status"],
        "reason_unknown": solver["reason_unknown"],
        "backend": solver["backend"],
        "backend_version": solver["backend_version"],
        "elapsed_seconds": solver["elapsed_seconds"],
        "trusted_for_proof": False,
    }


def build() -> dict[str, object]:
    source, support_sha = _support_record()
    incidence = audit_record35_incidence()
    facets = incidence.facet_complement_masks
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        facets,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("record-35 completion failed its full link audit")
    if global_audit.face_vector != TARGET_F_VECTOR:
        raise AssertionError("record-35 face vector changed")
    expected_nonfaces = {
        2: tuple(
            sorted(sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES)
        ),
        5: tuple(
            sorted(sum(1 << vertex for vertex in face) for face in TARGET_MISSING_FIVES)
        ),
    }
    if global_audit.minimal_nonfaces_by_size != expected_nonfaces:
        raise AssertionError("record-35 face table has unexpected minimal nonfaces")

    base = audit_record35_base_realization()
    multiwedge = audit_record35_multiwedge()
    if not base.verified or not multiwedge.verified:
        raise AssertionError("record-35 rational multiwedge proof failed")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record35_rational_multiwedge_polytopality",
        "claim": (
            "H9/14 survivor record 35 is the nerve of an ordinary rational "
            "convex simple 9-polytope with fourteen facets."
        ),
        "theorem_status": "certified_explicit_ordinary_polytopal",
        "proof_backend": (
            "complete 2^14 blocker replay, exact integer determinants, and "
            "an exact inverse-simplicial-multiwedge audit; no SAT/SMT backend"
        ),
        "source": {
            "support_path": "data/enumerations/h9_14_nerve_support.json",
            "support_sha256": support_sha,
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
        },
        "proof_source": {
            "path": "src/coxeter_polytopes/combinatorics/h9_14_record35.py",
            "sha256": _file_sha256(PROOF_SOURCE),
        },
        "incidence": incidence.to_record(),
        "global_selection_audit": global_audit.to_record(),
        "integral_base_realization": base.to_record(),
        "inverse_multiwedge": multiwedge.to_record(),
        "construction": (
            "Collapse clone blocks {3,6,8} and {5,7,9,10}. The quotient is "
            "the boundary of the explicitly realized integral simplicial "
            "4-polytope. Replace its two carrier vertices by three and four "
            "clones. The five simplicial wedges give a simplicial "
            "9-polytope boundary with exactly the target minimal nonfaces; "
            "duality gives the desired simple 9-polytope."
        ),
        "proof_lemma": {
            "statement": (
                "A simplicial wedge of a simplicial convex-polytope boundary "
                "is again a simplicial convex-polytope boundary. Iterating "
                "realizes every positive integral multiwedge. Dually this is "
                "the facet-wedge operation on simple convex polytopes."
            ),
            "reference": (
                "S. Choi and H. Park, Wedge Operations and Torus Symmetries, "
                "Tohoku Math. J. 68 (2016), Proposition 2.2(3) and Lemma 2.3"
            ),
            "url": (
                "https://www.jstage.jst.go.jp/article/tmj/68/1/"
                "68_1458248864/_pdf"
            ),
        },
        "algebro_geometric_corollary": {
            "claim": (
                "There is a projective Q-factorial toric 9-fold whose toric "
                "boundary-stratum dual complex is the record-35 nerve."
            ),
            "construction": (
                "Translate an interior rational point of the integral base "
                "4-polytope to the origin and polarize. Its polar is a "
                "rational simple 4-polytope. Five rational facet wedges give "
                "a rational simple 9-polytope; its rational simplicial normal "
                "fan defines the claimed projective Q-factorial toric variety."
            ),
            "caveat": (
                "No unimodularity or Delzant audit is supplied, so smoothness "
                "and an SNC toric boundary are not claimed."
            ),
        },
        "backend_corroboration": _backend_corroboration(),
        "updated_theorem_partition": {
            "independently_certified_impossible": 120,
            "certified_explicit_ordinary_polytopal": 6,
            "theorem_undecided": 16,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "No Coxeter labels are supplied.",
            "No Lorentzian Coxeter Gram matrix is supplied.",
            "Ordinary polytopality does not imply hyperbolic Coxeter realizability.",
            "The archived bounded Z3 timeout contributes no logical step.",
            "This does not produce a new compact hyperbolic Coxeter polytope.",
        ],
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != (
        "h9_14_record35_rational_multiwedge_polytopality"
    ):
        raise ValueError("wrong record-35 certificate kind")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-35 certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-35 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-35 certificate is not frozen")
    return {
        "facet_count": len(rebuilt["incidence"]["facet_complement_masks"]),
        "face_vector": rebuilt["global_selection_audit"]["face_vector"],
        "base_facet_count": len(
            rebuilt["integral_base_realization"]["support_facets"]
        ),
        "minimum_facet_side_absolute_determinant": rebuilt[
            "integral_base_realization"
        ]["minimum_facet_side_absolute_determinant"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    del argv
    manifest = build()
    OUTPUT.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"output": str(OUTPUT), **verify(manifest)}, indent=2))


if __name__ == "__main__":
    main()
