#!/usr/bin/env python3
"""Freeze and verify the ordinary-polytopal completion of H9/14 record 32."""

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
from coxeter_polytopes.combinatorics.h9_14_record32 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    TARGET_RECORD_INDEX,
    audit_record32_incidence,
    audit_record32_multiwedge,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record32_polytopality.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "94a65380b2d56d4787a3de28d5cca5dede0540cc2ec0ff23b2a3d65b4cfbd318"
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
            "convex_multiwedge_audit",
            "gale_chirotope_signs",
            "gale_audit",
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


def _gale_audit_record(
    facets: tuple[int, ...], signs: Sequence[int]
) -> dict[str, object]:
    audit = audit_uniform_rank_four_gale(facets, signs)
    if not audit.passes:
        raise AssertionError("archived record-32 Gale signs failed independent audit")
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
            "Independent oriented-matroid cross-check. Ordinary convex "
            "polytopality follows separately from the multiwedge construction."
        ),
    }


def build(signs: Sequence[int] | None = None) -> dict[str, object]:
    source, support_sha = _support_record()
    incidence = audit_record32_incidence(TARGET_GRAPH6)
    facets = incidence.facet_complement_masks
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        facets,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("record-32 completion failed the full global audit")
    if global_audit.face_vector != TARGET_F_VECTOR:
        raise AssertionError("record-32 face vector changed")
    expected_nonfaces = {
        2: tuple(sorted(sum(1 << vertex for vertex in face) for face in TARGET_MISSING_EDGES)),
        4: tuple(sorted(sum(1 << vertex for vertex in face) for face in TARGET_MISSING_FOURS)),
        5: tuple(sorted(sum(1 << vertex for vertex in face) for face in TARGET_MISSING_FIVES)),
    }
    if global_audit.minimal_nonfaces_by_size != expected_nonfaces:
        raise AssertionError("generated face table has unexpected minimal nonfaces")
    multiwedge = audit_record32_multiwedge()
    if not multiwedge.verified:
        raise AssertionError("record-32 convex multiwedge audit failed")

    if signs is None:
        result = solve_uniform_rank_four_gale(facets, timeout_ms=300_000)
        if result.status != "sat" or len(result.signs) != EXPECTED_CHIROTOPE_SIGN_COUNT:
            raise RuntimeError(
                f"record-32 Gale search did not finish SAT: {result.status} "
                f"({result.reason_unknown})"
            )
        normalized_signs = result.signs
        search: dict[str, object] = {
            "backend": result.backend,
            "backend_version": result.backend_version,
            "status": result.status,
            "elapsed_seconds": result.elapsed_seconds,
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
    gale_audit = _gale_audit_record(facets, normalized_signs)

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record32_ordinary_polytopal_completion",
        "claim": (
            "H9/14 survivor record 32 has an explicit 110-facet completion "
            "which is the boundary of a convex simplicial 9-polytope: two "
            "simplicial wedges over a pentagonal bipyramid, joined with the "
            "boundary of a four-simplex."
        ),
        "nonclaims": [
            "The construction supplies no Coxeter labeling.",
            "No Lorentzian Gram realization is supplied.",
            "The archived Gale signs are a cross-check, not the proof of real realizability.",
        ],
        "source": {
            "support_path": "data/enumerations/h9_14_nerve_support.json",
            "support_sha256": support_sha,
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
            "search_provenance": (
                "A finite forced-shadow/blocker census found the displayed "
                "a=2 clutter; all claims replay from the explicit clutter."
            ),
        },
        "incidence": incidence.to_record(),
        "global_selection_audit": global_audit.to_record(),
        "convex_multiwedge_audit": multiwedge.to_record(),
        "gale_search": search,
        "gale_chirotope_signs": list(normalized_signs),
        "gale_sign_sha256": sha256(
            bytes(1 if sign > 0 else 0 for sign in normalized_signs)
        ).hexdigest(),
        "gale_audit": gale_audit,
    }
    manifest["certificate_sha256"] = _critical_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_record32_ordinary_polytopal_completion":
        raise ValueError("wrong record-32 certificate kind")
    signs = manifest.get("gale_chirotope_signs")
    if not isinstance(signs, list):
        raise ValueError("record-32 certificate lacks Gale signs")
    rebuilt = build(signs)
    for key in (
        "claim",
        "nonclaims",
        "source",
        "incidence",
        "global_selection_audit",
        "convex_multiwedge_audit",
        "gale_chirotope_signs",
        "gale_sign_sha256",
        "gale_audit",
        "certificate_sha256",
    ):
        if manifest.get(key) != rebuilt.get(key):
            raise ValueError(f"record-32 certificate field does not verify: {key}")
    digest = _critical_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-32 critical digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-32 certificate digest is not frozen")
    return {
        "facet_count": len(rebuilt["incidence"]["facet_complement_masks"]),
        "face_vector": rebuilt["global_selection_audit"]["face_vector"],
        "ordinary_polytopal": rebuilt["convex_multiwedge_audit"]["ordinary_polytopal"],
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
        OUTPUT.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    print(json.dumps({"output": str(OUTPUT), **summary}, indent=2))


if __name__ == "__main__":
    main()
