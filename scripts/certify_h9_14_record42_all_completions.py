#!/usr/bin/env python3
"""Freeze the record-42 whole-shadow compact-Coxeter no-go theorem."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record42_all_completions import (
    PRODUCT_MISSING_FIVE,
    PRODUCT_MISSING_FOUR,
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record42_all_completions,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest
from scripts.independent_h9_14_record42_all_completions_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay as independent_replay,
)


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
SELECTED_LABEL = ROOT / "data/enumerations/h9_14_selected_nerve_coxeter_no_go.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record42_all_completions.json"

PROOF_PATHS = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record42_all_completions.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record46.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record92.py",
    ROOT / "src/coxeter_polytopes/combinatorics/facet_completion.py",
    ROOT / "src/coxeter_polytopes/combinatorics/simplicial_homology.py",
    ROOT / "src/coxeter_polytopes/combinatorics/local_link.py",
    ROOT / "src/coxeter_polytopes/combinatorics/missing_edges.py",
    ROOT / "src/coxeter_polytopes/combinatorics/nerve_support.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_selected_nerve_coxeter_no_go.py",
    ROOT / "scripts/certify_h9_14_selected_nerve_coxeter_no_go.py",
    ROOT / "scripts/independent_h9_14_selected_nerve_coxeter_no_go_audit.py",
    ROOT / "scripts/independent_h9_14_record42_all_completions_audit.py",
)

EXPECTED_SUPPORT_FILE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_SELECTED_LABEL_FILE_SHA256 = (
    "028bbc8e266145e59fa1ee9f9796a6e2c241f4a4f09637bd73f98bced9908629"
)
EXPECTED_SELECTED_LABEL_CERTIFICATE_SHA256 = (
    "c2dfcbfc03dc5dbae34938ab32f3a1c1111fd729d7844fd76d04dfc6cfa44de8"
)
EXPECTED_RECORD42_SELECTED_PROJECTION_SHA256 = (
    "9238c596b14426d120c25d8041c770b6dcce1dbf2526c2e2ea43c04894f66598"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "4590a0f0365dd26abd416cf4f73b909028a2863913f7c6d179a2a3e248122e9d"
)

EXPECTED_SELECTED_LABEL_KIND = (
    "h9_14_seven_selected_nerve_compact_coxeter_no_go_reaudited"
)
EXPECTED_RECORD42_SUPPORT_SYMMETRY_ORDERS = [24, 120]
EXPECTED_SCOPE_REAUDIT_STATUS = (
    "hostile_false_alarm_rejected_by_principal_block_proof"
)


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical_digest(document: Mapping[str, object]) -> str:
    # Normalize integer dictionary keys through the same JSON representation
    # that is written to disk.  Without this round trip, ``sort_keys=True``
    # orders integer histogram keys numerically before writing but string keys
    # lexicographically after loading, so an otherwise unchanged archive would
    # fail certificate closure.
    payload = json.loads(json.dumps(document, ensure_ascii=True))
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode()
    ).hexdigest()


def _projection_digest(value: object) -> str:
    return sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode()
    ).hexdigest()


def _path_record(path: Path) -> dict[str, str]:
    return {"path": str(path.relative_to(ROOT)), "sha256": _file_sha256(path)}


def _require_frozen_digest(expected: str, observed: str, message: str) -> None:
    if expected not in ("TO_BE_FROZEN", observed):
        raise ValueError(message)


def _support_record() -> tuple[dict[str, object], str]:
    raw = SUPPORT.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != EXPECTED_SUPPORT_FILE_SHA256:
        raise ValueError("record-42 support archive changed")
    document = json.loads(raw)
    verify_support_manifest(document)
    survivors = [row for row in document["records"] if row["survives_second_layer"]]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("record-42 survivor index changed")
    return record, digest


def _selected_product_label_chain() -> dict[str, object]:
    raw = SELECTED_LABEL.read_bytes()
    file_digest = sha256(raw).hexdigest()
    _require_frozen_digest(
        EXPECTED_SELECTED_LABEL_FILE_SHA256,
        file_digest,
        "selected-nerve label archive changed",
    )
    document = json.loads(raw)
    selected_certificate = document.get("certificate_sha256")
    if (
        not isinstance(selected_certificate, str)
        or _canonical_digest(document) != selected_certificate
    ):
        raise ValueError("selected-nerve label certificate closure failed")
    _require_frozen_digest(
        EXPECTED_SELECTED_LABEL_CERTIFICATE_SHA256,
        selected_certificate,
        "selected-nerve label certificate changed",
    )
    scope_reaudit = document.get("scope_reaudit")
    if (
        document.get("kind") != EXPECTED_SELECTED_LABEL_KIND
        or not isinstance(scope_reaudit, Mapping)
        or scope_reaudit.get("status") != EXPECTED_SCOPE_REAUDIT_STATUS
        or TARGET_RECORD_INDEX not in scope_reaudit.get("reaudited_indices", [])
        or scope_reaudit.get("normalizing_support_symmetry_orders", {}).get(
            str(TARGET_RECORD_INDEX)
        )
        != EXPECTED_RECORD42_SUPPORT_SYMMETRY_ORDERS
        or "principal submatrix" not in scope_reaudit.get(
            "principal_block_lemma", ""
        )
        or "inertia" not in scope_reaudit.get("principal_block_lemma", "")
    ):
        raise ValueError("selected-nerve principal-block scope re-audit changed")
    record = next(
        row
        for row in document["six_nerve_census"]
        if row["record_index"] == TARGET_RECORD_INDEX
    )
    projection_digest = _projection_digest(record)
    gate = record["support_gates"][0]
    if (
        gate["left_support"] != list(PRODUCT_MISSING_FOUR)
        or gate["right_support"] != list(PRODUCT_MISSING_FIVE)
        or gate["single_edge_test_count"] != 900
        or gate["single_edge_passing_count"] != 2
        or record.get("normalizing_support_symmetry_orders")
        != EXPECTED_RECORD42_SUPPORT_SYMMETRY_ORDERS
        or record["raw_two_edge_trial_count"] != 8_550
        or record["monotonicity_residual_two_edge_trial_count"] != 1
        or record["two_edge_passing_count"] != 0
        or record["stronger_single_edge_survivors"] != {"4": 0, "5": 0}
        or record["exact_index_two_witnesses"][0]["inertia"] != [7, 2, 0]
        or record["unrestricted_compact_coxeter_no_go"] is not True
    ):
        raise ValueError("selected record-42 unrestricted label theorem changed")
    _require_frozen_digest(
        EXPECTED_RECORD42_SELECTED_PROJECTION_SHA256,
        projection_digest,
        "selected record-42 projection changed",
    )
    return {
        "path": str(SELECTED_LABEL.relative_to(ROOT)),
        "file_sha256": file_digest,
        "certificate_sha256": document["certificate_sha256"],
        "record42_projection_sha256": projection_digest,
        "aggregate_kind": document["kind"],
        "scope_reaudit_status": scope_reaudit["status"],
        "principal_block_lemma": scope_reaudit["principal_block_lemma"],
        "product_missing_four": list(PRODUCT_MISSING_FOUR),
        "product_missing_five": list(PRODUCT_MISSING_FIVE),
        "normalizing_support_symmetry_orders": (
            record["normalizing_support_symmetry_orders"]
        ),
        "rank_nine_vertex_count": record["rank_nine_vertex_count"],
        "single_m3_cross_edge_trials": gate["single_edge_test_count"],
        "single_m3_cross_edge_survivors": gate["single_edge_passing_count"],
        "m4_and_m5_single_edge_survivors": record[
            "stronger_single_edge_survivors"
        ],
        "raw_two_edge_trials": record["raw_two_edge_trial_count"],
        "monotonicity_residual_two_edge_trials": record[
            "monotonicity_residual_two_edge_trial_count"
        ],
        "two_edge_survivors": record["two_edge_passing_count"],
        "exceptional_block_inertia": record["exact_index_two_witnesses"][0][
            "inertia"
        ],
        "unrestricted_compact_coxeter_no_go": True,
        "scope_upgrade": (
            "the selected theorem is selected-nerve-only; the new a=1 census "
            "proves that every surviving a=1 completion is in this exact S9 "
            "orbit, while the re-audited fixed principal block makes outside "
            "entries irrelevant by interlacing"
        ),
    }


def build() -> dict[str, object]:
    support_record, support_digest = _support_record()
    serialized = json.dumps(support_record, sort_keys=True, separators=(",", ":"))
    production = audit_record42_all_completions(serialized)
    if not production.verified:
        raise AssertionError("record-42 production census failed")
    product_chain = _selected_product_label_chain()
    hostile = independent_replay()
    if EXPECTED_SUMMARY_SHA256 not in (
        "TO_BE_FROZEN",
        hostile["summary_sha256"],
    ):
        raise AssertionError("record-42 hostile replay digest changed")
    if (
        hostile["theorem"][
            "whole_record42_graph_shadow_compact_coxeter_impossible"
        ]
        is not True
        or hostile["theorem"]["a0_exact_target_families_excluded"] != 15
        or hostile["a0_census"]["exact_target_family_count"] != 15
        or hostile["a0_census"]["ridge_valid_orbit_count"] != 1
        or hostile["positive_branches"]["a2_strict_count"] != 0
    ):
        raise AssertionError("record-42 hostile theorem chain is incomplete")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record42_whole_shadow_compact_coxeter_no_go",
        "claim": (
            "Every incidence completion of H9/14 graph shadow record 42 "
            "compatible with the compact Coxeter missing-face reductions, "
            "Dehn--Sommerville schedules, strict local inequalities, and the "
            "facet/blocker conditions fails the unrestricted compact Coxeter "
            "label tests. Hence the whole record-42 shadow is impossible in "
            "the compact case."
        ),
        "scope": {
            "dimension": 9,
            "facet_count": 14,
            "record_index": TARGET_RECORD_INDEX,
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "all_five_dehn_sommerville_schedules_covered": True,
            "finite_label_alphabet": [2, 3, 4, 5],
            "compact_only": True,
            "finite_volume_noncompact_chamber_excluded": False,
        },
        "nonclaims": [
            "No ideal-vertex or parabolic cusp census is performed.",
            "No finite-volume noncompact chamber is excluded.",
            "No graph shadow other than record 42 is addressed here.",
            "No arbitrary ordinary completion with a minimal nonface above the compact Coxeter Lanner bound is classified.",
        ],
        "proof_dependencies": [
            "Burcroff's no-missing-triples theorem for compact H9/14 nerves",
            "Burcroff Lemma 10.10: finite meeting labels are exactly 2,3,4,5 in this compact target",
            "every compact Coxeter missing face is Lanner of order at most five",
            "Dehn--Sommerville equations for the simplicial polar boundary",
            "strict generalized-lower-bound and Macaulay inequalities in every vertex link",
            "facet complements are the blocker of the full minimal-nonface clutter",
            "Perron--Frobenius monotonicity for unassigned finite non-right labels",
            "full S4 x S5 support normalization of all 108 x 420 labelled a=1 Lanner placements",
            "exact negative-index-two obstruction for the two exceptional a=1 cross edges",
            "Cauchy interlacing for the fixed a=1 nine-node principal Gram block",
        ],
        "support_source": {
            "path": str(SUPPORT.relative_to(ROOT)),
            "file_sha256": support_digest,
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
        },
        "proof_sources": [_path_record(path) for path in PROOF_PATHS],
        "whole_shadow_production_audit": production.to_record(),
        "a1_unrestricted_product_label_chain": product_chain,
        "independent_no_project_import_replay": {
            "path": "scripts/independent_h9_14_record42_all_completions_audit.py",
            "source_sha256": _file_sha256(
                ROOT / "scripts/independent_h9_14_record42_all_completions_audit.py"
            ),
            "summary_sha256": hostile["summary_sha256"],
            "a0_label_stage_transcript_sha256": hostile[
                "a0_label_stage_transcript_sha256"
            ],
            "a0_exact_target_families_excluded": hostile["theorem"][
                "a0_exact_target_families_excluded"
            ],
            "whole_shadow_compact_coxeter_excluded": True,
        },
        "theorem": {
            "record_index": TARGET_RECORD_INDEX,
            "a0_exact_target_families_excluded_by_complete_label_gate": 15,
            "a1_incidence_orbits_excluded_by_reaudited_principal_block_gate": 1,
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    manifest["certificate_sha256"] = _canonical_digest(manifest)
    return manifest


def verify(
    document: Mapping[str, object], *, recompute: bool = False
) -> dict[str, object]:
    if document.get("kind") != "h9_14_record42_whole_shadow_compact_coxeter_no_go":
        raise ValueError("unexpected record-42 certificate kind")
    if _file_sha256(SUPPORT) != EXPECTED_SUPPORT_FILE_SHA256:
        raise ValueError("record-42 support source drifted")
    _require_frozen_digest(
        EXPECTED_SELECTED_LABEL_FILE_SHA256,
        _file_sha256(SELECTED_LABEL),
        "record-42 selected label dependency drifted",
    )
    if document.get("a1_unrestricted_product_label_chain") != (
        _selected_product_label_chain()
    ):
        raise ValueError("record-42 principal-block dependency changed")
    if document.get("proof_sources") != [_path_record(path) for path in PROOF_PATHS]:
        raise ValueError("record-42 proof source provenance changed")
    digest = _canonical_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("record-42 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 not in ("TO_BE_FROZEN", digest):
        raise ValueError("record-42 certificate is not frozen")
    theorem = document.get("theorem")
    scope = document.get("scope")
    if (
        not isinstance(theorem, Mapping)
        or theorem.get("whole_graph_shadow_compact_coxeter_excluded") is not True
        or theorem.get("a0_exact_target_families_excluded_by_complete_label_gate") != 15
        or not isinstance(scope, Mapping)
        or scope.get("compact_only") is not True
        or scope.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("record-42 theorem scope changed")
    if recompute and json.loads(json.dumps(build(), sort_keys=True)) != json.loads(
        json.dumps(document, sort_keys=True)
    ):
        raise ValueError("record-42 full theorem did not replay")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "a0_exact_target_families_excluded": 15,
        "a1_incidence_orbits_excluded": 1,
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args(argv)
    document = (
        json.loads(OUTPUT.read_text(encoding="utf-8"))
        if args.verify_only
        else build()
    )
    if args.write:
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps(verify(document, recompute=args.recompute), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
