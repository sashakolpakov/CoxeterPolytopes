#!/usr/bin/env python3
"""Freeze the record-43 whole-shadow compact-Coxeter no-go theorem."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record43 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_RECORD_INDEX,
)
from coxeter_polytopes.combinatorics.h9_14_record43_all_completions import (
    SECOND_TARGET_MISSING_FIVES,
    audit_record43_all_completions,
)
from coxeter_polytopes.combinatorics.h9_14_record43_second_nerve_coxeter_no_go import (
    audit_record43_second_nerve_labels,
)
from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.nerve_support import (
    verify_support_manifest,
)
from scripts.independent_h9_14_record43_all_completions_audit import (
    EXPECTED_SUMMARY_SHA256 as EXPECTED_HOSTILE_SUMMARY_SHA256,
    replay as independent_replay,
)


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
SELECTED_COMPLETION = ROOT / "data/enumerations/h9_14_record43_selection.json"
SELECTED_LABEL_NO_GO = (
    ROOT / "data/enumerations/h9_14_selected_nerve_coxeter_no_go.json"
)
PROOF_SOURCE = (
    ROOT
    / "src/coxeter_polytopes/combinatorics/h9_14_record43_all_completions.py"
)
SECOND_LABEL_SOURCE = (
    ROOT
    / "src/coxeter_polytopes/combinatorics/"
    "h9_14_record43_second_nerve_coxeter_no_go.py"
)
RECORD43_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record43.py"
)
FACET_AUDIT_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/facet_completion.py"
)
HOMOLOGY_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/simplicial_homology.py"
)
LOCAL_LINK_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/local_link.py"
)
MISSING_EDGES_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/missing_edges.py"
)
RECORD92_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record92.py"
)
NERVE_SUPPORT_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/nerve_support.py"
)
SELECTED_LABEL_SOURCE = (
    ROOT
    / "src/coxeter_polytopes/combinatorics/"
    "h9_14_selected_nerve_coxeter_no_go.py"
)
SELECTED_LABEL_CERTIFIER = (
    ROOT / "scripts/certify_h9_14_selected_nerve_coxeter_no_go.py"
)
SELECTED_LABEL_HOSTILE = (
    ROOT / "scripts/independent_h9_14_selected_nerve_coxeter_no_go_audit.py"
)
INDEPENDENT_SOURCE = (
    ROOT / "scripts/independent_h9_14_record43_all_completions_audit.py"
)
OUTPUT = (
    ROOT / "data/enumerations/h9_14_record43_all_completions.json"
)

EXPECTED_SUPPORT_FILE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_SELECTED_COMPLETION_FILE_SHA256 = (
    "2c6463598ab4ae2073b54fc8a9345ab099acba0492c786ba163d0ad5d29d3c66"
)
EXPECTED_SELECTED_COMPLETION_CERTIFICATE_SHA256 = (
    "af521ed4e7476c95945df9b770f36a343dcf3fc79dabed557aaab0fd97a671a3"
)
EXPECTED_SELECTED_LABEL_FILE_SHA256 = (
    "028bbc8e266145e59fa1ee9f9796a6e2c241f4a4f09637bd73f98bced9908629"
)
EXPECTED_SELECTED_LABEL_CERTIFICATE_SHA256 = (
    "c2dfcbfc03dc5dbae34938ab32f3a1c1111fd729d7844fd76d04dfc6cfa44de8"
)
EXPECTED_SELECTED_LABEL_KIND = (
    "h9_14_seven_selected_nerve_compact_coxeter_no_go_reaudited"
)
EXPECTED_RECORD43_PROJECTION_SHA256 = (
    "6437c7a8881def45ba1b5abb620b1c0574f1e831f34806b303336a6dd6a6f66b"
)
EXPECTED_RECORD43_SUPPORT_SYMMETRY_ORDERS = [120]
EXPECTED_SCOPE_REAUDIT_STATUS = (
    "hostile_false_alarm_rejected_by_principal_block_proof"
)
EXPECTED_RECORD43_SOURCE_SHA256 = (
    "49623375ad2f9a3aa19865e740719daceb93d399db1bd9e38ec21b1f9e4d2ed5"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "ca83707489c65f23711849417a28bf548cb5d456a80b9f075e0e642cec785530"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _projection_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _path_record(path: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": _file_sha256(path),
    }


def _support_record() -> tuple[dict[str, object], str]:
    raw = SUPPORT.read_bytes()
    file_digest = sha256(raw).hexdigest()
    if file_digest != EXPECTED_SUPPORT_FILE_SHA256:
        raise ValueError("record-43 support archive file hash changed")
    source = json.loads(raw)
    verify_support_manifest(source)
    survivors = [
        row for row in source["records"] if row["survives_second_layer"]
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("second-layer record-43 graph changed")
    return record, file_digest


def _selected_completion_chain(
    selected_missing_fives: Sequence[Sequence[int]],
) -> dict[str, object]:
    raw = SELECTED_COMPLETION.read_bytes()
    file_digest = sha256(raw).hexdigest()
    if file_digest != EXPECTED_SELECTED_COMPLETION_FILE_SHA256:
        raise ValueError("selected record-43 completion archive changed")
    document = json.loads(raw)
    if document.get("certificate_sha256") != (
        EXPECTED_SELECTED_COMPLETION_CERTIFICATE_SHA256
    ):
        raise ValueError("selected record-43 certificate digest changed")
    incidence = document.get("incidence")
    global_audit = document.get("global_selection_audit")
    if not isinstance(incidence, Mapping) or not isinstance(global_audit, Mapping):
        raise ValueError("selected record-43 certificate is malformed")
    if incidence.get("missing_fives") != [list(row) for row in selected_missing_fives]:
        raise ValueError("first completion does not equal the selected nerve")
    if incidence.get("missing_edges") != [list(row) for row in TARGET_MISSING_EDGES]:
        raise ValueError("selected record-43 missing edges changed")
    required_global = {
        "facet_count": 90,
        "face_vector": list(TARGET_F_VECTOR),
        "all_ridge_degrees_two": True,
        "dual_graph_connected": True,
        "all_links_have_spherical_euler_characteristic": True,
        "all_links_are_mod2_homology_spheres": True,
        "passes_selection_model": True,
        "passes_connected_homology_audit": True,
    }
    for key, expected in required_global.items():
        if global_audit.get(key) != expected:
            raise ValueError(f"selected record-43 global audit changed: {key}")
    projection = {
        "certificate_sha256": document["certificate_sha256"],
        "incidence": incidence,
        "global_selection_audit": global_audit,
    }
    return {
        "path": str(SELECTED_COMPLETION.relative_to(ROOT)),
        "file_sha256": file_digest,
        "certificate_sha256": document["certificate_sha256"],
        "projection_sha256": _projection_digest(projection),
        "facet_complement_sha256": incidence["facet_complement_sha256"],
        "packed_face_table_sha256": incidence["packed_face_table_sha256"],
        "facet_count": global_audit["facet_count"],
        "face_vector": global_audit["face_vector"],
        "passes_connected_homology_audit": global_audit[
            "passes_connected_homology_audit"
        ],
        "first_census_row_equals_selected_nerve": True,
    }


def _selected_label_chain() -> dict[str, object]:
    raw = SELECTED_LABEL_NO_GO.read_bytes()
    file_digest = sha256(raw).hexdigest()
    if file_digest != EXPECTED_SELECTED_LABEL_FILE_SHA256:
        raise ValueError("selected-nerve compact no-go archive file hash changed")
    document = json.loads(raw)
    if (
        document.get("certificate_sha256")
        != EXPECTED_SELECTED_LABEL_CERTIFICATE_SHA256
        or _canonical_digest(document)
        != EXPECTED_SELECTED_LABEL_CERTIFICATE_SHA256
    ):
        raise ValueError("selected-nerve compact no-go archive closure failed")
    scope_reaudit = document.get("scope_reaudit")
    if (
        document.get("kind") != EXPECTED_SELECTED_LABEL_KIND
        or not isinstance(scope_reaudit, Mapping)
        or scope_reaudit.get("status") != EXPECTED_SCOPE_REAUDIT_STATUS
        or TARGET_RECORD_INDEX not in scope_reaudit.get("reaudited_indices", [])
        or scope_reaudit.get("normalizing_support_symmetry_orders", {}).get(
            str(TARGET_RECORD_INDEX)
        )
        != EXPECTED_RECORD43_SUPPORT_SYMMETRY_ORDERS
        or "principal submatrix" not in scope_reaudit.get(
            "principal_block_lemma", ""
        )
        or "inertia" not in scope_reaudit.get("principal_block_lemma", "")
    ):
        raise ValueError("selected-nerve principal-block scope re-audit changed")
    rows = document.get("six_nerve_census")
    if not isinstance(rows, list):
        raise ValueError("selected-nerve compact no-go rows are missing")
    record = next(
        (row for row in rows if row.get("record_index") == TARGET_RECORD_INDEX),
        None,
    )
    if not isinstance(record, Mapping):
        raise ValueError("selected-nerve compact no-go lacks record 43")
    if (
        record.get("rank_nine_vertex_count") != 90
        or record.get("selected_nerve_only") is not True
        or record.get("unrestricted_compact_coxeter_no_go") is not True
        or record.get("normalizing_support_symmetry_orders")
        != EXPECTED_RECORD43_SUPPORT_SYMMETRY_ORDERS
        or record.get("two_edge_passing_count") != 0
        or record.get("raw_two_edge_trial_count") != 1_890_000
        or record.get("monotonicity_residual_two_edge_trial_count") != 0
        or record.get("stronger_single_edge_survivors") != {"4": 0, "5": 0}
    ):
        raise ValueError("record-43 unrestricted label obstruction changed")
    support_gates = record.get("support_gates")
    if not isinstance(support_gates, list) or len(support_gates) != 3:
        raise ValueError("record-43 support gates changed")
    if any(
        gate.get("left_support") != list(TARGET_MISSING_FIVES[0])
        for gate in support_gates
    ) or {
        tuple(gate.get("right_support", ())) for gate in support_gates
    } != set(TARGET_MISSING_FIVES[1:]):
        raise ValueError("record-43 label supports do not match the selected nerve")
    sources = document.get("sources")
    if not isinstance(sources, Mapping) or not isinstance(sources.get("code"), list):
        raise ValueError("selected-nerve source closure is malformed")
    record43_source = {
        "path": "src/coxeter_polytopes/combinatorics/h9_14_record43.py",
        "sha256": EXPECTED_RECORD43_SOURCE_SHA256,
    }
    if record43_source not in sources["code"]:
        raise ValueError("selected-nerve archive lost record-43 source provenance")
    if _file_sha256(RECORD43_SOURCE) != EXPECTED_RECORD43_SOURCE_SHA256:
        raise ValueError("record-43 source file changed")
    witnesses = record.get("exact_index_two_witnesses")
    if not isinstance(witnesses, list) or [row.get("inertia") for row in witnesses] != [
        [8, 2, 0],
        [8, 2, 0],
    ]:
        raise ValueError("record-43 exceptional-block inertia changed")
    projection_digest = _projection_digest(record)
    if projection_digest != EXPECTED_RECORD43_PROJECTION_SHA256:
        raise ValueError("record-43 selected-nerve projection changed")
    return {
        "path": str(SELECTED_LABEL_NO_GO.relative_to(ROOT)),
        "file_sha256": file_digest,
        "certificate_sha256": document["certificate_sha256"],
        "record43_projection_sha256": projection_digest,
        "aggregate_kind": document["kind"],
        "scope_reaudit_status": scope_reaudit["status"],
        "principal_block_lemma": scope_reaudit["principal_block_lemma"],
        "normalizing_support_symmetry_orders": record[
            "normalizing_support_symmetry_orders"
        ],
        "rank_nine_vertex_count": record["rank_nine_vertex_count"],
        "selected_nerve_only_in_chained_certificate": record[
            "selected_nerve_only"
        ],
        "support_gate_count": len(support_gates),
        "single_edge_test_count_per_gate": [
            row["single_edge_test_count"] for row in support_gates
        ],
        "left_support": support_gates[0]["left_support"],
        "right_supports": [row["right_support"] for row in support_gates],
        "raw_two_edge_trial_count": record["raw_two_edge_trial_count"],
        "two_edge_passing_count": record["two_edge_passing_count"],
        "exceptional_block_inertias": [
            row["inertia"] for row in record["exact_index_two_witnesses"]
        ],
        "unrestricted_compact_coxeter_no_go": record[
            "unrestricted_compact_coxeter_no_go"
        ],
        "scope_upgrade": (
            "the re-audited certificate is selected-nerve-only; its fixed "
            "principal-block obstruction is unrestricted by interlacing, and "
            "the exhaustive "
            "census identifies it as the first of exactly two residual "
            "incidence orbits compatible with the compact missing-face, "
            "Dehn--Sommerville, and link necessities"
        ),
    }


def _minimal_transversals(hyperedges: Sequence[int]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, 1 << 14):
        if not all(candidate & hyperedge for hyperedge in hyperedges):
            continue
        if all(
            not all((candidate ^ (1 << vertex)) & hyperedge for hyperedge in hyperedges)
            for vertex in range(14)
            if candidate & (1 << vertex)
        ):
            answers.append(candidate)
    return tuple(answers)


def _second_completion_chain(
    support_record: Mapping[str, object],
) -> dict[str, object]:
    edge_masks = tuple(
        sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES
    )
    five_masks = tuple(
        sum(1 << vertex for vertex in face)
        for face in SECOND_TARGET_MISSING_FIVES
    )
    blocker = _minimal_transversals((*edge_masks, *five_masks))
    if len(blocker) != 90 or any(value.bit_count() != 5 for value in blocker):
        raise AssertionError("record-43 second target blocker changed")
    audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        support_record["ridge_support_core_masks"],  # type: ignore[arg-type]
        support_record["dehn_sommerville_schedule"],  # type: ignore[arg-type]
        blocker,
    )
    expected_nonfaces = {
        2: tuple(sorted(edge_masks)),
        5: tuple(sorted(five_masks)),
    }
    if (
        not audit.passes_connected_homology_audit
        or audit.face_vector != TARGET_F_VECTOR
        or audit.minimal_nonfaces_by_size != expected_nonfaces
    ):
        raise AssertionError("record-43 second target failed its full incidence audit")
    facet_digest = sha256(
        b"".join(value.to_bytes(2, "little") for value in blocker)
    ).hexdigest()
    if facet_digest != (
        "7afc91cca45c5d5796350082aa06d63e0df593734b2fb3fa200c615b6a20faa9"
    ):
        raise AssertionError("record-43 second target facet digest changed")
    return {
        "missing_fives": [list(face) for face in SECOND_TARGET_MISSING_FIVES],
        "facet_complement_count": len(blocker),
        "facet_complement_sha256": facet_digest,
        "global_selection_audit": audit.to_record(),
        "passes_connected_homology_audit": True,
    }


def _second_label_chain() -> dict[str, object]:
    audit = audit_record43_second_nerve_labels()
    if not audit.unrestricted_compact_coxeter_no_go:
        raise AssertionError("record-43 second target label obstruction changed")
    return {
        "source": _path_record(SECOND_LABEL_SOURCE),
        "audit": audit.to_record(),
        "unrestricted_compact_coxeter_no_go": True,
    }


def _proof_sources() -> list[dict[str, str]]:
    return [
        _path_record(path)
        for path in (
            PROOF_SOURCE,
            RECORD43_SOURCE,
            FACET_AUDIT_SOURCE,
            HOMOLOGY_SOURCE,
            LOCAL_LINK_SOURCE,
            MISSING_EDGES_SOURCE,
            RECORD92_SOURCE,
            NERVE_SUPPORT_SOURCE,
            SELECTED_LABEL_SOURCE,
            SELECTED_LABEL_CERTIFIER,
            SELECTED_LABEL_HOSTILE,
            SECOND_LABEL_SOURCE,
            INDEPENDENT_SOURCE,
        )
    ]


def build() -> dict[str, object]:
    source, support_digest = _support_record()
    schedules_json = json.dumps(
        source["dehn_sommerville_schedule"],
        sort_keys=True,
        separators=(",", ":"),
    )
    census = audit_record43_all_completions(schedules_json)
    if not census.verified or census.exact_residual_orbit_count != 2:
        raise AssertionError("record-43 whole-shadow incidence census failed")
    target_states = tuple(
        tuple(tuple(row) for row in state)
        for state in census.a0.to_record()["target_missing_five_states"]
    )
    if target_states != (TARGET_MISSING_FIVES, SECOND_TARGET_MISSING_FIVES):
        raise AssertionError("record-43 two completion states changed")
    first_completion_chain = _selected_completion_chain(target_states[0])
    second_completion_chain = _second_completion_chain(source)
    first_label_chain = _selected_label_chain()
    second_label_chain = _second_label_chain()
    hostile = independent_replay()
    if hostile.get("summary_sha256") != EXPECTED_HOSTILE_SUMMARY_SHA256:
        raise AssertionError("record-43 independent replay digest changed")
    hostile_theorem = hostile.get("theorem")
    hostile_a0 = hostile.get("a0_census")
    if (
        not isinstance(hostile_theorem, Mapping)
        or hostile_theorem.get(
            "whole_record43_graph_shadow_compact_coxeter_impossible"
        )
        is not True
        or hostile_theorem.get("incidence_completion_orbit_count") != 2
        or hostile_theorem.get("both_incidence_targets_compact_label_impossible")
        is not True
        or not isinstance(hostile_a0, Mapping)
        or hostile_a0.get("enumeration_path_count") != 1073
        or hostile_a0.get("exact_mask_clutter_count") != 887
        or hostile_a0.get("s4_s8_canonical_state_count") != 394
        or hostile_a0.get("scheduled_target_enumeration_path_count") != 2
        or hostile_a0.get("scheduled_target_exact_mask_count") != 2
        or hostile_a0.get("scheduled_target_s4_s8_canonical_count") != 2
        or hostile_a0.get("scheduled_target_full_graph_orbit_count") != 2
    ):
        raise AssertionError("record-43 independent replay did not close the theorem")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record43_whole_shadow_compact_coxeter_no_go",
        "claim": (
            "Under the compact H9/14 missing-face reductions, all four "
            "Dehn--Sommerville schedules for record 43 are exhaustive. The "
            "a=1,2,3 branches are impossible and the a=0 schedule has exactly "
            "two residual incidence orbits. The first is excluded by the "
            "frozen, principal-block-re-audited selected-nerve theorem and "
            "the second by a complete "
            "distinguished-support S4 x S4 finite-label and exact-inertia "
            "gate. Hence the whole graph shadow is not compact Coxeter."
        ),
        "scope": {
            "record_index": TARGET_RECORD_INDEX,
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "all_four_dehn_sommerville_schedules_covered": True,
            "exact_residual_incidence_orbits": 2,
            "compact_only": True,
            "finite_volume_noncompact_chamber_excluded": False,
            "global_lorentzian_gram_values_reenumerated_here": False,
        },
        "nonclaims": [
            "No finite-volume noncompact chamber is excluded.",
            "No ideal-vertex or maximal-parabolic census is performed.",
            "No other H9/14 graph shadow is excluded by this certificate.",
            "Arbitrary ordinary completions with minimal nonfaces above the compact Coxeter order bound are not classified.",
            "The first selected label obstruction is chained as a frozen theorem; the second label obstruction is recomputed here.",
        ],
        "proof_dependencies": [
            "Burcroff's no-missing-triples theorem for compact H9/14 nerves",
            "every compact Coxeter missing face is Lanner of order at most five",
            "Dehn--Sommerville equations for the simplicial polar boundary",
            "strict generalized-lower-bound inequalities in every vertex link",
            "Macaulay inequalities in every vertex link",
            "facet complements are the blocker of the full minimal-nonface clutter",
            "Perron--Frobenius monotonicity, the principal-submatrix interlacing lemma, full order-120 support normalization, and exact inertia in the frozen unrestricted selected-nerve label theorem",
            "the verified S4 x S4 automorphism subgroup and complete distinguished-support placement orbit census for the second nerve",
        ],
        "support_source": {
            "path": str(SUPPORT.relative_to(ROOT)),
            "file_sha256": support_digest,
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
        },
        "proof_sources": _proof_sources(),
        "whole_shadow_incidence_census": census.to_record(),
        "first_completion_chain": first_completion_chain,
        "second_completion_chain": second_completion_chain,
        "first_unrestricted_compact_label_chain": first_label_chain,
        "second_unrestricted_compact_label_chain": second_label_chain,
        "independent_hostile_replay": {
            "path": str(INDEPENDENT_SOURCE.relative_to(ROOT)),
            "source_sha256": _file_sha256(INDEPENDENT_SOURCE),
            "summary_sha256": hostile["summary_sha256"],
            "a0_enumeration_path_count": hostile_a0[
                "enumeration_path_count"
            ],
            "a0_exact_mask_clutter_count": hostile_a0[
                "exact_mask_clutter_count"
            ],
            "a0_s4_s8_canonical_state_count": hostile_a0[
                "s4_s8_canonical_state_count"
            ],
            "a0_s4_s8_canonical_facet_count_histogram": hostile_a0[
                "s4_s8_canonical_facet_count_histogram"
            ],
            "whole_shadow_compact_coxeter_excluded": hostile_theorem[
                "whole_record43_graph_shadow_compact_coxeter_impossible"
            ],
        },
        "theorem": {
            "dimension": 9,
            "facet_count": 14,
            "record_index": TARGET_RECORD_INDEX,
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    manifest["certificate_sha256"] = _canonical_digest(manifest)
    return manifest


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != (
        "h9_14_record43_whole_shadow_compact_coxeter_no_go"
    ):
        raise ValueError("wrong record-43 all-completion certificate kind")
    digest = _canonical_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("record-43 all-completion certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 not in ("TO_BE_FROZEN", digest):
        raise ValueError("record-43 all-completion certificate is not frozen")
    rebuilt = build()
    if dict(document) != rebuilt:
        raise ValueError("record-43 all-completion certificate does not replay")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "exact_residual_incidence_orbits": rebuilt["scope"][
            "exact_residual_incidence_orbits"
        ],
        "whole_graph_shadow_compact_coxeter_excluded": rebuilt["theorem"][
            "whole_graph_shadow_compact_coxeter_excluded"
        ],
        "finite_volume_noncompact_chamber_excluded": rebuilt["theorem"][
            "finite_volume_noncompact_chamber_excluded"
        ],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", action="store_true")
    args = parser.parse_args(argv)
    if args.write_certificate:
        document = build()
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps({"output": str(OUTPUT), **verify(document)}, indent=2))


if __name__ == "__main__":
    main()
