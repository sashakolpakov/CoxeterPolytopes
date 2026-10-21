#!/usr/bin/env python3
"""Build or verify the all-371 compact-Coxeter no-go for record 12."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record12_all_coxeter_labels import (
    audit_record12_all_coxeter_labels,
)
from scripts.independent_h9_14_record12_all_coxeter_audit import (
    compact_summary,
    replay as independent_replay,
)
from scripts.independent_h9_14_record12_coxeter_probe import (
    EXPECTED_ALL_SUMMARY_SHA256 as EXPECTED_SECONDARY_SUMMARY_SHA256,
    replay as secondary_independent_replay,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/enumerations/h9_14_record12_residual.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record12_all_coxeter_labels.json"
SOURCE_PATHS = (
    "src/coxeter_polytopes/combinatorics/h9_14_record12_all_coxeter_labels.py",
    "src/coxeter_polytopes/combinatorics/h9_14_record12_residual.py",
    "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    "scripts/independent_h9_14_record12_all_coxeter_audit.py",
    "scripts/independent_h9_14_record12_coxeter_probe.py",
)
INPUT_PATHS = ("data/enumerations/h9_14_record12_residual.json",)
EXPECTED_HOSTILE_SUMMARY_SHA256 = (
    "f1b14eaab1be6cc4ac7d1eb49cd79ed1401d6f14cd9ab8ee81efd18b8b863eb3"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "8e06523196d30b344a0b7cec5771d968d75ab253e5c4a9ac493bdf0dec21d3d1"
)


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _path_records(paths: Sequence[str]) -> list[dict[str, str]]:
    return [
        {"path": relative, "sha256": _file_sha256(ROOT / relative)}
        for relative in paths
    ]


def certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def _compare_replays(
    production_orbits: Sequence[Mapping[str, object]],
    independent_orbits: Sequence[Mapping[str, object]],
) -> None:
    if len(production_orbits) != 371 or len(independent_orbits) != 371:
        raise AssertionError("record-12 production/hostile orbit count changed")
    for production, independent in zip(
        production_orbits, independent_orbits, strict=True
    ):
        if (
            production["orbit_index"] != independent["orbit_index"]
            or production["orbit_signature"] != independent["orbit_signature"]
            or production["rank_nine_vertex_count"]
            != independent["rank_nine_vertex_count"]
            or production["minimal_nonface_histogram"]
            != independent["minimal_nonface_histogram"]
            or production["decisive_support_indices"]
            != independent["decisive_support_indices"]
            or production["incremental_lanner_and_actual_vertex_gate"]
            != independent["stages"]
            or not production["compact_coxeter_no_go"]
            or not production["verified"]
        ):
            raise AssertionError(
                f"record-12 production/hostile orbit {production['orbit_index']} differs"
            )


def _secondary_compact_summary(
    replayed: Mapping[str, object],
) -> dict[str, object]:
    orbits = replayed["orbits"]
    stage_histogram = Counter(len(orbit["stages"]) for orbit in orbits)
    stages = [stage for orbit in orbits for stage in orbit["stages"]]
    return {
        "method": replayed["method"],
        "input_file_sha256": replayed["input_file_sha256"],
        "input_certificate_sha256": replayed["input_certificate_sha256"],
        "lanner_tables": replayed["lanner_tables"],
        "orbit_count": replayed["orbit_count"],
        "schedule_histogram": replayed["schedule_histogram"],
        "support_count_histogram": {
            str(key): value for key, value in sorted(stage_histogram.items())
        },
        "stage_count": len(stages),
        "maximum_transient_state_count": max(
            stage["elliptic_actual_vertex_survivor_count"] for stage in stages
        ),
        "all_371_orbits_excluded": replayed["all_selected_orbits_excluded"],
        "full_replay_summary_sha256": replayed["summary_sha256"],
    }


def _compare_secondary_replay(
    production_orbits: Sequence[Mapping[str, object]],
    secondary_orbits: Sequence[Mapping[str, object]],
) -> None:
    if len(production_orbits) != 371 or len(secondary_orbits) != 371:
        raise AssertionError("record-12 secondary replay orbit count changed")
    for production, secondary in zip(
        production_orbits, secondary_orbits, strict=True
    ):
        signature = production["orbit_signature"]
        final_stage = secondary["stages"][-1]
        if (
            production["orbit_index"] != secondary["orbit_index"]
            or signature["x_matrix"] != secondary["x_matrix"]
            or signature["z_matrix"] != secondary["z_matrix"]
            or signature["missing_four_count"]
            != secondary["missing_four_count"]
            or signature["selected_complement_sha256"]
            != secondary["selected_complement_sha256"]
            or production["rank_nine_vertex_count"]
            != secondary["actual_rank_nine_face_count"]
            or not secondary["compact_coxeter_no_go"]
            or final_stage["elliptic_actual_vertex_survivor_count"] != 0
        ):
            raise AssertionError(
                f"record-12 secondary orbit {production['orbit_index']} differs"
            )


def generate() -> dict[str, object]:
    audit = audit_record12_all_coxeter_labels(INPUT)
    if not audit.verified or not audit.whole_shadow_compact_coxeter_no_go:
        raise AssertionError("record-12 whole-shadow compact no-go failed")
    theorem = audit.to_record()
    hostile = independent_replay(INPUT)
    if not hostile["all_371_orbits_excluded"]:
        raise AssertionError("record-12 independent all-orbit replay failed")
    if (
        EXPECTED_HOSTILE_SUMMARY_SHA256 != "TO_BE_FROZEN"
        and hostile["summary_sha256"] != EXPECTED_HOSTILE_SUMMARY_SHA256
    ):
        raise AssertionError("record-12 independent all-orbit replay changed")
    _compare_replays(theorem["orbit_gates"], hostile["orbit_gates"])
    secondary = secondary_independent_replay(True)
    if (
        not secondary["all_selected_orbits_excluded"]
        or secondary["summary_sha256"] != EXPECTED_SECONDARY_SUMMARY_SHA256
    ):
        raise AssertionError("record-12 secondary all-orbit replay changed")
    _compare_secondary_replay(theorem["orbit_gates"], secondary["orbits"])

    document: dict[str, object] = {
        "schema_version": 1,
        "kind": "h9_14_record12_all_orbit_compact_coxeter_no_go",
        "claim": (
            "None of the exactly 371 corrected residual incidence-complex "
            "orbits for compact H9/14 record 12 is a compact hyperbolic "
            "Coxeter nerve; hence the whole graph shadow is compact-Coxeter "
            "impossible."
        ),
        "theorem": theorem,
        "proof": {
            "compact_inputs": [
                "The corrected whole-shadow residual theorem leaves exactly the 371 audited incidence-complex orbits and no others.",
                "Every order-four/five minimal nonface of a compact Coxeter nerve supports a classical Lanner diagram.",
                "Every actual rank-nine vertex diagram of a compact Coxeter 9-polytope is elliptic.",
                "The complete classical Lanner tables contain 108 and 420 labelled placements in orders four and five.",
            ],
            "global_finite_label_bound": (
                "Burcroff, European Journal of Combinatorics 120 (2024), "
                "Lemma 10.10, proves that a compact hyperbolic Coxeter "
                "9-polytope with 14 facets has no multi-multiple edge. With "
                "Definition 4.1 this makes m in {2,3,4,5} the exhaustive "
                "finite dihedral-label alphabet for every meeting facet pair."
            ),
            "finite_search": (
                "For each orbit, reconstruct every minimal nonface and every "
                "rank-nine face from the archived matrices. Starting with "
                "support zero, repeatedly choose the remaining support with "
                "maximum overlap with the current cluster, then minimum new "
                "vertices, then minimum archived index. Join its complete "
                "Lanner placement table and test every actual-vertex "
                "restriction against the complete finite classification in "
                "the globally reduced label alphabet. Every orbit reaches "
                "zero states."
            ),
            "perron_frobenius_step": (
                "On an actual-vertex restriction write the normalized Gram "
                "matrix as G=I-A, with A symmetric and entrywise nonnegative. "
                "G is positive definite exactly when rho(A)<1. Replacing an "
                "unforced right angle by any allowed finite non-right angle "
                "gives A'>=A entrywise, hence rho(A')>=rho(A). A sparse "
                "restriction already rejected as non-elliptic cannot be "
                "repaired by any omitted finite label, so incremental "
                "pruning is lossless."
            ),
            "why_no_global_gram_search_is_needed": (
                "Every possible compact Coxeter labeling already fails a "
                "necessary positive-definite principal vertex condition; "
                "continuous dotted entries are never reached."
            ),
        },
        "independent_hostile_replay": compact_summary(hostile),
        "secondary_independent_replay": _secondary_compact_summary(secondary),
        "scope": {
            "record_index": 12,
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "selected_nerve_only": False,
            "all_371_residual_incidence_orbits_covered": True,
            "compact_only": True,
            "finite_volume_noncompact_chamber_excluded": False,
            "pl_sphericity_of_abstract_residuals_claimed": False,
            "convex_polytopality_of_abstract_residuals_claimed": False,
            "global_lorentzian_gram_values_enumerated": False,
        },
        "sources": _path_records(SOURCE_PATHS),
        "inputs": _path_records(INPUT_PATHS),
    }
    document["certificate_sha256"] = certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record12_all_orbit_compact_coxeter_no_go":
        raise ValueError("unexpected record-12 all-orbit certificate kind")
    if document.get("sources") != _path_records(SOURCE_PATHS):
        raise ValueError("record-12 all-orbit proof source closure changed")
    if document.get("inputs") != _path_records(INPUT_PATHS):
        raise ValueError("record-12 all-orbit proof input closure changed")
    if document.get("certificate_sha256") != certificate_digest(document):
        raise ValueError("record-12 all-orbit certificate digest changed")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and document.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-12 all-orbit certificate is not frozen")
    rebuilt = generate()
    if document != rebuilt:
        raise ValueError("record-12 all-orbit certificate does not replay")
    theorem = document["theorem"]
    if (
        not theorem["verified"]
        or not theorem["whole_shadow_compact_coxeter_no_go"]
        or theorem["orbits_excluded"] != 371
        or not document["independent_hostile_replay"][
            "all_371_orbits_excluded"
        ]
        or not document["secondary_independent_replay"][
            "all_371_orbits_excluded"
        ]
    ):
        raise ValueError("record-12 all-orbit theorem disappeared")
    return {
        "record_index": document["scope"]["record_index"],
        "residual_orbits": len(theorem["orbit_gates"]),
        "support_count_histogram": theorem["support_count_histogram"],
        "whole_graph_shadow_compact_coxeter_excluded": document["scope"][
            "whole_graph_shadow_compact_coxeter_excluded"
        ],
        "hostile_summary_sha256": document["independent_hostile_replay"][
            "full_replay_summary_sha256"
        ],
        "secondary_summary_sha256": document["secondary_independent_replay"][
            "full_replay_summary_sha256"
        ],
        "certificate_sha256": document["certificate_sha256"],
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    document = (
        json.loads(OUTPUT.read_text(encoding="utf-8"))
        if args.verify
        else generate()
    )
    if args.write:
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
