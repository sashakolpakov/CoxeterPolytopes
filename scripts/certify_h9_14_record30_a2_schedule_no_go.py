#!/usr/bin/env python3
"""Freeze the complete compact Record-30 a=2 schedule no-go theorem."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import permutations
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_h9_14_record30_positive_scout import (
    ARCHIVE as POSITIVE_ARCHIVE,
    verify as verify_positive,
)
from scripts.explore_h9_14_record30_a2_orbit10_cross_supersets import (
    verify as verify_orbit10,
)
from scripts.independent_h9_14_record30_a2_all_target_labels import (
    replay as replay_all_labels,
)
from scripts.independent_h9_14_record30_a2_exceptional_whole_nerves import (
    replay as replay_exceptional,
)


PRODUCTION_SOURCE = (
    ROOT / "scripts/explore_h9_14_record30_a2_orbit14_nauty.cpp"
)
INDEPENDENT_INCIDENCE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record30_a2_all_orbits_explicit.cpp"
)
LABEL_ENGINE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record30_a2_orbit14_labels.py"
)
ALL_LABEL_SOURCE = (
    ROOT / "scripts/independent_h9_14_record30_a2_all_target_labels.py"
)
EXCEPTIONAL_SOURCE = (
    ROOT / "scripts/independent_h9_14_record30_a2_exceptional_whole_nerves.py"
)
SELECTED_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_selected_nerve_coxeter_no_go.json"
)
ORBIT10_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record30_a2_orbit10_whole_nerve_no_go.json"
)
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_a2_schedule_no_go.json"

PROOF_PATHS = (
    PRODUCTION_SOURCE,
    INDEPENDENT_INCIDENCE_SOURCE,
    LABEL_ENGINE_SOURCE,
    ALL_LABEL_SOURCE,
    EXCEPTIONAL_SOURCE,
    ROOT / "scripts/certify_h9_14_record30_positive_scout.py",
    ROOT / "scripts/certify_h9_14_selected_nerve_coxeter_no_go.py",
    ROOT / "scripts/explore_h9_14_record30_a2_orbit10_cross_supersets.py",
)

EXPECTED_SOURCE_SHA256: dict[str, str] = {
    "scripts/explore_h9_14_record30_a2_orbit14_nauty.cpp": (
        "1fc448ee77a5d2eab0deb478db8e02302e3cc638b121ffcef0e0a1c5194183dd"
    ),
    "scripts/independent_h9_14_record30_a2_all_orbits_explicit.cpp": (
        "502056425673e681056aecc398de0e29587ee7bd160d8f8e30f8ad7754307a57"
    ),
    "scripts/independent_h9_14_record30_a2_orbit14_labels.py": (
        "d5b51d5d6dea0f645df308332d7e4a9b062e5c2aaca7ae835475a3e1cb662c81"
    ),
    "scripts/independent_h9_14_record30_a2_all_target_labels.py": (
        "93d06d49e7f9dc366dd115df1ed3ee63ddbc5ddd2100a7a18f59969f7674b27e"
    ),
    "scripts/independent_h9_14_record30_a2_exceptional_whole_nerves.py": (
        "841651f45ff877a580f2a6dbd024b44dedb2677437d82c617d708a2fb8d0da43"
    ),
    "scripts/certify_h9_14_record30_positive_scout.py": (
        "dfa7f24230126ae3bd07b289152c38738418f4902757001ce79d54dd521592d5"
    ),
    "scripts/certify_h9_14_selected_nerve_coxeter_no_go.py": (
        "62d0f779125e03f0cce8d190867caeedc68118ff950254b21b0226204da9ad11"
    ),
    "scripts/explore_h9_14_record30_a2_orbit10_cross_supersets.py": (
        "27218bf0071ab488e85701715f94f4408d775ed045bddc31f9cf88333d2e014f"
    ),
}
EXPECTED_POSITIVE_FILE_SHA256 = (
    "841ebf50748747a8dc96fb4ce337570f29f1f9455f2206609e18d31187eb3771"
)
EXPECTED_SELECTED_FILE_SHA256 = (
    "028bbc8e266145e59fa1ee9f9796a6e2c241f4a4f09637bd73f98bced9908629"
)
EXPECTED_SELECTED_CERTIFICATE_SHA256 = (
    "c2dfcbfc03dc5dbae34938ab32f3a1c1111fd729d7844fd76d04dfc6cfa44de8"
)
EXPECTED_ORBIT10_FILE_SHA256 = (
    "8696dd652caed5a4a0aceba8283bc27e91fe8abd8fa3048ce556a9252d9026ff"
)
EXPECTED_ORBIT10_CERTIFICATE_SHA256 = (
    "879c927f7bedb2d6fdb59a0eaccab7f936ad688e67444b8e74dddb750f8253b3"
)
EXPECTED_PRODUCTION_SHA256 = (
    "46e613f5d2102c9aba0701888e062d57228b20bc0f90f3d2db971ec837120a9c"
)
EXPECTED_INDEPENDENT_INCIDENCE_SHA256 = (
    "974cfc3fe3b3023efd44387f463bb18ab16caf4745328aed38db849b54eadf29"
)
EXPECTED_LABEL_SUMMARY_SHA256 = (
    "34c1fe75cea741d7e56afc093f2b9b29f438cecd79dd39e9a35982e52213b39e"
)
EXPECTED_EXCEPTIONAL_SUMMARY_SHA256 = (
    "453eb9ef4b5e7578e965fd37ea375b434f30a2cab873213a399e23dc0829ca78"
)
EXPECTED_SELECTED_RECORD30_PROJECTION_SHA256 = (
    "2c4aa79bef0d968189ea505cc4bfdfb8cefb194f31b5cbae6cb347866cfff9f7"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "3ac5e4ea4ddd624aef17dec45f64d403165b736bab015d088e2b8d78f43a5adf"
)


def _canonical_sha(value: object) -> str:
    normalized = json.loads(json.dumps(value, ensure_ascii=True))
    return sha256(
        json.dumps(
            normalized,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode()
    ).hexdigest()


def _certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return _canonical_sha(payload)


def _file_sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _source_record(path: Path) -> dict[str, str]:
    return {"path": str(path.relative_to(ROOT)), "sha256": _file_sha(path)}


def _require(expected: str, observed: str, message: str) -> None:
    if expected not in ("TO_BE_FROZEN", observed):
        raise ValueError(message)


def _positive_branch() -> tuple[dict[str, object], str]:
    file_digest = _file_sha(POSITIVE_ARCHIVE)
    if file_digest != EXPECTED_POSITIVE_FILE_SHA256:
        raise ValueError("Record-30 positive-A4 archive changed")
    document = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    verify_positive(document)
    branch = next(
        row for row in document["a1_a4"]["branches"] if row["a"] == 2
    )
    if len(branch["orbits"]) != 21:
        raise ValueError("Record-30 a=2 positive-A4 quotient changed")
    return branch, file_digest


def _representatives(branch: Mapping[str, object]) -> list[tuple[int, int]]:
    answer = []
    for index, orbit in enumerate(branch["orbits"]):  # type: ignore[index]
        representative = orbit["representative"]
        if len(representative) != 2:
            raise ValueError(f"a=2 orbit {index} is not an A4 pair")
        masks = tuple(
            sum(1 << int(vertex) for vertex in face) for face in representative
        )
        if any(mask.bit_count() != 4 for mask in masks):
            raise ValueError("an upstream A4 mask changed")
        answer.append(tuple(sorted(masks)))
    return answer


def _compiler() -> str:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("a C++20 compiler is required for a fresh replay")
    return compiler


def _nauty_paths() -> tuple[Path, Path]:
    candidates = []
    if os.environ.get("NAUTY_PREFIX"):
        candidates.append(Path(os.environ["NAUTY_PREFIX"]))
    candidates.extend(
        (Path("/opt/homebrew/opt/nauty"), Path("/usr/local/opt/nauty"))
    )
    for prefix in candidates:
        header = prefix / "include/nauty/nauty.h"
        library = prefix / "lib/libnauty.a"
        if header.is_file() and library.is_file():
            return header.parent, library
    raise RuntimeError(
        "nauty headers/static library were not found; set NAUTY_PREFIX"
    )


def _compile(
    source: Path, executable: Path, *, nauty: bool = False
) -> None:
    command = [
        _compiler(),
        "-std=c++20",
        "-O3",
        "-Wall",
        "-Wextra",
        "-Werror",
        "-pedantic",
    ]
    if nauty:
        include, library = _nauty_paths()
        command.extend((f"-I{include}", str(source), str(library)))
    else:
        command.append(str(source))
    command.extend(("-o", str(executable)))
    subprocess.run(command, cwd=ROOT, check=True)


def _run_json(
    command: Sequence[str], *, standard_input: str | None = None
) -> dict[str, object]:
    completed = subprocess.run(
        list(command),
        cwd=ROOT,
        input=standard_input,
        text=True,
        capture_output=True,
        check=True,
    )
    parsed = json.loads(completed.stdout)
    if not isinstance(parsed, dict):
        raise ValueError(f"{command[0]} did not emit a JSON object")
    return parsed


def _production_replay(
    representatives: Sequence[tuple[int, int]],
) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="record30-a2-nauty-") as directory:
        executable = Path(directory) / "production"
        _compile(PRODUCTION_SOURCE, executable, nauty=True)
        rows = [
            _run_json(
                (
                    str(executable),
                    str(index),
                    str(fours[0]),
                    str(fours[1]),
                )
            )
            for index, fours in enumerate(representatives)
        ]
    result = {
        "kind": "h9_14_record30_a2_all_orbits_nauty_incidence_replay",
        "compact_only": True,
        "positive_A4_orbit_count": len(rows),
        "pure_128_facet_target_orbit_count": sum(
            int(row["pure_128_facet_target_orbit_count"]) for row in rows
        ),
        "orbit_audits": rows,
        "scope": (
            "Exact nauty colored-incidence quotient, cover intersection, and "
            "literal blocker replay for all 21 compact Record-30 a=2 "
            "positive-A4 representatives. No finite-volume inference."
        ),
    }
    return result


def _independent_incidence_replay(
    representatives: Sequence[tuple[int, int]],
) -> dict[str, object]:
    lines = [str(len(representatives))]
    lines.extend(
        " ".join(map(str, (index, *fours)))
        for index, fours in enumerate(representatives)
    )
    with tempfile.TemporaryDirectory(prefix="record30-a2-explicit-") as directory:
        executable = Path(directory) / "independent"
        _compile(INDEPENDENT_INCIDENCE_SOURCE, executable)
        return _run_json((str(executable),), standard_input="\n".join(lines) + "\n")


def _colored_incidence_code(
    fours: Sequence[int], additions: Sequence[int]
) -> tuple[int, ...]:
    edges = ((1 << 0) | (1 << 11), (1 << 1) | (1 << 12), (1 << 2) | (1 << 13))
    best = None
    for edge_order in permutations(range(3)):
        for four_order in permutations(range(2)):
            for addition_order in permutations(range(len(additions))):
                signatures = []
                for vertex in range(14):
                    signature = 0
                    for position, index in enumerate(edge_order):
                        signature |= bool(edges[index] >> vertex & 1) << position
                    for position, index in enumerate(four_order):
                        signature |= bool(fours[index] >> vertex & 1) << (3 + position)
                    for position, index in enumerate(addition_order):
                        signature |= bool(additions[index] >> vertex & 1) << (5 + position)
                    signatures.append(signature)
                code = tuple(sorted(signatures))
                if best is None or code < best:
                    best = code
    if best is None:
        raise AssertionError("colored incidence canonicalization failed")
    return best


COMMON_KEYS = (
    "positive_A4_orbit_index",
    "fixed_A4_masks",
    "stabilizer_order",
    "forced_shadow_count",
    "additional_A5_count",
    "candidate_count",
    "raw_unordered_addition_count",
    "base_small_transversal_count",
    "initial_order_five_cover_count",
    "search_node_count",
    "canonical_state_counts",
    "facet_lower_bound_prune_count",
    "uncoverable_prune_count",
    "scheduled_incidence_orbit_count",
    "scheduled_blocker_count_histogram",
    "pure_128_facet_target_orbit_count",
)


def _incidence_comparison(
    production: Mapping[str, object], independent: Mapping[str, object]
) -> dict[str, object]:
    production_rows = production.get("orbit_audits")
    independent_rows = independent.get("orbit_audits")
    if (
        not isinstance(production_rows, list)
        or not isinstance(independent_rows, list)
        or len(production_rows) != 21
        or len(independent_rows) != 21
    ):
        raise ValueError("incidence replay coverage changed")
    rows = []
    all_keys = []
    for index, (raw_production, raw_independent) in enumerate(
        zip(production_rows, independent_rows, strict=True)
    ):
        if not isinstance(raw_production, Mapping) or not isinstance(
            raw_independent, Mapping
        ):
            raise ValueError("incidence orbit row is malformed")
        left = {key: raw_production[key] for key in COMMON_KEYS}
        right = {key: raw_independent[key] for key in COMMON_KEYS}
        if left != right:
            raise ValueError(f"production/independent mismatch in orbit {index}")
        fours = [int(value) for value in raw_production["fixed_A4_masks"]]
        production_keys = sorted(
            _colored_incidence_code(fours, [int(face) for face in family])
            for family in raw_production["pure_target_representatives"]
        )
        independent_keys = sorted(
            _colored_incidence_code(fours, [int(face) for face in family])
            for family in raw_independent["pure_target_representatives"]
        )
        if production_keys != independent_keys:
            raise ValueError(f"target orbit sets disagree in orbit {index}")
        all_keys.extend((index, list(key)) for key in production_keys)
        rows.append(
            {
                "positive_A4_orbit_index": index,
                "target_orbit_count": len(production_keys),
                "target_orbit_key_sha256": _canonical_sha(production_keys),
                "stabilizer_order": int(raw_independent["stabilizer_order"]),
                "literal_target_count_with_fixed_A4_pair": int(
                    raw_independent["literal_target_count_with_fixed_A4_pair"]
                ),
                "target_orbit_size_histogram": raw_independent[
                    "target_orbit_size_histogram"
                ],
            }
        )
    if len(all_keys) != 77:
        raise AssertionError("complete target-orbit census changed")
    return {
        "production_and_independent_common_counters_identical": True,
        "production_and_independent_target_orbit_sets_identical": True,
        "positive_A4_orbit_count": 21,
        "target_orbit_count": len(all_keys),
        "all_target_orbit_key_sha256": _canonical_sha(all_keys),
        "orbit_rows": rows,
    }


def _selected_record30_dependency() -> dict[str, object]:
    file_digest = _file_sha(SELECTED_ARCHIVE)
    if file_digest != EXPECTED_SELECTED_FILE_SHA256:
        raise ValueError("selected-nerve archive changed")
    document = json.loads(SELECTED_ARCHIVE.read_text(encoding="utf-8"))
    certificate = document.get("certificate_sha256")
    if certificate != EXPECTED_SELECTED_CERTIFICATE_SHA256:
        raise ValueError("selected-nerve certificate changed")
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    if _canonical_sha(payload) != certificate:
        raise ValueError("selected-nerve certificate closure failed")
    record = next(
        row for row in document["six_nerve_census"] if row["record_index"] == 30
    )
    projection_digest = _canonical_sha(record)
    _require(
        EXPECTED_SELECTED_RECORD30_PROJECTION_SHA256,
        projection_digest,
        "selected Record-30 projection changed",
    )
    gate = record["support_gates"][0]
    if (
        gate["left_support"] != [3, 4, 5, 6]
        or gate["right_support"] != [7, 8, 9, 10]
        or gate["sparse_state_count"] != 81
        or gate["single_edge_passing_count"] != 0
        or record["unrestricted_compact_coxeter_no_go"] is not True
    ):
        raise ValueError("selected Record-30 product theorem changed")
    return {
        "path": str(SELECTED_ARCHIVE.relative_to(ROOT)),
        "file_sha256": file_digest,
        "certificate_sha256": certificate,
        "record30_projection_sha256": projection_digest,
        "fixed_A4_masks": [120, 1920],
        "additional_A5_masks": [],
        "normalizing_support_symmetry_orders": record[
            "normalizing_support_symmetry_orders"
        ],
        "single_m3_cross_edge_trial_count": gate["single_edge_test_count"],
        "single_m3_cross_edge_passing_count": gate["single_edge_passing_count"],
        "unrestricted_compact_coxeter_no_go": True,
    }


def _orbit10_dependency() -> dict[str, object]:
    file_digest = _file_sha(ORBIT10_ARCHIVE)
    if file_digest != EXPECTED_ORBIT10_FILE_SHA256:
        raise ValueError("orbit-10 whole-nerve archive changed")
    document = json.loads(ORBIT10_ARCHIVE.read_text(encoding="utf-8"))
    summary = verify_orbit10(document)
    if summary["certificate_sha256"] != EXPECTED_ORBIT10_CERTIFICATE_SHA256:
        raise ValueError("orbit-10 whole-nerve certificate changed")
    return {
        "path": str(ORBIT10_ARCHIVE.relative_to(ROOT)),
        "file_sha256": file_digest,
        "certificate_sha256": summary["certificate_sha256"],
        "fixed_A4_masks": [57, 2496],
        "additional_A5_masks": [1592, 1984],
        "sparse_state_count": summary["sparse_state_count"],
        "first_empty_connected_shell": summary[
            "first_empty_connected_shell"
        ],
        "whole_nerve_compact_coxeter_impossible": True,
        "exact_metric_upgrade_in_this_certificate": True,
    }


def build() -> dict[str, object]:
    branch, positive_digest = _positive_branch()
    representatives = _representatives(branch)
    production = _production_replay(representatives)
    independent = _independent_incidence_replay(representatives)
    comparison = _incidence_comparison(production, independent)
    labels = replay_all_labels(independent)
    exceptional = replay_exceptional()
    selected = _selected_record30_dependency()
    orbit10 = _orbit10_dependency()

    production_digest = _canonical_sha(production)
    independent_digest = _canonical_sha(independent)
    _require(
        EXPECTED_PRODUCTION_SHA256,
        production_digest,
        "production incidence transcript changed",
    )
    _require(
        EXPECTED_INDEPENDENT_INCIDENCE_SHA256,
        independent_digest,
        "independent incidence transcript changed",
    )
    _require(
        EXPECTED_LABEL_SUMMARY_SHA256,
        str(labels["summary_sha256"]),
        "all-target label transcript changed",
    )
    _require(
        EXPECTED_EXCEPTIONAL_SUMMARY_SHA256,
        str(exceptional["summary_sha256"]),
        "exceptional exact-metric transcript changed",
    )
    if (
        comparison["target_orbit_count"] != 77
        or labels["sparse_label_no_go_target_count"] != 74
        or labels["sparse_label_survivor_count"] != 3
        or labels["sparse_label_survivors"]
        != [
            {
                "positive_A4_orbit_index": 1,
                "orbit_target_index": 0,
                "terminal_sparse_state_count": 11664,
            },
            {
                "positive_A4_orbit_index": 7,
                "orbit_target_index": 1,
                "terminal_sparse_state_count": 648,
            },
            {
                "positive_A4_orbit_index": 10,
                "orbit_target_index": 0,
                "terminal_sparse_state_count": 648,
            },
        ]
        or any(
            case["zero_exact_field_norm_count"] != 0
            or case["zero_second_prime_residue_count"] != 0
            or case["whole_nerve_compact_coxeter_impossible"] is not True
            for case in exceptional["cases"]
        )
    ):
        raise AssertionError("Record-30 a=2 theorem chain is incomplete")

    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record30_a2_compact_schedule_no_go",
        "claim": (
            "Every completion of H9/14 graph-shadow Record 30 in the compact "
            "Dehn--Sommerville schedule a=2 fails an exact necessary compact "
            "Coxeter condition. The complete 21-orbit incidence quotient has "
            "77 pure 128-facet targets; 74 fail the unrestricted sparse "
            "Lanner/actual-vertex gate, and the remaining three fail complete "
            "whole-nerve label or exact Gram-rank obstructions."
        ),
        "scope": {
            "dimension": 9,
            "facet_count": 14,
            "record_index": 30,
            "excluded_schedule": 2,
            "a2_compact_schedule_excluded": True,
            "whole_record30_shadow_excluded": False,
            "finite_label_alphabet": [2, 3, 4, 5],
            "compact_only": True,
            "finite_volume_noncompact_chamber_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "nonclaims": [
            "Record-30 compact schedules a=0 and a=1 remain outside this theorem.",
            "No ideal-vertex or parabolic cusp census is performed.",
            "No finite-volume noncompact chamber is excluded.",
            "No Lorentzian Gram realization or new hyperbolic example is produced.",
        ],
        "proof_dependencies": [
            "Burcroff's compact H9/14 no-missing-triples and finite-label reductions",
            "Dehn--Sommerville schedule and strict local link/Macaulay inequalities",
            "facet complements equal the literal blocker of the minimal-nonface clutter",
            "complete order-four/order-five Lanner placement tables (108 and 420 placements)",
            "actual rank-nine vertex diagrams must be elliptic",
            "Perron--Frobenius monotonicity for omitted finite non-right labels",
            "negative-index-one connectivity of the nonzero Gram graph",
            "an H9 chamber Gram matrix has rank at most ten",
        ],
        "source_provenance": [_source_record(path) for path in PROOF_PATHS],
        "upstream_positive_A4_quotient": {
            "path": str(POSITIVE_ARCHIVE.relative_to(ROOT)),
            "file_sha256": positive_digest,
            "a2_positive_A4_orbit_count": 21,
            "representative_A4_masks": [list(row) for row in representatives],
        },
        "production_nauty_incidence_replay": production,
        "production_nauty_transcript_sha256": production_digest,
        "independent_explicit_incidence_replay": independent,
        "independent_explicit_incidence_transcript_sha256": independent_digest,
        "incidence_cross_audit": comparison,
        "independent_all_target_sparse_label_replay": labels,
        "surviving_target_dependencies": {
            "orbit1_target0_selected_product_theorem": selected,
            "orbit10_target0_prior_whole_nerve_theorem": orbit10,
        },
        "independent_exceptional_exact_whole_nerve_replay": exceptional,
        "theorem": {
            "positive_A4_orbit_count": 21,
            "pure_128_facet_target_orbit_count": 77,
            "sparse_label_no_go_target_orbit_count": 74,
            "selected_product_no_go_target_orbit_count": 1,
            "exact_principal_block_no_go_target_orbit_count": 2,
            "a2_compact_schedule_excluded": True,
            "whole_record30_shadow_excluded": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
    }
    document["certificate_sha256"] = _certificate_digest(document)
    return document


def verify(
    document: Mapping[str, object], *, recompute: bool = False
) -> dict[str, object]:
    if document.get("kind") != "h9_14_record30_a2_compact_schedule_no_go":
        raise ValueError("unexpected Record-30 a=2 certificate kind")
    if _file_sha(POSITIVE_ARCHIVE) != EXPECTED_POSITIVE_FILE_SHA256:
        raise ValueError("Record-30 positive dependency drifted")
    if _file_sha(SELECTED_ARCHIVE) != EXPECTED_SELECTED_FILE_SHA256:
        raise ValueError("selected-nerve dependency drifted")
    if _file_sha(ORBIT10_ARCHIVE) != EXPECTED_ORBIT10_FILE_SHA256:
        raise ValueError("orbit-10 dependency drifted")
    source_records = [_source_record(path) for path in PROOF_PATHS]
    if document.get("source_provenance") != source_records:
        raise ValueError("Record-30 a=2 proof source provenance changed")
    if EXPECTED_SOURCE_SHA256 and EXPECTED_SOURCE_SHA256 != {
        row["path"]: row["sha256"] for row in source_records
    }:
        raise ValueError("Record-30 a=2 frozen source hashes changed")
    digest = _certificate_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("Record-30 a=2 certificate closure failed")
    _require(
        EXPECTED_CERTIFICATE_SHA256,
        digest,
        "Record-30 a=2 certificate is not frozen",
    )
    theorem = document.get("theorem")
    scope = document.get("scope")
    if (
        not isinstance(theorem, Mapping)
        or theorem.get("pure_128_facet_target_orbit_count") != 77
        or theorem.get("a2_compact_schedule_excluded") is not True
        or theorem.get("whole_record30_shadow_excluded") is not False
        or not isinstance(scope, Mapping)
        or scope.get("compact_only") is not True
        or scope.get("finite_volume_noncompact_chamber_excluded") is not False
        or scope.get("new_hyperbolic_coxeter_polytope") is not False
    ):
        raise ValueError("Record-30 a=2 theorem scope changed")
    if recompute and json.loads(json.dumps(build(), sort_keys=True)) != json.loads(
        json.dumps(document, sort_keys=True)
    ):
        raise ValueError("Record-30 a=2 full theorem did not replay")
    return {
        "record_index": 30,
        "excluded_compact_schedule": 2,
        "positive_A4_orbit_count": 21,
        "pure_target_orbit_count": 77,
        "sparse_label_no_go_target_count": 74,
        "exceptional_target_count": 3,
        "a2_compact_schedule_excluded": True,
        "whole_record30_shadow_excluded": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args(argv)
    document = (
        json.loads(ARCHIVE.read_text(encoding="utf-8"))
        if args.verify_only
        else build()
    )
    if args.write:
        ARCHIVE.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(verify(document, recompute=args.recompute), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
