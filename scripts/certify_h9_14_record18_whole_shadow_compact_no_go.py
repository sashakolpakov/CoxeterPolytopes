#!/usr/bin/env python3
"""Freeze the complete compact Coxeter exclusion for H9/14 record 18.

The certificate composes independently replayed branch audits.  It deliberately
makes no assertion about ideal vertices or finite-volume noncompact chambers.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from types import ModuleType
from typing import Mapping


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/h9_14_record18_whole_shadow_compact_no_go.json"
)

POSITIVE_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record18_positive_independent_audit.json"
)
A3_ARCHIVE = ROOT / "data/enumerations/h9_14_record18_a3_independent_audit.json"
A0_ARCHIVE = ROOT / "data/enumerations/h9_14_record18_a0_independent_audit.json"
A0_DOWNSTREAM = (
    ROOT / "data/enumerations/h9_14_record18_a0_downstream_independent_audit.json"
)
HIGH_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record18_high_addition_independent_audit.json"
)
A1_RAW_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record18_a1_ridge_raw_families.json"
)
A1_DOWNSTREAM = (
    ROOT / "data/enumerations/h9_14_record18_a1_downstream_independent_audit.json"
)
A2_RAW_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record18_a2_add4_ridge_raw_families.json"
)
A2_DOWNSTREAM = (
    ROOT / "data/enumerations/h9_14_record18_a2_downstream_independent_audit.json"
)

POSITIVE_CERTIFIER = (
    ROOT / "scripts/certify_independent_h9_14_record18_positive_audit.py"
)
A3_CERTIFIER = ROOT / "scripts/certify_independent_h9_14_record18_a3_audit.py"
A1_RAW_CERTIFIER = (
    ROOT / "scripts/certify_h9_14_record18_a1_ridge_raw_families.py"
)
A2_RAW_CERTIFIER = (
    ROOT / "scripts/certify_h9_14_record18_a2_add4_ridge_raw_families.py"
)

EVIDENCE_FILES = (
    POSITIVE_ARCHIVE,
    A3_ARCHIVE,
    A0_ARCHIVE,
    A0_DOWNSTREAM,
    HIGH_ARCHIVE,
    A1_RAW_ARCHIVE,
    A1_DOWNSTREAM,
    A2_RAW_ARCHIVE,
    A2_DOWNSTREAM,
    ROOT / "scripts/independent_h9_14_record18_positive_orbit_audit.cpp",
    ROOT / "scripts/independent_h9_14_record18_a4_exception_audit.cpp",
    POSITIVE_CERTIFIER,
    ROOT / "scripts/independent_h9_14_record18_addition_cover_audit.cpp",
    A3_CERTIFIER,
    ROOT / "scripts/independent_h9_14_record18_a0_balanced_cut_audit.cpp",
    ROOT / "scripts/certify_independent_h9_14_record18_a0_audit.py",
    ROOT / "scripts/independent_h9_14_record18_a0_downstream_audit.py",
    ROOT / "scripts/independent_h9_14_record18_a2_add5_audit.cpp",
    ROOT / "scripts/certify_independent_h9_14_record18_high_addition_audit.py",
    ROOT / "scripts/explore_h9_14_record18_ridge_aware_additions.cpp",
    A1_RAW_CERTIFIER,
    ROOT / "scripts/independent_h9_14_record18_a1_downstream_audit.py",
    A2_RAW_CERTIFIER,
    ROOT / "scripts/independent_h9_14_record18_a2_downstream_audit.py",
    Path(__file__).resolve(),
)

EXPECTED_SCHEDULES = (
    (0, 8, 110, 42),
    (1, 14, 108, 36),
    (2, 20, 106, 30),
    (3, 26, 104, 24),
    (4, 32, 102, 18),
    (5, 38, 100, 12),
    (6, 44, 98, 6),
)
EXPECTED_A2_CASES = {
    (53, 197),
    (53, 54),
    (53, 86),
    (53, 198),
    (23, 101),
    (53, 58),
    (53, 90),
    (53, 202),
    (23, 57),
    (23, 105),
    (15, 53),
    (23, 27),
    (23, 43),
    (23, 45),
}


def _load_module(name: str, path: Path) -> ModuleType:
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load verifier {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _load(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"malformed evidence document {path}")
    return value


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_files() -> list[dict[str, str]]:
    return [
        {
            "path": str(path.relative_to(ROOT)),
            "sha256": _file_sha256(path),
        }
        for path in EVIDENCE_FILES
    ]


def _canonical_digest(document: Mapping[str, object]) -> str:
    critical = dict(document)
    critical.pop("certificate_sha256", None)
    payload = json.dumps(
        critical, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return sha256(payload).hexdigest()


def _all_terminals_zero(document: Mapping[str, object]) -> bool:
    classes = document.get("classes")
    return isinstance(classes, list) and all(
        isinstance(entry, Mapping)
        and isinstance(entry.get("coxeter_labels"), Mapping)
        and entry["coxeter_labels"].get("terminal_state_count") == 0
        for entry in classes
    )


def verify_evidence() -> dict[str, object]:
    positive_verifier = _load_module(
        "record18_positive_aggregate_verifier", POSITIVE_CERTIFIER
    )
    a3_verifier = _load_module("record18_a3_aggregate_verifier", A3_CERTIFIER)
    a1_raw_verifier = _load_module(
        "record18_a1_raw_aggregate_verifier", A1_RAW_CERTIFIER
    )
    a2_raw_verifier = _load_module(
        "record18_a2_raw_aggregate_verifier", A2_RAW_CERTIFIER
    )

    positive = _load(POSITIVE_ARCHIVE)
    positive_summary = positive_verifier.verify(positive)
    if tuple(map(tuple, positive["schedules"])) != EXPECTED_SCHEDULES:
        raise AssertionError("record-18 schedule table changed")
    # The independently derived global value is g4=7-a.  It is zero at the
    # first omitted schedule and strictly decreases thereafter.
    global_g4_intercept, global_g4_slope = 7, -1
    if global_g4_intercept + 7 * global_g4_slope != 0 or global_g4_slope >= 0:
        raise AssertionError("a>=7 did not fail the strict global g4 bound")
    if any(total >= 14 for a, _, _, total in EXPECTED_SCHEDULES if a in (5, 6)):
        raise AssertionError("a=5,6 did not fail the summed strict link-g4 bound")
    if (
        positive_summary["positive_orbit_count"] != 70
        or positive_summary["a4_orbit_count"] != 13
        or positive_summary["a4_schedule_incidence_empty"] is not True
    ):
        raise AssertionError("positive/a4 evidence chain is incomplete")

    a3 = _load(A3_ARCHIVE)
    a3_summary = a3_verifier.verify(a3)
    if (
        a3_summary["input_positive_orbit_count"] != 40
        or a3_summary["scheduled_pure_104_facet_target_count"] != 0
        or a3_summary["a3_schedule_incidence_empty"] is not True
    ):
        raise AssertionError("a3 evidence chain is incomplete")

    a0 = _load(A0_ARCHIVE)
    if (
        a0.get("literal_family_sets_equal") is not True
        or a0.get("canonical_normalized_multisets_equal") is not True
        or a0.get("independent_survivor_count") != 196
        or a0.get("production_survivor_count") != 196
        or a0.get("canonical_orbit_count") != 3
    ):
        raise AssertionError("a0 exact family evidence is incomplete")
    a0_downstream = _load(A0_DOWNSTREAM)
    if (
        a0_downstream.get("a0_compact_branch_excluded") is not True
        or a0_downstream.get("canonical_class_count") != 3
        or a0_downstream.get("class_multiplicities") != [16, 36, 144]
        or a0_downstream.get("all_classes_pass_connected_mod2_homology_audit")
        is not True
        or a0_downstream.get("all_classes_have_unrestricted_compact_coxeter_no_go")
        is not True
        or not _all_terminals_zero(a0_downstream)
    ):
        raise AssertionError("a0 downstream evidence is incomplete")

    high = _load(HIGH_ARCHIVE)
    cases = {
        tuple(case["representative_A4_masks"]): case for case in high["cases"]
    }
    a1_cases = {key: value for key, value in cases.items() if len(key) == 1}
    a2_cases = {key: value for key, value in cases.items() if len(key) == 2}
    if set(a1_cases) != {(15,), (23,), (53,)} or set(a2_cases) != EXPECTED_A2_CASES:
        raise AssertionError("high-addition orbit partition changed")
    if (
        a1_cases[(15,)]["canonical_target_count"] != 0
        or a1_cases[(23,)]["canonical_target_count"] != 0
        or a1_cases[(53,)]["canonical_target_count"] != 8
        or a1_cases[(53,)]["literal_target_orbit_weight_sum"] != 1080
    ):
        raise AssertionError("a1 high-addition transcript changed")

    a1_raw = _load(A1_RAW_ARCHIVE)
    a1_raw_verifier.verify_record(a1_raw)
    a1_downstream = _load(A1_DOWNSTREAM)
    if (
        a1_downstream.get("a1_compact_branch_excluded") is not True
        or a1_downstream.get("raw_target_weight") != 1080
        or a1_downstream.get("canonical_class_count") != 8
        or a1_downstream.get("all_classes_pass_connected_mod2_homology_audit")
        is not True
        or a1_downstream.get("all_classes_have_unrestricted_compact_coxeter_no_go")
        is not True
        or not _all_terminals_zero(a1_downstream)
    ):
        raise AssertionError("a1 downstream evidence is incomplete")

    a2_raw = _load(A2_RAW_ARCHIVE)
    a2_raw_verifier.verify_record(a2_raw)
    a2_downstream = _load(A2_DOWNSTREAM)
    if (
        a2_downstream.get("a2_compact_branch_excluded") is not True
        or a2_downstream.get("a2_representative_count") != 14
        or a2_downstream.get("a2_add4_target_raw_weight") != 8
        or a2_downstream.get("a2_add4_literal_family_count_audited") != 8
        or a2_downstream.get("all_target_classes_pass_connected_mod2_homology_audit")
        is not True
        or a2_downstream.get("all_target_classes_have_unrestricted_compact_coxeter_no_go")
        is not True
        or not _all_terminals_zero(a2_downstream)
    ):
        raise AssertionError("a2 downstream evidence is incomplete")

    add5 = [case for case in a2_cases.values() if case["additional_A5_count"] == 5]
    add4 = [case for case in a2_cases.values() if case["additional_A5_count"] == 4]
    if (
        len(add5) != 2
        or any(case["canonical_target_count"] != 0 for case in add5)
        or len(add4) != 12
        or sum(case["literal_target_orbit_weight_sum"] for case in add4) != 8
    ):
        raise AssertionError("a2 incidence partition changed")

    return {
        "schedule_rows": [list(row) for row in EXPECTED_SCHEDULES],
        "a_ge_7_global_g4_excluded": True,
        "a5_a6_summed_link_g4_excluded": True,
        "a4_positive_orbit_count": 13,
        "a4_incidence_target_count": 0,
        "a3_positive_orbit_count": 40,
        "a3_incidence_target_count": 0,
        "a2_positive_orbit_count": 14,
        "a2_add5_zero_target_orbit_count": 2,
        "a2_add4_literal_target_count": 8,
        "a2_label_survivor_count": 0,
        "a1_positive_orbit_count": 3,
        "a1_literal_target_count": 1080,
        "a1_target_orbit_count": 8,
        "a1_label_survivor_count": 0,
        "a0_literal_target_count": 196,
        "a0_target_orbit_count": 3,
        "a0_label_survivor_count": 0,
        "all_nonnegative_a_schedules_closed": True,
    }


def fresh_record() -> dict[str, object]:
    branch_summary = verify_evidence()
    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record18_whole_shadow_compact_coxeter_exclusion",
        "claim": (
            "No compact hyperbolic Coxeter 9-polytope with fourteen facets "
            "has record-18 missing-edge graph M???????????o?K?_."
        ),
        "proof_scope": (
            "Whole graph shadow under the compact H9/14 missing-face/Lanner "
            "reductions, not one selected nerve."
        ),
        "branch_summary": branch_summary,
        "closed_schedules": [0, 1, 2, 3, 4, 5, 6],
        "compact_record18_excluded": True,
        "whole_shadow_compact_coxeter_theorem": True,
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_example": False,
        "nonclaims": [
            "No finite-volume noncompact record-18 chamber is excluded.",
            "No ideal-vertex or maximal-parabolic Vinberg endpoint audit is run.",
            "No hyperbolic Coxeter polytope is constructed.",
            "This is one graph-shadow exclusion, not the full H9/14 theorem.",
        ],
        "source_files": _source_files(),
    }
    document["certificate_sha256"] = _canonical_digest(document)
    return document


def verify_record(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != (
        "h9_14_record18_whole_shadow_compact_coxeter_exclusion"
    ):
        raise ValueError("wrong record-18 aggregate certificate kind")
    if document.get("certificate_sha256") != _canonical_digest(document):
        raise ValueError("record-18 aggregate certificate payload changed")
    if document.get("source_files") != _source_files():
        raise ValueError("record-18 aggregate evidence files changed")
    current = verify_evidence()
    if document.get("branch_summary") != current:
        raise ValueError("record-18 aggregate branch summary changed")
    if (
        document.get("closed_schedules") != list(range(7))
        or document.get("compact_record18_excluded") is not True
        or document.get("whole_shadow_compact_coxeter_theorem") is not True
        or document.get("compact_only") is not True
        or document.get("finite_volume_noncompact_checked") is not False
        or document.get("finite_volume_noncompact_chamber_excluded") is not False
        or document.get("new_hyperbolic_example") is not False
    ):
        raise ValueError("record-18 aggregate scope or conclusion changed")
    return {
        "compact_record18_excluded": True,
        "closed_schedules": list(range(7)),
        "finite_volume_noncompact_checked": False,
        "new_hyperbolic_example": False,
        "certificate_sha256": document["certificate_sha256"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.fresh or args.write:
        document = fresh_record()
    else:
        document = _load(OUTPUT)
    summary = verify_record(document)
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps({"output": str(OUTPUT), **summary}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
