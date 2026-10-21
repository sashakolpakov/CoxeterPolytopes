#!/usr/bin/env python3
"""Certify the raw-Stage-0 to PB-Atlas-Stage-A interface.

The external input of the composed numerical compiler is only an integral
unimodular Lorentzian Gram matrix and primitive oriented full-span walls.
Raw Stage 0 derives the canonical class, physical scales, and the Coble
boundary.  This module converts each accepted index-two branch to the closed
PB-Atlas physical-branch schema and checks it with the independent Stage-A
parser.

This is deliberately only an interface theorem.  It does not assert that the
numerical classes are effective curves, that a surface exists, or that its
Mori cone is complete.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_raw_adjunction_log_stage0 import (
    FIXTURES as RAW_FIXTURES,
    solve_raw_input,
)
from scripts.pb_atlas_contraction_core import (
    INPUT_SCHEMA as PB_INPUT_SCHEMA,
    parse_physical_branch,
)

PB_FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_a_inputs.json"
ARCHIVE = ROOT / "data/enumerations/raw_to_pb_atlas_stage0_a_bridge.json"
SCHEMA = "coxeter.raw_to_pb_atlas_stage0_a_bridge.v1"

CASE_TO_ALIAS = {
    "coble_type_ii_h10_30": "H10/30",
    "coble_type_i_h11_26": "H11/26",
    "coble_type_vi_char3_h14_60": "H14/60",
    "coble_type_mii_char3_h17_48": "H17/48",
}


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def candidate_to_pb_input(
    raw_input: Mapping[str, Any], candidate: Mapping[str, Any]
) -> dict[str, Any]:
    """Convert one accepted index-two Stage-0 candidate to Stage-A input."""

    if not candidate.get("accepted"):
        raise ValueError("the Stage-0 candidate is not accepted")
    if candidate.get("branch") != "coble_index_two_log_boundary":
        raise ValueError("only the index-two Coble branch enters PB-Atlas")
    boundary = [int(index) for index in candidate["boundary_indices_zero_based_inferred"]]
    ramification = [int(value) for value in candidate["ramification_indices_inferred"]]
    if not boundary or any(value != 2 for value in ramification):
        raise ValueError("the inferred branch is not a nonempty index-two boundary")
    scales = [int(value) for value in candidate["physical_scales"]]
    weights = [scales[index] for index in boundary]
    payload = {
        "schema": PB_INPUT_SCHEMA,
        "picard_gram": [list(row) for row in raw_input["picard_gram"]],
        "wall_vectors": [list(row) for row in raw_input["wall_vectors"]],
        "boundary_indices_zero_based": boundary,
        "boundary_weights": weights,
        "canonical_class": [int(value) for value in candidate["canonical_class"]],
    }
    # The parser independently checks unimodularity, full span, acute
    # orientation, characteristic parity, K^2, and sum q_i r_i = -2K.
    return parse_physical_branch(payload).input_payload


def build_certificate() -> dict[str, Any]:
    raw_fixtures = _load(RAW_FIXTURES)
    pb_fixtures = _load(PB_FIXTURES)
    pb_by_alias = {
        case["calibration_alias"]: case["input"] for case in pb_fixtures["cases"]
    }
    rows: list[dict[str, Any]] = []
    for case in raw_fixtures["cases"]:
        case_id = case["case_id"]
        raw_input = case["raw_input"]
        stage0 = solve_raw_input(raw_input)
        checks: dict[str, bool] = {
            "raw_input_has_only_gram_and_walls": set(raw_input)
            == {"schema", "picard_gram", "wall_vectors"},
            "stage0_passes": stage0["status"] == "PASS",
        }
        row: dict[str, Any] = {
            "case_id": case_id,
            "stage0_accepted_candidate_count": stage0["accepted_candidate_count"],
        }
        if case_id in CASE_TO_ALIAS:
            alias = CASE_TO_ALIAS[case_id]
            accepted = stage0["accepted_candidates"]
            checks["unique_stage0_branch"] = len(accepted) == 1
            if len(accepted) != 1:
                raise ArithmeticError(f"{case_id} lost its unique Stage-0 branch")
            payload = candidate_to_pb_input(raw_input, accepted[0])
            expected = pb_by_alias[alias]
            checks.update(
                {
                    "stage0_branch_is_index_two": accepted[0]["branch"]
                    == "coble_index_two_log_boundary",
                    "derived_payload_equals_frozen_stage_a_input": payload == expected,
                    "derived_boundary_weights_are_physical_scales": payload[
                        "boundary_weights"
                    ]
                    == [
                        accepted[0]["physical_scales"][index]
                        for index in accepted[0]["boundary_indices_zero_based_inferred"]
                    ],
                    "stage_a_parser_accepts_derived_payload": parse_physical_branch(
                        payload
                    ).input_payload
                    == payload,
                }
            )
            row.update(
                {
                    "routing": "PB_ATLAS_NUMERICAL_STAGE_A",
                    "calibration_alias": alias,
                    "derived_pb_atlas_input": payload,
                }
            )
        elif case_id == "toric_h6_9":
            accepted = stage0["accepted_candidates"]
            checks.update(
                {
                    "unique_weak_branch": len(accepted) == 1
                    and accepted[0]["branch"] == "weak_minus_one_minus_two",
                    "no_boundary_was_invented": len(
                        accepted[0]["boundary_indices_zero_based_inferred"]
                    )
                    == 0,
                }
            )
            row["routing"] = "DIRECT_WEAK_ANTICANONICAL_BRANCH"
        elif case_id == "blowup_F4_maximal_triangle_negative_control":
            checks.update(
                {
                    "maximal_triangle_has_no_accepted_branch": stage0[
                        "accepted_candidate_count"
                    ]
                    == 0,
                    "maximal_triangle_is_not_routed_to_stage_a": not stage0[
                        "accepted_candidates"
                    ],
                }
            )
            row["routing"] = "REJECTED_AT_LOG_EQUATION"
        else:
            raise ArithmeticError(f"unrouted raw fixture {case_id}")
        row["checks"] = checks
        row["all_checks_pass"] = all(checks.values())
        rows.append(row)

    result = {
        "schema": SCHEMA,
        "status": "PASS_RAW_STAGE0_TO_PB_STAGE_A",
        "source_sha256": {
            str(RAW_FIXTURES.relative_to(ROOT)): _sha256(RAW_FIXTURES),
            str(PB_FIXTURES.relative_to(ROOT)): _sha256(PB_FIXTURES),
            "scripts/certify_raw_adjunction_log_stage0.py": _sha256(
                ROOT / "scripts/certify_raw_adjunction_log_stage0.py"
            ),
            "scripts/pb_atlas_contraction_core.py": _sha256(
                ROOT / "scripts/pb_atlas_contraction_core.py"
            ),
        },
        "case_count": len(rows),
        "derived_stage_a_case_count": sum(
            row["routing"] == "PB_ATLAS_NUMERICAL_STAGE_A" for row in rows
        ),
        "direct_weak_case_count": sum(
            row["routing"] == "DIRECT_WEAK_ANTICANONICAL_BRANCH" for row in rows
        ),
        "rejected_case_count": sum(
            row["routing"] == "REJECTED_AT_LOG_EQUATION" for row in rows
        ),
        "cases": rows,
        "scope": {
            "proved": [
                "the composed numerical interface starts from Gram plus oriented walls only",
                "accepted index-two Stage-0 branches determine every Stage-A input field",
                "the four frozen Coble Stage-A payloads are reproduced exactly",
            ],
            "not_proved": [
                "surface existence or nonempty incidence strata",
                "wall effectivity, irreducibility, or complete Mori generation",
                "same-polytope return or essential surjectivity",
            ],
        },
    }
    result["all_checks_pass"] = all(row["all_checks_pass"] for row in rows)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_certificate()
    if args.write:
        ARCHIVE.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.check:
        frozen = _load(ARCHIVE)
        if result != frozen:
            raise SystemExit("raw-to-PB bridge archive mismatch")
    print(
        json.dumps(
            {
                "status": result["status"],
                "all_checks_pass": result["all_checks_pass"],
                "derived_stage_a_case_count": result["derived_stage_a_case_count"],
                "direct_weak_case_count": result["direct_weak_case_count"],
                "rejected_case_count": result["rejected_case_count"],
            },
            sort_keys=True,
        )
    )
    return 0 if result["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
