#!/usr/bin/env python3
"""Replay the unconditional maximal arithmetic chamber theorem interfaces."""

from __future__ import annotations

import argparse
import hashlib
import json
from math import comb, factorial
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = (
    ROOT
    / "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness.md"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness.json"
)
FINITE_PAIR_REPORT = (
    ROOT / "results/projected_k3_maximal_abstract_chamber_finiteness.md"
)
RATIONAL_CARRIER_REPORT = ROOT / "results/rational_surface_eigencarrier_universality.md"
NIKULIN_AUDIT = ROOT / "results/nikulin_conceptual_dependency_hostile_audit.md"


DIMENSION_MAX = 29
RANK_MAX = DIMENSION_MAX + 1
FIELD_DEGREE_MAX = 25
FRAME_CONSTANT = 62
CONJUGATE_HEIGHT = FRAME_CONSTANT**RANK_MAX
BLOWUP_MAX = 4 * FIELD_DEGREE_MAX * RANK_MAX - 3


DEPENDENCY_SHA256 = {
    "finite_pair_report": (
        "59dd458acf14d18021f6c5583096da13f39c83115fc136d6ba91bcc850ce3603"
    ),
    "rational_carrier_report": (
        "57e607e7cd6646c5827b3d4d0bff86d7bf8294845e29fe95d4e8b1478b854333"
    ),
    "absw_primary_pdf": (
        "974ac7c5141ca1c45b515f145215ad115aba867e42551a9b103233ea52930ee4"
    ),
    "vinberg_1971_primary_pdf": (
        "16a4f497688e6f968b4c63d81a6b049b3f096f1cd4f43e458b2882a5c44ded38"
    ),
    "vinberg_1984_dimension_primary_pdf": (
        "fd25f8f892eed27e4303a8f3252fd384d8954587a6c1c67f1efed2c38aabe44f"
    ),
    "nikulin_2007_primary_pdf": (
        "9a820fef0b88158b283057f7d97acc89130610a9a4aacc6c53a8215e23ce874a"
    ),
    "nikulin_1982_primary_pdf": (
        "696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d"
    ),
    "nikulin_2011_primary_pdf": (
        "716e6e7333efbfaa22060b60dc97e8303363a4914409ccdd56ef7a90c923c7af"
    ),
    "belolipetsky_linowitz_arxiv_v2_pdf": (
        "026e0f635f71e214807e0b505bf7297673eb9e0d47a6a91b7bd6dec134eab97f"
    ),
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def simple_directed_cycle_slot_bound(vertices: int) -> int:
    """Number of ordered simple cycles, modulo cyclic rotation."""

    return sum(
        factorial(vertices) // factorial(vertices - length) // length
        for length in range(2, vertices + 1)
    )


def coefficient_bounds(degree: int, height: int) -> list[int]:
    """Vieta bounds for a monic polynomial with bounded roots."""

    return [comb(degree, index) * height**index for index in range(1, degree + 1)]


def validate_repository_dependencies() -> None:
    actual = {
        "finite_pair_report": sha256_file(FINITE_PAIR_REPORT),
        "rational_carrier_report": sha256_file(RATIONAL_CARRIER_REPORT),
    }
    for key, value in actual.items():
        if value != DEPENDENCY_SHA256[key]:
            raise AssertionError(f"dependency hash changed: {key}")


def validate_report_scope() -> None:
    text = REPORT.read_text(encoding="utf-8")
    required = (
        "Theorem 1.1 (unconditional maximal-chamber finiteness)",
        "2\\le n\\le29",
        "[k(P):\\mathbf Q]\\le25",
        "unconditional global finiteness theorem outright",
        "b_C^2=\\prod",
        "62^{30}",
        "Theorem 2.1.3(2)",
        "Theorem 5.2.1",
        "W(P)=R(L)",
        "=\\boxed{2997}",
        "Any weak-action upgrade is a",
    )
    for needle in required:
        if needle not in text:
            raise AssertionError(f"report scope marker missing: {needle}")


def build_certificate() -> dict[str, Any]:
    validate_repository_dependencies()
    validate_report_scope()

    if RANK_MAX != 30:
        raise AssertionError("Vinberg rank ceiling changed")
    if BLOWUP_MAX != 2997:
        raise AssertionError("universal rational-surface blow-up bound changed")

    bounds = coefficient_bounds(FIELD_DEGREE_MAX, CONJUGATE_HEIGHT)
    if len(bounds) != FIELD_DEGREE_MAX or any(value <= 0 for value in bounds):
        raise AssertionError("Vieta coefficient box is not finite")
    if bounds[0] != FIELD_DEGREE_MAX * CONJUGATE_HEIGHT:
        raise AssertionError("first Vieta bound changed")

    cycle_slots = simple_directed_cycle_slot_bound(RANK_MAX)
    if cycle_slots <= comb(RANK_MAX, 2):
        raise AssertionError("cycle slot bound lost longer cycles")

    return {
        "schema": (
            "coxeter.unconditional-maximal-arithmetic-reflection-chamber-"
            "finiteness.v1"
        ),
        "terminal_theorem": {
            "minimum_dimension": 2,
            "maximum_dimension": DIMENSION_MAX,
            "maximum_native_rank": RANK_MAX,
            "maximum_ground_field_degree": FIELD_DEGREE_MAX,
            "finite_ground_field_set": True,
            "finite_native_pair_set": True,
            "reflection_maximality_required_for_chamber_count": True,
            "finite_unmarked_normalized_chambers": True,
            "AG_carrier_hypothesis": False,
        },
        "two_terminal_routes": {
            "direct": (
                "ABSW 2008 main theorem (or Nikulin 2007 Theorem 2.2)"
            ),
            "finite_pair_factorization": [
                "Vinberg n<30",
                "Belolipetsky-Linowitz degree<=25",
                "Nikulin bounded connected root frame A=62",
                "field inheritance plus bounded integral cycles",
                "finite algebraic-integer coefficient box",
                "Nikulin fixed-(k,m) reflective-lattice finiteness",
                "reflection maximality W=R(L)",
                "one complete-mirror chamber orbit",
            ],
            "uses_ABSW_finiteness_conclusion_in_factorization": False,
        },
        "finite_field_box": {
            "frame_constant": FRAME_CONSTANT,
            "maximum_cycle_length": RANK_MAX,
            "uniform_conjugate_height": str(CONJUGATE_HEIGHT),
            "maximum_element_degree": FIELD_DEGREE_MAX,
            "simple_directed_cycle_slot_bound": cycle_slots,
            "Vieta_coefficient_bounds": [str(value) for value in bounds],
            "finite": True,
        },
        "weak_rational_surface_corollary": {
            "formula": "4*[k:Q]*(n+1)-3",
            "maximum_blowups": BLOWUP_MAX,
            "requires_reflection_maximality": False,
            "effective_or_nef_return_claimed": False,
            "reflection_action_return_claimed": False,
        },
        "dependency_sha256": DEPENDENCY_SHA256,
        "scope": {
            "published_inputs_imported": True,
            "project_native_reproof_of_transition_enumeration_claimed": False,
            "new_AG_only_proof_claimed": False,
            "unconditional_global_finiteness_claimed": True,
            "external_primary_hashes_are_audit_time_pins": True,
            "external_primary_pdf_bytes_required_for_replay": False,
        },
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    wrapper = {"payload": payload, "payload_sha256": digest(payload)}
    return json.dumps(wrapper, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.write:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(encoded, encoding="utf-8")
    if args.check_output:
        if not ARCHIVE.exists() or ARCHIVE.read_text(encoding="utf-8") != encoded:
            raise SystemExit("frozen archive differs from exact replay")
    if not args.write and not args.check_output:
        print(encoded, end="")


if __name__ == "__main__":
    main()
