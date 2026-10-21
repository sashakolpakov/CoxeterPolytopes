#!/usr/bin/env python3
"""Certify the unconditional NS/Chow encoding plus maximal-chamber theorem."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.md"
ARCHIVE = ROOT / "data/enumerations/unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.json"

DIMENSION_MAX = 29
FIELD_DEGREE_MAX = 25
TRACE_RANK_MAX = FIELD_DEGREE_MAX * (DIMENSION_MAX + 1)
BLOWUP_MAX = 4 * TRACE_RANK_MAX - 3

DEPENDENCY_SHA256 = {
    "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness.md":
        "6eee91512840287ce47a195a317fa4bbfd19678a95b52869ec0a45775d61e67b",
    "scripts/certify_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py":
        "c32c7969f900093db388718f3ac769f669d42e43b6c2a95c76297eb34899f541",
    "data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness.json":
        "2bce48c0a7156caea3f4b4890411d770450912b4d4b2d3827ce93cc7e9730acc",
    "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.md":
        "1f1592ef207580b765d56bdca9461429456ef9fa46eb180b0ae65d69dc68ef3c",
    "scripts/audit_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py":
        "e6fc5f7543093b6167e6d46612b64c7247007cc16b06f081e14260a283193400",
    "data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.json":
        "87ec814ae05d29b04ac068296c4f995144264fed353fa02f607eb22841cff734",
    "results/universal_weak_chow_action_carrier_and_leverage_no_go.md":
        "df38b26e0cffd8a391834be922869153af5badf41a632e74a48b472b332e101a",
    "scripts/certify_universal_weak_chow_action_carrier.py":
        "d3a19cfe7069b811721a77a624952d92ab5bc7c5574d9e332d2ff8c7f0205eba",
    "data/enumerations/universal_weak_chow_action_carrier.json":
        "d75897bbf8d80a662a7b9d3aa3ae729052494e20633d835566db70870bb0fc9e",
    "results/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.md":
        "91345590b506ef7c7d4287be5a4d94c7563ba7613167ffbe34272e6b7fc422c1",
    "scripts/audit_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile.py":
        "98765eaae927d9265cf69591f02c31d1d0adc71653936d2199f53eb75b16a347",
    "data/enumerations/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.json":
        "a9c19a9f2dce89c03182c9dc7ad8248a0ab75e32e5e58cd1bca6fac7683764f4",
    "results/unconditional_ag_absw_universality_leverage_no_go.md":
        "fcc3acb8bee4751bef083cd00a179f9b542602561ba1bf81d89f374dd5523150",
    "results/projected_k3_maximal_abstract_chamber_finiteness.md":
        "59dd458acf14d18021f6c5583096da13f39c83115fc136d6ba91bcc850ce3603",
    "results/projected_k3_universality_rank_saturation_no_go.md":
        "7ed5176c0c0a36e036c0e12dfa6acf1ee15d90d5ec7143ca5df59ebc3b7143d2",
    "results/borcherds_h21_outer_routing_counterexample.md":
        "fa899977f418d6fa93f7d5270e067c28a648c5e55f2636cb2b362e381534821d",
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def validate_dependencies() -> None:
    for relative, expected in DEPENDENCY_SHA256.items():
        path = ROOT / relative
        if not path.is_file():
            raise AssertionError(f"dependency missing: {relative}")
        actual = sha256_file(path)
        if actual != expected:
            raise AssertionError(
                f"dependency hash changed: {relative}: {actual} != {expected}"
            )


def validate_report_scope() -> None:
    text = REPORT.read_text(encoding="utf-8")
    required = (
        "unconditional theorem, with provenance separated",
        "Theorem 1.1 (uniform weak NS/Chow carrier)",
        "standard admissible Vinberg lattice-form model",
        "b\\le2997",
        "projected wall normals",
        "raw ambient trace",
        "Theorem 1.2 (unconditional maximal-chamber finiteness)",
        "only finitely many unmarked normalized chamber",
        "two independent classical endings",
        "does not invoke the ABSW finiteness",
        "finite Vieta boxes",
        "The first does **not** prove the second",
        "not an AG-only proof of finiteness",
        "No non-tautological effective/nef/mixed-action carrier",
    )
    for marker in required:
        if marker not in text:
            raise AssertionError(f"scope marker missing: {marker}")


def build_certificate() -> dict[str, Any]:
    validate_dependencies()
    validate_report_scope()
    if TRACE_RANK_MAX != 750:
        raise AssertionError("trace-rank arithmetic changed")
    if BLOWUP_MAX != 2997:
        raise AssertionError("blow-up arithmetic changed")

    return {
        "schema": "coxeter.unconditional-ns-chow-and-maximal-chamber.v1",
        "unconditional_weak_carrier": {
            "input": (
                "every finite-sided finite-volume arithmetic Coxeter chamber "
                "in the standard admissible Vinberg lattice-form model"
            ),
            "requires_reflection_maximality": False,
            "maximum_dimension_imported": DIMENSION_MAX,
            "maximum_field_degree_imported": FIELD_DEGREE_MAX,
            "maximum_trace_rank": TRACE_RANK_MAX,
            "maximum_blowups": BLOWUP_MAX,
            "full_trace_NS_injection": True,
            "faithful_rational_Chow_reflection_action": True,
            "exact_semialgebraic_eigenslice_return": True,
            "effective_or_nef_return": False,
        },
        "unconditional_maximal_chambers": {
            "input": "reflection-maximal finite-volume arithmetic Coxeter chambers",
            "finite_native_pair_set": True,
            "finite_native_similarity_hosts": True,
            "finite_unmarked_normalized_chambers": True,
            "direct_published_routes": ["ABSW 2008", "Nikulin 2007 Theorem 2.2"],
            "ABSW_free_factorization": True,
            "AG_carrier_hypothesis": False,
        },
        "provenance_boundary": {
            "finiteness_derived_from_weak_carrier": False,
            "classical_dimension_and_field_inputs_imported": True,
            "global_AG_only_finiteness_claimed": False,
            "weak_carrier_universal_because_numerical": True,
            "universal_effective_nef_or_mixed_action_compiler_proved": False,
            "current_complex_high_leverage_stack_universal": False,
        },
        "dependency_sha256": DEPENDENCY_SHA256,
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


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
        if not ARCHIVE.is_file() or ARCHIVE.read_text(encoding="utf-8") != encoded:
            raise SystemExit("frozen archive differs from exact replay")
    if not args.write and not args.check_output:
        print(encoded, end="")


if __name__ == "__main__":
    main()
