#!/usr/bin/env python3
"""Exact dependency and numerical replay for the weak-log transfer theorem."""

from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

try:
    from scripts.certify_pfirst_weak_log_semiampleness import (
        build_certificate as build_semiampleness_certificate,
        verify_certificate as verify_semiampleness_certificate,
    )
    from scripts.certify_reflection_maximal_weak_log_exact_nef_return import (
        payload as build_weak_log_calibration,
    )
except ModuleNotFoundError:  # Direct ``python scripts/this_file.py`` replay.
    from certify_pfirst_weak_log_semiampleness import (
        build_certificate as build_semiampleness_certificate,
        verify_certificate as verify_semiampleness_certificate,
    )
    from certify_reflection_maximal_weak_log_exact_nef_return import (
        payload as build_weak_log_calibration,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE = (
    ROOT / "data/enumerations/pfirst_weak_log_complete_carrier_transfer.json"
)
SCHEMA = "coxeter.pfirst.weak_log_complete_carrier_transfer.v1"

SOURCE_PATHS = {
    "composition_theorem": "results/pfirst_weak_log_complete_carrier_transfer.md",
    "semiampleness_theorem": "results/pfirst_weak_log_semiampleness_certificate.md",
    "weak_log_exact_return": "results/reflection_maximal_weak_log_exact_nef_return.md",
    "weak_log_exact_return_hostile": "results/reflection_maximal_weak_log_exact_nef_return_hostile_audit.md",
    "orbifold_rank_theorem": "results/rational_wrong_log_orbifold_rank_gate.md",
    "orbifold_rank_hostile": "results/rational_wrong_log_orbifold_rank_gate_hostile_audit.md",
    "maximal_finiteness_theorem": "results/rational_wrong_log_reflection_maximal_finiteness.md",
    "maximal_finiteness_hostile": "results/rational_wrong_log_maximal_finiteness_hostile_audit.md",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_digest(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def source_pins() -> dict[str, dict[str, str]]:
    return {
        key: {"path": path, "sha256": sha256(ROOT / path)}
        for key, path in SOURCE_PATHS.items()
    }


def build_certificate() -> dict[str, Any]:
    semiample = build_semiampleness_certificate()
    assert verify_semiampleness_certificate(semiample)
    weak_log = build_weak_log_calibration()

    semi_f4 = semiample["calibration"]["class_calculation"]
    weak_f4 = weak_log
    assert weak_f4["surface"].startswith("Bl_p(F_4)")
    semi_d = [Fraction(str(value)) for value in semi_f4["D"]]
    weak_d = [Fraction(value[0], value[1]) for value in weak_f4["weak_log_class"]]
    assert semi_d == weak_d
    assert semi_f4["mD"] == [int(2 * value) for value in weak_d]
    assert semi_f4["D_square"] == Fraction(*weak_f4["weak_log_square"])
    f4_exact_return = (
        weak_f4["area_over_pi"] == [1, 6]
        and weak_f4["proper_overgroup_area_upper_over_pi"] == [1, 12]
        and weak_f4["minimum_noncompact_coxeter_area_over_pi"] == [1, 6]
        and "return the complete nef chamber" in weak_f4["conclusion"]
    )

    raw: dict[str, Any] = {
        "schema": SCHEMA,
        "status": "proved carrier-relative composition; raw essential surjectivity open",
        "sources": source_pins(),
        "certificate_gates": {
            "smooth_complex_rational_surface": True,
            "disjoint_smooth_rational_minus4_boundary": True,
            "exact_effective_irreducible_oriented_facets": True,
            "finite_sided_finite_volume_chamber": True,
            "adapted_host_integrality": True,
            "relative_or_stronger_reflection_maximality": True,
            "complete_affine_cover_and_line_bundle_transitions": True,
            "compatible_global_sections": True,
            "localized_unit_identities_on_every_chart": True,
        },
        "logical_chain": [
            "localized unit identities imply O(mD) is globally generated",
            "global generation implies D is nef",
            "nef weak-log class classifies every hidden negative curve",
            "adapted-host maximality implies W(P)=R(M)",
            "ample orientation implies P=P(Nef(S)) including ideal closure",
            "boundary contraction and orbifold c2 imply rho<=26",
            "rational full-nef facets imply minimal Vinberg field Q",
            "finite adapted hosts plus one complete mirror orbit imply finitely many chambers",
        ],
        "rank_elimination": {
            "nef_square_inequality": "10-rho+b>=0",
            "orbifold_c2_inequality": "4*rho+8-7*b>=0",
            "nonnegative_linear_combination": (
                "7*(10-rho+b)+(4*rho+8-7*b)=78-3*rho>=0"
            ),
            "maximum_picard_rank": 26,
            "maximum_hyperbolic_dimension": 25,
            "empty_boundary_maximum_picard_rank": 10,
            "empty_boundary_maximum_hyperbolic_dimension": 9,
        },
        "field_and_finiteness": {
            "minimal_vinberg_field": "Q",
            "adapted_host_interval": "L(2) subset M subset (1/2)L",
            "finitely_many_hosts_in_rank_at_most_26": True,
            "relative_maximality_gives_full_reflection_group": True,
            "complete_mirror_arrangement_has_one_chamber_orbit": True,
        },
        "f4_composed_calibration": {
            "surface": weak_f4["surface"],
            "weak_log_class": weak_f4["weak_log_class"],
            "weak_log_square": weak_f4["weak_log_square"],
            "twice_weak_log_class": semi_f4["mD"],
            "semiampleness_provider": semiample["calibration"]["name"],
            "affine_chart_count": len(
                semiample["calibration"]["toric_divisor_2D"]["local_generation"]
            ),
            "global_section_count": len(
                semiample["calibration"]["toric_divisor_2D"][
                    "returned_global_sections"
                ]
            ),
            "coxeter_orders": weak_f4["coxeter_orders"],
            "area_over_pi": weak_f4["area_over_pi"],
            "exact_full_nef_return": f4_exact_return,
        },
        "scope": {
            "fixed_certificate_verifier_and_transfer": True,
            "composition_producer_is_dependency_and_calibration_replay": True,
            "composition_producer_accepts_arbitrary_carrier_packet": False,
            "complete_Mori_generator_list_assumed": False,
            "nefness_inferred_from_proposed_walls": False,
            "arbitrary_PB_Atlas_adapter_implemented": False,
            "raw_surface_free_essential_surjectivity_proved": False,
            "nonrational_eigen_nef_branch_exhausted": False,
            "global_ABSW_replacement_proved": False,
        },
    }
    raw["content_sha256_without_this_field"] = canonical_digest(raw)
    return raw


def verify_certificate(payload: Mapping[str, Any]) -> bool:
    if payload.get("schema") != SCHEMA:
        raise ValueError("wrong transfer-certificate schema")
    without_digest = deepcopy(dict(payload))
    archived_digest = without_digest.pop("content_sha256_without_this_field", None)
    if archived_digest != canonical_digest(without_digest):
        raise ArithmeticError("transfer-certificate content digest mismatch")

    for key, expected_path in SOURCE_PATHS.items():
        record = payload["sources"][key]
        if record["path"] != expected_path or record["sha256"] != sha256(
            ROOT / expected_path
        ):
            raise ArithmeticError(f"stale source pin: {key}")

    gates = payload["certificate_gates"]
    if not all(gates.values()):
        raise ValueError("a complete carrier certificate gate is false")

    rank = payload["rank_elimination"]
    if rank["maximum_picard_rank"] != 26 or rank[
        "maximum_hyperbolic_dimension"
    ] != 25:
        raise ArithmeticError("wrong weak-log rank ceiling")
    # The exact nonnegative combination used in the proof is
    # 7*(10-rho+b)+(4*rho+8-7*b)=78-3*rho.
    for rho in range(0, 80):
        for b in range(0, 80):
            left_a = 10 - rho + b
            left_b = 4 * rho + 8 - 7 * b
            if 7 * left_a + left_b != 78 - 3 * rho:
                raise ArithmeticError("rank-elimination identity failed")
            if left_a >= 0 and left_b >= 0 and rho > 26:
                raise ArithmeticError("rank ceiling has a counterexample")
    if rank["empty_boundary_maximum_picard_rank"] != 10:
        raise ArithmeticError("wrong empty-boundary rank ceiling")

    scope = payload["scope"]
    forbidden_promotions = (
        "composition_producer_accepts_arbitrary_carrier_packet",
        "arbitrary_PB_Atlas_adapter_implemented",
        "raw_surface_free_essential_surjectivity_proved",
        "nonrational_eigen_nef_branch_exhausted",
        "global_ABSW_replacement_proved",
    )
    if any(scope[key] for key in forbidden_promotions):
        raise ValueError("carrier-relative theorem was promoted beyond its scope")
    if not scope["composition_producer_is_dependency_and_calibration_replay"]:
        raise ValueError("composition producer role is misstated")
    if scope["complete_Mori_generator_list_assumed"]:
        raise ValueError("composition theorem must not assume a complete Mori list")
    if scope["nefness_inferred_from_proposed_walls"]:
        raise ValueError("circular nefness inference is forbidden")

    semiample = build_semiampleness_certificate()
    if not verify_semiampleness_certificate(semiample):
        raise ArithmeticError("upstream semiampleness certificate failed")
    weak_log = build_weak_log_calibration()
    calibration = payload["f4_composed_calibration"]
    if calibration["weak_log_class"] != weak_log["weak_log_class"]:
        raise ArithmeticError("F4 weak-log class mismatch")
    if calibration["weak_log_square"] != [8, 1]:
        raise ArithmeticError("F4 weak-log square mismatch")
    if calibration["affine_chart_count"] != 5 or calibration[
        "global_section_count"
    ] != 4:
        raise ArithmeticError("F4 cover/section count mismatch")
    if not calibration["exact_full_nef_return"]:
        raise ValueError("F4 exact-return flag is false")
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    payload = build_certificate()
    assert verify_certificate(payload)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.check.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive differs from exact replay")
    elif args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
