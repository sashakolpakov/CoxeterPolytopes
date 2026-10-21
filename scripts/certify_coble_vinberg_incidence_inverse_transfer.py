#!/usr/bin/env python3
"""Freeze the two high-dimensional calibrations of the Coble--Vinberg
incidence inverse-transfer theorem.

The general theorem is mathematical: it applies a finite all-contraction
incidence compiler to any raw datum in ``D_CobVin``.  This executable does
not pretend to enumerate that a-priori class.  It independently reloads the
four exact archives which certify that the H14/60 and H17/48 raw inputs lie
in the denominator, that their incidence branches are nonempty precisely in
characteristic three, and that their regenerated chambers are literally the
same finite-volume chambers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/coble_vinberg_incidence_inverse_transfer.json"

INPUTS = {
    "type_vi_raw": ROOT
    / "data/enumerations/polytope_first_coble_type_vi_h14.json",
    "blowdown_bridge": ROOT
    / "data/enumerations/coble_polytope_first_blowdown_incidence_bridge.json",
    "type_vi_full": ROOT
    / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json",
    "type_mii_full": ROOT
    / "data/enumerations/coble_type_mii_full_ns_h17.json",
    "type_mii_hostile": ROOT
    / "data/enumerations/coble_type_mii_h17_hostile_audit.json",
}

EXPECTED_SHA256 = {
    "type_vi_raw": "b1dfafffd87ff281d73c3517f18033ca7d89852e9a3e0dd39a035f890ddf26fb",
    "blowdown_bridge": "08d40b33637fe4d332dc2e62a368bee401c0eda385baa3f823dde3285a15dfdd",
    "type_vi_full": "8247ad85f7970d3bf1f08fe15517b70376409c17cc2784b8f14b129aaecc2ede",
    "type_mii_full": "6ed6fba936c362e02793375fa3adfea0279d312c853a8bb4bab30c23b4b8c79c",
    "type_mii_hostile": "124a299740c6a82c7ab543ce9236c1ad38d6e988e637a0d173cbe2203756c25f",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(name: str) -> dict[str, Any]:
    path = INPUTS[name]
    digest = _sha256(path)
    expected = EXPECTED_SHA256[name]
    if digest != expected:
        raise ArithmeticError(
            f"pinned input changed for {name}: expected {expected}, got {digest}"
        )
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("all_checks_pass") is not True:
        raise ArithmeticError(f"upstream certificate is not green: {path}")
    return payload


def build_certificate() -> dict[str, Any]:
    vi_raw = _load("type_vi_raw")
    bridge = _load("blowdown_bridge")
    vi = _load("type_vi_full")
    mii = _load("type_mii_full")
    mii_hostile = _load("type_mii_hostile")

    vi_taxonomy = vi_raw["wall_taxonomy"]["counts"]
    vi_cusps = vi_raw["cusp_lift_exhaustion"]
    vi_return = vi_raw["exact_same_P_return"]
    vi_cm_cone = vi["Coble_Mukai_chamber"]["exact_cone"]
    vi_vinberg = vi["literal_vinberg"]

    mii_curves = mii["curve_exhaustion"]
    mii_cm_cone = mii["Coble_Mukai_chamber"]["exact_double_description"]
    mii_ag = mii["AG_native_stratum"]
    mii_volume = mii["finite_volume"]

    vi_bridge = bridge["type_VI_plane_compiler"]
    mii_bridge = bridge["type_MII_raw_compiler"]

    checks = {
        "all_five_upstream_archives_green": all(
            payload["all_checks_pass"]
            for payload in (vi_raw, bridge, vi, mii, mii_hostile)
        ),
        "surface_free_denominator_has_two_distinct_positive_ranks": (
            vi_return["raw_mirror_count"] == 60
            and mii_curves["complete_negative_curve_count"] == 48
            and vi_raw["raw_arithmetic_lattice"]["rank"] == 15
            and len(mii["full_NS"]["basis"]) == 18
        ),
        "type_vi_boundary_and_root_stock": (
            vi_taxonomy
            == {
                "boundary_minus4": 5,
                "crossing_minus1_curve": 10,
                "cusp_lift_minus1_curve": 35,
                "minus2_curve": 10,
            }
            and vi["curve_exhaustion"]["negative_curve_count"] == 60
        ),
        "type_vi_CM_is_causal_and_vinberg_complete": (
            vi_cm_cone["hyperideal_vertex_count"] == 0
            and vi_cm_cone["ideal_vertex_count"] == 57
            and vi_vinberg["bad_edge_count"] == 0
            and vi_vinberg["bad_parabolic_neighborhood_count"] == 0
        ),
        "type_vi_all_cusp_lifts_regenerate_raw_walls": (
            vi_cusps["primitive_CM_cusps"] == 57
            and vi_cusps["integral_cusp_boundary_pairs"] == 35
            and vi_cusps[
                "all_and_only_integral_lifts_are_raw_minus_one_walls"
            ]
            and vi_return["oriented_mirror_sets_equal"]
        ),
        "type_vi_incidence_branch_exists_only_in_char3": (
            vi_bridge["boundary_preserving_P2_blowdowns"]["count"] == 5
            and vi_bridge["residual_root_incidence"][
                "allowed_characteristics_not_two"
            ]
            == [3]
        ),
        "type_mii_boundary_and_root_stock": (
            mii_curves["boundary_curves"] == 8
            and mii_curves["minus_one_roots_meeting_two_boundaries"] == 16
            and mii_curves["actual_minus_two_curves"] == 24
            and mii_curves["complete_negative_curve_count"] == 48
        ),
        "type_mii_CM_is_causal_and_cusp_complete": (
            mii_cm_cone["hyperideal_vertex_count"] == 0
            and mii_cm_cone["ideal_vertex_count"] == 202
            and mii_curves["total_single_boundary_lift_tests"] == 1616
            and mii_curves["integral_single_boundary_lifts"] == 0
        ),
        "type_mii_all_branch_grid_incidence_only_in_char3": (
            mii_bridge["boundary_frame_search"][
                "q_over_two_integral_on_P_generated_module"
            ]
            == 1
            and mii_bridge["P_compiled_blowdown"][
                "all_24_permutations_occur"
            ]
            and mii_bridge["incidence_nonemptiness_iff"][
                "allowed_characteristics_not_two"
            ]
            == [3]
        ),
        "both_returns_are_finite_volume_noncompact": (
            vi_return["finite_volume_noncompact"]
            and mii_volume["finite_volume"]
            and not mii_volume["compact"]
            and vi_return["full_H14_cone"]["hyperideal_vertex_count"] == 0
            and mii_hostile["literal_full_cone"]["finite_vertex_count"]
            == 190808
            and mii_hostile["literal_full_cone"]["ideal_vertex_count"] == 5244
            and mii_hostile["literal_full_cone"]["hyperideal_vertex_count"] == 0
        ),
        "type_mii_canonical_cover_is_supersingular_sigma_one": (
            mii_ag["covering_K3"]
            == "supersingular Kum(E x E), Artin invariant 1"
            and mii_ag["characteristic"] == 3
        ),
    }
    if not all(checks.values()):
        failed = [name for name, value in checks.items() if not value]
        raise ArithmeticError(f"Coble--Vinberg transfer checks failed: {failed}")

    return {
        "schema": "coble-vinberg-incidence-inverse-transfer-v1",
        "status_date": "2026-08-09",
        "direction": "POLYTOPE_TO_SURFACE_TO_EXACT_SAME_POLYTOPE",
        "surface_free_denominator": {
            "name": "D_CobVin",
            "definition": [
                "an oriented finite-volume integral Coxeter wall datum with a strict integral interior vector h",
                "all numerical Coble boundary frames and odd-unimodular de-scalings are exhaustively enumerated",
                "for each branch K=-sum(B_i)/2 is characteristic, rank(N)=10+b, and every wall has the boundary/adjunction alphabet",
                "the numerical CM roots reconstructed from (-2) candidates and two-boundary (-1) candidates pass a literal Vinberg finite-index certificate",
                "the remaining numerical walls are exactly all integral primitive CM-cusp/single-boundary lifts",
            ],
            "not_a_hypothesis": [
                "a surface or an effective curve",
                "a Kondo type or one of the nine classified diagrams",
                "a characteristic",
                "a K3 lattice, crystalline period, or involution",
                "equality with a nef cone",
                "the bound b in {1,2,5,8}",
            ],
            "known_distinct_positive_inputs": 2,
        },
        "finite_all_branch_compiler": {
            "name": "C_CobVin(P,h)",
            "branches": (
                "For every numerical boundary frame, enumerate every sequence of input "
                "classes which becomes a smooth (-1)-class under the exact blowdown "
                "update, ending at P2 or a Hirzebruch surface. Rank drops at every step."
            ),
            "incidence_space": (
                "Reverse each sequence as a universal iterated point blow-up; relative "
                "Hilbert schemes parameterize one smooth geometrically irreducible "
                "curve in every numerical wall class. Impose the marked incidences, "
                "h ample, disjoint reduced boundaries summing to -2K, h0(-K)=0, "
                "and h0(-2K)=1. This is a finite union of finite-type locally closed "
                "stacks over Spec Z[1/2]."
            ),
            "termination": (
                "There are finitely many frames, finitely many input classes, and at "
                "most rank(N) contractions per branch. Hilbert/flattening and "
                "saturated-ideal algorithms decide geometric nonemptiness and its "
                "constructible characteristic image."
            ),
            "necessity": (
                "On an exact carrier, the strict transform on the original surface "
                "of every curve in a minimal-model contraction has negative square; "
                "strong-Coble negativity and exact nefness put it among the input "
                "facets, so one enumerated branch contains every exact carrier."
            ),
        },
        "promotion_and_return_theorem": {
            "intrinsic_iff": (
                "For every (P,h) in D_CobVin and every algebraically closed field of "
                "characteristic not two, P has an exact marked strong-Coble full-nef "
                "carrier taking h to an ample class if and only if C_CobVin(P,h) has "
                "a geometric point over that field."
            ),
            "root_exhaustion": (
                "A compiler point makes every numerical CM root an effective "
                "irreducible root: a listed (-2) curve is irreducible, while the "
                "negative-definite support of 2E+B_i/2+B_j/2 forces a listed "
                "two-boundary (-1) root to have no proper effective root below it. "
                "Its raw Vinberg certificate makes W(Gamma) finite "
                "index in O(CM); hence W(S) is finite index, Aut(S) is finite, and "
                "Kondo's Proposition Vinberg2 says Gamma is the complete effective "
                "irreducible root set."
            ),
            "minus_one_exhaustion": (
                "Adjunction gives total boundary incidence two. A hidden (1,1) "
                "(-1)-curve gives a hidden effective irreducible CM root. A hidden "
                "(2) curve gives the primitive nef isotropic class f=E+B_i/2, hence "
                "an ideal extreme ray. Primitivity follows from the Halphen index-at-"
                "most-two theorem and the primitive fiber multiplicity vector (2,1) "
                "on 2E+B_i. The denominator's complete integral cusp-lift "
                "regeneration already lists it."
            ),
            "same_P": (
                "The unique anti-bicanonical boundary, complete root stock, and all "
                "single-boundary cusp lifts exhaust every negative curve. Therefore "
                "the marked nef cone is literally the oriented input chamber P."
            ),
            "canonical_K3_corollary": (
                "The double cover branched along the reconstructed reduced -2K "
                "divisor is a smooth K3 with its deck involution. In characteristic "
                "three its crystalline period and period-preserving involution exist "
                "geometrically; no abstract crystalline Torelli guess is required."
            ),
        },
        "classification_pushback": {
            "published_input": (
                "Kondo's exhaustive classification of finite-automorphism strong "
                "Coble surfaces in characteristic not two"
            ),
            "theorem_on_every_compiler_positive_input_not_every_denominator_input": {
                "boundary_counts": [1, 2, 5, 8],
                "hyperbolic_dimensions": [10, 11, 14, 17],
                "no_members_in_dimensions": "12, 13, 15, 16, or above 17",
                "H14_and_H17_characteristic": 3,
            },
            "positive_high_dimensional_calibrations": [
                {
                    "raw_id": "H14/60",
                    "boundary_count": 5,
                    "CM_root_count": 20,
                    "CM_cusp_count": 57,
                    "single_boundary_lift_count": 35,
                    "incidence_characteristics_not_two": [3],
                    "exact_return": "27845 finite + 1027 ideal + 0 hyperideal rays",
                },
                {
                    "raw_id": "H17/48",
                    "boundary_count": 8,
                    "CM_root_count": 40,
                    "CM_cusp_count": 202,
                    "single_boundary_lift_count": 0,
                    "incidence_characteristics_not_two": [3],
                    "exact_return": "190808 finite + 5244 ideal + 0 hyperideal rays",
                },
            ],
        },
        "five_part_score": {
            "C_surface_free_class": "pass",
            "Car_finite_exhaustive_carrier_enumeration": "pass",
            "Car_realized_branch_for_every_P_in_C": (
                "fail: geometric nonemptiness is not forced on all of D_CobVin"
            ),
            "Car_overall_under_strict_contract": "fail",
            "Return_exact_same_P": "pass for every geometric compiler branch",
            "AG_classification_applies": "pass after a geometric branch exists",
            "Pushback_new_theorem_on_all_inputs": (
                "fail: only the compiler-positive/carrier locus is classified; raw "
                "D_CobVin conditions do not yet force incidence-stack nonemptiness"
            ),
            "rating": (
                "carrier-relative L3 exact-return criterion; strict programme L2-L3, "
                "not L4, because Car/essential-surjectivity fails only at class-wide "
                "incidence nonemptiness"
            ),
        },
        "scope": {
            "new_project_theorem": (
                "the all-branch incidence iff, uniform Kondo--Vinberg/cusp hidden-wall "
                "closure, and exact same-P return on every compiler-positive branch"
            ),
            "exact_L4_missing_gate": (
                "Find a surface-free P-only subdenominator containing more than a "
                "prelisted answer set on which the raw numerical conditions force at "
                "least one incidence stack to be nonempty. Without that theorem, "
                "Kondo pushback classifies only R_CobVin={P:C_CobVin(P) nonempty}."
            ),
            "compactness": (
                "Every geometric branch is noncompact: pull back a ruling from a "
                "minimal Hirzebruch model (or from the penultimate F1 before a P2 "
                "model) to obtain a primitive nef isotropic ray. Finite volume is a "
                "separate raw Vinberg/cone hypothesis and is never inferred from "
                "finite sidedness or from the existence of a cusp."
            ),
            "published_ingredients": [
                "Vinberg's finite-index criterion",
                "Dolgachev--Kondo finite-index implies finite Aut",
                "Kondo Proposition Vinberg2",
                "Kondo's primitive-isotropic/fibration and Halphen results",
                "Kondo's nine-class theorem",
            ],
            "not_claimed": [
                "a new abstract H14 or H17 Coxeter polytope",
                "a compact-polytope theorem",
                "a theorem for raw wall systems outside D_CobVin",
                "a standalone lattice-only crystalline Torelli compiler",
                "a generic executable enumerator for all unsupplied wall matrices",
            ],
        },
        "pinned_inputs": {
            name: {
                "path": str(path.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256[name],
            }
            for name, path in INPUTS.items()
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(
            encoding="utf-8"
        ) != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "all_checks_pass": result["all_checks_pass"],
                "denominator": result["surface_free_denominator"]["name"],
                "rating": result["five_part_score"]["rating"],
                "positive_inputs": [
                    row["raw_id"]
                    for row in result["classification_pushback"][
                        "positive_high_dimensional_calibrations"
                    ]
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
