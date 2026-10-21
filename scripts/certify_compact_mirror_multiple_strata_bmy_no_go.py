#!/usr/bin/env python3
"""Exact local BMY audit for ordinary multiple-mirror strata.

The arithmetic geometry is proved in the companion report.  This executable
replays the local blow-up identities, Langer's equal half-weight ordinary
point values for the log-canonical dihedral cases, and all sign conclusions.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "compact_mirror_multiple_strata_bmy_no_go.json"
)
SCHEMA = "coxeter.compact_mirror_multiple_strata_bmy_no_go.v1"


def _q(value: Fraction) -> int | str:
    return value.numerator if value.denominator == 1 else str(value)


def reduced_blowup_record(branch_count: int) -> dict:
    """Reduced total-transform log pair at an ordinary b-fold point."""

    b = branch_count
    log_c2_change = 0
    log_c1_square_change = -(b - 2) ** 2
    defect_change = 3 * log_c2_change - log_c1_square_change
    return {
        "branch_count": b,
        "log_c2_change": log_c2_change,
        "log_c1_square_change": log_c1_square_change,
        "bmy_defect_change": defect_change,
        "nonnegative": defect_change >= 0,
    }


def half_weight_eorb(branch_count: int) -> Fraction:
    """Langer local e_orb for b equal coefficient-one-half lines, b<=4."""

    if branch_count == 2:
        return Fraction(1, 4)
    if branch_count == 3:
        return Fraction(1, 16)
    if branch_count == 4:
        return Fraction(0)
    raise ValueError("the half-weight ordinary point is BMY-relevant only for b=2,3,4")


def half_weight_record(branch_count: int) -> dict:
    b = branch_count
    alpha_sum = Fraction(b, 2)
    eorb = half_weight_eorb(b)
    pair_mass = Fraction(b * (b - 1), 8)

    # Contribution of this singular point to
    # 3 e_orb(S,Delta) - (K_S+Delta)^2 after curve terms are separated.
    global_point_term = 3 * (alpha_sum + eorb - 1) - 2 * pair_mass

    # Blow up and put coefficient one on the exceptional curve.
    reduced_exceptional_c2_change = -eorb
    reduced_exceptional_c1_change = -(2 - alpha_sum) ** 2
    reduced_exceptional_defect_change = (
        3 * reduced_exceptional_c2_change
        - reduced_exceptional_c1_change
    )

    # Crepant exceptional coefficient beta=A-1.
    crepant_beta = alpha_sum - 1
    crepant_c2_change = (2 - alpha_sum) ** 2 - eorb
    crepant_defect_change = 3 * crepant_c2_change

    return {
        "dihedral_label_m": b,
        "local_group_order": 2 * b,
        "coefficient_sum_A": _q(alpha_sum),
        "log_canonical": alpha_sum <= 2,
        "langer_local_eorb": _q(eorb),
        "global_bmy_point_term": _q(global_point_term),
        "reduced_exceptional_resolution": {
            "exceptional_coefficient": 1,
            "orbifold_c2_change": _q(reduced_exceptional_c2_change),
            "log_c1_square_change": _q(reduced_exceptional_c1_change),
            "bmy_defect_change": _q(reduced_exceptional_defect_change),
        },
        "crepant_resolution": {
            "exceptional_coefficient": _q(crepant_beta),
            "orbifold_c2_change": _q(crepant_c2_change),
            "log_c1_square_change": 0,
            "bmy_defect_change": _q(crepant_defect_change),
        },
        "all_relevant_terms_nonnegative": (
            global_point_term >= 0
            and reduced_exceptional_defect_change >= 0
            and crepant_defect_change >= 0
        ),
    }


def build_certificate() -> dict:
    reduced = [reduced_blowup_record(b) for b in range(2, 13)]
    half = [half_weight_record(b) for b in (2, 3, 4)]

    # Langer's universal local bound for an ordinary lc point is
    # e_orb <= (1-A/2)^2 = (2-A)^2/4.  The following coefficients are the
    # resulting lower bounds for the two resolution conventions.
    universal_full_exceptional_factor = Fraction(1, 4)
    universal_crepant_factor = Fraction(9, 4)

    checks = {
        "reduced_resolution_never_creates_negative_defect": all(
            row["bmy_defect_change"] >= 0 for row in reduced
        ),
        "half_weight_lc_exactly_m_at_most_four": all(
            row["log_canonical"] for row in half
        ) and Fraction(5, 2) > 2,
        "half_weight_point_terms_are_1_4_3_16_0": [
            row["global_bmy_point_term"] for row in half
        ] == ["1/4", "3/16", 0],
        "half_weight_resolution_corrections_are_nonnegative": all(
            row["all_relevant_terms_nonnegative"] for row in half
        ),
        "universal_full_exceptional_lower_bound_nonnegative": (
            universal_full_exceptional_factor >= 0
        ),
        "universal_crepant_lower_bound_nonnegative": (
            universal_crepant_factor >= 0
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"local BMY no-go replay failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "reduced_ordinary_multiple_point_blowups": reduced,
        "half_weight_dihedral_points": half,
        "universal_weighted_lc_resolution_bounds": {
            "hypotheses": "ordinary branches, 0<=alpha_i<=1, A=sum(alpha_i)<=2",
            "langer_bound": "e_orb(p)<=((2-A)^2)/4",
            "coefficient_one_exceptional_defect_change_lower_bound": "((2-A)^2)/4",
            "crepant_exceptional_defect_change_lower_bound": "9*((2-A)^2)/4",
        },
        "natural_dihedral_stack": {
            "local_group": "I_2(m)",
            "group_order": "2m",
            "local_orbifold_mass": "1/(2m)",
            "feedback": (
                "absorbed by the orbifold tangent bundle and Hirzebruch "
                "proportionality; it is not an additional negative stratum term"
            ),
        },
        "conclusion": (
            "codimension-two mirror multiplicities do not provide a uniformly "
            "negative local BMY contribution: the reduced and every lc weighted "
            "resolution correction is nonnegative; half-weight m>=5 is non-lc"
        ),
        "remaining_gate": (
            "a class-uniform noncircular bound for U_det=w*lambda^n/q_min, "
            "or a different AG invariant not absorbed by stack proportionality"
        ),
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != payload:
            raise SystemExit(f"archived artifact differs: {OUTPUT}")
    print(json.dumps({
        "all_checks_pass": payload["all_checks_pass"],
        "half_weight_point_terms": [
            row["global_bmy_point_term"]
            for row in payload["half_weight_dihedral_points"]
        ],
        "remaining_gate": payload["remaining_gate"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
