#!/usr/bin/env python3
"""Freeze the compact determinant-line L4 scope obstruction.

This certificate contains no experimental search.  It records the exact
restriction-of-scalars dimension calculation which prevents the classical
Ma/Wang reflective-divisor theorems from receiving a compact arithmetic
facet reflection when the ground field has degree greater than one.

The sign convention is the one used by the project: the determinant
extension has signatures ``(2,n)`` at one real place and ``(0,n+2)`` at all
other real places.  Its trace quadratic space therefore has signature
``(2, d*(n+2)-2)`` over Q.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/automorphic_pushback_l4_hostile_audit.json"
)
SCHEMA = "coxeter.automorphic_pushback_l4_hostile_audit.v1"


def _canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def trace_shadow_row(n: int, degree: int) -> dict[str, Any]:
    if n < 1:
        raise ValueError("hyperbolic dimension must be positive")
    if degree < 1:
        raise ValueError("field degree must be positive")

    field_rank = n + 2
    trace_rank = degree * field_rank
    trace_negative_index = trace_rank - 2
    domain_codimension = trace_negative_index - n
    group_dimension_res = degree * field_rank * (field_rank - 1) // 2
    group_dimension_full_trace = trace_rank * (trace_rank - 1) // 2
    group_dimension_difference = (
        group_dimension_full_trace - group_dimension_res
    )
    expected_difference = degree * (degree - 1) * field_rank**2 // 2

    # An F-reflection has a one-dimensional (-1)-eigenspace over F.
    # After restriction of scalars that eigenspace has Q-dimension d.
    minus_eigenspace_q_dimension = degree
    trace_fixed_locus_codimension = degree
    facet_image_codimension_lower_bound = domain_codimension + 1

    return {
        "hyperbolic_dimension_n": n,
        "field_degree_d": degree,
        "determinant_extension_F_rank": field_rank,
        "trace_Q_rank": trace_rank,
        "trace_signature": [2, trace_negative_index],
        "classical_type_IV_dimension_N": trace_negative_index,
        "Hilbert_type_IV_dimension": n,
        "Hilbert_domain_codimension_in_trace_domain": domain_codimension,
        "facet_image_codimension_lower_bound": (
            facet_image_codimension_lower_bound
        ),
        "Res_SO_Q_dimension": group_dimension_res,
        "full_trace_SO_Q_dimension": group_dimension_full_trace,
        "group_dimension_difference": group_dimension_difference,
        "group_dimension_difference_formula": expected_difference,
        "facet_reflection_minus_eigenspace_Q_dimension": (
            minus_eigenspace_q_dimension
        ),
        "facet_reflection_is_Q_pseudoreflection": degree == 1,
        "trace_fixed_locus_codimension": trace_fixed_locus_codimension,
        "trace_fixed_locus_is_divisor": degree == 1,
        "restriction_group_is_proper": degree > 1,
    }


def _counterfactual_wang_status(trace_dimension: int) -> str:
    """Apply only Wang's *classical* high-dimensional exclusion numerically.

    This is deliberately labelled counterfactual: a compact F-reflection is
    not a Q reflection on the trace lattice, and a Hilbert form need not
    extend to the full trace orthogonal group.
    """

    if trace_dimension == 21:
        return "classically_excluded"
    if trace_dimension >= 23 and trace_dimension != 26:
        return "classically_excluded"
    if trace_dimension == 26:
        return "classical_exception_dimension_only"
    return "not_excluded_by_Wang_high_dimension_theorem"


def analyze() -> dict[str, Any]:
    sample_pairs = [
        (4, 2),
        (8, 2),
        (9, 2),
        (10, 2),
        (11, 2),
        (12, 2),
        (13, 2),
        (4, 3),
        (6, 3),
        (7, 3),
        (4, 4),
        (5, 4),
        (6, 4),
        (4, 5),
        (29, 2),
    ]
    rows = []
    for n, degree in sample_pairs:
        row = trace_shadow_row(n, degree)
        row["counterfactual_Wang_trace_status"] = _counterfactual_wang_status(
            row["classical_type_IV_dimension_N"]
        )
        rows.append(row)

    compact_range_checks = {
        "Meyer_forces_degree_gt_one_for_compact_n_ge_4": True,
        "all_sample_compact_rows_have_proper_restriction_group": all(
            row["restriction_group_is_proper"] for row in rows
        ),
        "all_sample_facet_reflections_fail_Q_pseudoreflection": all(
            not row["facet_reflection_is_Q_pseudoreflection"] for row in rows
        ),
        "all_sample_trace_fixed_loci_have_codimension_at_least_two": all(
            row["trace_fixed_locus_codimension"] >= 2 for row in rows
        ),
        "all_group_dimension_formulas_match": all(
            row["group_dimension_difference"]
            == row["group_dimension_difference_formula"]
            for row in rows
        ),
    }

    payload: dict[str, Any] = {
        "schema_version": SCHEMA,
        "direction": (
            "POLYTOPE -> canonical determinant-line Hilbert carrier -> "
            "independent AG -> theorem on the same POLYTOPE"
        ),
        "scope": {
            "polytope_branch": "compact arithmetic Coxeter polytopes",
            "dimension_range": "n>=4",
            "carrier_species": (
                "canonical determinant line plus every reflection-stable "
                "isotropic discriminant gluing"
            ),
            "noncompact_branch_separated": True,
            "claim_about_other_surface_or_hyperkaehler_carriers": False,
        },
        "trace_shadow_barrier": {
            "theorem": (
                "For compact n>=4 the minimal ground field has d>1.  "
                "The trace type-IV domain has dimension N=d(n+2)-2, the "
                "native Hilbert domain has codimension (d-1)(n+2), and an "
                "F-facet reflection has a d-dimensional Q (-1)-eigenspace. "
                "It is therefore not a Q pseudoreflection and its trace "
                "fixed locus has codimension d, not one."
            ),
            "consequence": (
                "The native F-reflection does not become a branch "
                "pseudoreflection for the full rational trace orthogonal "
                "group.  A coincident Q-root hyperplane would be a separate "
                "root-by-root condition, not a functorial consequence.  The "
                "restriction-of-scalars arithmetic group is contained in a "
                "proper Q-subgroup and has infinite index in the full trace "
                "orthogonal arithmetic group."
            ),
            "what_remains_possible": (
                "A classical trace-lattice form may restrict to the Hilbert "
                "subdomain, but a native Hilbert form is not forced to "
                "extend and cannot be normed over an infinite coset set."
            ),
            "sample_rows": rows,
        },
        "classical_theorem_scope": {
            "Ma_2018": {
                "input": (
                    "an integral Z-lattice of signature (2,N), with an "
                    "arithmetic group finite-index in its rational "
                    "orthogonal group"
                ),
                "usable_on_native_compact_Hilbert_carrier": False,
                "bounded_slope_finiteness_requires_reflective_form": True,
            },
            "Wang_2025_2026": {
                "input": (
                    "an even Z-lattice of signature (N,2); the strongest "
                    "root-system theorems often additionally use 2U+L or U+K"
                ),
                "usable_on_native_compact_Hilbert_carrier": False,
                "trace_numerology_is_not_a_transfer": True,
            },
            "Gritsenko_Hulek_Kodaira": {
                "input": (
                    "a genuinely strongly reflective divisor on the same "
                    "orthogonal modular quotient"
                ),
                "strong_reflective_support_forced_by_P": False,
                "Kodaira_label_alone_bounds_P": False,
            },
        },
        "candidate_class_hostile_audit": [
            {
                "candidate": "all_compact_marked_arithmetic_P_n_ge_4",
                "surface_free_raw_definition": True,
                "known_non_singleton": True,
                "exact_pure_support_forced": False,
                "native_Hilbert_classification_available": False,
                "passes_L4": False,
            },
            {
                "candidate": "compact_complete_root_chambers",
                "surface_free_raw_definition": True,
                "same_open_chamber_return_if_no_hidden_cuts": True,
                "exact_pure_support_forced": False,
                "native_Hilbert_classification_available": False,
                "passes_L4": False,
            },
            {
                "candidate": "vanishing_Whittaker_obstruction_space",
                "surface_free_raw_definition": True,
                "independent_Coxeter_definition": False,
                "rejection": (
                    "Although computable from the lattice, this defines the "
                    "class by a modular carrier gate equivalent to universal "
                    "Borcherds principality; even granting it, no compact "
                    "Hilbert-reflective classification is available."
                ),
                "passes_L4": False,
            },
            {
                "candidate": "bounded_slope_or_reflective_form_exists",
                "surface_free_raw_definition": False,
                "rejection": "circular by the flagship acceptance test",
                "passes_L4": False,
            },
            {
                "candidate": "trace_lattice_in_Ma_or_Wang_list",
                "surface_free_raw_definition": True,
                "rejection": (
                    "membership does not turn an F-reflection into a Q "
                    "reflection, does not force restriction to have exactly "
                    "the P mirrors, and builds the desired bound into a "
                    "classified finite list"
                ),
                "passes_L4": False,
            },
        ],
        "native_Hilbert_exact_support_gates": {
            "fixed_Borcherds_principal_part": (
                "finite cusp-pairing linear algebra, hence decidable for a "
                "fixed enumerated support"
            ),
            "forced_by_compact_root_reflectivity": False,
            "Bugaenko_control": (
                "the unique determinant extension has an obstructed pure "
                "norm-two column and obstructed complete first compensator "
                "shell; this is not a no-go for every arbitrary reflective form"
            ),
            "arbitrary_automorphic_exact_support_decidable": False,
            "reason": (
                "known converse theorems require a hyperbolic split and do "
                "not cover the anisotropic compact Hilbert lattice"
            ),
            "compensated_divisor_counts_as_pure_support": False,
        },
        "same_polytope_return_gate": {
            "required": [
                "enumerate every determinant-line discriminant gluing",
                "prove every global divisor component is reflective",
                "prove no zero-auxiliary mirror cuts the open chamber",
                "prove the original P is the complete chamber, not a chosen subchamber",
                "retain Coxeter labels and wall normalizations",
            ],
            "finite_carrier_lattices_alone_imply_finite_polytopes": False,
        },
        "minimum_completion_package": [
            {
                "name": "automatic_bounded_purity",
                "statement": (
                    "For an independently defined non-singleton compact raw "
                    "Coxeter class C, every canonical determinant gluing has "
                    "a pure globally reflective divisor returning exactly P, "
                    "with multiplicity/slope bounded independently of the "
                    "facet count; or all alternative gluings are excluded."
                ),
                "known": False,
            },
            {
                "name": "Hilbert_reflective_finiteness",
                "statement": (
                    "For anisotropic Res_{F/Q} SO(2,n) quotients, pure "
                    "reflective forms of that bounded multiplicity/slope "
                    "occur in only finitely many effectively enumerable "
                    "F-lattices and dimensions."
                ),
                "known": False,
            },
        ],
        "verdict": {
            "strict_noncircular_compact_L4_subclass_found": False,
            "new_polytope_bound_claimed": False,
            "new_nonexistence_theorem_for_polytopes_claimed": False,
            "result_level": (
                "theorem-shaped trace/reflection obstruction plus exact "
                "two-lemma completion target"
            ),
            "rating": "L3_CARRIER_PLUS_L4_SCOPE_OBSTRUCTION_NOT_L4",
        },
        "primary_sources": [
            "https://arxiv.org/abs/0908.3076",
            "https://arxiv.org/abs/1701.03225",
            "https://arxiv.org/abs/2301.11536",
            "https://arxiv.org/abs/2301.12606",
            "https://arxiv.org/abs/1202.3361",
            "https://arxiv.org/abs/2512.20701",
        ],
        "checks": compact_range_checks,
    }
    payload["all_checks_pass"] = all(compact_range_checks.values())
    digest_payload = dict(payload)
    payload["certificate_sha256"] = _canonical_sha256(digest_payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    payload = analyze()
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if archived != payload:
            raise SystemExit("archived certificate differs from exact replay")
        print(payload["certificate_sha256"])
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(payload["certificate_sha256"])


if __name__ == "__main__":
    main()
