#!/usr/bin/env python3
"""Hostile operational audit of the prime-binomial Coble/K3-rank pilot.

The audit freezes the actual repository state: four row-specific calibrations
exist, but there is no shared arbitrary-input prime-binomial recognizer.  The
K3 arithmetic is recorded only as a conditional payoff, never as an earned
pre/post-AG delta.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/coble_prime_binomial_k3_rank_pilot.json"
CANDIDATE_JSON = ROOT / "data/enumerations/coble_lattice_determinantal_candidate.json"
CANDIDATE_SCRIPT = ROOT / "scripts/certify_coble_lattice_determinantal_candidate.py"
GENERIC_JSON = ROOT / "data/enumerations/coble_finite_geometry_raw_recognizer.json"
GENERIC_SCRIPT = ROOT / "scripts/certify_coble_finite_geometry_raw_recognizer.py"


def k3_rank(boundary_count: int) -> int:
    if boundary_count < 1:
        raise ValueError("a strong Coble boundary is nonempty")
    return 10 + boundary_count


def polytope_dimension(boundary_count: int) -> int:
    return k3_rank(boundary_count) - 1


def ceiling(characteristic: int) -> int:
    return 20 if characteristic == 0 else 22


def _cli_options(tree: ast.AST) -> list[str]:
    return sorted(
        {
            value.value
            for node in ast.walk(tree)
            for value in ast.walk(node)
            if isinstance(value, ast.Constant)
            and isinstance(value.value, str)
            and value.value.startswith("--")
        }
    )


def _top_level_path_constants(tree: ast.Module) -> list[str]:
    wanted = {"H10", "H11", "CI4", "HV"}
    return sorted(
        target.id
        for node in tree.body
        if isinstance(node, ast.Assign)
        for target in node.targets
        if isinstance(target, ast.Name) and target.id in wanted
    )


def _build_certificate_arity(tree: ast.Module) -> int:
    function = next(
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "build_certificate"
    )
    return len(function.args.posonlyargs) + len(function.args.args)


def build_certificate() -> dict:
    candidate = json.loads(CANDIDATE_JSON.read_text(encoding="utf-8"))
    generic = json.loads(GENERIC_JSON.read_text(encoding="utf-8"))
    candidate_source = CANDIDATE_SCRIPT.read_text(encoding="utf-8")
    generic_source = GENERIC_SCRIPT.read_text(encoding="utf-8")
    candidate_tree = ast.parse(candidate_source)
    generic_tree = ast.parse(generic_source)

    controls = candidate["calibrations"]
    expected = {
        "H10_30": (0, 1, 10),
        "H11_26": (0, 2, 11),
        "H14_60": (3, 5, 14),
        "H17_48": (3, 8, 17),
    }
    rows = []
    for key, (characteristic, boundary_count, dimension) in expected.items():
        source = controls[key]
        rows.append(
            {
                "id": key,
                "support": source["support"],
                "boundary_count": boundary_count,
                "dimension": dimension,
                "cover_embedded_rank": k3_rank(boundary_count),
                "within_conditional_k3_ceiling": k3_rank(boundary_count)
                <= ceiling(characteristic),
                "row_specific_program": source["raw_program"],
            }
        )

    generic_positives = [
        case["recognition"]
        for case in generic["case_results"]
        if case["expected_accepted"]
    ]
    generic_boundary_counts = sorted(
        {case["input_summary"]["boundary_count"] for case in generic_positives}
    )
    generic_dimensions = sorted(
        {case["input_summary"]["hyperbolic_dimension"] for case in generic_positives}
    )

    operational = {
        "candidate_cli_options": _cli_options(candidate_tree),
        "candidate_has_arbitrary_input_option": "--input"
        in _cli_options(candidate_tree),
        "candidate_hardwired_source_constants": _top_level_path_constants(
            candidate_tree
        ),
        "candidate_build_certificate_arity": _build_certificate_arity(candidate_tree),
        "candidate_hardcodes_accepted_component": "accepted_component" in candidate_source,
        "candidate_hardcodes_discarded_component": "discarded_component"
        in candidate_source,
        "generic_recognizer_has_input_option": "--input" in _cli_options(generic_tree),
        "generic_recognizer_boundary_counts": generic_boundary_counts,
        "generic_recognizer_dimensions": generic_dimensions,
    }

    counterexamples = [
        {
            "row": "H10_30",
            "component": "c=1",
            "prime_binomial": True,
            "saturated_exponent_lattice": True,
            "failure": controls["H10_30"]["discarded_component"],
        },
        {
            "row": "H11_26",
            "component": "a=-b",
            "prime_binomial": True,
            "saturated_exponent_lattice": True,
            "failure": controls["H11_26"]["discarded_component"],
        },
    ]

    checks = {
        "four_separately_routed_controls_exist": len(rows) == 4,
        "candidate_has_no_arbitrary_input_option": not operational[
            "candidate_has_arbitrary_input_option"
        ],
        "candidate_reads_four_hardwired_sources": operational[
            "candidate_hardwired_source_constants"
        ]
        == ["CI4", "H10", "H11", "HV"],
        "candidate_builder_accepts_no_raw_datum": operational[
            "candidate_build_certificate_arity"
        ]
        == 0,
        "good_and_bad_factors_are_hardcoded": operational[
            "candidate_hardcodes_accepted_component"
        ]
        and operational["candidate_hardcodes_discarded_component"],
        "prime_binomial_alone_has_two_exact_false_positives": all(
            row["prime_binomial"] and row["saturated_exponent_lattice"]
            for row in counterexamples
        ),
        "actual_generic_recognizer_is_thin_below_k3_ceiling": max(
            generic_dimensions
        )
        <= 17
        and max(generic_boundary_counts) <= 8,
        "no_recognized_pre_AG_b_value_exceeds_complex_K3_bound": max(
            generic_boundary_counts
        )
        <= 10,
        "PB_Atlas_is_open": True,
        "strict_L4_not_claimed": True,
    }

    return {
        "schema": "coxeter.coble_prime_binomial_k3_rank_l4_pilot.v2",
        "status_date": "2026-08-10",
        "status": "FAIL_NO_SHARED_RAW_RECOGNIZER; PB-Atlas is the one exact missing lemma",
        "verdict": {
            "four_rows_share_prime_binomial_outputs_after_separate_routing": True,
            "four_rows_pass_one_invariant_raw_prime_binomial_recognizer": False,
            "reason": (
                "the candidate executable has no raw-input path, reads four named source "
                "archives, and hardcodes accepted/discarded factors; the only generic raw "
                "recognizer covers a different rigid five-core class"
            ),
            "strict_L4": False,
            "pre_AG_allowed_b_values_wider_than_K3_bound_proved": False,
        },
        "operational_audit": operational,
        "prime_binomial_false_positives": counterexamples,
        "separately_routed_controls": rows,
        "only_current_generic_raw_class": {
            "name": "C_FG",
            "boundary_counts": generic_boundary_counts,
            "dimensions": generic_dimensions,
            "K3_rank_pushback_removes_recognized_members": 0,
            "nonthin_delta": False,
        },
        "conditional_payoff_after_PB_Atlas": {
            "complex": {"max_boundary_count": 10, "max_dimension": 19},
            "positive_characteristic_not_2": {
                "max_boundary_count": 12,
                "max_dimension": 21,
            },
            "uses_classification": False,
            "earned_today": False,
        },
        "missing_lemma": {
            "id": "PB-Atlas",
            "statement": (
                "from an arbitrary raw Gram/wall/contraction branch, canonically generate all "
                "incidence, collision, Jacobian, irreducibility and Hilbert equations and prove "
                "that every primitive-binomial regular open is exactly a nonempty atlas of the "
                "required strong-Coble incidence stratum"
            ),
            "proved": False,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def certify(certificate: dict) -> None:
    if not certificate["all_checks_pass"]:
        raise ArithmeticError("a hostile pilot audit check failed")
    if certificate["verdict"]["strict_L4"]:
        raise ArithmeticError("the absent shared recognizer was promoted to L4")
    if certificate["verdict"][
        "pre_AG_allowed_b_values_wider_than_K3_bound_proved"
    ]:
        raise ArithmeticError("a synthetic b-range was counted as recognized survivors")
    if certificate["missing_lemma"]["proved"]:
        raise ArithmeticError("PB-Atlas must remain explicitly open")


def main() -> int:
    certificate = build_certificate()
    certify(certificate)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if OUTPUT.exists() and OUTPUT.read_text(encoding="utf-8") != rendered:
        raise SystemExit(f"stale output: run with --write to refresh {OUTPUT}")
    print("PASS: hostile audit confirms no shared raw recognizer and no strict L4 theorem")
    return 0


if __name__ == "__main__":
    import sys

    if sys.argv[1:] == ["--write"]:
        certificate = build_certificate()
        certify(certificate)
        OUTPUT.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        raise SystemExit(0)
    raise SystemExit(main())
