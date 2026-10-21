#!/usr/bin/env python3
"""Replay the frozen fixed-data/outer-routing synthesis.

This is intentionally lightweight.  It pins the six theorem/audit reports
and checks that every designated synthesis document retains the three exact
scope boundaries: classification-importing finite-pair closure, the single
X3 counterexample refuting universality over the all-finite-volume class
without settling the maximal-only version, and the absolutely maximal H21
mixed-carrier obstruction.
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
    / "data/enumerations/ag_to_absw_fixed_data_outer_routing_synthesis.json"
)
SCHEMA = "coxeter.ag_to_absw_fixed_data_outer_routing_synthesis.v2"

SOURCES = {
    "finite_pair_owner": (
        "results/projected_k3_maximal_abstract_chamber_finiteness.md",
        "59dd458acf14d18021f6c5583096da13f39c83115fc136d6ba91bcc850ce3603",
    ),
    "finite_pair_audit": (
        "results/projected_k3_maximal_abstract_chamber_finiteness_independent_hostile_audit.md",
        "edc663c0e7f35880123149d81d832812d5f28837c837df16eaa0f1daf725916a",
    ),
    "x3_owner": (
        "results/projected_k3_universality_rank_saturation_no_go.md",
        "7ed5176c0c0a36e036c0e12dfa6acf1ee15d90d5ec7143ca5df59ebc3b7143d2",
    ),
    "x3_audit": (
        "results/projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.md",
        "89296f04bac6b0b4ed68d9b86936a08ede515ba5bf9c847abeaadbc628926071",
    ),
    "h21_owner": (
        "results/borcherds_h21_outer_routing_counterexample.md",
        "fa899977f418d6fa93f7d5270e067c28a648c5e55f2636cb2b362e381534821d",
    ),
    "h21_audit": (
        "results/borcherds_h21_outer_routing_counterexample_independent_hostile_audit.md",
        "8d705a619e9c1ffff3ad69765e90bb1a1a2a5d8713015ecca79382e8a1ccbc86",
    ),
}

DOCUMENTS = (
    "GOAL.md",
    "RESEARCH_PLAN.md",
    "results/canonical_chamber_four_property_completion_audit.md",
    "results/inverse_carrier_transfer_ledger.md",
    "results/nikulin_to_ag_translation_ledger.md",
    "results/ag_to_absw_fixed_data_closure_and_outer_routing_obstruction.md",
)

REQUIRED_LINKS = (
    "projected_k3_maximal_abstract_chamber_finiteness.md",
    "projected_k3_universality_rank_saturation_no_go.md",
    "borcherds_h21_outer_routing_counterexample.md",
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_certificate() -> dict[str, Any]:
    source_sha256: dict[str, str] = {}
    source_pins_exact = True
    for name, (relative, expected) in SOURCES.items():
        actual = _sha(ROOT / relative)
        source_sha256[name] = actual
        source_pins_exact &= actual == expected

    texts = {
        relative: (ROOT / relative).read_text(encoding="utf-8")
        for relative in DOCUMENTS
    }
    document_sha256 = {
        relative: _sha(ROOT / relative) for relative in DOCUMENTS
    }
    links_in_every_document = all(
        all(link in text for link in REQUIRED_LINKS) for text in texts.values()
    )
    aggregate = "\n".join(texts.values())
    forbidden_field_shorthands = (
        "every exact finite-volume full-trace K3 carrier",
        "every exact finite-volume K3 branch",
        "after exact finite-volume K3 return",
        "exact finite-volume primitive full-trace K3 carrier",
        "exact compact full-trace K3 carrier",
        "finite-volume exact-K3 field implication",
        "every finite-volume full-trace K3 branch",
        "closes every exact finite-volume carrier",
        "exact finite-volume field theorem",
        "Once the full trace carrier is present",
        "dimension and field implications are now proved once",
        "Nikulin/Vinberg are hostile benchmarks only",
    )
    checks = {
        "six_frozen_source_hashes_exact": source_pins_exact,
        "three_result_links_in_every_synthesis_document": links_in_every_document,
        "finite_pair_step_explicitly_imports_classification": (
            "classification-import" in aggregate
            and "Nikulin's 1982" in aggregate
            and "finite projected ambient hosts" in aggregate
        ),
        "x3_refutes_all_fv_universality_but_not_each_input_or_maximal_only": (
            "all-finite-volume" in aggregate
            and "one counterexample" in aggregate
            and "not proved absolutely" in aggregate
            and "mixed-divisor K3 realization" in aggregate
        ),
        "field_scope_retains_exact_admissible_full_trace_order_module": (
            "exact admissible finite-volume full restriction-of-scalars"
            in aggregate
            and "trace/order-module complex-K3" in aggregate
        ),
        "field_scope_has_no_local_shorthand_regression": not any(
            phrase in aggregate for phrase in forbidden_field_shorthands
        ),
        "finite_pair_count_retains_reflection_maximality": (
            "finite native pairs and reflection maximality" in aggregate
            and "maximal finite-volume arithmetic chambers" in aggregate
        ),
        "h21_scope_is_absolute_and_mixed_positive": (
            "absolutely reflection-maximal" in aggregate
            and "rank-$22$ Borcherds" in aggregate
            and "mixed divisor/action" in aggregate
            and "characteristic-$2$" in aggregate
        ),
        "global_absw_and_ag_only_claims_are_denied": (
            "not an AG-only" in aggregate
            and "not a new global ABSW proof" in aggregate
            and "does not prove global ABSW" in aggregate
        ),
        "no_h21_pending_marker": (
            "H21 independent audit pending" not in aggregate
            and "borcherds_h21_outer_routing_counterexample_independent_hostile_audit.md"
            in aggregate
        ),
    }
    scope = {
        "finite_pair_to_chambers": "proved for maximal finite-volume arithmetic chambers in dimension at least two, conditional on a finite native pair set",
        "classification_provenance": "imports Nikulin 1982 fixed-(k,m) reflective-lattice similarity finiteness; not AG-only",
        "projected_ambient_host_finiteness_needed_for_abstract_count": False,
        "x3_no_go_scope": "one X3 input refutes the universal assertion over all finite-volume exact-module projected-root complex-K3 inputs; it does not say every input fails",
        "x3_maximal_only_universality_settled": False,
        "x3_mixed_divisor_k3_positive": True,
        "h21_absolute_reflection_maximality": True,
        "h21_current_complex_stack_empty": True,
        "h21_characteristic_two_mixed_curve_cremona_positive": True,
        "universal_complex_projected_root_k3": False,
        "global_absw_proved": False,
        "future_carrier": "mixed divisor/action or another genuinely broader carrier",
    }
    return {
        "schema": SCHEMA,
        "source_sha256": source_sha256,
        "document_sha256": document_sha256,
        "checks": checks,
        "scope": scope,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_certificate()
    if args.check:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != result:
            raise ArithmeticError("synthesis archive does not match live inputs")
        if not result["all_checks_pass"]:
            raise ArithmeticError("a synthesis scope check failed")
    if args.write:
        OUTPUT.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
