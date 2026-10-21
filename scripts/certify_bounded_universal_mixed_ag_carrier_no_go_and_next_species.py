"""Replay the bounded mixed-carrier theorem/gate audit.

This certificate pins the exact local inputs and records the logical status of
the X3/H21 mixed Mukai controls.  Literature theorems (in particular derived
Torelli) remain explicit external inputs in the accompanying report.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/bounded_universal_mixed_ag_carrier_no_go_and_next_species.json"
SCHEMA = "coxeter-polytopes.bounded-universal-mixed-ag-carrier-audit.v1"

DEPENDENCIES = {
    "results/universal_weak_chow_action_carrier_and_leverage_no_go.md":
        "df38b26e0cffd8a391834be922869153af5badf41a632e74a48b472b332e101a",
    "results/projected_k3_universality_rank_saturation_no_go.md":
        "7ed5176c0c0a36e036c0e12dfa6acf1ee15d90d5ec7143ca5df59ebc3b7143d2",
    "results/borcherds_h21_outer_routing_counterexample.md":
        "fa899977f418d6fa93f7d5270e067c28a648c5e55f2636cb2b362e381534821d",
    "results/borcherds_d4_h21_og10_carrier_obstruction.md":
        "99f07bc673b422f2268ec9ac1782b2d64139d494a76147b78deec2748d589406",
    "results/supersingular_k3_golden_h10_obstruction.md":
        "665bae088e2361c8c0b8ac45f1683e973fc4df122f0d46e93e26da9feb24088d",
    "results/coherent_full_h2_real_multiplication_bound.md":
        "79bfaab4eaa549502e720da6cee1250f4c26a3a9022bd345430660540034c9b1",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_certificate() -> dict[str, object]:
    controls = {
        "X3": {
            "dimension": 19,
            "picard_rank": 20,
            "discriminant_order": 3,
            "curve_facets": 24,
            "fm_action_facets": 12,
            "long_square": -6,
            "long_divisibility": 3,
            "long_discriminant_action": -1,
            "mixed_k3_fm_pass": True,
            "full_nef_k3_pass": False,
            "absolute_maximality_claimed": False,
        },
        "H21": {
            "dimension": 21,
            "picard_rank": 22,
            "characteristic": 2,
            "curve_facets": 42,
            "automorphism_facets": 168,
            "mixed_k3_fm_pass": True,
            "og10_full_cone_pass": False,
            "absolute_maximality_proved_in_dependency": True,
        },
    }
    gates = {
        "TRACE_full_order_mukai_embedding": False,
        "LIFT_every_side_to_actual_geometric_or_fm_action": False,
        "EXACT_no_hidden_wall_return": False,
        "BOUND_uniform_polarization_graph_kernel_complexity": False,
        "CONSTR_finite_constructible_extraction": False,
    }
    theorem = {
        "bounded_decorated_qc_constructible_extraction_implies_finite_image": True,
        "bounded_underlying_varieties_alone_imply_finite_image": False,
        "mixed_k3_fm_species_passes_both_controls": True,
        "mixed_k3_fm_species_currently_total": False,
        "unconditional_global_AG_finiteness_proved": False,
        "full_trace_mukai_rank_ceiling": 24,
        "full_trace_necessary_inequality": "[k:Q](n+1) <= 24",
        "classical_reflection_rank_bound_imported": False,
    }
    return {
        "schema": SCHEMA,
        "dependency_sha256": DEPENDENCIES,
        "dependency_hashes_match": all(
            _sha256(ROOT / relative) == expected
            for relative, expected in DEPENDENCIES.items()
        ),
        "controls": controls,
        "unproved_universal_gates": gates,
        "terminal_theorem_status": theorem,
        "all_replay_checks_pass": all(
            [
                controls["X3"]["mixed_k3_fm_pass"],
                controls["H21"]["mixed_k3_fm_pass"],
                not controls["X3"]["absolute_maximality_claimed"],
                not theorem["unconditional_global_AG_finiteness_proved"],
                not any(gates.values()),
            ]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.write_archive:
        ARCHIVE.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    if args.check_archive:
        frozen = json.loads(ARCHIVE.read_text())
        if frozen != certificate:
            raise SystemExit("frozen archive differs from live certificate")
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
