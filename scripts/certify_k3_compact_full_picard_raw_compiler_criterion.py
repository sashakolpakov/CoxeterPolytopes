#!/usr/bin/env python3
"""Replay finite interfaces in the compact full-Picard K3 iff compiler."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_compact_full_picard_raw_compiler_criterion.md"
ARCHIVE = (
    ROOT
    / "data/enumerations/k3_compact_full_picard_raw_compiler_criterion.json"
)


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quadratic_sign(a: int, b: int) -> int:
    """Return sign(a + b*sqrt(5)) exactly."""

    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    comparison = a * a - 5 * b * b
    if comparison == 0:
        raise AssertionError("sqrt(5) would be rational")
    if a > 0:
        return 1 if comparison > 0 else -1
    return 1 if comparison < 0 else -1


def golden_energy_shell() -> dict[str, Any]:
    """Exact trace -2, totally-negative shell in O_K^vee for Q(sqrt(5)).

    Write z=(u+v*phi)/sqrt(5).  Its trace is v.  For v=-2, negativity
    at the two embeddings is equivalent to

        1-sqrt(5) < u < 1+sqrt(5),

    hence u is one of -1,0,1,2,3.
    """

    v = -2
    coordinates: list[list[int]] = []
    for u in range(-1, 4):
        tau0_sign = quadratic_sign(2 * u + v, v)
        tau1_sign = -quadratic_sign(2 * u + v, -v)
        if tau0_sign != -1 or tau1_sign != -1:
            raise AssertionError("claimed golden energy point is not negative")
        coordinates.append([u, v])

    lower_tau1_sign = -quadratic_sign(2 * (-2) + v, -v)
    upper_tau0_sign = quadratic_sign(2 * 4 + v, v)
    if lower_tau1_sign != 1 or upper_tau0_sign != 1:
        raise AssertionError("golden energy-shell integer endpoints changed")

    return {
        "field": "Q(sqrt(5))",
        "order": "Z[(1+sqrt(5))/2]",
        "codifferent_coordinates": "(u+v*phi)/sqrt(5)",
        "trace_coordinate": "v",
        "required_trace": -2,
        "integer_u_lower": -1,
        "integer_u_upper": 3,
        "coordinates": coordinates,
        "cardinality": len(coordinates),
        "lower_excluded_u": -2,
        "upper_excluded_u": 4,
    }


def hidden_root_ellipsoid_toy() -> dict[str, Any]:
    """Replay the explicit compact-window root bound on diag(2,-2,-2)."""

    # p=(1,0,0), cosh(R)=2, hence 2*cosh(2R)=14.
    ellipsoid_bound = 14
    coordinate_bound = 2  # 2*(a^2+b^2+c^2) <= 14.
    signed_roots: list[list[int]] = []
    for a, b, c in itertools.product(
        range(-coordinate_bound, coordinate_bound + 1), repeat=3
    ):
        lorentz_square = 2 * a * a - 2 * b * b - 2 * c * c
        q_p = 2 * (a * a + b * b + c * c)
        if lorentz_square == -2 and q_p <= ellipsoid_bound:
            signed_roots.append([a, b, c])

    if len(signed_roots) != 12:
        raise AssertionError("toy hidden-root window changed")
    if any(abs(row[0]) > 1 for row in signed_roots):
        raise AssertionError("a root in the window misses the radius bound")

    # The next root shell has |a|=2, Q_p=18, and sinh distance=2>sqrt(3).
    next_shell: list[list[int]] = []
    for a, b, c in itertools.product(range(-3, 4), repeat=3):
        lorentz_square = 2 * a * a - 2 * b * b - 2 * c * c
        q_p = 2 * (a * a + b * b + c * c)
        if lorentz_square == -2 and q_p == 18:
            next_shell.append([a, b, c])
    if len(next_shell) != 16 or any(abs(row[0]) != 2 for row in next_shell):
        raise AssertionError("toy first excluded root shell changed")

    return {
        "lattice_gram": [[2, 0, 0], [0, -2, 0], [0, 0, -2]],
        "center": [1, 0, 0],
        "cosh_radius": 2,
        "sinh_radius_squared": 3,
        "positive_definite_Qp_gram": [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
        "ellipsoid_bound_2cosh_2R": ellipsoid_bound,
        "coordinate_search_bound": coordinate_bound,
        "signed_root_count": len(signed_roots),
        "unsigned_mirror_count": len(signed_roots) // 2,
        "signed_roots": sorted(signed_roots),
        "first_excluded_Qp": 18,
        "first_excluded_signed_root_count": len(next_shell),
    }


def nonrational_rank_window() -> list[dict[str, int | bool]]:
    rows: list[dict[str, int | bool]] = []
    for module_rank in range(3, 11):
        ceiling = 20 // module_rank
        rows.append(
            {
                "module_rank": module_rank,
                "hyperbolic_dimension": module_rank - 1,
                "field_degree_ceiling": ceiling,
                "nonrational_possible": ceiling >= 2,
            }
        )
    return rows


def build_certificate() -> dict[str, Any]:
    energy = golden_energy_shell()
    hidden = hidden_root_ellipsoid_toy()
    ranks = nonrational_rank_window()
    if energy["cardinality"] != 5:
        raise AssertionError("golden physical energy shell changed")
    if hidden["unsigned_mirror_count"] != 6:
        raise AssertionError("toy hidden mirror count changed")
    if max(row["field_degree_ceiling"] for row in ranks) != 6:
        raise AssertionError("K3 field-degree ceiling changed")

    return {
        "schema": "coxeter.k3-compact-full-picard-raw-compiler-criterion.v1",
        "report_sha256": sha256_file(REPORT),
        "input_scope": {
            "pinned_exact_algebraic_raw_realization": True,
            "compact": True,
            "totally_real_nonrational": True,
            "surface_or_integral_trace_module_in_input": False,
            "abstract_diagram_without_exact_dotted_data_sufficient": False,
        },
        "rank_gate": {
            "K3_Picard_ceiling": 20,
            "necessary_inequality": "[k:Q]*(n+1)<=20",
            "rows": ranks,
        },
        "physical_normalization_gate": {
            "cycle_order_generated_from_raw_scaling_invariants": True,
            "cycle_order_contained_in_every_successful_coefficient_order": True,
            "finite_overorders_between_cycle_and_maximal_order": True,
            "diagonal_shell": (
                "u in O0^vee; all embeddings u<0; Tr(u)=-2"
            ),
            "place_energies_sum": 2,
            "common_square_class_test": (
                "(u_i/h_i)/(u_1/h_1) is a square in k"
            ),
            "all_facet_lines_require_same_module_minus_two_roots": True,
            "golden_exact_shell": energy,
        },
        "module_gate": {
            "root_frame": "F0=sum_i O0*delta_i",
            "finite_interval": "F0 subset M subset F0#",
            "finite_discriminant_subgroup_enumeration": True,
            "multiplier_order": "O_M={a in k:aM subset M}",
            "all_orders_O0_subset_O_subset_O_M_enumerated": True,
            "projectivity_over_O0_required": False,
            "projectivity_over_multiplier_order_required": False,
            "filters": [
                "integral",
                "even",
                "O-stable for an enumerated overorder O",
                "projective over that enumerated order O",
                "all raw facet field reflections preserve M",
                "every raw facet line contains a trace square -2 root",
            ],
        },
        "embedding_gate": {
            "criterion": (
                "even T of signature (2,20-rk(N)) with q_T isometric to -q_N"
            ),
            "necessary_and_sufficient": True,
            "nonstable_local_discriminant_conditions_required": True,
            "stable_length_sieve_alone_is_exact": False,
        },
        "hidden_root_gate": {
            "all_place_energy_bound": 2,
            "zero_distinguished_projection_is_terminal_rejection": True,
            "Qp": "2*B(p,z)^2/B(p,p)-B(z,z)",
            "Qp_positive_definite": True,
            "relevant_root_bound": "Qp(delta)<=2*cosh(2R)",
            "exact_lattice_enumeration_finite": True,
            "toy": hidden,
        },
        "construction": {
            "generic_period_gives_exact_NS": True,
            "Weyl_framing_gives_exact_eigen_nef_section": True,
            "selected_p_lies_in_E_and_is_ample": True,
            "selected_projective_p_lies_in_P_interior": True,
            "return_is_closure_of_projectivized_Amp_intersection": True,
            "K3_root_effectivity_gives_physical_curves": True,
            "full_Picard_makes_physical_roots_in_trace": True,
            "canonical_class": 0,
        },
        "maximality_scope": {
            "ordinary_K3_root_reflection_preserves_eigenspace_in_general": False,
            "ordinary_root_Q_minus_dimension": 1,
            "native_root_Q_minus_dimension": "[k:Q]",
            "native_reflection_integrality_automatic_for_every_K3_root": False,
            "compact_hidden_root_test_removed_by_raw_maximality_alone": False,
            "same_host_native_closure_plus_relative_maximality_is_shortcut": True,
            "zero_projection_test_still_required_under_shortcut": True,
        },
        "cusped_comparison": {
            "exact_finite_volume_full_trace_forces_k_Q": True,
            "bounded_closed_star_physical_basis_constant": 13,
            "bounded_spanning_theorem_already_proved": True,
            "finite_rational_Picard_host_census": True,
            "full_host_maximality_excludes_hidden_root_cuts": True,
            "separate_from_compact_nonrational_compiler": True,
        },
        "verdict": {
            "finite_candidate_list": True,
            "necessary_and_sufficient_for_pinned_realization": True,
            "compiler_may_return_EMPTY": True,
            "universal_nonemptiness_for_every_raw_maximal_P": False,
            "Property_4_essential_surjectivity_proved": False,
            "first_nonautomatic_gate": "PHYS",
            "later_exact_gates": ["EMB", "ROOT"],
        },
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-output", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.write_output:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(encoded, encoding="utf-8")
    if args.check_output:
        if not ARCHIVE.exists() or ARCHIVE.read_text(encoding="utf-8") != encoded:
            raise SystemExit("archived compact K3 compiler certificate differs")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
