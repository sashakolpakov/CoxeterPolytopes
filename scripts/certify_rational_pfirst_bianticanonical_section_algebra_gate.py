#!/usr/bin/env python3
"""Exact finite replay for the bianticanonical Property-4 section-algebra gate.

The geometric existence and Bezout arguments are proved in the companion
report.  This certifier checks the marked I_(1,10) arithmetic, the two formal
fat-point systems, their intersection budgets, and the identical proper-point
proximity/contraction skeleton on the two geometric branches.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/enumerations/rational_pfirst_bianticanonical_section_algebra_gate.json"
)
SCHEMA = "rational_pfirst_bianticanonical_section_algebra_gate/v1"
WINGER_ARCHIVE = ROOT / "data/enumerations/a3_h3_plane_arrangement_coble_factories.json"
UPSTREAM_PINS = {
    "results/h3_winger_arrangement_coble_natural_cone.md": (
        "dafca20b2843488a501d7ac21bc618f380d9d660456864c068cbc0afd428e035"
    ),
    "scripts/certify_a3_h3_plane_arrangement_coble_factories.py": (
        "572e4dbbe246d75876a493b87067e62922f1a0cff65d396b1e30affcf7b5d31a"
    ),
    "data/enumerations/a3_h3_plane_arrangement_coble_factories.json": (
        "01f6649f696d8ed6064bc9f3e567e1c9a2d930bc0bee9c9c9720f8d9d6b43fb6"
    ),
    "tests/test_a3_h3_plane_arrangement_coble_factories.py": (
        "5243c2d67128c72e5a8bd713c8edb6a1a9bdfddf345db9df2035a836be96e3de"
    ),
    "results/reflection_maximal_coble_exact_nef_return.md": (
        "cfd6e175165244b5a4df6e94d888fbe7803d10b0ae7223dfd0c111a663eb67be"
    ),
    "results/reflection_maximal_coble_exact_nef_return_hostile_audit.md": (
        "9da4a73e731322bc7b0654224431ceb9b5f2aefac70bd0941ca5139fdb65e5e6"
    ),
}


def dot(x: tuple[int, ...], y: tuple[int, ...]) -> int:
    if len(x) != 11 or len(y) != 11:
        raise ValueError("I_(1,10) vectors must have length eleven")
    return x[0] * y[0] - sum(a * b for a, b in zip(x[1:], y[1:]))


def plane_monomial_count(degree: int) -> int:
    if isinstance(degree, bool) or not isinstance(degree, int) or degree < 0:
        raise ValueError("degree must be a nonnegative integer")
    return math.comb(degree + 2, 2)


def ordinary_jet_count(multiplicity: int) -> int:
    if (
        isinstance(multiplicity, bool)
        or not isinstance(multiplicity, int)
        or multiplicity < 0
    ):
        raise ValueError("multiplicity must be a nonnegative integer")
    return math.comb(multiplicity + 1, 2)


def sha256_json(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build() -> dict[str, Any]:
    live_upstream = {
        relative: sha256_file(ROOT / relative) for relative in UPSTREAM_PINS
    }
    if live_upstream != UPSTREAM_PINS:
        raise ArithmeticError("the frozen Winger existence dependency changed")
    winger = json.loads(WINGER_ARCHIVE.read_text())["H3_Winger_surface"]
    winger_projection = {
        "field": winger["field"],
        "center_count": len(winger["center_coordinates"]),
        "all_ten_nodes": winger["all_ten_prescribed_singularities_ordinary_nodes"],
        "absolutely_irreducible": winger["absolute_irreducibility_argument"][
            "absolutely_irreducible"
        ],
        "normalization_genus": winger["absolute_irreducibility_argument"][
            "normalization_genus"
        ],
        "cubic_rank": winger["cubic_rank"],
        "sextic_rank": winger["sextic_rank"],
        "strong_coble": winger["surface_conclusion"]["strong_coble"],
    }
    if winger_projection != {
        "field": "Q(sqrt(5)), characteristic zero",
        "center_count": 10,
        "all_ten_nodes": True,
        "absolutely_irreducible": True,
        "normalization_genus": 0,
        "cubic_rank": 10,
        "sextic_rank": 27,
        "strong_coble": True,
    }:
        raise ArithmeticError("the Winger ten-node projection changed")

    # Coordinates are (H,E_1,...,E_10) in the convention
    # x.y=x_0*y_0-sum x_i*y_i.
    canonical = (-3,) + (1,) * 10
    anticanonical = (3,) + (-1,) * 10
    bianticanonical = tuple(2 * value for value in anticanonical)

    lattice = {
        "lattice": "I_(1,10)",
        "rank": 11,
        "K": list(canonical),
        "F_minus_K": list(anticanonical),
        "B_minus_2K": list(bianticanonical),
        "K_square": dot(canonical, canonical),
        "F_square": dot(anticanonical, anticanonical),
        "B_square": dot(bianticanonical, bianticanonical),
        "K_dot_F": dot(canonical, anticanonical),
        "K_dot_B": dot(canonical, bianticanonical),
        "F_arithmetic_genus": 1
        + (dot(anticanonical, anticanonical) + dot(canonical, anticanonical))
        // 2,
        "B_arithmetic_genus": 1
        + (dot(bianticanonical, bianticanonical) + dot(canonical, bianticanonical))
        // 2,
    }
    if lattice != {
        "lattice": "I_(1,10)",
        "rank": 11,
        "K": [-3] + [1] * 10,
        "F_minus_K": [3] + [-1] * 10,
        "B_minus_2K": [6] + [-2] * 10,
        "K_square": -1,
        "F_square": -1,
        "B_square": -4,
        "K_dot_F": 1,
        "K_dot_B": 2,
        "F_arithmetic_genus": 1,
        "B_arithmetic_genus": 0,
    }:
        raise ArithmeticError("marked I_(1,10) arithmetic changed")

    formal_systems = {
        "ordered_centers": 10,
        "all_centers_proper_and_distinct": True,
        "proximity_parent_lists": [[] for _ in range(10)],
        "proximity_matrix": [
            [1 if i == j else 0 for j in range(10)] for i in range(10)
        ],
        "labelled_exceptional_intersection_matrix": [
            [-1 if i == j else 0 for j in range(10)] for i in range(10)
        ],
        "labelled_exceptional_contraction_subsets": 2**10,
        "labelled_exceptional_complete_contraction_orders": math.factorial(10),
        "minus_K": {
            "plane_degree": 3,
            "multiplicities": [1] * 10,
            "matrix_shape": [
                10 * ordinary_jet_count(1),
                plane_monomial_count(3),
            ],
        },
        "minus_2K": {
            "plane_degree": 6,
            "multiplicities": [2] * 10,
            "matrix_shape": [
                10 * ordinary_jet_count(2),
                plane_monomial_count(6),
            ],
        },
    }
    if formal_systems["minus_K"]["matrix_shape"] != [10, 10]:
        raise ArithmeticError("cubic formal matrix shape changed")
    if formal_systems["minus_2K"]["matrix_shape"] != [30, 28]:
        raise ArithmeticError("sextic formal matrix shape changed")

    bezout = {
        "smooth_cubic_branch": {
            "cubic_against_cubic": {
                "global_intersection": 3 * 3,
                "forced_local_lower_bound": 10 * 1,
            },
            "cubic_against_double_point_sextic": {
                "global_intersection": 3 * 6,
                "forced_local_lower_bound": 10 * 2,
            },
            "conclusion": (
                "every cubic through the centers is C and every sextic double "
                "at them is C squared"
            ),
        },
        "ten_node_sextic_branch": {
            "nodal_sextic_against_cubic": {
                "global_intersection": 6 * 3,
                "forced_local_lower_bound": 10 * 2 * 1,
            },
            "nodal_sextic_against_double_point_sextic": {
                "global_intersection": 6 * 6,
                "forced_local_lower_bound": 10 * 2 * 2,
            },
            "conclusion": (
                "no cubic passes through all nodes and the nodal sextic is the "
                "unique sextic double at them"
            ),
        },
    }
    for branch in bezout.values():
        for key, row in branch.items():
            if key == "conclusion":
                continue
            if row["forced_local_lower_bound"] <= row["global_intersection"]:
                raise ArithmeticError("Bezout strict inequality changed")

    geometric_branches = {
        "smooth_cubic_ten_points": {
            "minus_K_matrix_rank": 9,
            "h0_minus_K": 1,
            "minus_2K_matrix_rank": 27,
            "h0_minus_2K": 1,
            "degree_two_generator": "square_of_degree_one_generator",
            "unique_bianticanonical_divisor": "2F_nonreduced",
            "physical_boundary_gate": False,
        },
        "irreducible_ten_node_rational_sextic": {
            "minus_K_matrix_rank": 10,
            "h0_minus_K": 0,
            "minus_2K_matrix_rank": 27,
            "h0_minus_2K": 1,
            "degree_two_generator": "divisor_prime_generator",
            "unique_bianticanonical_divisor": "B_reduced_irreducible_smooth_rational",
            "physical_boundary_gate": True,
        },
    }
    for row in geometric_branches.values():
        if row["h0_minus_K"] != 10 - row["minus_K_matrix_rank"]:
            raise ArithmeticError("cubic rank/kernel count changed")
        if row["h0_minus_2K"] != 28 - row["minus_2K_matrix_rank"]:
            raise ArithmeticError("sextic rank/kernel count changed")

    criterion = {
        "one_boundary_rank11_equivalence": (
            "B is a reduced irreducible member of |-2K| iff the degree -K "
            "piece is zero and the one-dimensional degree -2K piece is "
            "generated by a divisor-prime section"
        ),
        "exact_return_handoff": (
            "under the listed finite-volume, adapted-host, boundary-sum, and "
            "relative-reflection-maximal hypotheses, divisor-prime sections "
            "for every proposed facet on one oriented marked rational surface "
            "imply literal full-nef return"
        ),
        "exact_return_handoff_hypotheses": [
            "finite_sided_finite_volume_coxeter_chamber",
            "distinct_irredundant_physical_negative_facets",
            "pairwise_orthogonal_minus4_boundaries_summing_to_minus2K",
            "minus1_minus2_ordinary_adjunction_rows",
            "all_side_reflections_preserve_adapted_host",
            "no_proper_discrete_reflection_overgroup_in_Oplus_adapted_host",
            "one_oriented_marked_smooth_projective_rational_surface",
        ],
        "not_determined_by": [
            "marked_NS_and_K",
            "proper_point_proximity_DAG",
            "labelled_exceptional_contraction_boolean_fan",
            "degree_and_multiplicity_vectors",
            "formal_jet_matrix_descriptors",
        ],
        "first_required_geometric_data": [
            "actual_evaluation_matrix_rank",
            "kernel_section_multiplication",
            "divisor_primality_or_absolute_factorization",
        ],
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "status": "PROVED_SCOPED_PROPERTY4_NO_GO_AND_COX_PRIME_CRITERION",
        "lattice": lattice,
        "common_formal_frame": formal_systems,
        "bezout_replay": bezout,
        "geometric_branches": geometric_branches,
        "explicit_ten_node_existence_dependency": {
            "pins": UPSTREAM_PINS,
            "projected_certificate": winger_projection,
        },
        "criterion": criterion,
        "scope": {
            "base_field": "C",
            "same_full_nef_cone_claimed": False,
            "same_raw_coxeter_chamber_claimed": False,
            "global_property4_claimed": False,
            "universal_incidence_nonemptiness_claimed": False,
            "finite_generation_of_Cox_ring_required": False,
        },
    }
    payload["content_sha256"] = sha256_json(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()
    payload = build()
    if args.check_archive:
        frozen = json.loads(ARCHIVE.read_text())
        if payload != frozen:
            raise ArithmeticError("frozen archive changed")
    if args.print_json or not args.check_archive:
        print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
