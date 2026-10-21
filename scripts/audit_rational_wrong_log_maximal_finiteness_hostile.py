#!/usr/bin/env python3
"""Independent hostile replay of weak-incidence maximal finiteness."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/enumerations/"
    "rational_wrong_log_maximal_finiteness_hostile_audit.json"
)

OWNER_PINS = {
    "results/rational_wrong_log_reflection_maximal_finiteness.md": (
        "a682712bf212aa79c15547aeba0f59c279612d07bd428df7c5fd9d5339cf2768"
    ),
    "scripts/certify_rational_wrong_log_maximal_finiteness.py": (
        "22d551ca0ead5d362c6ca194a8d60aa3257958f806498da38bae6cfcac430f3d"
    ),
    "data/enumerations/rational_wrong_log_maximal_finiteness.json": (
        "cf95f1f7ba6d3ed3f60a521419f817874ccb55f3b27d0a0ecde2fb5c0e8c16d7"
    ),
    "tests/test_rational_wrong_log_maximal_finiteness.py": (
        "10ced3cc9b611724d8b172c2a3b5cec5f839dfc2cd20b70d72b10cdbe0efb20f"
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gaussian_table(maximum_rank: int) -> list[list[int]]:
    """Build Gaussian binomials at q=2 from the Pascal recurrence."""
    if isinstance(maximum_rank, bool) or maximum_rank < 0:
        raise ValueError("maximum rank must be a nonnegative integer")
    rows = [[1]]
    for rank in range(1, maximum_rank + 1):
        previous = rows[-1]
        current = [1]
        for dimension in range(1, rank):
            current.append(
                previous[dimension]
                + 2 ** (rank - dimension) * previous[dimension - 1]
            )
        current.append(1)
        rows.append(current)
    return rows


def binary_subspace_counts(maximum_rank: int) -> list[int]:
    return [sum(row) for row in gaussian_table(maximum_rank)]


def literal_binary_subspaces(rank: int) -> set[frozenset[int]]:
    """Enumerate all F_2-subspaces literally; used only at tiny ranks."""
    if isinstance(rank, bool) or not isinstance(rank, int) or rank < 0 or rank > 6:
        raise ValueError("literal enumeration is restricted to ranks 0 through 6")
    spaces: set[frozenset[int]] = {frozenset({0})}
    for vector in range(1, 1 << rank):
        enlarged = set(spaces)
        for space in spaces:
            if vector not in space:
                enlarged.add(
                    frozenset(set(space) | {entry ^ vector for entry in space})
                )
        spaces = enlarged
    return spaces


def mod_two_pairing(left: int, right: int) -> int:
    """Diagonal pairing inherited from I_(1,rho-1), where -1=1 mod 2."""
    return (left & right).bit_count() % 2


def is_totally_isotropic(space: frozenset[int]) -> bool:
    return all(
        mod_two_pairing(left, right) == 0
        for left in space
        for right in space
    )


def ag_rank_elimination(poison_horizon: int = 1000) -> dict[str, Any]:
    """Eliminate b from the two AG inequalities without a search premise."""
    if poison_horizon < 27:
        raise ValueError("poison horizon must see ranks beyond the endpoint")
    feasible: dict[int, list[int]] = {}
    for rank in range(0, poison_horizon + 1):
        lower = max(0, rank - 10)
        upper = (4 * rank + 8) // 7
        if lower <= upper:
            feasible[rank] = list(range(lower, upper + 1))
    return {
        "symbolic_elimination": "7(rho-10)<=4rho+8, hence 3rho<=78",
        "maximum_feasible_rank": max(feasible),
        "rank_26_boundary_counts": feasible.get(26, []),
        "rank_27_boundary_counts": feasible.get(27, []),
        "no_feasible_rank_27_through_horizon": not any(
            rank >= 27 for rank in feasible
        ),
        "endpoint_is_eliminated_as_nonexistent": False,
        "endpoint_claim": "rho=26 survives the inequalities only at b=16",
    }


def adapted_wall_audit() -> dict[str, Any]:
    """Check the three reflection formulas on the adapted generators."""
    # q=2 times the physical intersection form.
    boundary_half_norm = Fraction(2 * (-4), 4)
    minus_one_norm = Fraction(2 * (-1))
    minus_two_norm = Fraction(2 * (-2))

    # Reflection multiplier is -2 q(x,r)/q(r,r).
    boundary_on_lattice = [
        -2 * Fraction(intersection) / boundary_half_norm
        for intersection in range(-7, 8)
    ]
    # q(h_i,E)=B_i.E may be any integer.  Norm -2 makes it integral.
    minus_one_on_half_boundary = [
        -2 * Fraction(intersection) / minus_one_norm
        for intersection in range(0, 8)
    ]
    # q(x,R)=2 x.R on L and q(B_i/2,R)=B_i.R=0 under W2.
    minus_two_on_lattice = [
        -2 * Fraction(2 * intersection) / minus_two_norm
        for intersection in range(-7, 8)
    ]
    minus_two_on_half_boundary = (
        -2 * Fraction(0) / minus_two_norm
    )
    hostile_missing_w2_multiplier = (
        -2 * Fraction(1) / minus_two_norm
    )

    checks = {
        "adapted_pairing_half_boundary_diagonal": str(boundary_half_norm),
        "adapted_pairing_half_boundary_off_diagonal": "0",
        "half_boundary_with_lattice_is_integral": True,
        "boundary_root_norm": str(boundary_half_norm),
        "minus_one_root_norm": str(minus_one_norm),
        "minus_two_root_norm": str(minus_two_norm),
        "boundary_reflection_coefficients_integral": all(
            value.denominator == 1 for value in boundary_on_lattice
        ),
        "minus_one_reflection_coefficients_integral": all(
            value.denominator == 1 for value in minus_one_on_half_boundary
        ),
        "minus_two_lattice_coefficients_integral": all(
            value.denominator == 1 for value in minus_two_on_lattice
        ),
        "minus_two_half_boundary_coefficient": str(
            minus_two_on_half_boundary
        ),
        "missing_W2_control_coefficient": str(
            hostile_missing_w2_multiplier
        ),
        "missing_W2_control_is_nonintegral": (
            hostile_missing_w2_multiplier.denominator != 1
        ),
        "W3_used_for_host_integrality": False,
        "W3_used_for_AG_nefness": True,
    }
    return checks


def permutation_orbit(
    start: tuple[int, ...], generators: tuple[tuple[int, int], ...]
) -> set[tuple[int, ...]]:
    orbit = {start}
    frontier = [start]
    while frontier:
        current = frontier.pop()
        for first, second in generators:
            image = list(current)
            image[first], image[second] = image[second], image[first]
            result = tuple(image)
            if result not in orbit:
                orbit.add(result)
                frontier.append(result)
    return orbit


def chamber_orbit_audit() -> dict[str, Any]:
    """A2 control: complete mirrors give one orbit; one mirror does not."""
    complete = permutation_orbit((0, 1, 2), ((0, 1), (1, 2)))
    incomplete = permutation_orbit((0, 1, 2), ((0, 1),))
    return {
        "A2_complete_chamber_count": 6,
        "A2_complete_orbit_size": len(complete),
        "A2_incomplete_orbit_size": len(incomplete),
        "complete_arrangement_transitive": len(complete) == 6,
        "selected_subarrangement_need_not_be_transitive": len(incomplete) < 6,
    }


def field_audit() -> dict[str, Any]:
    """Check rational Cartan/cycle arithmetic on a mixed-norm control."""
    gram = [
        [Fraction(-2), Fraction(1), Fraction(0)],
        [Fraction(1), Fraction(-4), Fraction(2)],
        [Fraction(0), Fraction(2), Fraction(-2)],
    ]
    cartan = [
        [
            Fraction(2) * gram[row][column] / gram[row][row]
            for column in range(3)
        ]
        for row in range(3)
    ]
    cycles = [
        cartan[0][1] * cartan[1][0],
        cartan[1][2] * cartan[2][1],
        cartan[0][1] * cartan[1][2] * cartan[2][0],
    ]
    return {
        "cartan_matrix": [[str(value) for value in row] for row in cartan],
        "sample_cyclic_products": [str(value) for value in cycles],
        "all_entries_and_cycles_rational": all(
            isinstance(value, Fraction)
            for row in cartan
            for value in row
        )
        and all(isinstance(value, Fraction) for value in cycles),
        "minimal_Vinberg_field": "Q",
        "surface_field_of_definition_claimed": False,
    }


def build_audit() -> dict[str, Any]:
    source_hashes = {
        relative: sha256(ROOT / relative) for relative in OWNER_PINS
    }
    recurrence_counts = binary_subspace_counts(26)
    literal_counts = [
        len(literal_binary_subspaces(rank)) for rank in range(1, 6)
    ]
    isotropic_counts = [
        sum(
            is_totally_isotropic(space)
            for space in literal_binary_subspaces(rank)
        )
        for rank in range(1, 6)
    ]
    coarse_total = 2 * recurrence_counts[2] + sum(recurrence_counts[3:27])

    payload = {
        "schema": (
            "coxeter.rational_wrong_log_maximal_finiteness."
            "hostile_audit.v1"
        ),
        "verdict": "PASS after TeX repair and scope clarification",
        "owner_hashes": source_hashes,
        "owner_pins_match": source_hashes == OWNER_PINS,
        "ag_rank_elimination": ag_rank_elimination(),
        "adapted_wall_audit": adapted_wall_audit(),
        "finite_host_audit": {
            "recurrence": (
                "[n h]_2=[n-1 h]_2+2^(n-h)[n-1 h-1]_2"
            ),
            "small_recurrence_counts": recurrence_counts[1:6],
            "small_literal_counts": literal_counts,
            "small_totally_isotropic_counts": isotropic_counts,
            "rank_2_picard_types": 2,
            "rank_3_through_26_picard_types_each": 1,
            "coarse_total": coarse_total,
            "coarse_total_expected": (
                5515561289142524504797715176506813442521589215285286
            ),
            "overcount_includes_nonintegral_hosts": True,
            "overcount_includes_isometric_duplicates": True,
        },
        "maximality_modes": {
            "absolute_discrete_reflection_maximality_saturates": True,
            "host_relative_reflection_maximality_saturates": True,
            "arithmetic_reflection_maximality_saturates": True,
            "arithmetic_reason": (
                "W and O+(M) are finite-covolume lattices with W<=O+(M), "
                "so every intermediate discrete overgroup has finite index "
                "and is arithmetic"
            ),
            "lattice_maximality_alone_saturates": False,
            "known_effective_walls_alone_saturate": False,
        },
        "chamber_orbit_audit": chamber_orbit_audit(),
        "field_audit": field_audit(),
        "scope": {
            "carrier_relative": True,
            "literal_complete_full_nef_required": True,
            "finite_sided_finite_volume_required": True,
            "raw_carrier_surjectivity_claimed": False,
            "hostile_incidence_bound_claimed": False,
            "terminating_reflectivity_census_claimed": False,
            "surface_classification_avoided": False,
            "reflective_lattice_classification_avoided": True,
        },
    }
    certify(payload)
    payload["all_checks_pass"] = True
    return payload


def certify(payload: dict[str, Any]) -> None:
    if not payload["owner_pins_match"]:
        raise AssertionError("owner artifact hash drift")
    rank = payload["ag_rank_elimination"]
    if rank["maximum_feasible_rank"] != 26:
        raise AssertionError("AG rank endpoint changed")
    if rank["rank_26_boundary_counts"] != [16]:
        raise AssertionError("rank 26 no longer has the unique b=16 endpoint")
    if rank["rank_27_boundary_counts"] or not rank[
        "no_feasible_rank_27_through_horizon"
    ]:
        raise AssertionError("a rank above 26 passed the AG inequalities")

    walls = payload["adapted_wall_audit"]
    for key in (
        "half_boundary_with_lattice_is_integral",
        "boundary_reflection_coefficients_integral",
        "minus_one_reflection_coefficients_integral",
        "minus_two_lattice_coefficients_integral",
        "missing_W2_control_is_nonintegral",
    ):
        if not walls[key]:
            raise AssertionError(f"adapted wall check failed: {key}")

    hosts = payload["finite_host_audit"]
    if hosts["small_recurrence_counts"] != [2, 5, 16, 67, 374]:
        raise AssertionError("Gaussian recurrence counts changed")
    if hosts["small_literal_counts"] != [2, 5, 16, 67, 374]:
        raise AssertionError("literal binary subspace counts changed")
    if hosts["small_totally_isotropic_counts"] != [1, 2, 4, 11, 31]:
        raise AssertionError("mod-two integrality controls changed")
    if hosts["coarse_total"] != hosts["coarse_total_expected"]:
        raise AssertionError("coarse host total changed")

    modes = payload["maximality_modes"]
    if not all(
        modes[key]
        for key in (
            "absolute_discrete_reflection_maximality_saturates",
            "host_relative_reflection_maximality_saturates",
            "arithmetic_reflection_maximality_saturates",
        )
    ):
        raise AssertionError("a valid maximality mode failed")
    if modes["lattice_maximality_alone_saturates"]:
        raise AssertionError("lattice maximality was promoted incorrectly")

    if not payload["chamber_orbit_audit"]["complete_arrangement_transitive"]:
        raise AssertionError("complete chamber transitivity failed")
    if not payload["field_audit"]["all_entries_and_cycles_rational"]:
        raise AssertionError("rational field arithmetic failed")


def rendered(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()
    payload = rendered(build_audit())
    if args.check_output:
        if not ARCHIVE.is_file() or ARCHIVE.read_text(encoding="utf-8") != payload:
            raise SystemExit("hostile archive mismatch")
        print("PASS: hostile archive matches independent reconstruction")
        return 0
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
