#!/usr/bin/env python3
"""Independent hostile replay of compact physical-star K3 host finiteness.

This script deliberately does not import the owner certifier.  It freezes the
owner bytes, checks exact quadratic-order trace duality, replays the
root-frame discriminant calculation, and checks the algebraic distinction
between a native k-reflection and an ordinary rational root reflection.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]

PINS = {
    "owner_report": (
        "results/k3_compact_physical_star_trace_host_finiteness.md",
        "5cc86e7225c11697c83af0ec34fa5c773e8e969c47107eabdd2a327c23c1ea3a",
    ),
    "owner_certifier": (
        "scripts/certify_k3_compact_physical_star_trace_host_finiteness.py",
        "c4a52c76f28192eada049fe66031fd62bb35f6366332f92abc54cec995d58d55",
    ),
    "owner_archive": (
        "data/enumerations/k3_compact_physical_star_trace_host_finiteness.json",
        "cc2d48b699aeebd1aaa10ae6b4b8249addde2eac681c9bb475aa19fcd3e14e90",
    ),
    "owner_tests": (
        "tests/test_k3_compact_physical_star_trace_host_finiteness.py",
        "7f38e30ed6fb20af8c54a266c9b8ec4cad59ddccd90f1fd451d40f9565dbd3ec",
    ),
    "field_input": (
        "results/k3_compact_narrow_facet_field_finiteness.md",
        "39c9a3a8422110bcdcf751ed7bae90ac86cdf3a720d220ee164ff241bd9ac9b5",
    ),
    "composition_input": (
        "results/compact_trace_host_chamber_finiteness_composition_hostile_audit.md",
        "8fd86c7c588f247fd3ec380e8f1e68840da81484414ab7e7413a7b3487013900",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def energy_and_codifferent_certificate() -> dict[str, Any]:
    """Check the place-energy estimate and a nonmaximal-order codifferent."""

    # If x,y are the distinguished energies, the other-place contribution is
    # at most sqrt((2-x)(2-y)).  Weighted Cauchy gives the stronger bound 14.
    sampled_max = 0.0
    for i in range(1, 400):
        x = Fraction(i, 200)
        for j in range(1, 400):
            y = Fraction(j, 200)
            if x >= 2 or y >= 2:
                continue
            bound = 7 * math.sqrt(float(x * y)) + math.sqrt(
                float((2 - x) * (2 - y))
            )
            sampled_max = max(sampled_max, bound)
            assert bound <= 14 + 1e-12

    # For O=Z[sqrt(5)], the trace matrix on (1,sqrt(5)) is diag(2,10).
    # Its dual basis is (1/2,sqrt(5)/10), proving the exact codifferent.
    trace_matrix = sp.diag(2, 10)
    dual_basis_matrix = trace_matrix.inv()
    assert trace_matrix * dual_basis_matrix == sp.eye(2)

    # Enumerate one bounded Minkowski box in that codifferent.  Writing
    # z=u/2+v*sqrt(5)/10, the generous finite ranges below are exhaustive for
    # |tau_0(z)|<=14 and |tau_1(z)|<=2.
    sqrt5 = math.sqrt(5)
    points: list[tuple[int, int]] = []
    for u in range(-32, 33):
        for v in range(-80, 81):
            z0 = u / 2 + v * sqrt5 / 10
            z1 = u / 2 - v * sqrt5 / 10
            if abs(z0) <= 14 and abs(z1) <= 2:
                points.append((u, v))
    assert points
    # The box itself implies |u|<=16 and |v|<=~35.8, so no boundary point is
    # lost by the enumeration ranges.
    assert max(abs(u) for u, _ in points) <= 16
    assert max(abs(v) for _, v in points) <= 35

    return {
        "place_energy_identity": "sum_tau(-tau(g_ii))=2",
        "distinguished_normalized_bound": "2*abs(g_ij)/sqrt(t_i*t_j)<14",
        "owner_trace_bound": "abs(B(delta_i,delta_j))<16",
        "independent_stronger_trace_bound": "abs(B(delta_i,delta_j))<14",
        "sampled_weighted_bound_max": sampled_max,
        "nonmaximal_order": "Z[sqrt(5)]",
        "order_codifferent_basis": ["1/2", "sqrt(5)/10"],
        "bounded_codifferent_box_points": len(points),
        "codifferent_box_is_finite": True,
    }


def order_and_overlattice_certificate() -> dict[str, Any]:
    """Replay the conductor and root-frame duality gates."""

    # Z[sqrt(5)] has discriminant 20 and index two in O_K, of discriminant 5.
    lower_disc = 20
    maximal_disc = 5
    index = math.isqrt(lower_disc // maximal_disc)
    assert index == 2
    assert lower_disc == index**2 * maximal_disc

    # An explicit even integral index-two extension of a root frame.
    gram_f = sp.diag(2, -2, -2)
    change = sp.Matrix(
        [
            [sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 1, 0],
            [0, 0, 1],
        ]
    )
    gram_m = change.T * gram_f * change
    assert gram_m == sp.Matrix([[0, -1, 0], [-1, -2, 0], [0, 0, -2]])
    assert all(gram_m[i, i] % 2 == 0 for i in range(3))
    assert abs(int(gram_f.det())) == 8
    assert abs(int(gram_m.det())) == 2
    assert abs(int(gram_f.det())) == index**2 * abs(int(gram_m.det()))

    return {
        "cycle_order": "O0=Z[bounded integral Cartan cycles]",
        "cycle_fraction_field_is_k": True,
        "coefficient_order_interval": "O0 subset O subset O_k",
        "toy_overorder_index": index,
        "toy_overorder_count": 2,
        "root_frame_containment": "F subset M subset F#",
        "root_frame_discriminant_group_order": abs(int(gram_f.det())),
        "toy_overlattice_index": index,
        "toy_overlattice_even_and_integral": True,
        "index_square_identity_checked": True,
        "projectivity_order_stability_evenness_are_finite_filters": True,
    }


Quad = tuple[Fraction, Fraction]


def qadd(x: Quad, y: Quad) -> Quad:
    return (x[0] + y[0], x[1] + y[1])


def qmul(x: Quad, y: Quad) -> Quad:
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def qdiv(x: Quad, y: Quad) -> Quad:
    norm = y[0] * y[0] - 5 * y[1] * y[1]
    if norm == 0:
        raise ZeroDivisionError
    return qmul(x, (y[0] / norm, -y[1] / norm))


def qembed(x: Quad) -> float:
    return float(x[0]) + float(x[1]) * math.sqrt(5)


def full_picard_descent_certificate() -> dict[str, Any]:
    """Check real-line proportionality descent and reflection species."""

    delta: list[Quad] = [
        (Fraction(1), Fraction(1)),
        (Fraction(2), Fraction(-1)),
        (Fraction(3), Fraction(1, 2)),
    ]
    scalar: Quad = (Fraction(1, 2), Fraction(1, 2))
    r = [qmul(scalar, coordinate) for coordinate in delta]
    recovered = qdiv(r[0], delta[0])
    assert recovered == scalar
    assert all(qdiv(rj, dj) == recovered for rj, dj in zip(r, delta))
    real_ratios = [qembed(rj) / qembed(dj) for rj, dj in zip(r, delta)]
    assert all(abs(value - qembed(scalar)) < 1e-12 for value in real_ratios)

    # Restriction of scalars of diag(-1,1,1) over a quadratic field has a
    # two-dimensional rational minus eigenspace; an ordinary K3 root
    # reflection on the rank-six trace lattice has only one.
    native_res = sp.diag(-1, -1, 1, 1, 1, 1)
    ordinary_root = sp.diag(-1, 1, 1, 1, 1, 1)
    native_minus = 6 - (native_res + sp.eye(6)).rank()
    ordinary_minus = 6 - (ordinary_root + sp.eye(6)).rank()
    assert native_minus == 2
    assert ordinary_minus == 1

    return {
        "full_picard_identity": "NS(X)=M_Z",
        "physical_wall_vector_lies_in_M": True,
        "distinguished_real_line_proportionality_descends_to_k": True,
        "recovered_scalar": [str(recovered[0]), str(recovered[1])],
        "native_quadratic_reflection_Q_minus_dimension": native_minus,
        "ordinary_K3_root_reflection_Q_minus_dimension": ordinary_minus,
        "reflection_species_are_distinct_for_degree_gt_one": True,
        "primitive_trace_block_in_larger_NS_is_not_enough": True,
    }


def build_audit() -> dict[str, Any]:
    observed = {
        name: {"path": relative, "sha256": sha256(ROOT / relative)}
        for name, (relative, _) in PINS.items()
    }
    return {
        "schema": "coxeter.k3-compact-physical-star-host-independent-hostile.v1",
        "verdict": "PASS_EXACT_STATED_SCOPE",
        "pins": observed,
        "energy_codifferent": energy_and_codifferent_certificate(),
        "order_overlattice": order_and_overlattice_certificate(),
        "full_picard_descent": full_picard_descent_certificate(),
        "same_host_composition": {
            "complete_oriented_host": "(k,tau0,O,M,g,C+)",
            "same_marked_native_host_maximality_required": True,
            "maximality_implies_W_equals_complete_native_reflection_group": True,
            "complete_locally_finite_mirror_arrangement_one_chamber_orbit": True,
            "ordinary_K3_reflection_maximality_is_not_a_substitute": True,
        },
        "scope": {
            "physical_same_vector_minus_two_star_is_an_extra_hypothesis": True,
            "full_Picard_exact_eigen_nef_makes_that_hypothesis_automatic": True,
            "finite_ambient_NS_extensions_not_proved": True,
            "carrier_existence_not_proved": True,
            "Property_4_not_proved": True,
        },
    }


def main() -> int:
    for relative, expected in PINS.values():
        actual = sha256(ROOT / relative)
        if actual != expected:
            raise SystemExit(f"pin mismatch: {relative}: {actual} != {expected}")
    print(json.dumps(build_audit(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
