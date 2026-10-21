from __future__ import annotations

import json

from coxeter_polytopes.exact.silver import Silver
from scripts.certify_silver_k3_726_independent_replay import (
    ARCHIVE,
    _determinant_nonzero_mod_prime,
    _finite_34,
    verify_archive,
)


def test_frozen_silver_726_independent_replay() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))

    assert verify_archive(document) == {
        "elliptic_rank_eight_count": 25_764,
        "elliptic_rank_nine_count": 3_604,
        "exact_inertia_check_count": 29_368,
        "finite_incidence_sum": 32_436,
        "connected_affine_count": 0,
        "maximum_affine_order_checked": 11,
        "finite_volume_open_edge_count": 16_760,
        "new_hyperbolic_example": False,
        "eventual_infinite_volume_proved": False,
    }
    enumeration = document["enumeration"]
    assert enumeration["elliptic_subsets_by_order"] == {
        "1": 726,
        "2": 9_098,
        "3": 42_057,
        "4": 103_017,
        "5": 151_756,
        "6": 140_105,
        "7": 79_842,
        "8": 25_764,
        "9": 3_604,
    }
    assert enumeration["modular_nonzero_determinant_witness_count"] == 39_043
    assert enumeration["exact_inertia_affine_fallback_count"] == 0
    assert enumeration["spectral_diagnostics"] == {
        "available": True,
        "closest_rejected_minimum_eigenvalue": -0.193993277713,
        "method": (
            "numpy.linalg.eigvalsh on each candidate extension of order at least three"
        ),
        "proof_role": (
            "none; finite type and affine determinant/inertia decisions are exact"
        ),
        "sign_mismatch_count": 0,
        "smallest_accepted_minimum_eigenvalue": 0.010956209263,
    }


def test_independent_finite_and_modular_affine_primitives() -> None:
    finite_a3 = tuple(
        tuple(
            Silver(2 if left == right else -1 if abs(left - right) == 1 else 0)
            for right in range(3)
        )
        for left in range(3)
    )
    affine_a2 = tuple(
        tuple(Silver(2 if left == right else -1) for right in range(3))
        for left in range(3)
    )

    assert _finite_34(finite_a3, (0, 1, 2)) is True
    assert _finite_34(affine_a2, (0, 1, 2)) is False
    assert _determinant_nonzero_mod_prime(affine_a2, (0, 1, 2), 5) is False
    assert _determinant_nonzero_mod_prime(affine_a2, (0, 1, 2), 13) is False
