from __future__ import annotations

import pytest

from scripts.certify_silver_k3_norm2_infinite_stabilizer import verify_frozen
from scripts.independent_silver_k3_norm2_infinite_stabilizer import (
    EXPECTED_SUMMARY_SHA256,
    verify_frozen as verify_hostile_frozen,
)


EXPECTED_CERTIFICATE_SHA256 = (
    "9b6fb3a41affed6e6b43b967c89d124daeab8ba5f4d0fe05a0c7ff6ccf2b9b12"
)


@pytest.fixture(scope="module")
def production_document() -> dict[str, object]:
    return verify_frozen()


def test_production_infinite_stabilizer_certificate_is_frozen(
    production_document: dict[str, object],
) -> None:
    document = production_document

    assert document["verified"]
    assert all(document["checks"].values())
    assert document["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert document["lattice"]["field_rank"] == 10
    assert document["lattice"]["field_determinant"] == {"a": "-2", "b": "-2"}
    assert document["lattice"]["identity_inertia_positive_negative_zero"] == [9, 1, 0]
    assert document["lattice"]["conjugate_inertia_positive_negative_zero"] == [10, 0, 0]
    assert document["height_zero"]["total_field_norm_two_roots"] == 82
    assert document["height_zero"]["reflection_closure_matches_exhaustion"]


def test_both_complete_chamber_symmetries_and_product_are_exact(
    production_document: dict[str, object],
) -> None:
    document = production_document

    seed = document["seed_C5_reversal"]
    swap = document["wall_4_18_swap"]
    product = document["infinite_order_product"]
    assert seed["matrix_sha256"] == (
        "bf8b56c5d8faeec170f2bb863a3526a1d525353a9c278fd3ee661ec7f5282a37"
    )
    assert swap["matrix_sha256"] == (
        "c754c65f8ea3f64d76163bdb9d077c292292329ecdf05b0e56f59990c9e419cb"
    )
    assert product["matrix_sha256"] == (
        "f8672080597bbb10daf78bc893ffe4e17a763cb0e04df18d114ddce22c22db43"
    )
    assert seed["future_cone_preserved"]
    assert swap["future_cone_preserved"]
    assert seed["separator_hypothesis"].endswith("including equality")
    assert swap["separator_hypothesis"].endswith("including equality")
    assert swap["boundary_minimizer_zero_wall_indices"] == [0, 1, 2, 3, 5, 6, 7]
    assert product["trace"] == {"a": "10", "b": "4"}
    assert product["characteristic_polynomial_factorization"] == (
        "(t-1)^7*(t+1)*(t^2-(4+4sqrt(2))*t+1)"
    )


def test_theorem_boundary_marks_noncomputational_inputs(
    production_document: dict[str, object],
) -> None:
    document = production_document
    conclusion = document["arithmetic_conclusion"]

    assert conclusion["anisotropic_over_Q_sqrt2_derived_from_conjugate_definiteness"]
    assert conclusion["complete_chamber_stabilizer"] == "infinite"
    assert conclusion["reflection_group_index"] == "infinite"
    assert conclusion["complete_chamber_hyperbolic_volume"] == "infinite"
    assert len(conclusion["citation_dependent_inputs_not_reproved_by_this_certificate"]) == 4
    assert "does not computationally reprove" in conclusion["computational_certificate_scope"]


def test_standard_library_hostile_replay_is_frozen() -> None:
    replay = verify_hostile_frozen()

    assert replay["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert replay["field_form_identity_signature"] == [9, 1, 0]
    assert replay["field_form_conjugate_signature"] == [10, 0, 0]
    assert replay["field_anisotropic_over_Qsqrt2_derived"]
    assert replay["height_zero_root_count"] == 82
    assert replay["height_zero_closure_exhaustive"]
    assert replay["seed_preserves_future_component"]
    assert replay["swap_preserves_future_component"]
    assert replay["separator_hypothesis_includes_equality"]
    assert replay["complete_chamber_symmetries"]
    assert replay["infinite_order_product"]
