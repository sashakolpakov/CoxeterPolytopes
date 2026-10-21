from fractions import Fraction

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.geometric.silver_k3_carrier import (
    audit_silver_c5_saturation,
    audit_silver_k3_embedding,
    silver_carrier_model_audits,
    silver_finite_completion_audits,
    silver_k3_carrier_theorem,
)


def test_c5_has_one_nonzero_even_sqrt2_stable_index_two_class() -> None:
    seed = audit_silver_c5_saturation()

    assert seed.raw_trace_determinant == -4
    assert seed.raw_trace_inertia == Inertia(1, 9, 0)
    assert seed.raw_smith_invariants == (2, 2)
    assert [
        (item.bits, item.square, item.is_even, item.sqrt2_stable)
        for item in seed.order_two_classes
    ] == [
        ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0), Fraction(0), True, True),
        ((0, 0, 0, 1, 0, 1, 0, 0, 0, 1), Fraction(-2), True, True),
        ((0, 1, 1, 0, 1, 0, 0, 1, 1, 0), Fraction(-1), False, False),
        ((0, 1, 1, 1, 1, 1, 0, 1, 1, 1), Fraction(0), True, False),
    ]
    assert [item.bits for item in seed.stable_nonzero_even_classes] == [
        (0, 0, 0, 1, 0, 1, 0, 0, 0, 1)
    ]


def test_c5_saturation_is_e8_plus_u_and_exposes_the_extra_root() -> None:
    seed = audit_silver_c5_saturation()

    assert seed.standard_trace_type == "E8(-1) + U"
    assert seed.saturated_trace_determinant == -1
    assert seed.saturated_trace_inertia == Inertia(1, 9, 0)
    assert seed.saturated_smith_invariants == ()
    assert seed.trace_isometry_determinant == 1
    assert seed.trace_isometry_verified is True
    assert seed.extra_root_square == -2
    assert seed.extra_root_sqrt2_pairing == -1
    assert seed.sqrt2_extra_root == (0, 0, 1, 0, 1, 0, 0, 0, 1, 0)
    assert seed.sqrt2_extra_root_square == -4
    assert seed.extra_root_field_norm == Silver(3, -1)
    assert seed.passes is True


def test_two_rank_five_silver_completions_have_explicit_trace_types() -> None:
    completions = silver_finite_completion_audits()

    assert [item.name for item in completions] == ["B5", "F4 + A1"]
    assert [item.trace_type for item in completions] == [
        "D10(-1)",
        "E8(-1) + A1(-1)^2",
    ]
    assert [item.field_determinant for item in completions] == [Silver(2), Silver(2)]
    assert all(item.trace_determinant == 4 for item in completions)
    assert all(item.trace_smith_invariants == (2, 2) for item in completions)
    assert all(item.coparity == 1 for item in completions)
    assert all(abs(item.trace_isometry_determinant) == 1 for item in completions)
    assert all(item.trace_isometry_verified for item in completions)
    assert all(item.passes for item in completions)


def test_both_completed_carriers_are_the_20_2_1_lattice() -> None:
    models = silver_carrier_model_audits()

    assert [item.trace_type for item in models] == [
        "U + E8(-1) + D10(-1)",
        "U + E8(-1)^2 + A1(-1)^2",
    ]
    assert all(item.field_determinant == Silver(-4, -4) for item in models)
    assert all(item.field_identity_inertia == Inertia(9, 1, 0) for item in models)
    assert all(item.field_conjugate_inertia == Inertia(10, 0, 0) for item in models)
    assert all(item.trace_determinant == -4 for item in models)
    assert all(item.trace_inertia == Inertia(1, 19, 0) for item in models)
    assert all(item.trace_smith_invariants == (2, 2) for item in models)
    assert all(item.two_elementary_invariants == (20, 2, 1) for item in models)
    assert all(item.sqrt2_action_self_adjoint for item in models)
    assert all(item.sqrt2_action_squares_to_two for item in models)
    assert all(item.passes for item in models)


def test_alexeev_aligned_carrier_embeds_primitively_with_t_a1_squared() -> None:
    embedding = audit_silver_k3_embedding()

    assert embedding.transcendental_gram == ((2, 0), (0, 2))
    assert embedding.discriminant_anti_isometry_squares == (
        (Fraction(-1, 2), Fraction(1, 2)),
        (Fraction(-1, 2), Fraction(1, 2)),
    )
    assert embedding.gluing_index == 4
    assert embedding.ambient_determinant == -1
    assert embedding.ambient_inertia == Inertia(3, 19, 0)
    assert embedding.ambient_even is True
    assert embedding.ambient_unimodular is True
    assert embedding.orthogonal_complement_verified is True
    assert embedding.primitive_embedding_verified is True
    assert embedding.passes is True


def test_silver_carrier_theorem_retains_the_wall_closure_boundary() -> None:
    theorem = silver_k3_carrier_theorem()

    assert theorem.alexeev_invariants == (20, 2, 1)
    assert theorem.alexeev_lattice_type == "U + E8(-1)^2 + A1(-1)^2"
    assert theorem.wall_closure_required is True
    assert theorem.passes is True
