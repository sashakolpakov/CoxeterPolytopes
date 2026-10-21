from fractions import Fraction

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_k3_vinberg import (
    admissible_silver_control_coordinates,
    audit_silver_k3_height_zero,
    audit_silver_saturation_wall,
    run_silver_k3_vinberg_prefix,
)


def test_silver_k3_height_zero_norm_two_system_is_exhaustive() -> None:
    audit = audit_silver_k3_height_zero()

    assert audit.passes
    assert audit.trace_short_vector_counts == (240, 240, 4)
    assert audit.field_norm_two_counts == (32, 48, 2)
    assert audit.total_field_norm_two_roots == 82
    assert audit.reflection_closure_matches_exhaustion


def test_silver_k3_control_scan_and_prefix_are_exact_through_two() -> None:
    assert admissible_silver_control_coordinates(2) == (
        Silver(1),
        SQRT2,
        Silver(2),
    )

    prefix = run_silver_k3_vinberg_prefix()
    assert prefix.passes
    assert prefix.control_norm == 1 - SQRT2
    assert prefix.field_identity_inertia == Inertia(9, 1, 0)
    assert prefix.field_conjugate_inertia == Inertia(10, 0, 0)
    assert tuple(len(level.alphabet.values) for level in prefix.levels) == (4, 5, 4)
    assert tuple(level.exact_weight_solutions for level in prefix.levels) == (1, 5, 7)
    assert tuple(level.nonintegral_reconstructions for level in prefix.levels) == (0, 3, 2)
    assert tuple(level.previous_wall_rejections for level in prefix.levels) == (0, 1, 4)
    assert tuple(len(level.roots) for level in prefix.levels) == (1, 1, 1)
    assert len(prefix.roots) == 12


def test_silver_k3_prefix_does_not_cross_the_vinberg_stopping_gate() -> None:
    prefix = run_silver_k3_vinberg_prefix()

    assert prefix.wall_identity_inertia == Inertia(9, 1, 2)
    assert prefix.wall_conjugate_inertia == Inertia(10, 0, 2)
    assert prefix.diagram.elliptic_rank_eight_count == 69
    assert prefix.diagram.elliptic_rank_nine_count == 10
    assert prefix.diagram.finite_endpoint_histogram == ((1, 48), (2, 21))
    assert prefix.diagram.connected_affine_subdiagrams == ()
    assert prefix.diagram.compact_stopping_available is False
    assert prefix.diagram.finite_volume_stopping_available is False
    assert prefix.stopping_status == "undecidable_incomplete_simple_roots"


def test_saturation_minus_two_divisor_is_not_an_integral_field_reflection() -> None:
    audit = audit_silver_saturation_wall()

    assert audit.passes
    assert audit.field_norm == Silver(3, -1)
    assert audit.trace_square == -2
    assert audit.reflection_coefficient_on_e0 == Silver(
        Fraction(-4, 7), Fraction(-6, 7)
    )
    assert audit.reflection_is_integral is False
