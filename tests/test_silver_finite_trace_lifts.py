from fractions import Fraction

import pytest

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import (
    SILVER_K3_TRACE_SCALAR,
    SQRT2,
    Silver,
)
from coxeter_polytopes.geometric.silver_finite import (
    forced_silver_core,
    silver_finite_trace_lift,
    silver_lanner_blocks,
    silver_lanner_census,
)


def test_silver_k3_trace_scalar_is_totally_negative_codifferent() -> None:
    alpha = SILVER_K3_TRACE_SCALAR

    assert alpha == Silver(Fraction(-1, 2), Fraction(-1, 4))
    assert alpha.trace() == -1
    assert alpha.sign(0) == alpha.sign(1) == -1
    assert 2 * SQRT2 * alpha == -1 - SQRT2
    assert (2 * SQRT2 * alpha).is_integral is True


@pytest.mark.parametrize(
    ("field_type", "trace_type", "field_det", "trace_det"),
    [
        ("I2(4)", "D4(-1)", Silver(2), 4),
        ("B3", "D6(-1)", Silver(2), 4),
        ("F4", "E8(-1)", Silver(1), 1),
    ],
)
def test_silver_finite_trace_lift_has_explicit_integral_isometry(
    field_type: str,
    trace_type: str,
    field_det: Silver,
    trace_det: int,
) -> None:
    lift = silver_finite_trace_lift(field_type)

    assert lift.trace_type == trace_type
    assert lift.field_determinant == field_det
    assert lift.trace_determinant == trace_det
    assert lift.field_identity_inertia == Inertia(lift.field_rank, 0, 0)
    assert lift.field_conjugate_inertia == Inertia(lift.field_rank, 0, 0)
    assert lift.trace_inertia == Inertia(0, lift.trace_rank, 0)
    assert abs(lift.isometry_determinant) == 1
    assert lift.isometry_verified is True
    assert lift.passes is True


def test_forced_silver_missing_face_cores_match_golden_trace_cores() -> None:
    order_four = forced_silver_core(4)
    order_five = forced_silver_core(5)

    assert (order_four.field_type, order_four.trace_type) == ("B3", "D6(-1)")
    assert (order_five.field_type, order_five.trace_type) == ("F4", "E8(-1)")
    assert order_four.orthogonally_splits is False
    assert order_five.orthogonally_splits is True

    with pytest.raises(ValueError):
        forced_silver_core(3)


def test_silver_lanner_census_is_exactly_the_one_four_cycles() -> None:
    census = silver_lanner_census()
    blocks = census.blocks

    assert census.passes is True
    assert [
        (
            stage.order,
            stage.graph_masks,
            stage.graph_prefilter_survivors,
            stage.decorated_candidates,
            stage.labelled_exact_survivors,
            stage.isomorphism_types,
        )
        for stage in census.stages
    ] == [
        (4, 64, 19, 176, 12, 1),
        (5, 1024, 137, 2384, 60, 1),
    ]
    assert [block.order for block in blocks] == [4, 5]
    assert [block.canonical_label_key.count(4) for block in blocks] == [1, 1]
    assert [block.field_determinant for block in blocks] == [
        Silver(-1, -2),
        Silver(-2, -2),
    ]
    assert [block.trace_determinant for block in blocks] == [-7, -4]
    assert [block.trace_type for block in blocks] == [
        "U + A6(-1)",
        "E8(-1) + U(2)",
    ]
    assert all(abs(block.trace_isometry_determinant) == 1 for block in blocks)
    assert all(block.trace_isometry_verified for block in blocks)
    assert all(block.passes for block in blocks)


def test_each_silver_lanner_forces_a_primitive_u2_fiber_pair() -> None:
    for block in silver_lanner_blocks():
        assert block.curve_inertia == Inertia(0, block.order - 1, 1)
        assert block.fiber_square == 0
        assert block.fiber_sqrt2_pairing == 2
        assert block.sqrt2_fiber_square == 0
        assert block.primitive_u2_gram == ((0, 2), (2, 0))
