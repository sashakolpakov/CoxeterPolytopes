import pytest

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.geometric.golden_finite import (
    forced_golden_core,
    golden_finite_trace_lift,
    h4_lanner_extensions,
)


@pytest.mark.parametrize(
    ("field_type", "trace_type", "field_det", "trace_det"),
    [
        ("H2", "A4(-1)", Golden(3, -1), 5),
        ("H3", "D6(-1)", Golden(4, -2), 4),
        ("H4", "E8(-1)", Golden(5, -3), 1),
    ],
)
def test_golden_finite_trace_lift_is_exact_integral_isometry(
    field_type: str,
    trace_type: str,
    field_det: Golden,
    trace_det: int,
) -> None:
    lift = golden_finite_trace_lift(field_type)

    assert lift.trace_type == trace_type
    assert lift.field_determinant == field_det
    assert lift.trace_determinant == trace_det
    assert lift.field_identity_inertia == Inertia(lift.field_rank, 0, 0)
    assert lift.field_conjugate_inertia == Inertia(lift.field_rank, 0, 0)
    assert lift.trace_inertia == Inertia(0, lift.trace_rank, 0)
    assert abs(lift.isometry_determinant) == 1
    assert lift.isometry_verified is True
    assert lift.passes is True


def test_h4_trace_core_is_the_unimodular_splitting_case() -> None:
    assert golden_finite_trace_lift(
        "H4"
    ).orthogonally_splits_every_integral_overlattice is True
    assert golden_finite_trace_lift(
        "H3"
    ).orthogonally_splits_every_integral_overlattice is False


def test_forced_missing_face_cores() -> None:
    order_four = forced_golden_core(4)
    order_five = forced_golden_core(5)

    assert (order_four.field_type, order_four.trace_type) == ("H3", "D6(-1)")
    assert (order_five.field_type, order_five.trace_type) == ("H4", "E8(-1)")
    assert order_four.orthogonally_splits is False
    assert order_five.orthogonally_splits is True

    with pytest.raises(ValueError):
        forced_golden_core(3)


def test_h4_has_exactly_three_golden_lanner_extensions() -> None:
    extensions = h4_lanner_extensions()

    assert [extension.attachment_labels for extension in extensions] == [
        (2, 2, 2, 3),
        (2, 2, 2, 5),
        (2, 2, 3, 2),
    ]
    assert [extension.field_determinant for extension in extensions] == [
        Golden(6, -4),
        Golden(8, -6),
        Golden(4, -4),
    ]
    assert [extension.complement_norm for extension in extensions] == [
        Golden(0, -2),
        Golden(-2, -6),
        Golden(-4, -8),
    ]
    assert [extension.trace_determinant for extension in extensions] == [
        -4,
        -20,
        -16,
    ]
    assert all(
        extension.trace_inertia == Inertia(1, 9, 0)
        for extension in extensions
    )
    assert [extension.visible_curve_type for extension in extensions] == [
        "A5(-1)",
        "A5(-1)",
        "D5(-1)",
    ]
