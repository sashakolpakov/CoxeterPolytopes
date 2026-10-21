from fractions import Fraction

from coxeter_polytopes.classification.pell_boost_nonreflectivity import (
    PellBoostCriterion,
)
from coxeter_polytopes.exact.silver import Silver


def test_lambda22_criterion_constants() -> None:
    criterion = PellBoostCriterion(
        coefficient=Silver(-1, 1),
        direction_norm=4,
        pell_A=Silver(5, 4),
        pell_U=Silver(4, 3),
        witness_scale=Fraction(1, 3),
    )
    assert criterion.audit()
    assert criterion.separator_height_bound == Silver(56, 40)
    assert criterion.control_square_bound == Silver(136, 96)
    assert criterion.witness_separator_coefficient == Silver(10, -7) / 6


def test_odd_unimodular_criterion_constants() -> None:
    criterion = PellBoostCriterion(
        coefficient=Silver(-1, 2),
        direction_norm=3,
        pell_A=Silver(16, 12),
        pell_U=Silver(7, 5),
        witness_scale=Fraction(43, 55),
    )
    assert criterion.audit()
    assert criterion.time_from_direction == Silver(7, 5)
    assert criterion.spatial_from_time == Silver(13, 9)
    assert criterion.direction_from_direction == Silver(5, 4)
    assert criterion.separator_direction_slope == Silver(8, -4)
    assert criterion.control_square_bound == Silver(297, 210)
    assert criterion.separator_height_bound == Silver(543, 384)
    assert criterion.witness_separator_coefficient == Silver(399, -282) / 55
