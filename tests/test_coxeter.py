from fractions import Fraction

import pytest

from coxeter_polytopes.exact.coxeter import classify_root_pair
from coxeter_polytopes.exact.golden import Golden, PHI


@pytest.mark.parametrize(
    ("norm_i", "inner", "norm_j", "label"),
    [
        (2, 0, 2, 2),
        (2, -1, 2, 3),
        (2, -1, 1, 4),
        (2, -PHI, 2, 5),
        (3, Fraction(-3, 2), 1, 6),
        (PHI + 2, -(PHI + 2), 4, 10),
    ],
)
def test_recognized_finite_labels(norm_i, inner, norm_j, label) -> None:
    relation = classify_root_pair(norm_i, inner, norm_j)
    assert relation.kind == "finite"
    assert relation.label == label


def test_parallel_ultraparallel_and_unrecognized() -> None:
    assert classify_root_pair(1, -1, 1).kind == "parallel"
    assert classify_root_pair(1, -2, 1).kind == "ultraparallel"
    relation = classify_root_pair(2, Golden(-1, -1), 2)
    assert relation.q == 2 + 3 * PHI
    assert relation.kind == "ultraparallel"
    unrecognized = classify_root_pair(8, -1, 1)
    assert unrecognized.kind == "finite_unrecognized"
    assert unrecognized.label is None


def test_invalid_roots_are_rejected() -> None:
    with pytest.raises(ValueError):
        classify_root_pair(0, 0, 1)
    with pytest.raises(ValueError):
        classify_root_pair(1, 1, 1)
