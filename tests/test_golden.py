from fractions import Fraction

import pytest

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO


def test_golden_relation_and_arithmetic() -> None:
    assert PHI * PHI == PHI + 1
    assert (2 + PHI) * (3 - PHI) == 5
    assert hash(Golden(5)) == hash(5)
    assert PHI.conjugate() == 1 - PHI
    assert PHI.trace() == 1
    assert PHI.norm() == -1
    assert PHI * PHI.inverse() == ONE
    assert PHI**-2 == (PHI + 1).inverse()


def test_exact_sign_in_both_embeddings() -> None:
    assert ZERO.sign() == 0
    assert PHI.sign(0) == 1
    assert PHI.sign(1) == -1
    assert (2 - PHI).sign(0) == 1
    assert (2 - PHI).sign(1) == 1
    assert (Golden(Fraction(-3, 2), 1)).sign(0) == 1
    assert (Golden(Fraction(-3, 2), 1)).sign(1) == -1
    assert Golden(-2, 1).sign(0) == -1
    assert Golden(-2, 1).sign(1) == -1
    with pytest.raises(ValueError):
        PHI.sign(2)


def test_ordering_uses_identity_embedding() -> None:
    assert PHI > 1
    assert 2 > PHI
    assert -PHI < 0


def test_json_round_trip_and_integrality() -> None:
    value = Golden(Fraction(3, 2), -1)
    assert value.to_json() == {"a": "3/2", "b": "-1"}
    assert Golden.from_json(value.to_json()) == value
    assert PHI.is_integral
    assert not value.is_integral
    with pytest.raises(ValueError):
        Golden.from_json({"a": "1", "b": "0", "extra": "0"})


def test_zero_division_and_bad_coercion() -> None:
    with pytest.raises(ZeroDivisionError):
        ONE / ZERO
    with pytest.raises(TypeError):
        Golden.coerce(1.5)  # type: ignore[arg-type]
