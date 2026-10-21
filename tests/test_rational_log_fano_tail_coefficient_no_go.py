import itertools
import json
from fractions import Fraction
from pathlib import Path

import pytest

from scripts.certify_rational_log_fano_tail_coefficient_no_go import (
    ARCHIVE,
    build_certificate,
    h1_nef_intersections,
    h2_nef_intersections,
    physical_relative_coefficients,
    tail_energy,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_log_fano_tail_coefficient_no_go.md"


def test_tail_energy_identity_and_sharp_zero_boundary():
    for length in range(1, 51):
        row = tail_energy(length)
        assert row["exceptional_coordinates"] == ["1"] * length
        assert row["coordinate_sum"] == str(length)
        assert row["square_loss"] == str(length)
        assert row["local_tail_intersections"] == ["0"] * (length - 1) + ["1"]


def test_cauchy_bound_and_unique_equality_on_exact_coefficient_grid():
    grid = (Fraction(0), Fraction(1, 2), Fraction(1))
    for length in range(1, 7):
        for values in itertools.product(grid, repeat=length):
            chain, endpoint = values[:-1], values[-1]
            row = tail_energy(length, chain, endpoint)
            coordinate_sum = Fraction(row["coordinate_sum"])
            loss = Fraction(row["square_loss"])
            assert coordinate_sum == length + endpoint
            assert loss >= coordinate_sum * coordinate_sum / length
            assert loss >= length
            if loss == length:
                assert values == (Fraction(0),) * length


def test_prime_basis_relative_coefficients_are_the_discrepancy_ladder():
    for length in range(1, 41):
        assert physical_relative_coefficients(length) == [
            Fraction(index) for index in range(1, length + 1)
        ]
        half = [Fraction(1, 2)] * (length - 1)
        assert physical_relative_coefficients(length, half, Fraction(1, 2)) == [
            Fraction(index) + Fraction(1, 2)
            for index in range(1, length + 1)
        ]


def test_h1_and_h2_nef_tests_force_full_hostile_coefficients():
    grid = [Fraction(index, 20) for index in range(21)]
    h1_feasible = []
    h2_feasible = []
    for first in grid:
        for second in grid:
            if max(h1_nef_intersections(first, second)) <= 0:
                h1_feasible.append((first, second))
            if max(h2_nef_intersections(first, second)) <= 0:
                h2_feasible.append((first, second))
    assert h1_feasible == [(Fraction(1), Fraction(1))]
    assert h2_feasible == [(Fraction(1), Fraction(1))]

    # Any distinct effective boundary meeting the hostile pair makes the
    # inequalities strictly worse, including at coefficient one.
    assert max(h1_nef_intersections(1, 1, Fraction(1, 7), 0)) > 0
    assert max(h2_nef_intersections(1, 1, 0, Fraction(1, 7))) > 0


def test_tail_input_rejection():
    for bad_length in (0, -1, True, Fraction(3, 2), "3"):
        with pytest.raises(ValueError):
            tail_energy(bad_length)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        tail_energy(3, [0])
    for bad_coefficient in (-1, 2, 0.5, True, "not-a-number"):
        with pytest.raises(ValueError):
            tail_energy(1, (), bad_coefficient)  # type: ignore[arg-type]


def test_archive_scope_and_report_controls():
    payload = build_certificate()
    assert json.loads(ARCHIVE.read_text(encoding="utf-8")) == payload
    assert payload["scope"] == {
        "method": "surface blowup, adjunction, canonical section, exact intersections",
        "absolute_tail_bound_claimed": False,
        "conditional_bound": "m<[-(K0+Delta0)]^2 if -(K+Delta) is nef and big",
        "minimal_hostile_H1_H2_weak_log_Fano_exists": False,
        "higher_parallel_dotted_H1_H2_claimed": False,
        "classification_used": False,
    }

    raw = REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)
    text = raw.decode("utf-8")
    for needle in (
        r"L_m^2\leq L_0^2-\frac{(m+b)^2}{m}\leq L_0^2-m",
        r"m<L_0^2",
        "a_1=a_2=1",
        "a=b=1",
        "not even\nnef",
        "crepant pullback would require",
        "becomes easier, not harder",
        "no effective boundary on the physical",
        "No Nikulin--Vinberg",
        "Higher, parallel, or dotted H1/H2 incidences are not\nclaimed",
    ):
        assert needle in text
