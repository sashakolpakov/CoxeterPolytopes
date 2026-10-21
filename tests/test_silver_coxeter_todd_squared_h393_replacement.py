import hashlib
import json
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.silver import SILVER_ONE, Silver
from coxeter_polytopes.search.silver_coxeter_todd_squared import (
    COEFFICIENT,
    FIRST_CUTTING_CONTROL,
    FIRST_CUTTING_D,
    FIRST_CUTTING_E,
    FIRST_CUTTING_HEIGHT,
    FIRST_CUTTING_NORM,
    FIRST_CUTTING_U,
    FIRST_CUTTING_V,
    FIRST_PACKET_HEIGHT,
    REPLACEMENT_WITNESS_SCALE,
    WITNESS_SCALE,
    cross_gram_maximum_against_first_cutter,
    dual_coordinate_coset_minimum,
    dual_direction_level_bound,
    enumerate_controls,
    first_cutting_root,
    first_cutting_root_audit,
    gated_gram_triples,
    k12_pair,
    lorentz_pair,
    replacement_hyperideal_witness,
    replacement_witness_norm,
    strict_hyperideal_witness,
    transformed_difference_gram,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_coxeter_todd_squared_h393_replacement.json"
FIRST_PACKET_ARTIFACT = (
    ROOT / "data/enumerations/silver_coxeter_todd_squared_first_packet.json"
)


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def control(height: Silver, norm: Silver):
    rows = tuple(
        row
        for row in enumerate_controls(FIRST_CUTTING_HEIGHT)
        if row.height == height and row.norm_class.norm == norm
    )
    assert len(rows) == 1
    return rows[0]


def test_complete_control_and_necessary_gram_chronology_through_h393() -> None:
    controls = enumerate_controls(FIRST_CUTTING_HEIGHT)
    assert len(controls) == 128
    assert len({row.height for row in controls}) == 87
    assert sum(row.height <= FIRST_PACKET_HEIGHT for row in controls) == 50
    assert sum(FIRST_PACKET_HEIGHT < row.height < FIRST_CUTTING_HEIGHT for row in controls) == 76
    assert sum(row.height == FIRST_CUTTING_HEIGHT for row in controls) == 2
    assert tuple(
        (row.height, row.norm_class.norm, len(gated_gram_triples(row)))
        for row in controls
        if gated_gram_triples(row)
    ) == (
        (FIRST_PACKET_HEIGHT, Silver(1), 1),
        (FIRST_PACKET_HEIGHT, Silver(2), 2),
        (Silver(356, 252), Silver(6, 2), 2),
        (Silver(362, 256), Silver(10, -4), 1),
        (FIRST_CUTTING_HEIGHT, Silver(1), 4),
        (FIRST_CUTTING_HEIGHT, Silver(2), 11),
    )


def test_exact_dual_coset_formula_and_first_cutting_direction() -> None:
    assert [dual_coordinate_coset_minimum(k) for k in range(4)] == [
        Fraction(0),
        Fraction(4),
        Fraction(4),
        Fraction(10),
    ]
    for k in range(-20, 21):
        assert dual_coordinate_coset_minimum(k + 4) - dual_coordinate_coset_minimum(k) == 6 * k + 12

    rows = tuple(
        dual_direction_level_bound(row)
        for row in enumerate_controls(FIRST_CUTTING_HEIGHT)
        if gated_gram_triples(row)
    )
    assert [row.directional_pair_count for row in rows] == [5, 10, 42, 19, 140, 715]
    assert [row.maximum_spatial_pairing for row in rows] == [
        Silver(6, 4),
        Silver(8, 6),
        Silver(36, 26),
        Silver(30, 22),
        Silver(15, 11),
        Silver(22, 16),
    ]
    old_margins = [
        WITNESS_SCALE * row.maximum_spatial_pairing - COEFFICIENT * row.control
        for row in rows
    ]
    assert [value.sign(0) for value in old_margins[2:]] == [-1, -1, -1, 1]


def test_explicit_primitive_reflecting_wall_cuts_old_witness() -> None:
    assert FIRST_CUTTING_HEIGHT == Silver(393, 278)
    assert FIRST_CUTTING_CONTROL == Silver(10, 7)
    assert FIRST_CUTTING_NORM == Silver(2)
    assert k12_pair(FIRST_CUTTING_D, FIRST_CUTTING_D) == 28
    assert k12_pair(FIRST_CUTTING_E, FIRST_CUTTING_E) == 84
    assert k12_pair(FIRST_CUTTING_D, FIRST_CUTTING_E) == 46
    assert (
        k12_pair(FIRST_CUTTING_U, FIRST_CUTTING_U),
        k12_pair(FIRST_CUTTING_V, FIRST_CUTTING_V),
        k12_pair(FIRST_CUTTING_U, FIRST_CUTTING_V),
    ) == (380, 204, 278)
    assert (FIRST_CUTTING_U[0], FIRST_CUTTING_V[0]) == (22, 16)
    assert first_cutting_root_audit()
    root = first_cutting_root()
    assert lorentz_pair(root, root) == 2
    assert lorentz_pair(strict_hyperideal_witness(), root) == Silver(
        Fraction(8, 5), Fraction(-1, 5)
    )
    assert lorentz_pair(strict_hyperideal_witness(), root).sign(0) > 0


def test_exact_cross_gram_bounds_make_the_cutter_accepted() -> None:
    expected_earlier = (
        ((82, 40, 57), 379, 671, 47_659, Silver(356, 252), Silver(0)),
        ((78, 42, 57), 285, 865, 41_815, Silver(356, 252), Silver(0)),
        ((166, 84, 118), 671, 1_439, 30_749, Silver(509, 360), Silver(-7, -4)),
    )
    index = 0
    for height, norm in (
        (Silver(356, 252), Silver(6, 2)),
        (Silver(362, 256), Silver(10, -4)),
    ):
        row = control(height, norm)
        for triple in gated_gram_triples(row):
            expected = expected_earlier[index]
            assert triple == expected[0]
            maximum = cross_gram_maximum_against_first_cutter(
                transformed_difference_gram(triple)
            )
            assert (
                maximum.first_column_count,
                maximum.second_column_count,
                maximum.feasible_cross_gram_count,
                maximum.maximum,
            ) == expected[1:5]
            assert (
                row.norm_class.odd_spatial_divisor * maximum.maximum
                - COEFFICIENT * FIRST_CUTTING_CONTROL * row.control
            ) == expected[5]
            index += 1

    q1 = control(FIRST_CUTTING_HEIGHT, SILVER_ONE)
    expected_q1 = (
        (109_347, Silver(554, 392), Silver(-2, -1)),
        (128_591, Silver(555, 392), Silver(-1, -1)),
        (121_857, Silver(555, 392), Silver(-1, -1)),
        (91_205, Silver(555, 392), Silver(-1, -1)),
    )
    for triple, expected in zip(gated_gram_triples(q1), expected_q1, strict=True):
        maximum = cross_gram_maximum_against_first_cutter(
            transformed_difference_gram(triple)
        )
        assert maximum.feasible_cross_gram_count == expected[0]
        assert maximum.maximum == expected[1]
        assert maximum.maximum - COEFFICIENT * FIRST_CUTTING_CONTROL * q1.control == expected[2]


def test_seven_over_four_replacement_is_strict_hyperideal_through_bound() -> None:
    assert REPLACEMENT_WITNESS_SCALE == Fraction(7, 4)
    assert replacement_witness_norm() == Silver(Fraction(205, 12), -12)
    assert replacement_witness_norm().sign(0) > 0
    replacement = replacement_hyperideal_witness()
    assert lorentz_pair(replacement, replacement) == replacement_witness_norm()
    assert lorentz_pair(replacement, first_cutting_root()) == Silver(Fraction(1, 2), -1)
    assert lorentz_pair(replacement, first_cutting_root()).sign(0) < 0
    rows = tuple(
        dual_direction_level_bound(row)
        for row in enumerate_controls(FIRST_CUTTING_HEIGHT)
        if gated_gram_triples(row)
    )
    margins = tuple(
        Silver(REPLACEMENT_WITNESS_SCALE) * row.maximum_spatial_pairing
        - COEFFICIENT * row.control
        for row in rows
    )
    assert all(value.sign(0) < 0 for value in margins)
    assert margins[-1] == Silver(Fraction(1, 2), -1)


def test_frozen_artifact_scope_dependency_and_conclusion() -> None:
    data = payload()
    assert data["host"] == "<-(12*sqrt(2)-13)> + K12 + K12"
    assert data["chronology"]["control_count"] == 128
    assert data["chronology"]["distinct_height_count"] == 87
    assert data["dual_coordinate_gate"]["formula"] == (
        "mu(k)=3*k^2/4 + rho(k mod 4), rho=(0,13/4,1,13/4)"
    )
    assert data["first_packet_dependency"]["sha256"] == hashlib.sha256(
        FIRST_PACKET_ARTIFACT.read_bytes()
    ).hexdigest()
    assert data["first_cutting_wall"]["accepted"] is True
    assert data["replacement_hyperideal_ray"]["strict_against_every_necessary_row_through_bound"] is True
    assert data["conclusion"] == {
        "compact_example": False,
        "finite_volume_example": False,
        "first_later_cutting_height": {
            "display": "393 + 278*sqrt(2)",
            "exact": {"a": "393", "b": "278"},
            "signs": [1, -1],
        },
        "old_witness_is_cut": True,
        "prefix_through_bound_has_infinite_volume": True,
        "replacement_is_strict_hyperideal": True,
        "vinberg_terminated": False,
    }
    assert "not Vinberg termination" in data["scope"]
