import json
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7, THETA7
from coxeter_polytopes.search.cubic7_leech_continuation import (
    DELTA7,
    P7_LEECH_SEED,
    SUPPORTED_LEECH_SEEDS,
    TWO_P7_LEECH_SEED,
    TWO_THETA_LEECH_SEED,
    UNIT_LEECH_SEED,
    cubic7_leech_control_levels,
    cubic7_leech_height_zero_audit,
    cubic7_leech_necessary_gram_types,
    cubic7_leech_prefix_obstruction,
    minimal_leech_inner_product_distribution,
    norm_six_direction_audit,
    two_theta_first_two_packets,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/cubic7_leech_h24_h1000_continuation.json"


def test_four_small_coefficient_directions_have_complete_local_norms_and_no_height_zero():
    assert all(seed.admissible for seed in SUPPORTED_LEECH_SEEDS)
    assert UNIT_LEECH_SEED.canonical_mirror_norms == (Cubic7(1), Cubic7(2))
    assert TWO_THETA_LEECH_SEED.canonical_mirror_norms == (Cubic7(2),)
    assert P7_LEECH_SEED.canonical_mirror_norms == (
        Cubic7(1),
        Cubic7(2),
        DELTA7,
        2 * DELTA7,
    )
    assert TWO_P7_LEECH_SEED.canonical_mirror_norms == (
        Cubic7(2),
        2 * DELTA7,
    )
    for seed in SUPPORTED_LEECH_SEEDS:
        assert len(cubic7_leech_height_zero_audit(seed)) == len(
            seed.canonical_mirror_norms
        )
        assert all(
            row.necessary_gram_type_count == 0
            for row in cubic7_leech_height_zero_audit(seed)
        )


def test_explicit_norm_six_direction_and_minimal_shell_distribution_are_exact():
    direction = norm_six_direction_audit()
    assert direction["square"] == 6
    assert minimal_leech_inner_product_distribution() == {
        -4: 1,
        -2: 4600,
        -1: 47104,
        0: 93150,
        1: 47104,
        2: 4600,
        4: 1,
    }


def test_two_theta_first_two_packets_are_complete_and_acute():
    packets = two_theta_first_two_packets()
    assert packets.first_height == Cubic7(5, 14, 6)
    assert packets.first_wall_count == 196_560
    assert packets.first_same_height_maximum == Cubic7(-4, -14, -6)
    assert packets.second_height == Cubic7(16, 45, 20)
    assert packets.minimal_pairs_per_fixed_vector == 47_104
    assert packets.second_wall_count == 9_258_762_240
    assert packets.first_second_maximum == Cubic7(-3, -7, -3)
    assert packets.second_same_height_cross_type_count == 494
    assert packets.second_same_height_maximum == Cubic7(-1)


def test_two_theta_control_box_has_eighteen_levels_and_expected_first_two_nonempty_types():
    levels = cubic7_leech_control_levels(TWO_THETA_LEECH_SEED, 1000)
    assert len(levels) == 18
    nonempty = tuple(
        (level, cubic7_leech_necessary_gram_types(level))
        for level in levels[:6]
        if cubic7_leech_necessary_gram_types(level)
    )
    assert tuple(level.height for level, _types in nonempty) == (
        Cubic7(5, 14, 6),
        Cubic7(16, 45, 20),
    )
    assert tuple(len(types) for _level, types in nonempty) == (1, 1)


def test_secondary_h300_prefix_obstructions_are_strict():
    unit = cubic7_leech_prefix_obstruction(
        UNIT_LEECH_SEED, 300, Fraction(7, 15)
    )
    p7 = cubic7_leech_prefix_obstruction(
        P7_LEECH_SEED, 300, Fraction(1, 2)
    )
    two_p7 = cubic7_leech_prefix_obstruction(
        TWO_P7_LEECH_SEED, 300, Fraction(3, 5)
    )
    assert (unit.control_level_count, len(unit.nonempty_necessary_levels)) == (37, 1)
    assert (p7.control_level_count, len(p7.nonempty_necessary_levels)) == (51, 0)
    assert (two_p7.control_level_count, len(two_p7.nonempty_necessary_levels)) == (13, 2)
    assert unit.strict_hyperideal and p7.strict_hyperideal and two_p7.strict_hyperideal


def test_frozen_h1000_artifact_records_the_strict_scope():
    payload = json.loads(ARTIFACT.read_text())
    assert payload["schema"] == "coxeter.cubic7_leech_h24_continuation.v1"
    assert payload["two_theta_actual_packets"]["first"]["wall_count"] == 196_560
    assert payload["two_theta_actual_packets"]["second"]["wall_count"] == 9_258_762_240
    deep = payload["two_theta_h1000"]
    assert deep["control_level_count"] == 18
    assert deep["nonempty_necessary_level_count"] == 9
    assert deep["total_necessary_gram_type_count"] == 47
    assert deep["witness"]["coordinate"] == "2/3"
    assert deep["witness"]["strict_hyperideal"] is True
    conclusion = payload["conclusion"]
    assert conclusion["new_compact_example"] is False
    assert conclusion["new_finite_volume_example"] is False
