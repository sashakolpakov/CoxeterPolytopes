import json
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_coxeter_todd_squared import (
    CUTTING_HEIGHT,
    OLD_WITNESS_SCALE,
    REPLACEMENT_WITNESS_SCALE,
    cutting_prefix_directional_maxima,
    cutting_root,
    cutting_root_audit,
    cubic7_k12_lorentz_pair,
    dual_witness,
    dual_witness_norm,
    replacement_prefix_pairing_bounds,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h33_cutter.json"


def test_explicit_cutting_root_is_norm_two_primitive_and_reflective():
    audit = cutting_root_audit()
    root = cutting_root()
    assert audit.height == CUTTING_HEIGHT == Cubic7(33, 94, 42)
    assert audit.control == Cubic7(1, 6, 3)
    assert audit.norm == Cubic7(2)
    assert audit.coefficient_gram == (
        (12, 30, 12),
        (30, 84, 36),
        (12, 36, 16),
    )
    assert audit.h_pairings == (4, 10, 4)
    assert cubic7_k12_lorentz_pair(root, root) == 2


def test_old_witness_is_cut_and_replacement_is_strict_hyperideal():
    audit = cutting_root_audit()
    root = cutting_root()
    assert OLD_WITNESS_SCALE == Fraction(3, 2)
    assert REPLACEMENT_WITNESS_SCALE == Fraction(7, 5)
    assert dual_witness_norm(OLD_WITNESS_SCALE) == Cubic7(3, -2)
    assert dual_witness_norm(OLD_WITNESS_SCALE).sign(0) > 0
    assert audit.old_witness_pairing == Cubic7(0, 1)
    assert cubic7_k12_lorentz_pair(dual_witness(OLD_WITNESS_SCALE), root) == Cubic7(0, 1)
    assert audit.old_witness_pairing.sign(0) > 0
    assert dual_witness_norm(REPLACEMENT_WITNESS_SCALE) == Cubic7(
        Fraction(196, 75), -2
    )
    assert dual_witness_norm(REPLACEMENT_WITNESS_SCALE).sign(0) > 0
    assert audit.replacement_witness_pairing == Cubic7(
        Fraction(-2, 5), 0, Fraction(-2, 5)
    )
    assert audit.replacement_witness_pairing.sign(0) < 0


def test_replacement_is_strict_against_both_h33_necessary_types():
    assert cutting_prefix_directional_maxima() == (
        Cubic7(0, 4, 2),
        Cubic7(2, 6, 3),
        Cubic7(2, 8, 4),
        Cubic7(3, 10, 5),
        Cubic7(4, 10, 4),
    )
    assert replacement_prefix_pairing_bounds() == (
        Cubic7(-2, Fraction(-2, 5), Fraction(4, 5)),
        Cubic7(Fraction(-6, 5), Fraction(-8, 5), Fraction(1, 5)),
        Cubic7(Fraction(-6, 5), Fraction(-4, 5), Fraction(-2, 5)),
        Cubic7(Fraction(-9, 5), 0, 1),
        Cubic7(Fraction(-2, 5), 0, Fraction(-2, 5)),
    )
    assert all(value.sign(0) < 0 for value in replacement_prefix_pairing_bounds())


def test_cutting_root_is_chronologically_accepted():
    audit = cutting_root_audit()
    assert audit.first_packet_maximum == Cubic7(0)
    assert audit.second_packet_maximum == Cubic7(0, -3, -2)
    assert audit.third_packet_maximum == Cubic7(-4, -8, -3)
    assert audit.other_cutting_type_schur_count == 1_225
    assert audit.other_cutting_type_maximum == Cubic7(-2, -8, -4)
    assert audit.same_height_schur_count == 8_173
    assert audit.same_height_minimum_gate_count == 6_960
    assert audit.same_height_maximum == Cubic7(0, -4, -2)
    assert all(
        value.sign(0) <= 0
        for value in (
            audit.first_packet_maximum,
            audit.second_packet_maximum,
            audit.third_packet_maximum,
            audit.other_cutting_type_maximum,
            audit.same_height_maximum,
        )
    )


def test_frozen_cutter_artifact_has_bounded_scope():
    payload = json.loads(ARTIFACT.read_text())
    assert payload["schema"] == "coxeter.cubic7_k12_squared_h33_cutter.v1"
    assert payload["chronology"]["nonempty_necessary_level_count"] == 4
    assert payload["chronology"]["necessary_gram_type_count"] == 5
    assert payload["first_cutting_wall"]["accepted"] is True
    assert payload["replacement_hyperideal_ray"][
        "strict_against_every_possible_wall_through_bound"
    ] is True
    assert payload["conclusion"]["prefix_through_bound_has_infinite_volume"] is True
    assert payload["conclusion"]["compact_example"] is False
    assert payload["conclusion"]["finite_volume_example"] is False
    assert payload["conclusion"]["vinberg_terminated"] is False
