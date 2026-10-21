import json
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_coxeter_todd_squared import (
    POST_H33_CUTOFF_HEIGHT,
    POST_H33_REPLACEMENT_SCALE,
    dual_witness_norm,
    post_h33_directional_audit,
    post_h33_threat_audit,
)

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h169.json"


def test_post_h33_frontier_and_first_realized_threat():
    rows = post_h33_directional_audit()
    assert [len(row["types"]) for row in rows] == [4, 1, 8, 13, 16, 19, 37]
    assert rows[-1]["level"].height == POST_H33_CUTOFF_HEIGHT
    threat = rows[-1]["types"][35]
    assert threat["actual_arg"] == (9, 23, 10)
    assert threat["actual_direction"] == Cubic7(9, 23, 10)
    assert threat["seven_fifths_bound"] == Cubic7(Fraction(3, 5), Fraction(1, 5))
    assert threat["shell_realization_count"] == 1_693_440
    audit = post_h33_threat_audit()
    assert audit["norm"] == Cubic7(2)
    assert audit["h33_cutter_pairing"] == Cubic7(-15, -46, -21)


def test_eleven_eighths_replacement_is_strict_through_cutoff():
    assert POST_H33_REPLACEMENT_SCALE == Fraction(11, 8)
    assert dual_witness_norm(POST_H33_REPLACEMENT_SCALE) == Cubic7(Fraction(121, 48), -2)
    rows = post_h33_directional_audit()
    assert all(item["eleven_eighths_bound"].sign(0) < 0 for row in rows for item in row["types"])
    assert rows[-1]["types"][35]["eleven_eighths_bound"] == Cubic7(Fraction(3, 8), Fraction(-3, 8), Fraction(-1, 4))


def test_h169_frozen_artifact_has_strict_bounded_scope():
    payload = json.loads(ARTIFACT.read_text())
    assert payload["first_realized_threat"]["status"] == "realized_threat_not_claimed_accepted"
    assert payload["replacement"]["strict_against_every_necessary_type_through_cutoff"] is True
    assert payload["conclusion"]["bounded_prefix_has_infinite_volume"] is True
    assert payload["conclusion"]["accepted_new_wall"] is False
    assert payload["conclusion"]["compact_example"] is False
    assert payload["conclusion"]["finite_volume_example"] is False
    assert payload["conclusion"]["vinberg_terminated"] is False
