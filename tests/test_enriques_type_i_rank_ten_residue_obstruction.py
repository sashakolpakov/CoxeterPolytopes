from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from scripts.certify_enriques_type_i_rank_ten_residue_obstruction import (
    CENSUS,
    OUTPUT,
    SHIMADA,
    analyze,
)
from scripts.certify_enriques_variable_lift_package import parse_gap_records


def _payload() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_archived_certificate_replays_exactly() -> None:
    assert analyze() == _payload()
    assert _payload()["all_checks_pass"] is True


def test_forced_root_uses_only_the_last_two_walls_and_half_sum_glue() -> None:
    census = json.loads(CENSUS.read_text(encoding="utf-8"))
    anti_gram = sp.Matrix(census["rank_ten_residue"]["canonical_gram"])
    walls = sp.Matrix(parse_gap_records(SHIMADA)["12A"]["walls"])

    assert all(int(value) % 2 == 0 for value in walls.row(10) - walls.row(11))
    assert anti_gram[10, 10] == anti_gram[11, 11] == -4
    assert anti_gram[10, 11] == 0
    for sign in (-1, 1):
        coefficient = sp.zeros(12, 1)
        coefficient[10] = sp.Rational(1, 2)
        coefficient[11] = sp.Rational(sign, 2)
        assert (coefficient.T * anti_gram * coefficient)[0] == -2


def test_raw_rootless_lattice_becomes_rootful_after_compatible_saturation() -> None:
    payload = _payload()
    raw = payload["raw_generated_anti_lattice"]
    saturated = payload["compatible_anti_saturation"]
    gates = payload["good_package_gates"]

    assert raw["isometry_type"] == "E8(-2) + <-4> + <-4>"
    assert raw["rootless"] is True
    assert raw["determinant_absolute"] == 4096
    assert saturated["index_over_raw_generated_lattice"] == 2
    assert saturated["isometry_type"] == "E8(-2) + A1(-1)^2"
    assert saturated["contains_minus_two_roots"] is True
    assert saturated["determinant_absolute"] == 1024
    assert gates["L2_rootlessness_of_saturated_anti_lattice"] == "FAILS"
    assert gates["L4_irreducibility"].startswith("not_reached")
