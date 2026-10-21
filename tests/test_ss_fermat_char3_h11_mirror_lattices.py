from __future__ import annotations

import json

from scripts.certify_ss_fermat_char3_h11_mirror_lattices import (
    OUTPUT,
    build_certificate,
)


def test_h11_mirror_lattice_archive_replays_exactly() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert build_certificate() == archived
    assert archived["all_checks_pass"] is True


def test_c2_h11_all_mirrors_are_integral_on_the_fixed_lattice() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))["c2_h11"]
    assert payload["facet_count"] == 206
    assert payload["all_reflections_integral_on_fixed_lattice"] is True
    assert payload["primitive_square_divisibility_integrality_histogram"] == {
        "-2,1,integral": 12,
        "-2,2,integral": 164,
        "-4,2,integral": 30,
    }


def test_v4_h11_minimal_stable_overlattice_is_index_four() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))["v4_h11"]
    assert payload["integral_reflection_count_on_fixed_lattice"] == 292
    assert payload["nonintegral_reflection_count_on_fixed_lattice"] == 36
    assert payload["primitive_square_divisibility_integrality_histogram"][
        "-24,6,nonintegral"
    ] == 36
    stable = payload["minimal_reflection_stable_overlattice"]
    assert stable["overlattice_index"] == 4
    assert stable["determinant_absolute"] == 48
    assert stable["smith_diagonal_absolute"] == [1] * 8 + [2] * 3 + [6]
    assert stable["all_328_reflections_integral_on_stable_lattice"] is True
    assert stable["primitive_square_divisibility_histogram"] == {
        "-2,1": 28,
        "-4,2": 24,
        "-6,3": 276,
    }
    assert payload["arithmeticity"] == {
        "field": "Q",
        "reason": (
            "all facet reflections preserve the integral Lorentzian stable "
            "overlattice and the Coxeter group has finite covolume"
        ),
        "arithmetic": True,
    }
