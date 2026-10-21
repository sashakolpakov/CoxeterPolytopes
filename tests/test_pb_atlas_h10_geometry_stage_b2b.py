from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
import sympy as sp

from scripts.certify_pb_atlas_h10_geometry_stage_b2b import (
    _singular_absolute_irreducible,
    build_certificate,
)
from scripts.pb_atlas_geometry_stage_b2b import _closed_substitution


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/pb_atlas_h10_geometry_stage_b2b.json"


@pytest.fixture(scope="module")
def certificate() -> dict:
    return build_certificate()


def test_generated_incidence_and_collision_component(certificate: dict) -> None:
    incidence = certificate["generated_incidence"]
    assert incidence["viable_chart_count"] == 1
    assert len(incidence["all_chart_records"]) == 8
    assert "(q8_1 - 1)*(q8_1 + 1)" in incidence["groebner_basis"]
    records = incidence["component_records"]
    assert sum(record["collision_free"] for record in records) == 1
    assert all(
        record["selected_chart_open_audit"][
            "every_declared_chart_open_factor_nonzero"
        ]
        for record in records
    )
    bad = next(record for record in records if not record["collision_free"])
    assert len(bad["proper_center_collision_pairs_zero_based"]) == 6


def test_strong_coble_and_all_physical_facets(certificate: dict) -> None:
    geometry = certificate["geometric_component"]
    assert geometry["minus_K"] == {
        "rows": 10,
        "columns": 10,
        "rank": 10,
        "kernel_dimension": 0,
    }
    assert geometry["minus_two_K"] == {
        "rows": 30,
        "columns": 28,
        "rank": 27,
        "kernel_dimension": 1,
    }
    assert geometry["strong_Coble"] is True
    assert geometry["physical_wall_count"] == 30
    assert geometry["constructed_wall_indices_partition_supplied_wall_list"] is True
    assert geometry["every_proposed_facet_effective_irreducible"] is True
    assert geometry["contracted_exceptional_facets"]["count"] == 10
    assert (
        geometry["contracted_exceptional_facets"]["wall_indices_zero_based"]
        == geometry["contracted_exceptional_facets"][
            "source_wall_indices_from_reverse_blowups_zero_based"
        ]
    )
    assert geometry["contracted_exceptional_facets"][
        "every_source_exceptional_matches_coefficientwise"
    ] is True
    assert all(
        row["coefficientwise_match"]
        for row in geometry["contracted_exceptional_facets"][
            "coefficientwise_source_exceptional_checks"
        ]
    )
    assert len(geometry["noncontracted_nonboundary_facets"]) == 19
    assert geometry["physical_boundary"]["physical_curve_is_twice_wall_root"] is True
    assert geometry["physical_boundary"]["wall_root_self_intersection"] == -1
    assert geometry["physical_boundary"]["physical_curve_self_intersection"] == -4


def test_every_generated_curve_has_exact_multiplicity_and_is_absolute(certificate: dict) -> None:
    geometry = certificate["geometric_component"]
    curves = [
        geometry["physical_boundary"],
        *geometry["noncontracted_nonboundary_facets"],
    ]
    assert all(curve["every_declared_multiplicity_is_exact"] for curve in curves)
    assert all(
        curve["absolute_factorization"]["absolutely_irreducible"]
        for curve in curves
    )
    assert all(
        curve["absolute_factorization"]["sole_nontrivial_factor_multiplicity"]
        == 1
        for curve in curves
    )


def test_dimension_bridge_is_not_same_p_circular(certificate: dict) -> None:
    bridge = certificate["dimension_bridge"]
    assert bridge["does_not_require_same_nef_cone_for_rank_bound"] is True
    assert bridge["current_calibration_picard_rank"] == 11
    assert bridge["current_calibration_hyperbolic_dimension"] == 10
    assert bridge["current_calibration_not_claimed_to_be_reflection_maximal_here"] is True
    assert bridge["current_compiler_base_field"] == "Q, with geometric conclusion over C"
    assert "does not construct" in bridge["odd_characteristic_analogue"]


def test_closed_coordinate_parser_rejects_inexact_missing_and_alias_inputs() -> None:
    chart = {"coordinate_variables": ["u", "v"]}
    assert _closed_substitution(chart, {"u": 1, "v": "3/2"}) == {
        sp.Symbol("u"): sp.Rational(1),
        sp.Symbol("v"): sp.Rational(3, 2),
    }
    with pytest.raises(ValueError, match="exact rational"):
        _closed_substitution(chart, {"u": 0.1, "v": 1})
    with pytest.raises(ValueError, match="missing"):
        _closed_substitution(chart, {"u": 1})
    with pytest.raises(ValueError, match="occurs more than once"):
        _closed_substitution(
            chart,
            {"u": 1, sp.Symbol("u"): 1, "v": 1},
        )


def test_absolute_factor_gate_rejects_a_nonreduced_power() -> None:
    with pytest.raises(ArithmeticError, match="not absolutely irreducible"):
        _singular_absolute_irreducible("(X + Y)^2", 2)


def test_frozen_archive_equals_exact_replay(certificate: dict) -> None:
    assert json.loads(ARCHIVE.read_text(encoding="utf-8")) == certificate
    for relative, expected in certificate["sources"].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        assert actual == expected
    assert certificate["all_checks_pass"] is True
