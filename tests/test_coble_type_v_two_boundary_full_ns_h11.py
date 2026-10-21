from __future__ import annotations

import json

import pytest

from scripts.certify_coble_type_v_two_boundary_full_ns_h11 import OUTPUT, build_certificate


@pytest.fixture(scope="module")
def certificate() -> dict[str, object]:
    return build_certificate()


def test_type_v_two_boundary_certificate_matches_archive(
    certificate: dict[str, object],
) -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert certificate == archived
    assert archived["all_checks_pass"] is True


def test_type_v_full_nef_h11_volume_is_literal(
    certificate: dict[str, object],
) -> None:
    assert certificate["direction"] == "SURFACE_TO_POLYTOPE_only"
    assert certificate["classification"] == (
        "finite_volume_noncompact_full_nef_H11_36_arithmetic_Coxeter_polytope"
    )
    cone = certificate["exact_double_description"]
    assert (cone["finite_vertex_count"], cone["ideal_vertex_count"], cone["hyperideal_vertex_count"]) == (
        6_356,
        432,
        0,
    )
    assert cone["all_facets_irredundant"] is True
    literal = certificate["literal_vinberg"]
    assert literal["elliptic_edge_count"] == 47_929
    assert literal["endpoint_type_histogram"] == {"FF": 25_372, "FI": 19_172, "II": 3_385}
    assert literal["bad_edge_count"] == literal["bad_parabolic_neighborhood_count"] == 0


def test_type_v_curve_lifts_and_arithmetic_are_exact(
    certificate: dict[str, object],
) -> None:
    curves = certificate["curve_exhaustion"]
    assert curves["negative_curve_types"] == {
        "boundary_minus4": 2,
        "minus2": 16,
        "minus1": 18,
    }
    assert curves["split_CM_cusps_zero_based"] == [4, 5, 6, 7, 8, 9, 10]
    assert len(curves["together_CM_cusps_zero_based"]) == 22
    lattice = certificate["full_NS_lattice"]
    assert lattice["unique_integral_index_two_overlattice"] is True
    assert lattice["NS_identification"] == "odd unimodular I_(1,11)"
    arithmetic = certificate["arithmetic"]
    assert arithmetic["primitive_square_divisibility_histogram"] == {
        "-4,2": 16,
        "-2,1": 6,
        "-2,2": 14,
    }
    assert arithmetic["all_36_reflections_integral_on_arithmetic_lattice"] is True
    assert arithmetic["walls_span_whole_arithmetic_lattice"] is True
    assert certificate["lineage"]["worldwide_literal_priority"] == "open"
