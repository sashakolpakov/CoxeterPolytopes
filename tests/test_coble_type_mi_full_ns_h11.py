from __future__ import annotations

import json

import pytest

from scripts.certify_coble_type_mi_full_ns_h11 import OUTPUT, build_certificate


@pytest.fixture(scope="module")
def certificate() -> dict[str, object]:
    return build_certificate()


def test_type_mi_certificate_matches_archive(certificate: dict[str, object]) -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert certificate == archived
    assert archived["all_checks_pass"] is True


def test_type_mi_full_nef_cone(certificate: dict[str, object]) -> None:
    chamber = certificate["full_nef_chamber"]
    assert chamber["hyperbolic_dimension"] == 11
    assert chamber["facet_count"] == 206
    assert chamber["norm_histogram"] == {"-4": 2, "-2": 30, "-1": 174}
    assert chamber["all_facets_irredundant"] is True
    assert chamber["exact_double_description"]["finite_vertex_count"] == 327_566
    assert chamber["exact_double_description"]["ideal_vertex_count"] == 18_354
    assert chamber["exact_double_description"]["hyperideal_vertex_count"] == 0
    assert chamber["literal_vinberg"]["bad_edge_count"] == 0
    assert chamber["literal_vinberg"]["invalid_full_zero_neighborhood_count"] == 0


def test_type_mi_exact_fermat_identity(certificate: dict[str, object]) -> None:
    bridge = certificate["exact_Fermat_polytope_identity"]
    assert bridge["all_42436_Gram_entries_match"] is True
    assert bridge["marked_isometry_determinant"] in (-1, 1)
    assert bridge["same_oriented_polytope"] is True
    assert certificate["coble_mukai"]["primitive_cusp_count"] == 247
    assert certificate["coble_mukai"]["cusp_lift_membership_histogram"] == {
        "none": 165,
        "both_boundaries": 82,
    }
