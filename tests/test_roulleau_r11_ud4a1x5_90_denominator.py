import json

from scripts.certify_roulleau_r11_ud4a1x5_90_denominator import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text())


def test_archive_recomputes_exactly() -> None:
    assert build_certificate() == _archive()


def test_exact_ninety_wall_lattice_reconstruction() -> None:
    payload = _archive()
    assert payload["all_checks_pass"]
    assert payload["lattice"]["picard_rank"] == 11
    assert payload["lattice"]["hyperbolic_dimension"] == 10
    assert payload["lattice"]["determinant"] == 128
    assert payload["lattice"]["smith_invariants"] == [1, 1, 1, 1] + [2] * 7
    assert payload["walls"]["count"] == 90
    assert len(payload["walls"]["coordinate_rows_in_A1_through_A11"]) == 90
    assert len(payload["walls"]["gram"]) == 90
    assert all(payload["walls"]["reconstruction_checks"].values())


def test_full_weighted_group_is_exactly_s8_times_c2() -> None:
    group = _archive()["weighted_wall_group"]
    assert group["exact_structure"] == "S8 x C2"
    assert group["exact_order"] == 80640
    assert group["upper_bound_certificate"]["upper_bound"] == 80640
    assert all(group["upper_bound_certificate"]["checks"].values())
    assert all(group["generator_certificate"]["presentation_checks"].values())
    assert all(group["exactness_checks"].values())
    assert [
        row["fixed_rank"] for row in group["generator_certificate"]["generators"]
    ] == [10, 10, 10, 10, 10, 10, 10, 3]


def test_vinberg_cusp_and_structural_finite_volume_are_separate() -> None:
    geometry = _archive()["hyperbolic_geometry"]
    assert geometry["finite_volume_certified"]
    assert geometry["compact"] is False
    cusp = geometry["vinberg_noncompactness_certificate"]
    assert cusp["wall_degree_histogram"] == {"0": 18, "1": 2, "2": 70}
    assert len(cusp["affine_components_one_based"]) == 9
    assert all(len(component) == 2 for component in cusp["affine_components_one_based"])
    assert all(cusp["checks"].values())


def test_no_nontrivial_fold_can_remain_in_h10() -> None:
    denominator = _archive()["proper_high_fold_denominator"]
    assert denominator["status"] == "closed_without_subgroup_enumeration"
    assert denominator["proper_fixed_rank_at_least_11_count"] == 0
    assert denominator["proper_fold_hyperbolic_dimension_at_least_10_count"] == 0
    assert all(denominator["checks"].values())
    lineage = _archive()["lineage"]
    assert lineage["new_abstract_polytope_count"] == 0
    assert lineage["new_higher_dimensional_fold_count"] == 0
