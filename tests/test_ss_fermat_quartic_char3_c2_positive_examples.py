from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/enumerations"


def load(name: str) -> dict:
    return json.loads((DATA / name).read_text(encoding="utf-8"))


def test_both_cyclic_casts_have_complete_exact_coxeter_reductions() -> None:
    census = load("ss_fermat_quartic_char3_c2_element_class_screen.json")
    assert census["schema"] == "ss_fermat_quartic_char3_c2_element_class_screen_v1"
    assert census["full_group_order"] == 26_127_360
    assert census["element_class_count"] == 61

    screen = load("ss_fermat_quartic_char3_c2_survivor_screen.json")
    cases = {case["hyperbolic_dimension"]: case for case in screen["casts"]}
    assert set(cases) == {10, 11}
    assert (
        cases[10]["candidate_ray_count"],
        cases[10]["facet_count"],
        cases[10]["redundant_ray_count"],
    ) == (2031, 196, 1835)
    assert (
        cases[11]["candidate_ray_count"],
        cases[11]["facet_count"],
        cases[11]["redundant_ray_count"],
    ) == (1892, 206, 1686)
    for case in cases.values():
        assert case["exact_complete"] is True
        assert case["unresolved_count"] == 0
        assert case["coxeter_pair_gate"]["all_pairs_coxeter"] is True
        assert case["coxeter_pair_gate"]["bad_pair_count"] == 0
        assert len(case["strict_facet_witnesses"]) == case["facet_count"]
        assert len(case["redundancy_certificates"]) == case["redundant_ray_count"]


def test_literal_vinberg_certificates_are_finite_volume_noncompact() -> None:
    expected = {
        10: (196, 1_427_720, 212_824, 10_560),
        11: (206, 2_688_957, 327_566, 18_354),
    }
    for dimension, counts in expected.items():
        certificate = load(
            f"ss_fermat_quartic_char3_c2_h{dimension}_vinberg_certificate.json"
        )
        facets, edges, finite, ideal = counts
        assert certificate["facet_count"] == facets
        assert certificate["elliptic_edge_rank"] == dimension - 1
        assert certificate["elliptic_edge_count"] == edges
        assert certificate["finite_vertex_count"] == finite
        assert certificate["ideal_vertex_count"] == ideal
        assert certificate["endpoint_count_histogram"] == {"2": edges}
        assert certificate["bad_edge_count"] == 0
        assert certificate["unrecognized_endpoint_edge_count"] == 0
        assert certificate["invalid_full_zero_neighborhood_count"] == 0
        assert certificate["finite_volume"] is True
        assert certificate["compact"] is False
        assert certificate["vinberg_decision"] == "finite_volume_noncompact"


def test_h10_independent_exact_cone_matches_literal_counts() -> None:
    cone = load("ss_fermat_quartic_char3_c2_h10_exact_cone_certificate.json")
    assert cone["ray_count"] == 223_384
    assert cone["finite_ray_count"] == 212_824
    assert cone["ideal_ray_count"] == 10_560
    assert cone["hyperideal_ray_count"] == 0
    assert cone["literal_vinberg_endpoint_counts_match"] is True
    assert cone["finite_volume_certified"] is True
    for dimension, facets in ((10, 196), (11, 206)):
        text = (DATA / f"ss_fermat_quartic_char3_c2_h{dimension}_exact_cone.ine").read_text()
        assert text.startswith("H-representation\nbegin\n")
        assert f"{facets} {dimension + 2} rational" in text


def test_h11_fixed_lattice_has_claimed_two_elementary_invariants() -> None:
    fixed = load("ss_fermat_quartic_char3_survivor_fixed_lattices.json")
    case = next(case for case in fixed["cases"] if case["id"] == "char3-C2-H11")
    assert case["rank"] == 12
    assert case["signature"] == [1, 11]
    assert case["even"] is True
    assert case["determinant_absolute"] == 1024
    assert case["smith_diagonal_absolute"] == [1, 1] + [2] * 10
