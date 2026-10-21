"""Regression tests for the direct full-NS rational n<=11 theorem."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest

from scripts.certify_rational_direct_full_ns_n11_dimension_theorem import (
    build_certificate,
    degree_one_root_records,
    exceptional_chain_coordinates,
    exceptional_prime_record,
    exceptional_remainder_records,
    f0_swap_record,
    f2_minus_four_low_pairing,
    f2_mixed_ray,
    f2_pair,
    f2_standard_walls,
    mixed_ray_records,
    plane_minus_four_low_pairing,
    plane_mixed_ray,
    plane_pair,
    plane_standard_walls,
    section_integrality_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_direct_full_ns_n11_dimension_theorem.md"
CERTIFIER = (
    ROOT / "scripts/certify_rational_direct_full_ns_n11_dimension_theorem.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/rational_direct_full_ns_n11_dimension_theorem.json"
)


def test_minimal_base_mirrors_are_exact_on_arbitrary_exceptional_rank() -> None:
    for rank in range(0, 51):
        row = f0_swap_record(rank)
        assert row["root_square"] == -2
        assert row["c_maps_to_d"] is True
        assert row["d_maps_to_c"] is True
        assert row["c_plus_d_is_on_mirror"] is True
        assert row["all_total_exceptional_coordinates_fixed"] is True

    for e in range(2, 20):
        for hits in range(0, 10):
            row = section_integrality_record(e, hits)
            assert row["negative_section_square"] == e + hits
            assert row["reflection_integral_on_full_NS"] is (
                e + hits in (1, 2)
            )
    assert section_integrality_record(2, 0)[
        "reflection_integral_on_full_NS"
    ] is True
    assert section_integrality_record(2, 1)[
        "reflection_integral_on_full_NS"
    ] is False


def test_every_possible_final_satellite_multiplicity_is_forbidden() -> None:
    # A satellite gives an older exceptional prime at least two children.
    # k=2 is square -3; k=3 is primitive square -4 but nonintegral; larger
    # k leave the stipulated physical alphabet.
    for children in range(2, 101):
        row = exceptional_prime_record(children)
        assert row["strict_prime_square"] == -(children + 1)
        assert row["strict_prime_primitive"] is True
        assert row["compatible_with_both_alphabet_and_integrality"] is False
    assert exceptional_prime_record(2)["square_in_physical_alphabet"] is False
    assert exceptional_prime_record(3)["square_in_physical_alphabet"] is True
    assert exceptional_prime_record(3)[
        "reflection_integral_on_unimodular_exceptional_lattice"
    ] is False


def test_degree_one_physical_roots_and_half_minus_four_exclusions() -> None:
    for length in range(3, 18):
        rows = degree_one_root_records(length)
        assert rows["plane_physical_minus_one_patterns"] == [
            [1, 1, *([0] * (length - 2))]
        ]
        assert rows["plane_physical_minus_two_patterns"] == [
            [1, 1, 1, *([0] * (length - 3))]
        ]
        assert rows["plane_half_minus_four_equations"]["solutions"] == []
        assert rows["F2_physical_minus_one_patterns"] == [
            [1, *([0] * (length - 1))]
        ]
        assert rows["F2_physical_minus_two_patterns"] == [
            [1, 1, *([0] * (length - 2))]
        ]
        assert rows["F2_half_minus_four_equations"]["solutions"] == []


def test_exceptional_remainders_are_outside_the_free_chain_cone() -> None:
    for m in range(3, 50):
        rows = exceptional_remainder_records(m)
        assert rows["plane_in_nonnegative_exceptional_cone"] is False
        assert rows["F2_in_nonnegative_exceptional_cone"] is False
        assert rows["plane_generator_coefficients"][2] == -1
        assert rows["F2_generator_coefficients"][1] == -1
    assert exceptional_chain_coordinates((1, 0, -1, 2)) == (1, 1, 0, 2)


def test_both_first_mixed_rays_are_spacelike_and_in_standard_wall_cones() -> None:
    for m in range(12, 61):
        ray = plane_mixed_ray(m)
        pairings = [plane_pair(ray, wall) for wall in plane_standard_walls(m)]
        assert plane_pair(ray, ray) == -4
        assert plane_pair((-3, *((1,) * m)), ray) == 4
        assert all(value >= 0 for value in pairings)
        assert pairings[1] == 1
        assert pairings[12] == 2
        assert sum(value != 0 for value in pairings) == 2

    for m in range(11, 61):
        ray = f2_mixed_ray(m)
        pairings = [f2_pair(ray, wall) for wall in f2_standard_walls(m)]
        assert f2_pair(ray, ray) == -4
        assert f2_pair((-2, -4, *((1,) * m)), ray) == 4
        assert all(value >= 0 for value in pairings)
        assert pairings[0] == 1
        assert pairings[12] == 2
        assert sum(value != 0 for value in pairings) == 2


def test_hidden_minus_four_low_pairing_branches_are_exhausted() -> None:
    plane = plane_minus_four_low_pairing()
    ruled = f2_minus_four_low_pairing()
    assert plane["solutions"] == 0
    assert ruled["solutions"] == 0
    assert [row["R"] for row in plane["rows"]] == [0, 1, 2, 3]
    assert [row["a"] for row in ruled["rows"]] == [0, 1, 2, 3]
    assert all(row["square_solutions"] == 0 for row in plane["rows"])
    assert all(row["square_solutions"] == 0 for row in ruled["rows"])


def test_certificate_scope_and_dimension_conclusion_are_literal() -> None:
    payload = build_certificate(40)
    assert payload["status"] == "PASS_DIRECT_FULL_NS_DIMENSION_BOUND"
    scope = payload["scope"]
    assert scope["host"] == "direct full Neron-Severi lattice only"
    assert scope["physical_negative_extremal_alphabet"] == (
        "smooth rational (-1),(-2),(-4)"
    )
    assert scope["square_zero_tangent_Mori_rays_allowed"] is True
    assert scope["every_facet_reflection_integral_on_full_NS"] is True
    assert scope["adapted_half_boundary_host_claimed"] is False
    assert scope["reflective_lattice_classification_or_census_used"] is False
    assert payload["minimal_base"]["survivors"] == ["P2", "F2"]
    assert payload["exceptional_geometry"][
        "all_successors_forced_free_on_newest_exceptional"
    ] is True
    assert payload["exact_frames"][
        "satellite_terminology_used_for_F2_second_centre"
    ] is False
    assert payload["conclusion"] == {
        "plane_blowups_at_most": 11,
        "F2_blowups_at_most": 10,
        "projective_nef_dimension": (
            "rho-1 even when nonnegative isotropic extreme rays occur"
        ),
        "hyperbolic_dimension_at_most": 11,
    }
    records = mixed_ray_records(12)
    assert records["plane"][0]["m"] == 12
    assert records["F2"][0]["m"] == 11


def test_report_closes_every_link_and_excludes_adapted_hosts() -> None:
    report = REPORT.read_text(encoding="utf-8")
    required = (
        "complete Mori cone",
        "every **negative** extremal ray",
        "Square-zero tangent Mori rays are allowed",
        "projectivization has dimension \\(\\rho(S)-1\\) even when",
        "every literal nef facet preserves",
        "finite-volume Coxeter",
        "no proper reflection overgroup",
        "Relative ampleness forces width one",
        "Facet integrality removes satellites",
        "free root module",
        "nonsimpliciality causes no gap",
        "Exactly one nonexceptional facet root",
        r"\sum_i n_i^2=2",
        r"\sum_i n_i=4",
        r"\sum_i n_i^2=1",
        r"\sum_i n_i=3",
        r"\alpha_0-(H-e_1-e_2)=-e_3",
        r"\delta_0\) itself is an actual physical facet",
        "fibre-following free infinitely-near point",
        "this point is not satellite",
        r"w_{\mathbf P^2}=7H-3e_1-2e_2-\cdots-2e_{12}",
        r"w_{\mathbb F_2}=4s+9f-2e_1-\cdots-2e_{11}",
        r"\rho(S)-1\le11",
        "does **not** cover an adapted",
    )
    for phrase in required:
        assert phrase in report
    assert "Nikulin" not in report
    assert "Vinberg" not in report


def test_archive_replays_exactly() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
    )


def test_certifier_is_standalone() -> None:
    source = CERTIFIER.read_text(encoding="utf-8")
    assert "certify_rational_single_chain_first_mixed_ray_bound" not in source
    assert "certify_rational_i1n_inverse_design" not in source
    assert "nikulin" not in source.lower()
    assert "vinberg" not in source.lower()


def test_invalid_bounds_are_rejected() -> None:
    with pytest.raises(ValueError):
        build_certificate(11)
    with pytest.raises(ValueError):
        f0_swap_record(-1)
    with pytest.raises(ValueError):
        section_integrality_record(-1, 0)
    with pytest.raises(ValueError):
        section_integrality_record(2, -1)
    with pytest.raises(ValueError):
        exceptional_prime_record(-1)
    with pytest.raises(ValueError):
        degree_one_root_records(2)
    with pytest.raises(ValueError):
        exceptional_remainder_records(2)
    with pytest.raises(ValueError):
        plane_mixed_ray(11)
    with pytest.raises(ValueError):
        f2_mixed_ray(10)
