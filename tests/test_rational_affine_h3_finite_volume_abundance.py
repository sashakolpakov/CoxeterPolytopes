"""Exact and scope checks for finite-volume affine-H3 abundance."""

from __future__ import annotations

from pathlib import Path

from scripts.certify_rational_affine_h3_finite_volume_abundance import certificate


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_affine_h3_finite_volume_abundance.md"


def test_affine_star_null_cycle_host_scaling_and_minimum_rank() -> None:
    cert = certificate()
    assert cert["physical_null_cycle"] == [4, 1, 1, 1, 1]
    assert cert["host_null_root"] == [2, 1, 1, 1, 1]
    assert cert["star_gram_rank"] == 4
    assert cert["minimum_lattice_picard_rank"] == 6
    assert cert["minimum_geometric_picard_rank"] == 11
    first = cert["rank_rows"][0]
    assert first == {
        "rho": 6,
        "cusp_dimension": 4,
        "affine_D4_factor_dimension": 4,
        "outside_factor_dimension": 0,
        "K_square": 4,
        "chi_minus_K": 5,
        "forbidden_by_effective_minus_K": True,
    }


def test_positive_root_relation_clears_half_boundary_denominators() -> None:
    relation = certificate()["sample_denominator_clearing"]
    assert relation["effective_curve_coefficients"] == [2, 1, 2]
    assert relation["F_coefficient"] == 3
    assert relation["root_to_curve_scales"] == ["1", "1/2", "1"]


def test_all_apparent_low_rank_endpoints_are_empty() -> None:
    rows = certificate()["rank_rows"]
    forbidden = [row for row in rows if row["forbidden_by_effective_minus_K"]]
    assert [row["rho"] for row in forbidden] == [6, 7, 8, 9, 10]
    assert all(row["chi_minus_K"] >= 1 for row in forbidden)
    first_possible = next(row for row in rows if not row["forbidden_by_effective_minus_K"])
    assert first_possible["rho"] == 11
    assert first_possible["outside_factor_dimension"] == 5


def test_nonreduced_picard_descent_forces_index_one_and_genus_three() -> None:
    descent = certificate()["picard_descent"]
    assert descent == {
        "reduced_components": 5,
        "reduced_nodes": 4,
        "dual_graph_first_betti": 0,
        "component_multidegrees_O_F_F": [0, 0, 0, 0, 0],
        "nilpotent_kernel_type": "complex_vector_group",
        "nilpotent_kernel_has_nontrivial_torsion": False,
        "O_F_F_power_trivial": True,
        "O_F_F_trivial": True,
        "F_base_point_free": True,
    }
    assert certificate()["fibration_case"] == {
        "index": 1,
        "genus": 3,
        "K_dot_general_fibre": 4,
        "minimal_core_rho_bound": 26,
    }


def test_scope_does_not_promote_abundance_to_tail_or_bare_cone_theorem() -> None:
    scope = certificate()["scope"]
    assert scope == {
        "bare_rational_polyhedrality_suffices": False,
        "literal_finite_volume_full_nef_suffices": True,
        "reflection_maximality_used": False,
        "vertical_tail_bound": False,
        "fibration_index_forced_one": True,
    }
    text = " ".join(REPORT.read_text().split())
    for phrase in (
        "no reflection maximality",
        "No claim about the bare variant is made here",
        "The nonreduced Picard group kills the index",
        "remaining affine obstruction is the global vertical-tail bound",
        "apparent rank-six through rank-ten endpoints are empty",
        "The former genus-two/divisibility-two alternative does not occur",
    ):
        assert phrase in text


def test_report_records_why_the_tvav_endpoint_is_vacuous() -> None:
    text = " ".join(REPORT.read_text().split())
    assert "Lemma 2.5, is correctly stated" in text
    assert "preceding Lemma 2.4 forces" in text
    assert "does not supply a nonvacuous affine-H3 endpoint" in text
