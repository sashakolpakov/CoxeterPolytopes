"""Independent hostile checks for finite-volume affine-H3 abundance."""

from __future__ import annotations

from pathlib import Path

from scripts.certify_rational_affine_h3_finite_volume_abundance_hostile import (
    certificate,
)


ROOT = Path(__file__).resolve().parents[1]
OWNER = ROOT / "results/rational_affine_h3_finite_volume_abundance.md"
AUDIT = ROOT / "results/rational_affine_h3_finite_volume_abundance_hostile_audit.md"


def test_hostile_replay_is_independent_and_preserves_the_abundance_core() -> None:
    cert = certificate()
    assert cert["verdict"] == "PASS_AFTER_REPAIR"
    assert cert["owner_core"] == {
        "finite_volume_implies_F_base_point_free": True,
        "reflection_maximality_used": False,
        "external_classification_used": False,
    }
    assert cert["physical_null_cycle"] == [4, 1, 1, 1, 1]
    assert cert["host_null_cycle"] == [2, 1, 1, 1, 1]


def test_termwise_orthogonality_and_half_root_clearing_are_exact() -> None:
    cert = certificate()
    assert cert["termwise_zero_profiles"] == [[0, 0, 0, 0, 0]]
    relation = cert["sample_positive_relation"]
    assert relation["coefficients"] == [1, 1, 1]
    assert relation["physical_root_scales"] == ["1", "1/2", "1"]
    assert relation["cleared_physical_coefficients"] == [2, 1, 2]


def test_rank_six_through_ten_are_empty_not_abundant_endpoints() -> None:
    rows = certificate()["empty_endpoint_rows"]
    assert [row["rho"] for row in rows] == [6, 7, 8, 9, 10]
    assert [row["chi_minus_K"] for row in rows] == [5, 4, 3, 2, 1]
    assert all(row["effective_minus_K"] for row in rows)
    assert all(row["minus_K_dot_F"] == -4 for row in rows)
    assert all(row["contradicts_nef_F"] for row in rows)
    assert certificate()["repairs"]["minimum_possible_picard_rank"] == 11


def test_nonreduced_picard_gate_forces_the_moving_index_to_one() -> None:
    cert = certificate()
    assert cert["supported_fibre_cycles"] == [
        [4, 1, 1, 1, 1],
        [8, 2, 2, 2, 2],
        [12, 3, 3, 3, 3],
    ]
    gate = cert["picard_gate"]
    assert gate["components"] == 5
    assert gate["nodes"] == 4
    assert gate["dual_graph_first_betti"] == 0
    assert gate["multidegrees_O_F_F"] == [0, 0, 0, 0, 0]
    assert gate["nilpotent_kernel_torsion_free"] is True
    assert gate["power_trivial_normal_bundle_implies_trivial"] is True
    assert gate["F_base_point_free"] is True
    assert cert["fibration_case"] == {
        "moving_index": 1,
        "general_fibre_class": "F",
        "genus": 3,
        "relative_minimal_rho_bound": 26,
    }
    assert cert["repairs"]["nonreduced_picard_forces_index_one"] is True
    assert cert["repairs"]["genus_two_branch"] is False
    assert cert["repairs"]["relative_minimal_core_bound_26"] is True


def test_fixed_exponent_two_lattice_does_not_control_null_divisibility() -> None:
    poison = certificate()["lattice_poison"]
    assert poison["ambient"] == "U_plus_D4_minus_one"
    assert poison["determinant_absolute"] == 4
    assert poison["discriminant_exponent"] == 2
    assert [row["N"] for row in poison["rows"]] == list(range(1, 13))
    assert all(row["same_affine_D4_gram"] for row in poison["rows"])
    assert poison["rows"][-1]["null_vector"] == [12, 0, 0, 0, 0, 0]
    assert certificate()["repairs"]["lattice_alone_forces_index_one"] is False


def test_owner_and_audit_have_the_repaired_scope() -> None:
    owner = " ".join(OWNER.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        r"\rho(S)\ge11.",
        r"\(|F|\) itself is base-point-free",
        r"2g(H)-2=K_S.F=4,\qquad g(H)=3.",
        r"\rho(S_{\min})\le26.",
        "does not supply a nonvacuous affine-H3 endpoint",
        "remaining affine obstruction is the global vertical-tail bound",
    ):
        assert phrase in owner
    for phrase in (
        "PASS AFTER REPAIR, with semiampleness upgraded to base-point-freeness",
        "The finite-volume cusp argument really does produce",
        "The first draft nevertheless contained two substantive downstream errors",
        "The fibration-index attack and its AG repair",
        "Lattice poison: why the AG step is essential",
        r"Therefore \(\rho(S_{\min})\le26\)",
    ):
        assert phrase in audit


def test_tvav_primary_source_is_pinned_and_used_only_as_a_vacuity_check() -> None:
    source = certificate()["source_audit"]
    assert source["arxiv"] == "0901.1094v2"
    assert len(source["source_archive_sha256"]) == 64
    assert len(source["main_tex_sha256"]) == 64
    assert source["lemma_2_4"].startswith("big -K forces")
    assert source["lemma_2_5"].startswith("nef non-big N")
    assert source["theorem_2_9"].startswith("big rational surfaces")
