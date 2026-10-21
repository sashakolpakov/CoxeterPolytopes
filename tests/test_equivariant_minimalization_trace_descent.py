from pathlib import Path

from scripts.certify_equivariant_minimalization_trace_descent import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/equivariant_minimalization_trace_descent_theorem.md"


def certificate():
    return build_certificate()


def test_certificate_passes():
    assert certificate()["all_checks_pass"]


def test_gamma16_is_the_reflection_maximal_control():
    payload = certificate()
    assert all(payload["maximality_checks"].values())
    assert payload["theorem_boundary"]["trace_lattice"] == (
        "pi^*(U + E8(-1) + E7(-1))"
    )


def test_exceptional_rank_is_the_free_blowup_parameter():
    rows = certificate()["sample_rows"]
    assert [row["blowups"] for row in rows] == [0, 1, 2, 5, 11]
    assert all(row["orthogonal_complement_rank"] == row["blowups"] for row in rows)
    assert all(row["canonical_square"] == -row["blowups"] for row in rows)


def test_every_extended_wall_reflection_fixes_the_canonical_vector():
    rows = certificate()["sample_rows"]
    assert all(
        row["checks"][
            "all_19_reflections_integral_isometric_involutions_fixing_K"
        ]
        for row in rows
    )


def test_exact_correspondence_formula_has_the_claimed_block_action():
    rows = certificate()["sample_rows"]
    assert all(
        row["checks"][
            "all_19_Gamma16_correspondence_formulas_equal_block_reflections"
        ]
        for row in rows
    )


def test_raw_wall_gram_does_not_change_under_blowup():
    rows = certificate()["sample_rows"]
    assert all(row["checks"]["trace_block_unchanged"] for row in rows)
    assert all(
        row["checks"]["trace_and_exceptional_blocks_orthogonal"] for row in rows
    )


def test_minimal_K3_and_K3_cover_rows_have_zero_defect():
    calibration = certificate()["D_pi_calibration"]
    rows = calibration["minimal_K3_and_K3_cover_zero_controls"]
    assert [row["id"] for row in rows] == [
        "wehler_h2_3",
        "bugaenko_h8_11",
        "coble_type_ii_h10_30",
        "coble_type_i_h11_26",
        "coble_type_vi_char3_h14_60",
        "coble_type_mii_char3_h17_48",
    ]
    assert all(row["D_pi_rank"] == 0 and row["D_pi_zero"] for row in rows)
    assert all(
        row["control_strength"] == "vacuous identity-contraction zero control"
        for row in rows
    )
    assert [row["complex_rank_20_scope"] for row in rows] == [
        True,
        True,
        True,
        True,
        False,
        False,
    ]


def test_rational_full_NS_controls_have_full_rank_nonzero_defect():
    calibration = certificate()["D_pi_calibration"]
    rows = calibration["rational_full_NS_nonzero_controls"]
    assert [row["D_pi_rank"] for row in rows] == [6, 10, 11, 14, 16]
    assert all(not row["D_pi_zero"] for row in rows)
    assert all(not row["canonical_operator_scope"] for row in rows)
    assert all(
        row["D_pi_shape"]
        == [row["trace_rank"], row["exceptional_frame_rank"]]
        for row in rows
    )
    assert all(
        len(row["D_pi_matrix"]) == row["D_pi_shape"][0]
        and all(
            len(matrix_row) == row["D_pi_shape"][1]
            for matrix_row in row["D_pi_matrix"]
        )
        for row in rows
    )
    assert all(calibration["blowdown_source_checks"].values())


def test_report_has_integral_descent_and_scope_guards():
    raw = REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)

    report = raw.decode("utf-8")
    assert r"T_{\mathbf Z}\subset T_{\mathbf Q}" in report
    assert r"\operatorname {Num}(Y)" in report
    assert r"E_\pi^\perp=\pi^*\operatorname {Num}(X)" in report
    assert "not merely the\nrational splitting (1.7)" in report
    assert r"D_\pi=J^tG_YC_\pi" in report
    assert "finite matrix gate **per certified contraction frame**" in report
    assert "not prove that all geometric contraction frames form a finite" in report
    assert "**defect-only negative controls**" in report
    assert "no statement about reflection maximality is inferred" in report
    assert "across all\nintrinsic correspondence carriers" in report
