import json
from pathlib import Path

import pytest

from scripts.audit_rational_h1_intersecting_minus4_log_core import (
    ARCHIVE,
    EXPECTED_LI_ARCHIVE_SHA256,
    EXPECTED_LI_TEX_SHA256,
    build_certificate,
    build_li_h1_geometry,
    free_tail,
)


ROOT = Path(__file__).resolve().parents[1]
AUDIT_REPORT = (
    ROOT / "results/rational_h1_intersecting_minus4_log_core_hostile_audit.md"
)


def test_core_identity_discrepancies_and_scope():
    payload = build_certificate()
    core = payload["core_arithmetic"]
    assert core["pair_matrix"] == [[-4, 2], [2, -4]]
    assert core["pair_determinant"] == "12"
    assert core["arithmetic_genus"] == 1
    assert core["adjoint_multidegrees"] == [0, 0]
    assert core["adjoint_square_formula"] == "14-rho_core"
    assert core["discrepancies"] == ["-1", "-1"]
    assert payload["verdict"] == "PASS_with_scope_warning"


def test_li_geometry_reconstructed_in_total_transform_coordinates():
    geometry = build_li_h1_geometry()
    assert geometry["picard_rank"] == 14
    assert len(geometry["physical_names"]) == 16
    assert geometry["physical_species_counts"] == {
        "B(-4)": 2,
        "R(-2)": 12,
        "E(-1)": 2,
    }
    assert geometry["host_root_norm_counts"] == {"-2": 4, "-4": 12}
    assert geometry["host_determinant_absolute"] == str(2**10)
    assert geometry["all_wall_reflections_integral"] is True
    assert geometry["all_wall_reflections_isometric"] is True
    assert geometry["unsupported_wall_pairs"] == []
    assert geometry["coxeter_order_pair_counts"] == {"2": 104, "3": 14, "4": 2}


def test_free_tail_is_saturated_Bm_not_a_missing_reflection_argument():
    for length in (1, 2, 3, 9, 32):
        row = free_tail(length)
        assert row["determinant"] == str((-2) ** length)
        assert row["generated_group"] == "W(B_m)=signed permutations"
        assert row["generated_group_order"] == str(
            2**length * __import__("math").factorial(length)
        )
    for bad in (0, -2, True, 1.5, "4"):
        with pytest.raises(ValueError):
            free_tail(bad)  # type: ignore[arg-type]


def test_source_hashes_archive_and_conditional_statements():
    payload = build_certificate()
    assert payload["li_source_provenance"]["archive_sha256"] == (
        EXPECTED_LI_ARCHIVE_SHA256
    )
    assert payload["li_source_provenance"]["main_tex_sha256"] == (
        EXPECTED_LI_TEX_SHA256
    )
    assert json.loads(ARCHIVE.read_text(encoding="utf-8")) == payload
    assert payload["conditional_consequences"] == {
        "total_marked_facet_degree": "m<=T and rho(S)<=14+T",
        "trace_descent": "d(n+1)<=rho(core)<=14",
        "rational_field": "n<=13",
        "neither_hypothesis_is_automatic": True,
    }


def test_audit_report_bytes_and_no_maximality_overclaim():
    raw = AUDIT_REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)
    text = raw.decode("utf-8")
    for needle in (
        "PASS, with one scope warning",
        "rank-14",
        "two distinct transverse points",
        "ordinary tangency",
        "2 + 12 + 2",
        "4 + 12",
        "signed-permutation",
        "does not prove an original-rank bound",
        "unbounded reflection-maximal family",
    ):
        assert needle in text
