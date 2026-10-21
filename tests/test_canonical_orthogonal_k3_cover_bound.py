import json
from pathlib import Path

import sympy as sp

from scripts.certify_canonical_orthogonal_k3_cover_bound import (
    LEDGER,
    carrier_ceiling,
    certify,
    certify_pullback_pairing,
    coble_dimension,
    coble_k3_boundary_ceiling,
    coble_k3_dimension_ceiling,
    coble_rank,
    dimension_ceiling,
    known_coble_boundary_ceiling,
    known_coble_dimension_ceiling,
    single_root_reflection,
    split_root_restriction_identity,
    trace_rank,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/canonical_orthogonal_k3_cover_dimension_theorem.md"


def test_full_certificate():
    certify()


def test_weakest_rank_arithmetic():
    for characteristic in (0, 3, 5):
        for degree in range(1, 8):
            for dimension in range(2, 20):
                rank = trace_rank(degree, dimension)
                assert rank == degree * (dimension + 1)
                if dimension <= dimension_ceiling(degree, characteristic):
                    assert rank <= carrier_ceiling(characteristic)


def test_strong_coble_cover_ceiling():
    for boundary_count in range(1, 13):
        assert coble_rank(boundary_count) == 10 + boundary_count
        assert coble_dimension(boundary_count) == 9 + boundary_count

    # The cover/Picard calculation alone is weaker in positive
    # characteristic: it permits b <= 12 numerically.  Do not confuse that
    # arithmetic ceiling with existence of a strong Coble surface.
    assert coble_k3_boundary_ceiling(0) == 10
    assert coble_k3_dimension_ceiling(0) == 19
    assert coble_k3_boundary_ceiling(3) == 12
    assert coble_k3_dimension_ceiling(3) == 21

    # The independent class-wide Coble theorem gives b <= 10 in every
    # characteristic != 2, so the best known nef-dimension ceiling is 19.
    for characteristic in (0, 3, 5):
        assert known_coble_boundary_ceiling(characteristic) == 10
        assert known_coble_dimension_ceiling(characteristic) == 19


def test_root_and_paired_root_lifts():
    assert single_root_reflection() == sp.Matrix([[-1]])
    assert split_root_restriction_identity()
    certify_pullback_pairing()


def test_action_levels_are_not_conflated():
    data = json.loads(LEDGER.read_text())
    rows = {row["id"]: row for row in data["rows"]}

    assert rows["wehler_h2_3"]["actual_automorphisms"]
    assert rows["bugaenko_h8_11"]["action_level"] == "M_coherent"

    for key in ("coble_type_ii_h10_30", "coble_type_i_h11_26"):
        assert rows[key]["action_level"] == "M_generatorwise"
        assert not rows[key]["coherent_coxeter_monodromy_representation_certified"]

    for key in (
        "coble_type_vi_char3_h14_60",
        "coble_type_mii_char3_h17_48",
    ):
        assert rows[key]["action_level"] == "C"
        assert not rows[key]["positive_characteristic_monodromy_certified"]


def test_seven_enriques_rows_are_one_sharp_rank_ten_class():
    data = json.loads(LEDGER.read_text())
    row = next(
        row for row in data["rows"] if row["id"] == "seven_complex_enriques_full_nef_h9"
    )
    assert row["row_count"] == 7
    assert row["trace_rank"] == row["n"] + 1 == 10
    assert row["facets"] == [12, 12, 20, 20, 20, 20, 20]


def test_report_has_clean_tex_and_explicit_scope_guards():
    raw = REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)

    report = raw.decode("utf-8")
    assert "(operatorname" not in report
    assert "Here **strong Coble** includes\nrationality" in report
    assert r"\jmath:\widetilde N\to\operatorname {NS}(X)" in report
    assert r"gives only \(b\leq12\) in positive characteristic" in report
    assert "best known\nclass-wide projectivized-nef ceiling is \\(19\\)" in report
