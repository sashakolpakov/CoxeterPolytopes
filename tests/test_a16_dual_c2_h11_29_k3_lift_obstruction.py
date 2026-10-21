import json
from pathlib import Path

import pytest

from scripts.certify_a16_dual_c2_h11_29_k3_lift_obstruction import (
    OUTPUT,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "a16_dual_c2_h11_29_k3_lift_obstruction.md"


@pytest.fixture(scope="module")
def rebuilt():
    return build_certificate()


@pytest.fixture(scope="module")
def archived():
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly(rebuilt, archived):
    assert rebuilt == archived
    assert archived["all_checks_pass"] is True
    assert archived["case_id"] == "A16_6_1_C01"


def test_generic_picard_finite_order_denominator_is_closed(archived):
    generic = archived["generic_picard_no_go"]
    assert generic["target_rank_B_minus_I_mod_two"] == 2
    assert [row["order"] for row in generic["cyclotomic_exhaustion"]] == [
        1,
        2,
        3,
        4,
        6,
        7,
        9,
        14,
        18,
    ]
    assert generic["order_four_determinant_obstruction"] == {
        "contradiction": True,
        "determinant_B_minus_I": 8,
        "forced_divisor": 16,
        "required_even_smith_factor_count": 4,
    }


def test_nearest_special_glue_is_integral_but_fails_kahler(archived):
    special = archived["explicit_special_order_four_near_miss"]
    assert special["all_64_discriminant_classes_checked"] is True
    assert special["k3_even_unimodular_signature"] == [3, 19]
    assert special["k3_action_order"] == 4
    assert special["picard_rank_signature_determinant"] == [18, [1, 17], -16]
    assert special["fixed_NS_is_exactly_archived_fixed_lattice"] is True
    assert special["fixed_NS_rank"] == 12
    assert special["obstruction_root_gram"] == [["-2", "0"], ["0", "-2"]]
    assert special["both_roots_are_negated"] is True
    assert special["kahler_gate_passes"] is False


def test_hyperbolic_handoff_remains_the_certified_h11_29(archived):
    section = archived["hyperbolic_section_handoff"]
    assert section == {
        "bad_edge_count": 0,
        "bad_parabolic_neighborhood_count": 0,
        "compact": False,
        "dimension": 11,
        "elliptic_edge_count": 28798,
        "facet_count": 29,
        "finite_vertex_count": 3876,
        "finite_volume": True,
        "ideal_vertex_count": 253,
    }


def test_report_keeps_the_special_period_scope_open():
    text = REPORT.read_text(encoding="utf-8")
    assert "no finite-order lift with `NS(X)=N`" in text
    assert "does not prove" in text
    assert "`(-2)`-classes" in text
    assert "literal nef cone" in text
