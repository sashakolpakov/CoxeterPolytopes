import json
from pathlib import Path

import pytest

from scripts.certify_rational_h2_minus4_minus2_genus_one_core import (
    ARCHIVE,
    build_certificate,
    reflection_matrix,
    tail_row,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_h2_minus4_minus2_genus_one_core.md"


def test_exact_h2_pair_and_sharp_core_bound():
    payload = build_certificate()
    pair = payload["pair_arithmetic"]
    assert pair["intersection_matrix"] == [[-4, 2], [2, -2]]
    assert pair["determinant"] == "4"
    assert pair["discrepancies"] == ["-1", "-1"]
    assert pair["arithmetic_genus_D"] == 1
    assert pair["adjoint_multidegrees_on_B_R"] == [0, 0]
    assert pair["adjoint_square_formula"] == "12-rho"
    equality = payload["sharp_equality_control"]
    assert equality["strict_transform_squares"] == [-4, -2]
    assert equality["picard_rank"] == 12
    assert equality["K_plus_B_plus_R"] == "0"


def test_free_tail_generators_are_the_Bm_signed_permutation_generators():
    for length in (1, 2, 3, 7, 19, 40):
        row = tail_row(length)
        assert row["root_squares_qM"] == [-4] * (length - 1) + [-2]
        assert row["root_gram_determinant"] == str((-2) ** length)
        assert row["simple_reflections_are_adjacent_swaps_then_terminal_sign"]


def test_tail_and_reflection_input_rejection():
    for bad in (0, -1, True, 1.5, "3"):
        with pytest.raises(ValueError):
            tail_row(bad)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        reflection_matrix([0, 0])


def test_archive_and_scope_controls():
    payload = build_certificate()
    assert json.loads(ARCHIVE.read_text(encoding="utf-8")) == payload
    assert payload["scope"]["core_only"] is True
    assert payload["scope"]["original_rank_bound_claimed"] is False
    assert payload["scope"]["reflection_maximal_counterfamily_claimed"] is False
    assert payload["conditional_total_marked_degree_bound"] == {
        "premise": (
            "every contracted tail component is a physical facet and "
            "T is the total integral ample-marking degree of all "
            "physical facets"
        ),
        "tail_length_bound": "m<=T",
        "picard_rank_bound": "rho(S)<=12+T",
        "hyperbolic_dimension_bound": "dim<=11+T",
        "per_facet_degree_bound_alone_is_sufficient": False,
    }
    assert payload["conditional_trace_bound"] == {
        "rank": "d(n+1)",
        "bound_if_the_trace_block_descends": "d(n+1)<=12",
        "rational_field_bound": "n<=11",
        "descent_is_automatic": False,
    }


def test_report_is_clean_and_states_the_exact_open_gate():
    raw = REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)
    text = raw.decode("utf-8")
    for needle in (
        r"\rho(S_0)\leq12",
        r"d(n+1)\leq12",
        r"m\leq T",
        r"\rho(S)\leq12+T",
        "does not bound",
        "bound on each individual facet degree",
        "signed permutation",
        "reflection-maximal counterfamily",
        "two distinct transverse points",
        "tangency",
    ):
        assert needle in text
