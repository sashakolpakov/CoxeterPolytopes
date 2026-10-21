from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import pytest

from scripts.certify_shimada_aut3singk3_sparse_fixed_section import (
    CASES,
    DEFAULT_SOURCE,
    _classify_pair,
    _iter_rank3_vectors,
    analyze_scaled,
    build_case,
    resource_estimate,
)


def _fixture() -> tuple[
    tuple[tuple[int, ...], ...],
    tuple[int, ...],
    list[tuple[int, ...]],
    tuple[tuple[int, ...], ...],
]:
    rank = 20
    gram = tuple(
        tuple((1 if i == 0 else -1) if i == j else 0 for j in range(rank))
        for i in range(rank)
    )
    ample = (1,) + (0,) * (rank - 1)
    first = (1, -8, -8) + (0,) * (rank - 3)
    second = (1, -8, 1) + (0,) * (rank - 3)
    moving_left = (1, 0, 0, 1, 0) + (0,) * (rank - 5)
    moving_right = (1, 0, 0, 0, 1) + (0,) * (rank - 5)
    duplicate_fixed_restriction = (2, 0, 0, 1, 1) + (0,) * (rank - 5)
    action_rows = [list(row) for row in tuple(
        tuple(int(i == j) for j in range(rank)) for i in range(rank)
    )]
    action_rows[3][3] = action_rows[4][4] = 0
    action_rows[3][4] = action_rows[4][3] = 1
    action = tuple(tuple(row) for row in action_rows)
    return (
        gram,
        ample,
        [first, second, moving_left, moving_right, duplicate_fixed_restriction],
        action,
    )


def test_rank_three_stream_parser_selects_outer_vectors(tmp_path: Path) -> None:
    source = tmp_path / "sample.txt"
    source.write_text(
        "junk[0]:=0; A[2]:=[[[1, -2, 3/2], [4, 5, 6]], "
        "[[7, 8, 9], [-1/3, 0, 2]]];\n",
        encoding="utf-8",
    )
    assert list(_iter_rank3_vectors(source, "A", 2)) == [
        (0, (Fraction(1), Fraction(-2), Fraction(3, 2))),
        (0, (Fraction(4), Fraction(5), Fraction(6))),
        (1, (Fraction(7), Fraction(8), Fraction(9))),
        (1, (Fraction(-1, 3), Fraction(0), Fraction(2))),
    ]


def test_exact_pair_classifier_covers_every_rational_crystallographic_role() -> None:
    assert _classify_pair(-2, -8, 0) == "2"
    assert _classify_pair(-4, -4, 2) == "3"
    assert _classify_pair(-2, -4, 2) == "4"
    assert _classify_pair(-12, -4, 6) == "6"
    assert _classify_pair(-4, -4, 4) == "infinity"
    assert _classify_pair(-4, -4, 5) == "ultraparallel"
    assert _classify_pair(-4, -6, 4) == "forbidden"
    assert _classify_pair(-4, -4, -1) == "forbidden"


def test_sparse_fixture_orbits_deduplicates_witnesses_and_forbidden_pair() -> None:
    gram, ample, walls, action = _fixture()
    result = analyze_scaled(
        gram=gram,
        ample=ample,
        walls=walls,
        action=action,
    )
    assert all(result["source_action_checks"].values())
    assert result["complete_source_wall_count"] == 5
    assert result["wall_orbit_count"] == 4
    assert result["wall_orbit_size_histogram"] == {"1": 3, "2": 1}
    assert result["unique_restriction_count"] == 3
    assert result["restriction_square_sign_histogram"] == {
        "negative": 2,
        "zero": 0,
        "positive": 1,
    }
    assert result["strict_facet_count"] == 2
    assert result["projection_witness_failure_count"] == 0
    assert result["pair_gate"]["coxeter_edge_counts"]["forbidden"] == 1
    forbidden = result["pair_gate"]["first_forbidden_pairs"][0]
    assert forbidden["inner_product"] == -55
    assert forbidden["squares"] == [-127, -64]
    assert len(result["first_forbidden_pair_strict_witnesses"]) == 2
    assert all(
        witness["strict"]
        and witness["target_slack"] == 0
        and witness["minimum_other_slack"] > 0
        for witness in result["first_forbidden_pair_strict_witnesses"]
    )
    assert result["classification"] == "non_Coxeter_fixed_section"
    assert result["volume_status"].startswith("not tested")


def test_sparse_fixture_rejects_nonisometric_or_nonample_fixing_actions() -> None:
    gram, ample, walls, action = _fixture()
    bad_gram = [list(row) for row in gram]
    bad_gram[3][3] = -2
    with pytest.raises(ArithmeticError, match="does not preserve GramS"):
        analyze_scaled(
            gram=bad_gram,
            ample=ample,
            walls=walls,
            action=action,
        )

    moved_ample = list(ample)
    moved_ample[3] = 1
    with pytest.raises(ArithmeticError, match="does not fix the ample class"):
        analyze_scaled(
            gram=gram,
            ample=moved_ample,
            walls=walls,
            action=action,
        )


def test_x0_resource_model_avoids_quadratic_storage() -> None:
    estimate = resource_estimate(CASES["x0"])
    assert estimate["source_wall_count"] == 3278
    assert estimate["maximum_pair_count"] == 3278 * 3277 // 2
    assert estimate["dense_source_gram_entry_count_avoided"] == 3278**2
    assert estimate["retained_vector_coordinate_slots_upper_bound"] == 3 * 3278 * 20
    assert estimate["conservative_incremental_peak_mib"] < 32
    assert estimate["hard_operational_incremental_budget_mib"] == 64


def test_frozen_x0_terminal_forbidden_pair() -> None:
    archive = (
        Path(__file__).resolve().parents[1]
        / "data/enumerations/shimada_aut3singk3_x0_symplectic_h11_sparse_screen.json"
    )
    payload = json.loads(archive.read_text(encoding="utf-8"))
    section = payload["fixed_section"]
    assert payload["case"]["key"] == "x0"
    assert payload["all_checks_pass"] is True
    assert section["wall_orbit_count"] == 1734
    assert section["strict_facet_count"] == 314
    assert section["classification"] == "non_Coxeter_fixed_section"
    assert section["pair_gate"]["exhaustive"] is False
    forbidden = section["pair_gate"]["first_forbidden_pairs"]
    assert len(forbidden) == 1
    assert forbidden[0]["facet_pair_one_based"] == [5, 37]
    assert forbidden[0]["squares"] == [-4, -6]
    assert forbidden[0]["inner_product"] == 4
    assert forbidden[0]["normalized_square_ratio"] == "2/3"
    witnesses = section["first_forbidden_pair_strict_witnesses"]
    assert len(witnesses) == 2
    assert all(
        witness["strict"]
        and witness["target_slack"] == 0
        and witness["minimum_other_slack"] > 0
        for witness in witnesses
    )


@pytest.mark.skipif(
    not (DEFAULT_SOURCE / "AutXa.txt").exists(),
    reason="pinned Shimada supplement is not present",
)
@pytest.mark.parametrize("case_key", ["x2-rational", "x2-symplectic", "x1-symplectic"])
def test_sparse_engine_replays_frozen_x1_x2(case_key: str) -> None:
    payload = build_case(CASES[case_key], DEFAULT_SOURCE)
    assert payload["all_checks_pass"] is True
    assert all(payload["frozen_replay_checks"].values())
    assert payload["algorithm"]["dense_pairing_matrix_materialized"] is False
