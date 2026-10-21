from __future__ import annotations

import json

from scripts.certify_alexeev_g0_split_lattice_census import (
    OUTPUT,
    SOURCE_SHA256,
    build,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_live_raw_gram_replay_matches_the_archive() -> None:
    assert build() == _archive()


def test_all_nine_complex_budgets_fail_and_five_are_absolute_gram_no_gos() -> None:
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert payload["raw_source"]["sha256"] == SOURCE_SHA256
    rows = {tuple(row["key"]): row for row in payload["rows"]}
    assert len(rows) == 9
    assert all(row["complex_carrier_excluded"] for row in rows.values())

    forest_expected = {
        (11, 11, 1): ([9, 1, 0], -1),
        (12, 10, 1): ([10, 1, 0], -8),
        (13, 9, 1): ([9, 1, 0], -16),
        (14, 8, 0): ([9, 1, 0], -16),
        (14, 8, 1): ([9, 1, 0], -64),
    }
    for key, (inertia, determinant) in forest_expected.items():
        row = rows[key]
        assert row["status"] == "no_positive_semidefinite_completion"
        assert row["representative_inertia_positive_negative_zero"] == inertia
        assert row["representative_determinant"] == determinant
        assert row["all_characteristic_split_system_excluded"] is True


def test_rank_minima_radical_parity_and_saturated_host_determinants() -> None:
    rows = {tuple(row["key"]): row for row in _archive()["rows"]}
    expected = {
        (15, 7, 1): (6, "E6", 21, 6),
        (16, 6, 1): (6, "A5 + A1", 22, 12),
        (18, 4, 0): (4, "A2 + A2", 22, 9),
        (20, 2, 1): (1, "A1", 21, 4),
    }
    for key, (minimum_rank, root_lattice, host_rank, determinant) in expected.items():
        row = rows[key]
        assert row["minimum_anti_invariant_rank"] == minimum_rank
        assert row["difference_root_lattice_H"] == root_lattice
        assert row["radical"]["every_relation_has_invariant_sum_in_2M"] is True
        assert row["checks"]["split_component_overlattice_integral"] is True
        assert row["checks"]["glue_rank_equals_difference_rank"] is True
        host = row["split_component_overlattice"]
        assert host["rank"] == host_rank
        assert host["absolute_determinant_K"] == determinant

    assert rows[(15, 7, 1)]["rank_minimality_certificate"][
        "positive_semidefinite_switching_classes"
    ] == 1
    assert rows[(18, 4, 0)]["rank_minimality_certificate"][
        "lower_bound_basis_determinant"
    ] == 9
    for key in ((16, 6, 1), (18, 4, 0), (20, 2, 1)):
        certificate = rows[key]["rank_minimality_certificate"]
        assert certificate["lower_bound_sign_support"]["forest"] is True
        assert certificate["lower_bound_sign_support"][
            "all_signings_switch_to_the_displayed_basis"
        ] is True
        assert certificate["all_candidate_sets_are_one_opposite_pair"] is True
        assert certificate[
            "opposite_choice_is_exactly_the_independent_sheet_swap"
        ] is True
        assert certificate["all_opposite_choices_remain_pairwise_compatible"] is True
        assert certificate["rank_minimal_switching_orbits"] == 1
    assert rows[(16, 6, 1)]["all_characteristic_split_system_excluded"] is True
    assert "characteristic 3" in _archive()["theorems"]["type_MII_host_forced"]


def test_scope_is_decorated_and_does_not_split_short_walls() -> None:
    payload = _archive()
    definition = payload["decorated_polytope_class"]["definition"]
    warning = payload["decorated_polytope_class"]["short_wall_warning"]
    assert "complete mixed-reflection chambers" in definition
    assert "long (-4)" in definition
    assert "short square-minus-two walls" in warning
    assert "not put through the forced split-root test" in warning
    assert "H10/22" in payload["scope"]["type_I_II_noncontradiction"]
    assert "H10/30" in payload["scope"]["type_I_II_noncontradiction"]
