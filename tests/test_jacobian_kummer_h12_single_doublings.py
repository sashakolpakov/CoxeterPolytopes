from __future__ import annotations

from functools import lru_cache
import json

from scripts.certify_jacobian_kummer_h12_single_doublings import OUTPUT, build


@lru_cache(maxsize=1)
def _payload():
    return build()


def test_archive_rebuilds_exactly() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True


def test_complete_doubling_wall_census() -> None:
    census = _payload()["complete_doubling_wall_census"]
    assert census["doubling_wall_count"] == 124
    assert census["non_doubling_wall_count"] == 48
    classes = census["single_double_facet_count_classes"]
    assert {key: len(value) for key, value in classes.items()} == {
        "296": 12,
        "317": 16,
        "326": 96,
    }


def test_three_exact_pairwise_distinct_h12_descendants() -> None:
    payload = _payload()
    rows = payload["representatives"]
    assert [row["facet_count"] for row in rows] == [296, 317, 326]
    assert all(
        row["all_exact_checks_pass"]
        and row["hyperbolic_dimension"] == 12
        and row["full_linear_rank"] == 13
        and row["finite_volume"]
        and not row["compact"]
        and row["integral_reflection_count"] == row["facet_count"]
        and row["strict_facet_witnesses"]["count"] == row["facet_count"]
        and len(row["boundary_records"]) == row["facet_count"]
        for row in rows
    )
    assert payload["decision"]["new_to_project_exact_polytope_count"] == 3
    assert payload["decision"]["infinite_H12_family"] is True
    redoubling = payload["infinite_redoubling_corollary"]
    assert redoubling["parallel_doubling_walls_one_based"] == [31, 36]
    assert redoubling["generator_squares_and_pairing"] == [-4, -4, 4]
    assert redoubling["all_translation_powers_distinct"] is True
    assert payload["scope"]["all_124_single_doubles_constructible"] is True
    assert payload["scope"]["all_124_labelled_isomorphism_types_classified"] is False
    assert payload["scope"]["infinite_parallel_redoubling_family_proved"] is True
    assert payload["scope"]["all_parameter_facet_formula_proved_here"] is False
