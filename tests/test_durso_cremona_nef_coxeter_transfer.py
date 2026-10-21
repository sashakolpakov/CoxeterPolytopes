from __future__ import annotations

import json

from scripts.certify_durso_cremona_nef_coxeter_transfer import (
    DEFAULT_OUTPUT,
    build_payload,
)


def test_exact_coxeter_dimensions_and_ray_splits() -> None:
    payload = build_payload()
    replay = payload["exact_arithmetic_replay"]
    assert replay["all_n_ge_10_coxeter_criterion"] == [10, 11, 13]
    expected = {
        "10": (17, 2, 0),
        "11": (26, 2, 0),
        "13": (44, 2, 0),
    }
    for dimension, counts in expected.items():
        row = replay["packages"][dimension]
        assert row["is_coxeter"]
        assert (
            row["ray_census"]["finite"],
            row["ray_census"]["ideal"],
            row["ray_census"]["hyperideal"],
        ) == counts
        assert row["checks"]["all_facets_irredundant"]
        assert row["checks"]["finite_volume_from_no_hyperideal_rays"]


def test_adjacent_noncoxeter_controls_are_retained() -> None:
    packages = build_payload()["exact_arithmetic_replay"]["packages"]
    assert not packages["12"]["is_coxeter"]
    assert packages["12"]["terminal_cosine_squared"] == "1/3"
    assert not packages["14"]["is_coxeter"]
    assert packages["14"]["terminal_cosine_squared"] == "1/5"
    assert packages["12"]["ray_census"]["hyperideal"] == 0
    assert packages["14"]["ray_census"]["hyperideal"] == 0


def test_h13_is_tumarkin_01_not_project_vinberg_p13() -> None:
    payload = build_payload()
    separation = payload["h13_lineage_separation"]
    assert separation["durso"]["name"] == "Tumarkin 13-tum04_13_01"
    assert separation["durso"]["ray_split"] == [44, 2, 0]
    assert separation["project_existing"]["ray_split"] == [43, 3, 0]
    assert separation["weighted_graphs_nonisomorphic"]
    h13_reflections = payload["exact_arithmetic_replay"]["packages"]["13"][
        "reflection_rows"
    ]
    assert not h13_reflections[-1]["integral_on_NS"]
    assert payload["claim_boundary"]["reflection_integral_on_NS_in_H13"] is False


def test_claim_boundary_and_archive_replay() -> None:
    payload = build_payload()
    assert payload["all_checks_pass"]
    assert payload["claim_boundary"]["new_abstract_polytopes"] is False
    assert payload["claim_boundary"]["full_nef_cone_claim"] is False
    assert payload["claim_boundary"]["strict_inverse_L4"] is False
    archived = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
    assert archived == payload
