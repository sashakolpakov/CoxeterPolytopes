from __future__ import annotations

import json

from scripts.independent_finite_volume_coxeter_dessin_cusp_hostile_audit import (
    OUTPUT,
    build_certificate,
)


def test_hostile_archive_matches_helper_free_builder() -> None:
    archived = json.loads(OUTPUT.read_text())
    assert archived == build_certificate()
    assert archived["all_checks_pass"]
    assert archived["verdict"] == (
        "PASS_AFTER_BETA_ZERO_REPAIR_WITH_HYPOTHESES_MADE_EXPLICIT"
    )


def test_beta_zero_information_loss_and_repair_are_literal() -> None:
    audit = build_certificate()["identity_audit"]
    counterexample = audit["beta_zero_information_loss_counterexample"]
    assert counterexample["x_square"] == "1"
    assert counterexample["x_pair_r"] == "-1"
    assert counterexample["boxed_norm_equality_holds"]
    assert not counterexample["mirror_equation_holds"]
    assert audit["exact_grid_trials"] == 300


def test_complete_toric_ray_and_cusp_replay() -> None:
    replay = build_certificate()["toric_h6_9_replay"]
    assert replay["active_subset_trials"] == 84
    assert replay["extreme_ray_count"] == 21
    assert replay["finite_ideal_hyperideal"] == [18, 3, 0]
    assert len(replay["cusps"]) == 3
    for cusp in replay["cusps"]:
        assert cusp["active_gram_rank"] == 5
        assert cusp["compact_rank_five_horosphere_cell"]
        assert len(cusp["vertical_wall_labels"]) == 7
        assert len(cusp["bounded_wall_labels"]) == 2
        assert cusp["bounded_height_cap_squares"] == ["2", "2"]
        assert sorted(
            (row["node_count"], row["rank"], row["positive_kernel_relation"])
            for row in cusp["affine_components"]
        ) == [(3, 2, [1, 2, 1]), (4, 3, [1, 1, 1, 1])]


def test_hostile_certifier_has_no_project_import() -> None:
    source = (
        __import__(
            "scripts.independent_finite_volume_coxeter_dessin_cusp_hostile_audit",
            fromlist=["ROOT"],
        )
        .__file__
    )
    text = open(source, encoding="utf-8").read()
    assert "from scripts" not in text
    assert "import scripts" not in text
