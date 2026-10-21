from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_roulleau_rank12_bh_k3_lifts import (
    OUTPUT,
    build_certificate,
)


def certificate() -> dict:
    return build_certificate()


def test_four_effective_realizations_of_three_polytopes() -> None:
    data = certificate()
    assert data["all_checks_pass"] is True
    assert data["effective_K3_group_realization_count"] == 4
    assert data["distinct_polytope_realization_count"] == 3
    assert [
        row["brandhorst_hofmann"]["id"] for row in data["realizations"]
    ] == ["0.4.4.1", "0.4.4.2", "0.4.4.6", "1.4.4.2"]


def test_full_integral_action_conjugacy_and_rootless_coinvariant() -> None:
    data = certificate()
    for row in data["realizations"]:
        conjugacy = row["exact_action_conjugacy"]
        assert conjugacy["full_NS_action_pairs_integrally_conjugate"] is True
        assert conjugacy["all_glue_pairs_carried"] is True
        assert conjugacy["algebraic_coinvariant_rootless"] is True
        assert conjugacy["fixed_discriminant_isometry"]["bijective"] is True
        assert (
            conjugacy["fixed_discriminant_isometry"][
                "quadratic_form_preserved_mod_2Z"
            ]
            is True
        )


def test_all_three_vinberg_polytopes_are_finite_volume_noncompact() -> None:
    data = certificate()
    census = {
        row["polytope_id"]: (
            row["coxeter_polytope"]["facets"],
            row["coxeter_polytope"]["finite_vertices"],
            row["coxeter_polytope"]["cusps"],
            row["coxeter_polytope"]["elliptic_edges"],
        )
        for row in data["realizations"]
    }
    assert census == {
        "R12-E8A1x2-C2-H10-12": (12, 17, 2, 99),
        "R12-D8A1x2-C2a-H10-15": (15, 42, 10, 294),
        "R12-D4x2A1x2-p67-H10-38": (38, 526, 141, 4518),
    }
    assert all(
        row["coxeter_polytope"]["classification"]
        == "finite_volume_noncompact"
        and row["coxeter_polytope"]["Vinberg_bad_audits"] == 0
        for row in data["realizations"]
    )


def test_full_group_enrichment_is_effective_v4_on_picard() -> None:
    full = certificate()["full_group_enrichment"]
    assert full["group"] == "C4 x C2"
    assert full["NS_image"] == "C2 x C2"
    assert full["common_fixed_rank_signature_determinant"] == [
        4,
        [1, 3, 0],
        -16,
    ]
    assert full["common_coinvariant_rank_signature_determinant"] == [
        8,
        [0, 8, 0],
        256,
    ]
    assert full["checks"]["common_coinvariant_rootless"] is True


def test_exhaustive_rank12_rank13_filter_is_recorded() -> None:
    scan = certificate()["finite_positive_filter"]
    assert (scan["pure_rows_scanned"], scan["full_rows_scanned"]) == (587, 3561)
    assert (scan["pure_rank_12_or_13_rows"], scan["full_rank_12_or_13_rows"]) == (
        53,
        92,
    )
    assert scan["pure_exact_action_hits"] == ["0.4.4.1", "0.4.4.2", "0.4.4.6"]
    assert scan["full_exact_action_hits"] == ["1.4.4.2"]
    assert scan["rank13_exact_action_hits"] == []


def test_archived_json_is_exactly_reproducible() -> None:
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert archived == certificate()
