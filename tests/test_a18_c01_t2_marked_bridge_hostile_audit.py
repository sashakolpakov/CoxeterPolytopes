from __future__ import annotations

from functools import lru_cache
import json

import sympy as sp

from scripts.audit_a18_c01_t2_marked_bridge_hostile import (
    CONSTRUCTION,
    FACTORY,
    HOSTILE,
    OUTPUT,
    PROMOTION,
    SELECTED_SHEET_MASK,
    audit,
)


def _matrix(rows: list[list[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(str(value)) for value in row] for row in rows])


@lru_cache(maxsize=1)
def _payload() -> dict:
    return audit()


def test_marked_bridge_archive_replays_exactly() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["case_id"] == "A18_4_0_C01"
    assert payload["sheet_order_search"]["searched_mask_count"] == 1024
    assert payload["sheet_order_search"]["exact_matching_masks"] == [
        192,
        440,
        583,
        831,
    ]
    assert payload["sheet_order_search"]["selected_mask"] == SELECTED_SHEET_MASK


def test_literal_form_action_and_ordered_24_root_equalities() -> None:
    payload = _payload()
    hostile = json.loads(HOSTILE.read_text(encoding="utf-8"))
    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))

    candidate = hostile["unique_candidate"]
    hostile_gram = _matrix(candidate["rank18"]["gram"])
    hostile_action = _matrix(candidate["rank18"]["involution"])
    hostile_roots = _matrix(candidate["fold_back"]["source_roots"])
    inventory = next(
        row
        for row in construction["inventory"]["examples"]
        if row["key"] == [18, 4, 0]
    )
    archive_gram = _matrix(inventory["ambient"]["gram"])
    wall_rows = _matrix(inventory["walls"]["coordinates_in_leading_root_basis"])
    factory_case = next(
        row for row in factory["constructions"] if row["id"] == "A18_4_0_C01"
    )
    archive_row_action = _matrix(
        factory_case["section"]["case"][
            "representative_nonidentity_primal_actions"
        ][0]
    )
    source_ids = payload["sheet_order_search"][
        "selected_ordered_factory_source_wall_indices_one_based"
    ]
    archive_roots = wall_rows[[index - 1 for index in source_ids], :].T
    marked_isometry = _matrix(payload["marked_isometry"]["matrix_P"])

    assert marked_isometry.det() == 1
    assert all(value.q == 1 for value in marked_isometry)
    assert all(value.q == 1 for value in marked_isometry.inv())
    assert marked_isometry.T * archive_gram * marked_isometry == hostile_gram
    assert marked_isometry * hostile_action == archive_row_action.T * marked_isometry
    assert marked_isometry * hostile_roots == archive_roots


def test_full_48_wall_action_and_fourteen_labelled_orbits() -> None:
    payload = _payload()
    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    inventory = next(
        row
        for row in construction["inventory"]["examples"]
        if row["key"] == [18, 4, 0]
    )
    walls = _matrix(inventory["walls"]["coordinates_in_leading_root_basis"])
    factory_case = next(
        row for row in factory["constructions"] if row["id"] == "A18_4_0_C01"
    )
    action = _matrix(
        factory_case["section"]["case"][
            "representative_nonidentity_primal_actions"
        ][0]
    )
    wall_lookup = {tuple(walls.row(index)): index + 1 for index in range(48)}
    permutation = [
        wall_lookup[tuple(walls.row(index) * action)] for index in range(48)
    ]
    assert permutation == payload["source_chamber"]["action_permutation_one_based"]
    assert sorted(permutation) == list(range(1, 49))
    assert all(permutation[permutation[index] - 1] == index + 1 for index in range(48))

    records = payload["active_marked_orbits"]
    assert len(records) == 14
    assert [row["factory_orbit_index"] for row in records] == [
        1,
        3,
        4,
        5,
        6,
        8,
        9,
        10,
        11,
        12,
        13,
        16,
        17,
        20,
    ]
    assert sum(len(row["factory_source_wall_indices_ordered"]) for row in records) == 24
    assert all(row["orbit_and_folded_sum_transport_exact"] for row in records)
    for row in records:
        ordered = row["factory_source_wall_indices_ordered"]
        expected_images = ordered if len(ordered) == 1 else list(reversed(ordered))
        assert row["archive_action_images_in_order"] == expected_images


def test_bridge_imports_only_datum_C_and_V_while_variable_U_stays_open() -> None:
    payload = _payload()
    promotion = json.loads(PROMOTION.read_text(encoding="utf-8"))
    promotion_case = next(
        row for row in promotion["cases"] if row["id"] == "A18_4_0_C01"
    )
    imported = payload["imported_existing_C01_conclusions"]
    assert imported["strong_promotion_action_sha256"] == promotion_case[
        "ambient_NS"
    ]["action_sha256"]
    assert imported["action_digest_matches_exact_bridge_target"] is True
    assert imported["C_complete_invariant_nef_section"] is True
    assert imported["V_finite_volume_noncompact"] is True

    ledger = payload["obligation_ledger"]
    assert ledger["IC0"].startswith("closed for the single pinned")
    assert ledger["R"].startswith("closed for this labelled C01 datum")
    assert ledger["C"].startswith("imported for this datum")
    assert ledger["V"].startswith("imported for this datum")
    assert ledger["U"] == "open for any nontrivial variable H9/14 class"
    assert ledger["class_wide_surface_softening_theorem"] is False
    assert payload["scope"]["new_polytope_claim"] is False
    assert payload["scope"]["class_wide_IC2_or_IC4_claim"] is False
