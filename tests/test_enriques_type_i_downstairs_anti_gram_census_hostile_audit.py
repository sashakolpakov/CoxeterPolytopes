from __future__ import annotations

from functools import lru_cache
import hashlib
import json

from scripts.audit_enriques_type_i_downstairs_anti_gram_census_hostile import (
    OUTPUT,
    analyze,
)


ARCHIVE_SHA256 = "9843f1cfb999fb513f65df08d1a125bd988817c9376156aab50cc085b6217f6f"


@lru_cache(maxsize=1)
def _payload() -> dict:
    return analyze()


def test_hostile_archive_replays_and_freezes_strengthened_theorem() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert hashlib.sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    assert payload["verdict"] == "PASS/FREEZE (strengthened: rank ten eliminated)"
    assert payload["all_checks_pass"] is True
    assert "unique rank-nine geometric orbit" in payload["final_theorem"]


def test_independent_switching_gauge_and_sturm_inertia_match() -> None:
    quotient = _payload()["independent_switching_census"]
    assert quotient["different_from_production_gauge"] is True
    assert quotient["raw_assignment_count"] == 27648
    assert quotient["orbit_count"] == 18
    assert quotient["orbit_size_histogram"] == {"512": 2, "1024": 6, "2048": 10}
    inertia = _payload()["independent_exact_inertia"]
    assert inertia["engine"] == (
        "factorized exact characteristic polynomial plus Sturm root counts"
    )
    assert inertia["raw_assignment_histogram"] == {
        "p0_n10_z2": 512,
        "p0_n9_z3": 1024,
        "p1_n10_z1": 9216,
        "p1_n11_z0": 512,
        "p1_n9_z2": 2048,
        "p2_n10_z0": 14336,
    }
    assert inertia["rank_9_equals_geometric_package"] is True


def test_all_sheet_choices_force_rank_ten_saturation_root() -> None:
    killer = _payload()["rank_ten_hostile_killer"]
    assert killer["generated_lattice_before_saturation"] == "E8(-2)+<-4>+<-4>"
    assert killer["generated_lattice_determinant_absolute"] == 4096
    assert killer["generated_lattice_rootless_by_norm_divisibility"] is True
    assert killer["orthogonal_wall_pair_one_based"] == [11, 12]
    assert killer["downstairs_root_sum"] == [8, 4, 8, 14, 12, 10, 8, 6, 4, 2]
    assert killer["square_for_all_sheet_signs"] == {
        "+1,+1": -2,
        "+1,-1": -2,
        "-1,+1": -2,
        "-1,-1": -2,
    }
    assert killer["e_lies_in_saturation_of_anti_span"] is True
    assert killer["e_is_anti_invariant"] is True
    assert killer["L2_rootlessness_fails"] is True


def test_rank_ten_is_eliminated_without_overclaiming_remaining_compiler() -> None:
    scope = _payload()["remaining_scope"]
    assert scope == {
        "new_polytope_or_surface_claimed": False,
        "rank_ten_branch_survives": False,
        "still_not_a_downstairs_construction_of_the_rank_nine_embedding_glue_period_L4_package": True,
    }
