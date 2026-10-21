from __future__ import annotations

import json

from scripts.explore_h9_14_record30_a2_orbit10_two_cross import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify,
)
from scripts.independent_h9_14_record30_a2_orbit10_two_cross import replay


def test_record30_a2_orbit10_complete_two_cross_no_go() -> None:
    frozen = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    checked = verify(frozen)
    assert checked == {
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
        "literal_labelled_trial_count": 3_265_920,
        "local_survivor_count": 36_288,
        "rank_at_least_eleven_count": 36_288,
        "terminal_sparse_state_count": 648,
        "unresolved_local_survivor_count": 0,
    }

    hostile = replay()
    assert hostile["agrees_with_frozen_production"] is True
    assert hostile["implementation"]["imports_production_enumerator"] is False
    assert hostile["implementation"]["imports_coxeter_label_module"] is False
    assert hostile["terminal_states_sha256"] == frozen["sparse_gate"]["terminal_states_sha256"]
    assert hostile["local_survivor_sha256"] == frozen["two_cross_shell"]["local_survivor_sha256"]
    assert hostile["zero_fixed_determinant_residue_count"] == 0
    assert hostile["all_local_survivors_rank_at_least_eleven"] is True


def test_record30_a2_orbit10_two_cross_scope_is_not_overstated() -> None:
    frozen = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    scope = frozen["scope"]
    assert scope["exactly_two_cross_component_nonright_meeting_edges"] is True
    assert scope["three_or_more_cross_edges_excluded"] is False
    assert scope["whole_orbit10_nerve_excluded"] is False
    assert scope["whole_record30_a2_schedule_excluded"] is False
    assert scope["finite_volume_noncompact_statement"] is False
    assert frozen["conclusion"]["new_hyperbolic_polytope"] is False
