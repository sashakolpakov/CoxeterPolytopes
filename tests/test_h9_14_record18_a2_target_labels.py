from __future__ import annotations

import json
from pathlib import Path

from scripts.explore_h9_14_record18_a2_target_labels import (
    audit_target,
    sparse_label_gate,
)
from scripts.explore_h9_14_record18_a2_add5_z3 import solve as solve_add5


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/enumerations/h9_14_record18_a2_idx3_additions.json"
CANONICAL_SUPPORTS = (
    (0, 2, 4, 5),
    (0, 3, 6, 7),
    (0, 4, 5, 6, 7),
    (1, 8, 9, 10, 11),
    (2, 4, 5, 8, 9),
    (3, 6, 7, 10, 11),
)


def test_record18_a2_idx3_residual_target_is_a_homology_sphere() -> None:
    document = json.loads(INPUT.read_text(encoding="utf-8"))
    fours = tuple(tuple(row) for row in document["representative_A4s"])
    additions = tuple(
        tuple(row) for row in document["scheduled_target_additions"][0]
    )
    audit = audit_target(fours, additions)
    assert audit["blocker_size_histogram"] == {"5": 106}
    assert audit["global_audit"]["ridge_degree_histogram"] == {"2": 477}
    assert audit["global_audit"]["passes_connected_homology_audit"] is True
    assert tuple(tuple(row) for row in audit["canonical_supports"]) == (
        CANONICAL_SUPPORTS
    )


def test_record18_a2_idx3_residual_orbit_has_no_compact_finite_labels() -> None:
    gate = sparse_label_gate(CANONICAL_SUPPORTS, 250_000)
    assert gate["actual_rank_nine_vertex_count"] == 106
    assert [stage["passing_state_count"] for stage in gate["stages"]] == [
        108,
        84,
        4,
        28,
        196,
        0,
    ]
    assert gate["state_cap_exceeded"] is False
    assert gate["terminal_sparse_state_count"] == 0
    assert gate["unrestricted_compact_label_no_go"] is True


def test_record18_both_a2_add5_orbits_fail_the_exact_ridge_gate() -> None:
    for representative in ((53, 54), (23, 27)):
        audit = solve_add5(representative)
        assert audit["forced_A5_shadow_count"] == 15
        assert audit["additional_A5_count"] == 5
        assert audit["candidate_count"] == 777
        assert audit["base_small_transversal_count"] == 56
        assert audit["base_order_five_transversal_count"] == 233
        assert audit["ridge_constraint_count"] == 630
        assert audit["ridge_and_scheduled_model_count"] == 0
        assert audit["scheduled_pure_target_count"] == 0
        assert audit["solver_terminal_status"] == "unsat"
