from __future__ import annotations

import json

from scripts.audit_alexeev_parked_rank10_h9_examples_lineage import OUTPUT, audit


def test_parked_rank10_h9_hostile_lineage_rebuilds() -> None:
    payload = audit()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert payload["decision"] == {
        "geometry": "five_exact_finite_volume_noncompact_H9_polytopes_replay",
        "literal_diagram_lineage": "one_known_ImHof_and_four_bounded_unmatched_worldwide_priority_open",
        "arithmetic_lineage": "zero_new_arithmetic_classes",
    }
    assert payload["coxiter_bounded_comparison"]["unique_exact_collision"][
        "candidate_to_ImHof_permutation_one_based"
    ] == [6, 11, 10, 9, 8, 1, 2, 3, 4, 5, 7, 12]
    assert payload["coxiter_bounded_comparison"]["other_four_exact_match_count"] == 0
    assert payload["prior_project_H9_comparison"]["total_exact_match_count"] == 0
    assert all(
        case["all_exact_geometry_checks_pass"] is True
        and case["hyperideal_vertex_count"] == 0
        and case["reflection_stable_overlattice_lineage"]["new_arithmetic_class"] is False
        for case in payload["cases"]
    )
