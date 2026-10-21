from __future__ import annotations

from collections import Counter
import json

from scripts.inventory_shimada_disc16_atlas_fixed_sections import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_archived_stabilizer_denominator_counts() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    summary = payload["summary"]
    assert summary["chamber_count"] == 4_539
    assert summary["nontrivial_stabilizer_chamber_count"] == 98
    assert summary["total_subgroup_count_in_nontrivial_chambers"] == 346
    assert summary["total_subgroup_conjugacy_class_count"] == 257
    assert summary["identity_ambient_class_count"] == 98
    assert summary["proper_subgroup_conjugacy_class_count"] == 159
    assert summary["proper_high_dimensional_subgroup_class_count"] == 107
    assert summary["proper_high_dimensional_class_histogram"] == {
        "H10": 4,
        "H11": 101,
        "H14": 2,
    }
    assert summary["chambers_with_a_proper_class_above_h9"] == 94


def test_identity_classes_are_not_misreported_as_proper_sections() -> None:
    payload = _payload()
    chambers = payload["chambers"]
    assert len(chambers) == 98
    assert len({chamber["chamber_index"] for chamber in chambers}) == 98

    classes = [row for chamber in chambers for row in chamber["classes"]]
    identities = [row for row in classes if row["order"] == 1]
    proper_high = [
        row for row in classes if row["proper_fixed_section_above_h9"]
    ]
    assert len(identities) == 98
    assert all(row["fixed_rank"] == 20 for row in identities)
    assert all(row["proper_subgroup"] is False for row in identities)
    assert all(
        row["order"] > 1 and row["fixed_rank"] < 20 for row in proper_high
    )
    assert len(proper_high) == 107
    assert Counter(row["hyperbolic_dimension"] for row in proper_high) == {
        10: 4,
        11: 101,
        14: 2,
    }


def test_scope_does_not_conflate_archived_and_full_graph_stabilizers() -> None:
    payload = _payload()
    scope = payload["scope"]
    assert scope["all_published_chambers"] is True
    assert scope["all_nontrivial_archived_Aut_D_stabilizers"] is True
    assert scope["all_subgroups_of_each_archived_Aut_D"] is True
    assert scope["full_weighted_wall_graph_automorphism_groups_computed"] is False
    assert scope["full_weighted_wall_graph_symmetry_inventory_is_separate"] is True
    assert scope["does_not_claim_archived_Aut_D_is_full_weighted_wall_symmetry"] is True
    assert scope["wall_angle_gate_performed"] is False
    assert scope["volume_gate_performed"] is False

    decision = payload["decision"]
    assert decision["archived_stabilizer_denominator_complete"] is True
    assert decision["proper_fixed_section_classes_above_h9"] == 107
    assert decision["requires_wall_angle_reduction_and_volume_audit"] is True
    assert decision["full_weighted_wall_graph_denominator_complete"] is False

