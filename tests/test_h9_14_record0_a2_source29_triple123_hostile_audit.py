from scripts.independent_h9_14_record0_a2_source29_triple123_hostile_audit import (
    EXPECTED_AUDIT_DIGEST,
    audit,
)


def test_independent_triple123_fibre_and_budget_reconstruction() -> None:
    result = audit()
    assert result["audit_digest"] == EXPECTED_AUDIT_DIGEST
    assert result["base_link"] == {
        "extension_counts_f_minus_one_through_f2": [1, 9, 36, 84],
        "h0_through_h3": [1, 3, 6, 10],
        "g1_through_g3": [2, 3, 4],
    }
    ground = result["ground_and_fibres"]
    assert ground["candidate_A5_count"] == 1600
    assert ground["allowed_six_fibre_count"] == 84
    assert ground["T_containing_candidate_A5_count"] == 36
    assert ground["killers_per_fibre_histogram"] == {3: 84}
    assert ground["noncontaining_killer_union_count"] == 252
    assert [row["outside_fibre_budget"] for row in result["strict_budget"]] == [4, 2, 1, 0]


def test_independent_burnside_and_residual_stabilizer_audit() -> None:
    symmetry = audit()["symmetry_and_branches"]
    assert symmetry["zero_support_orbits_by_d"] == {"d0": 1, "d1": 2, "d2": 5}
    assert symmetry["active_support_orbits_by_d"] == {"d0": 1, "d1": 5, "d2": 24}
    assert symmetry["fixed_mask_kills_123567"] == {
        "227": False,
        "229": False,
        "230": True,
    }
    assert symmetry["mask230_excluded_disjoint_orbits_per_affected_branch"] == 7
    assert [row["configuration_branch_count"] for row in symmetry["branch_rows"]] == [
        8, 30, 30, 30, 23, 23,
    ]
    assert symmetry["combined_configuration_branch_count"] == 144


def test_independent_near_seed_and_scope_audit() -> None:
    result = audit()
    near = result["near_seed"]
    assert near["T_containing_masks"] == [2094, 302]
    assert near["base_covered_fibre_count"] == 13
    assert [row["mask"] for row in near["new_outside_killers"]] == [342, 1578]
    assert near["actual_missing_six_fibre_count"] == 15
    assert near["final_g3"] == -1
    assert near["negative_g3_defect"] == 1
    assert result["scope"] == {
        "source29_status": "OPEN",
        "enumerated_objects": "pair-configuration orbits, not completed strict-link survivors",
        "configuration_rows_are_survivors": False,
        "outside_killer_supports_still_unresolved": True,
        "compact_incidence_or_hyperbolic_witness_produced": False,
    }
