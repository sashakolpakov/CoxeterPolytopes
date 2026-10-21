import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LINEAGE = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_lineage.json"
)


def _load() -> dict:
    return json.loads(LINEAGE.read_text(encoding="utf-8"))


def test_production_and_primary_sources_are_pinned() -> None:
    record = _load()
    assert record["all_checks_pass"] is True
    assert record["production_certificates"]["construction"]["sha256"] == (
        "a9087200b360c4bed75cd93c987bd7e7d493ff1ce852bb9c5ea4dc3490f1fe11"
    )
    assert record["production_certificates"]["vinberg"]["sha256"] == (
        "032c8343a6dd7f4904b9daecb074b8b0bbbf1a94f9431252884484c98d2e8456"
    )
    assert record["primary_sources"]["alexeev"]["all_cases_sha256"] == (
        "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
    )
    assert record["primary_sources"]["scharlau_kirschmer"][
        "list13_entry_one_based"
    ] == 10
    assert record["primary_sources"]["scharlau_kirschmer"][
        "list14_entry_one_based"
    ] == 7


def test_h13_is_a_finite_index_chamber_on_known_14_10_0_lattice() -> None:
    record = _load()["H13_50"]
    assert record["literal_polytope"]["finite_volume_noncompact"] is True
    assert record["fixed_lattice"]["nikulin_invariants"] == [14, 10, 0]
    assert record["fixed_lattice"]["doubled_dual_determinant"] == -16
    assert record["fixed_lattice"]["doubled_dual_discriminant_q_histogram"] == {
        "0": 6,
        "1": 10,
    }
    full = record["known_full_reflection_chamber"]
    assert full["parent_case"] == "(15,11,1)"
    assert full["facet_count"] == 16
    assert full["root_norm_histogram"] == {"-2": 4, "-4": 12}
    assert full["coxeter_edge_histogram"] == {"2": 104, "3": 14, "4": 2}
    assert full["finite_special_subset_counts_by_rank"] == [
        16, 120, 560, 1820, 4364, 7965, 11232,
        12272, 10300, 6488, 2920, 842, 111,
    ]
    assert full["maximal_parabolic_neighborhood_count"] == 6
    assert record["group_relation"]["index_value"] is None
    assert record["priority"]["worldwide_literal_chamber_priority"] == "open"
    assert all(record["checks"].values())


def test_h14_is_exact_alexeev_doubled_dual_rediscovery() -> None:
    record = _load()["H14_98"]
    assert record["literal_polytope"]["finite_volume_noncompact"] is True
    assert record["literal_polytope"]["finite_special_subset_counts_by_rank"] == [
        98, 1666, 14616, 83559, 342888, 1056573, 2490674,
        4501714, 6177500, 6299076, 4598006, 2249968, 646786, 75898,
    ]
    assert record["literal_polytope"]["ideal_vertex_count"] == 2114
    closure = record["reflection_stable_overlattice"]
    assert closure["integral_nonzero_order_two_overlattice_class_count"] == 127
    assert closure["W98_stable_class_count"] == 1
    assert closure["unique_stable_class_vector_mod_2"] == [0, 1, 1] + [0] * 12
    assert closure["overlattice_index"] == 2
    assert closure["determinant"] == 64
    assert closure["smith_diagonal_absolute"] == [1] * 9 + [2] * 6
    match = record["exact_Alexeev_match"]
    assert match["source_case"] == "(15,9,1)"
    assert len(set(match["candidate_to_alexeev_permutation_zero_based"])) == 98
    assert match["direct_finite_special_counts_ranks_1_2_3"] == [98, 1666, 14616]
    assert match["transformed_facet_gram_sha256"] == (
        "257c0444e7774232672520096ef5639074d83ef2e1a58d625ab2435007f2f330"
    )
    assert all(match["checks"].values())
    assert record["priority"]["new_hyperbolic_polytope_claim"] is False
    assert all(record["checks"].values())


def test_bounded_priority_boundary_is_not_worldwide_novelty() -> None:
    record = _load()
    primary = record["primary_sources"][
        "bounded_H13_literal_primary_source_comparison"
    ]
    assert primary["literal_H13_50_found"] is False
    assert len(primary["alexeev_explicit_matrix_inventory"]) == 18
    assert [
        row
        for row in primary["alexeev_explicit_matrix_inventory"]
        if row["facet_count"] == 50
    ] == [{
        "facet_count": 50,
        "hyperbolic_dimension": 19,
        "lattice_invariants_r_a_delta": [20, 2, 1],
    }]
    assert len(primary["scharlau_kirschmer_List13_inventory"]) == 17
    assert not any(
        row["facet_count"] == 50
        for row in primary["scharlau_kirschmer_List13_inventory"]
    )
    assert "not an audit of all papers" in primary["warning"]
    assert all(primary["checks"].values())
    bounded = record["bounded_CoxIter_comparison"]
    assert bounded["targets"]["H13_facets_50"][
        "absent_from_both_by_dimension_and_facet_count"
    ] is True
    assert bounded["targets"]["H14_facets_98"][
        "absent_from_both_by_dimension_and_facet_count"
    ] is True
    assert "bounded negative check" in bounded["warning"]
    assert record["decision"]["new_abstract_hyperbolic_examples_certified_by_this_lineage_audit"] == 0
    assert record["decision"]["literal_priority_candidates_still_alive"] == [
        "SS2-2A-H13-50"
    ]
