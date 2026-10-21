import json
from pathlib import Path

from scripts.certify_borcherds_char2_v4_h10_34 import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "borcherds_char2_v4_h10_34_certificate.json"
)


def _payload() -> dict:
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    return payload


def test_unique_high_noncyclic_subgroup_class_and_fixed_lattice() -> None:
    payload = _payload()
    assert payload["high_subgroup_filter"]["noncyclic_conclusion"] == (
        "only V4 with fusion 2a + 2c + 2c can have fixed rank at least 11"
    )
    proof = payload["subgroup_class_proof"]
    assert proof["compatible_2c_elements_in_centralizer"] == 24
    assert proof["distinct_V4_subgroups_containing_fixed_2a"] == 12
    assert proof["centralizer_orbits_on_those_V4_subgroups"] == 1
    assert proof["V4_full_group_conjugacy_class_size"] == 3780
    assert proof["geometric_membership_2a_2c_product"] == [True, False, False]
    assert payload["fixed_lattice"]["determinant"] == 1024
    assert payload["fixed_lattice"]["smith_diagonal_absolute"] == (
        [1] * 4 + [2] * 4 + [4] * 3
    )


def test_v4_h10_34_literal_vinberg_certificate() -> None:
    payload = _payload()
    chamber = payload["restricted_chamber"]
    assert chamber["negative_facet_count"] == 34
    assert chamber["coxeter_edge_counts"] == {
        "2": 331,
        "3": 17,
        "4": 30,
        "6": 0,
        "infinity": 75,
        "ultraparallel": 108,
    }
    vinberg = payload["vinberg"]
    assert vinberg["finite_special_subset_counts_by_rank"] == [
        34, 378, 2280, 8725, 22392, 38931, 44922, 32466, 12792, 1774,
    ]
    assert vinberg["ideal_vertex_count"] == 146
    assert vinberg["edge_endpoint_type_histogram"] == {
        "0,2": 1098,
        "1,1": 5648,
        "2,0": 6046,
    }
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0


def test_exact_normalized_gram_deduplication() -> None:
    payload = _payload()
    dedup = payload["exact_deduplication"]
    assert dedup["target_id"] == "A15-C2-H10-34"
    assert dedup["all_1156_normalized_squared_entries_match"] is True
    assert dedup["same_abstract_hyperbolic_Coxeter_polytope"] is True
    assert len(dedup["source_to_target_facet_permutation_one_based"]) == 34
    assert dedup["fixed_lattices_are_not_integrally_isometric"] == {
        "borcherds_V4_determinant": 1024,
        "borcherds_V4_smith": [1] * 4 + [2] * 4 + [4] * 3,
        "A15_C2_determinant": 2048,
        "A15_C2_smith": [1] * 2 + [2] * 7 + [4] * 2,
    }
    assert payload["decision"]["abstract_novelty"].startswith("not new")


def test_full_v4_certificate_rebuilds_exactly() -> None:
    assert build_certificate() == _payload()
