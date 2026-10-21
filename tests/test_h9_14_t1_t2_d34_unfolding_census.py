from functools import lru_cache
import json

from scripts.certify_h9_14_t1_t2_d34_unfolding_census import (
    OUTPUT,
    build_certificate,
)


@lru_cache(maxsize=1)
def payload() -> dict:
    return build_certificate()


def test_exact_archive_and_bounded_scope() -> None:
    result = payload()
    assert result == json.loads(OUTPUT.read_text())
    assert result["all_checks_pass"] is True
    assert result["bounded_primary_corpus"]["marking_count"] == 6
    assert result["bounded_primary_corpus"]["diagram_type_count"] == 4
    assert result["surface_free_D34_specification"]["retained_marking_count"] == 3
    assert result["surface_free_D34_specification"]["retained_diagram_type_count"] == 2
    assert result["new_polytope_claim"] is False


def test_complete_scaling_and_overlattice_inventory() -> None:
    cases = {row["type_id"]: row for row in payload()["compiler_inventory"]}
    assert set(cases) == {"D34-8-8", "D34-9-8"}
    assert sum(row["scaling_count"] for row in cases.values()) == 36
    assert sum(
        lane["candidate_count"]
        for row in cases.values()
        for lane in row["lanes"]
    ) == 28
    assert cases["D34-8-8"]["diagram_automorphism_order"] == 48
    assert cases["D34-9-8"]["diagram_automorphism_order"] == 8
    assert [lane["candidate_count"] for lane in cases["D34-8-8"]["lanes"]] == [16, 5]
    assert [lane["candidate_count"] for lane in cases["D34-9-8"]["lanes"]] == [2, 5]


def test_t1_new_unfolding_and_t2_replay() -> None:
    cases = {row["type_id"]: row for row in payload()["compiler_inventory"]}
    t2 = cases["D34-8-8"]["unique_survivor"]
    t1 = cases["D34-9-8"]["unique_survivor"]
    assert t1["fixed_lattice"]["isometry_class"] == "U+E8(-2)"
    assert t2["fixed_lattice"]["isometry_class"] == "U+E8(-2)"
    assert t1["rank18_unfolding"]["two_elementary_invariants"] == [18, 2, 0]
    assert t1["rank18_unfolding"]["signature_determinant"] == [1, 17, -4]
    assert t1["weak_K3_lattice_action_extension"]["complement"] == "U+U(2)"
    assert t1["roots_and_fold"]["roots_span_index_in_rank18_lattice"] == 2
    assert t1["roots_and_fold"]["off_diagonal_pairing_histogram"] == {"0": 242, "1": 34}
    assert t2["rank18_unfolding"]["two_elementary_invariants"] == [18, 4, 0]
    assert t2["weak_K3_lattice_action_extension"]["complement"] == "U(2)+U(2)"
    assert t2["roots_and_fold"]["roots_span_index_in_rank18_lattice"] == 1


def test_urcv_boundary_and_casewise_v() -> None:
    result = payload()
    assert result["URCV_boundary"]["U"].startswith("open:")
    assert "T1 remains open" in result["URCV_boundary"]["C"]
    for row in result["compiler_inventory"]:
        volume = row["casewise_V_import"]
        assert volume["classification"] == "finite_volume_noncompact"
        assert volume["bad_edge_count"] == 0
        assert volume["bad_parabolic_neighborhood_count"] == 0
