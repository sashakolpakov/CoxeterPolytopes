from functools import lru_cache
import json

from scripts.audit_h9_14_t1_t2_d34_unfolding_census_hostile import (
    OUTPUT,
    build_certificate,
)


@lru_cache(maxsize=1)
def payload() -> dict:
    return build_certificate()


def cases() -> dict[str, dict]:
    return {row["type_id"]: row for row in payload()["compiler_inventory"]}


def test_exact_archive_and_independent_method() -> None:
    result = payload()
    assert result == json.loads(OUTPUT.read_text())
    assert result["all_checks_pass"] is True
    method = result["method_independence"]
    assert method["production_code_imported"] is False
    assert method["Hermite_normal_form_used"] is False
    assert method["principal_facet_basis_assumed"] is False
    assert result["production_comparison"]["all_agree"] is True


def test_corpus_filter_scalings_and_all_overclosures() -> None:
    result = payload()
    assert result["bounded_primary_corpus"]["marking_count"] == 6
    assert result["bounded_primary_corpus"]["diagram_type_count"] == 4
    assert result["D34_filter"]["retained_marking_count"] == 3
    assert result["D34_filter"]["retained_diagram_type_count"] == 2
    assert sum(row["scaling_count"] for row in result["compiler_inventory"]) == 36
    assert sum(
        lane["candidate_count"]
        for row in result["compiler_inventory"]
        for lane in row["lanes"]
    ) == 28
    assert all(
        lane["smith_quotient"]["radical_saturated"]
        and lane["smith_quotient"]["labelled_gram_recovered"]
        for row in result["compiler_inventory"]
        for lane in row["lanes"]
    )


def test_t1_delicate_index_delta_integrality_and_fold() -> None:
    t1 = cases()["D34-9-8"]
    survivor = t1["unique_survivor"]
    fixed = survivor["fixed_lattice"]
    rank18 = survivor["rank18_unfolding"]
    roots = survivor["roots_and_fold"]
    k3 = survivor["weak_K3_lattice_action_extension"]
    lane = next(lane for lane in t1["lanes"] if lane["label"] == "short4_long10")
    candidate = next(row for row in lane["candidates"] if row["survives"])
    assert candidate["index_over_root_span"] == 2
    assert candidate["discriminant_delta"] == 0
    assert candidate["q_histogram"] == {"0": 136, "1": 120}
    assert fixed["two_elementary_invariants"] == [10, 8, 0]
    assert fixed["isometry_class"] == "U+E8(-2)"
    assert rank18["signature_determinant"] == [1, 17, -4]
    assert rank18["two_elementary_invariants"] == [18, 2, 0]
    assert rank18["q_histogram"] == {"0": 3, "1": 1}
    assert rank18["smith"] == [1] * 16 + [2, 2]
    assert roots["root_span_smith"] == [1] * 17 + [2]
    assert roots["roots_span_index_in_rank18_lattice"] == 2
    assert roots["off_diagonal_pairing_histogram"] == {"0": 242, "1": 34}
    assert roots["exact_labelled_fold_back"] is True
    assert all(survivor["checks"].values())
    assert k3["complement"] == "U+U(2)"
    assert k3["graph_glue_index"] == 4
    assert all(k3["checks"].values())


def test_t2_replay_and_primitive_k3_glue() -> None:
    survivor = cases()["D34-8-8"]["unique_survivor"]
    rank18 = survivor["rank18_unfolding"]
    roots = survivor["roots_and_fold"]
    k3 = survivor["weak_K3_lattice_action_extension"]
    assert survivor["fixed_lattice"]["two_elementary_invariants"] == [10, 8, 0]
    assert rank18["signature_determinant"] == [1, 17, -16]
    assert rank18["two_elementary_invariants"] == [18, 4, 0]
    assert rank18["q_histogram"] == {"0": 10, "1": 6}
    assert roots["root_span_smith"] == [1] * 18
    assert roots["roots_span_index_in_rank18_lattice"] == 1
    assert roots["off_diagonal_pairing_histogram"] == {"0": 244, "1": 32}
    assert k3["complement"] == "U(2)+U(2)"
    assert k3["graph_glue_index"] == 16
    assert all(k3["checks"].values())


def test_scope_stays_corpus_relative_and_U_open() -> None:
    boundary = payload()["scope_boundary"]
    assert boundary["U"].startswith("open:")
    assert "T1 cone completeness" in boundary["C"]
    assert "not re-audited" in boundary["V"]
    assert boundary["new_polytope"] is False
