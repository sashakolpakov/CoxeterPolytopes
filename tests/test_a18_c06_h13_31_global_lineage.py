from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/enumerations/a18_c06_h13_31_global_lineage.json"
REPORT = ROOT / "results/a18_c06_h13_31_global_lineage.md"
LEDGER = ROOT / "results/higher_dimensional_polytope_lineage_ledger.md"


def _payload() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_literal_h13_31_polytope_and_integral_fingerprint_are_frozen() -> None:
    payload = _payload()
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())
    target = payload["target"]
    assert target["id"] == "A18-C2c-H13-31"
    assert (target["hyperbolic_dimension"], target["facet_count"]) == (13, 31)
    assert target["volume_type"] == "finite-volume noncompact"
    assert (target["finite_vertex_count"], target["ideal_vertex_count"]) == (
        3433,
        238,
    )
    assert target["face_vector_f0_through_f13"] == [
        3671,
        28766,
        96280,
        194317,
        269179,
        270518,
        201943,
        112606,
        46543,
        13964,
        2925,
        399,
        31,
        1,
    ]
    assert target["coxeter_relation_histogram"] == {
        "0": 66,
        "2": 327,
        "3": 44,
        "4": 28,
    }
    assert target["coxeter_graph_automorphism_group_order"] == 48
    maximal = target["maximal_reflection_stable_overclosure"]
    assert maximal["index_over_fixed_lattice"] == 4
    assert maximal["determinant"] == -4
    assert maximal["primitive_target_mirror_type_histogram"] == {
        "norm_-2_div_1": 26,
        "norm_-4_div_2": 5,
    }


def test_Tumarkin_denominator_and_subgroup_word_handoff_are_exact() -> None:
    payload = _payload()
    denominator = payload["maximal_denominator"]
    assert all(denominator["checks"].values())
    assert denominator["primary_diagram"].endswith("13-tum04_13_01")
    assert denominator["primitive_radical_one_based_node_coefficients"] == [
        2,
        4,
        6,
        5,
        4,
        3,
        2,
        1,
        0,
        -1,
        -2,
        -3,
        -2,
        -1,
        3,
    ]
    assert denominator["scharlau_kirschmer_denominator"] == {
        "list": "List13",
        "entry_one_based": 4,
        "genus_symbol": "II_{13,1}(2^{-2})",
        "facets": 15,
        "finite_vertices": 44,
        "cusps": 2,
    }
    assert denominator["vinberg_audit"]["bad_edge_count"] == 0
    assert denominator["vinberg_audit"]["bad_parabolic_neighborhood_count"] == 0
    subgroup = payload["subgroup_presentation_handoff"]
    assert subgroup["exact_word_checks_pass"]
    assert len(subgroup["target_reflection_conjugates"]) == 31
    assert subgroup["coset_enumeration_experiment"]["result"].startswith(
        "OVERFLOW"
    )
    assert "neither an index lower bound" in subgroup[
        "coset_enumeration_experiment"
    ]["logical_inference"]


def test_Allcock_arbitrary_subtree_and_named_pure_double_lanes_are_closed() -> None:
    payload = _payload()
    audit = payload["Allcock_H13_subtree_gluing"]
    assert all(audit["checks"].values())
    assert audit["Q"]["Allcock_named_selected_gluing_costs"] == [17, 17, 18]
    assert [
        row["wall_zero_based"] for row in audit["Q"]["31_facet_tree_candidates"]
    ] == [1, 2]
    assert all(
        row["relation_histogram"]
        == {"0": 110, "2": 321, "3": 20, "4": 14}
        and not row["target_isomorphic"]
        for row in audit["Q"]["31_facet_tree_candidates"]
    )
    assert audit["F"]["31_facet_tree_candidates"] == []
    assert audit["F"]["two_colour_facet_formula"].endswith("never 31")
    pure = payload["named_pure_doubling_closures"]
    assert all(pure["checks"].values())
    assert all(row["target_absent"] for row in pure["records"])
    allcock = next(row for row in pure["records"] if row["seed"] == "Allcock_D6D6_F")
    assert allcock["isomorphism_classes_by_facet_count"] == {
        "21": 1, "24": 1, "25": 1, "27": 1, "30": 2, "31": 1
    }


def test_report_and_ledger_preserve_priority_and_algebraic_scope() -> None:
    report = REPORT.read_text(encoding="utf-8")
    ledger = LEDGER.read_text(encoding="utf-8")
    assert "strongest worldwide-priority-open" in report
    assert "worldwide priority theorem" in report.lower()
    assert "**not**" in report
    assert "order-four K3 lift theorem" in report
    assert "mixed-root reflection chamber" in report
    assert "A18-C2c-H13-31" in ledger
    assert "Allcock" in ledger
    assert "13-tum04_13_01" in ledger
