from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "data/enumerations/alexeev_118_dual_h11_29_construction.json"
VINBERG = ROOT / "data/enumerations/alexeev_118_dual_h11_29_vinberg.json"
LINEAGE = ROOT / "data/enumerations/alexeev_118_dual_h11_29_lineage.json"
STABLE = ROOT / "data/enumerations/alexeev_118_dual_h11_29_stable_overlattices.json"
REPORT = ROOT / "results/alexeev_118_dual_h11_29_positive_example.md"
LEDGER = ROOT / "results/higher_dimensional_polytope_lineage_ledger.md"


def _read(path: Path) -> dict:
    return json.loads(path.read_text())


def test_complete_two_stock_subgroup_denominator_and_strict_sections() -> None:
    payload = _read(CONSTRUCTION)
    assert payload["all_checks_pass"]
    assert payload["summary"] == {
        "both_strict_integral_coxeter_candidates": True,
        "candidate_dimensions": [11, 11],
        "candidate_facet_counts": [29, 29],
        "proper_high_class_counts": [1, 1],
        "stock_count": 2,
        "subgroup_conjugacy_class_counts": [56, 56],
        "weighted_group_orders": [720, 720],
    }
    assert {case["id"] for case in payload["cases"]} == {
        "A16_6_1_C01",
        "A16_12_1_C01",
    }
    for case in payload["cases"]:
        assert case["decision"] == {
            "classification": "coxeter_candidate_requires_volume_audit",
            "finite_volume_claimed": False,
        }
        assert case["exact_reduction"]["active_hyperbolic_facet_count"] == 29
        assert all(case["checks"].values())


def test_literal_vinberg_finite_volume_criterion_is_frozen() -> None:
    payload = _read(VINBERG)
    assert payload["all_checks_pass"]
    assert payload["summary"]["finite_volume_count"] == 2
    assert payload["summary"]["noncompact_count"] == 2
    assert payload["summary"]["shared_elliptic_edge_count"] == 28798
    assert payload["summary"]["shared_finite_vertex_count"] == 3876
    assert payload["summary"]["shared_cusp_count"] == 253
    for audit in payload["audits"]:
        assert audit["classification"] == "finite_volume_noncompact_Coxeter_polytope"
        assert audit["facet_count"] == 29
        assert audit["bad_edge_count"] == 0
        assert audit["bad_parabolic_neighborhood_count"] == 0
        assert all(audit["checks"].values())


def test_priority_scope_report_and_ledger_use_the_closed_lineage() -> None:
    lineage = _read(LINEAGE)
    stable = _read(STABLE)
    assert lineage["all_checks_pass"]
    assert all(lineage["checks"].values())
    integral = lineage["normalized_integral_lineage"]
    assert all(integral["checks"].values())
    assert integral["integral_isometry_identity"] == "A^T (G_12/2) A = M_6"
    assert integral["containing_full_reflection_chamber"] == {
        "census": "List11",
        "entry_one_based": 5,
        "genus_symbol": "I_{11,1}(2^{2})",
        "facet_count": 14,
        "finite_vertex_count": 56,
        "cusp_count": 5,
        "section_to_full_index_status": (
            "finite; numerical value open because both H11 Coxeter Euler "
            "characteristics vanish"
        ),
    }
    assert "3/4" in integral["A18_4_0_C03_noncollision"][
        "A18_C03_discriminant_q_mod_2_histogram"
    ]
    assert "3/4" not in integral["A18_4_0_C03_noncollision"][
        "A16_6_discriminant_q_mod_2_histogram"
    ]
    assert stable["all_checks_pass"]
    assert stable["summary"]["one_common_normalized_integral_reflection_system"]
    assert stable["summary"]["new_arithmetic_commensurability_class_count"] == 0
    report = REPORT.read_text()
    ledger = LEDGER.read_text()
    assert "newly extracted literal 29-wall chamber" in report
    assert "not a worldwide priority theorem" in report
    assert "A16-DUAL-C2-H11-29" in ledger
    assert "one normalized integral 29-reflection subgroup" in ledger
    for relative in (
        "../data/enumerations/alexeev_118_dual_h11_29_construction.json",
        "../data/enumerations/alexeev_118_dual_h11_29_vinberg.json",
        "../data/enumerations/alexeev_118_dual_h11_29_lineage.json",
        "../data/enumerations/alexeev_118_dual_h11_29_stable_overlattices.json",
    ):
        assert (REPORT.parent / relative).resolve().is_file()
