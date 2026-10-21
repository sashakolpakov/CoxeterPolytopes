from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a15_11_h10_21_x2_k3_carrier import OUTPUT, build_certificate


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text())
    assert rebuilt == archived
    return rebuilt


def test_exact_index_four_bridge_to_actual_x2_fixed_lattice(certificate: dict) -> None:
    bridge = certificate["lattice_bridge"]
    assert bridge["M"]["determinant"] == 8192
    assert bridge["N"]["determinant"] == 512
    assert bridge["N"]["isometry_class"] == "U+E8(-2)+<-2>"
    assert all(bridge["checks"].values())
    assert bridge["checks"]["P_congruence_to_exact_standard"] is True
    assert bridge["checks"]["Q_actual_X2_fixed_to_standard"] is True
    assert bridge["checks"]["explicit_M_index_four_in_actual_F"] is True
    assert bridge["checks"]["actual_fixed_gram_inside_full_NS"] is True


def test_all_21_on_F_but_only_11_on_full_x2_ns(certificate: dict) -> None:
    walls = certificate["walls"]
    assert walls["count"] == 21
    assert walls["primitive_F_type_histogram"] == {
        "-4,2": 10,
        "-2,1": 4,
        "-2,2": 7,
    }
    assert walls["full_NS_extension_histogram"] == {
        "-4,2,1,False": 10,
        "-2,1,1,True": 4,
        "-2,2,1,True": 7,
    }
    assert walls["all_21_integral_on_F"] is True
    assert all(all(row["checks"].values()) for row in walls["records"])
    assert sum(row["reflection_integral_on_full_NS"] for row in walls["records"]) == 11


def test_generic_full_ns_carrier_and_named_x2_scope(certificate: dict) -> None:
    generic = certificate["generic_F_polarized_K3_carrier"]
    named = certificate["named_X2_occurrence"]
    scope = certificate["scope"]
    assert generic["full_NS"] == "NS(Y)=F=U+E8(-2)+<-2>"
    assert generic["all_21_reflections_integral_on_full_NS"] is True
    assert generic["effective_minus_two_facets"] == 11
    assert generic["transported_strict_point_square"] == "159/2"
    assert generic["minimum_oriented_wall_pairing"] == "1"
    assert generic["oriented_wall_pairings"] == [
        "1", "1", "11", "1", "1", "1", "1", "1", "2", "5", "1",
        "1", "11", "14", "7", "7", "6", "11", "13", "15", "15",
    ]
    assert named["full_NS_rank"] == 20
    assert named["fixed_lattice_rank"] == 11
    assert named["full_NS_reflection_integral_count"] == 11
    assert named["full_NS_reflection_nonintegral_count"] == 10
    assert "not Shimada's 37-wall" in scope["not_fixed_nef_section"]
    assert "not asserted" in scope["generic_K3_Y"]


def test_two_exact_finite_volume_engines(certificate: dict) -> None:
    volume = certificate["finite_volume"]
    vinberg = volume["literal_vinberg"]
    cdd = volume["independent_GMP_cdd_double_description"]
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["finite_vertices"], vinberg["ideal_vertices"]) == (607, 50)
    assert vinberg["elliptic_edges"] == 4159
    assert all(vinberg["checks"].values())
    assert cdd["extreme_ray_count"] == 657
    assert cdd["causal_type_histogram"] == {"finite": 607, "ideal": 50}
    assert cdd["hyperideal_ray_count"] == 0
    assert cdd["ideal_facet_incidence_histogram"] == {"11": 18, "12": 32}


def test_reusable_criterion_and_bounded_novelty(certificate: dict) -> None:
    criterion = certificate["reusable_carrier_criterion"]
    lineage = certificate["arithmetic_lineage"]
    assert "very general N-polarized projective K3" in criterion["statement"]
    assert "separate" in criterion["full_nef_gate"]
    assert lineage["index_in_full_reflection_group"] == 6720
    assert lineage["priority_scope"] == "bounded evidence only; no worldwide-priority claim"
    assert all(row["upstream_all_checks_pass"] for row in certificate["sources"].values())
    assert certificate["all_checks_pass"] is True


def test_report_freezes_two_tiers_and_nonextension() -> None:
    report = Path(__file__).resolve().parents[1] / "results/a15_11_h10_21_x2_k3_carrier.md"
    text = " ".join(report.read_text().lower().split())
    assert "toot" in text
    assert "very general" in text
    assert "607 finite" in text and "50 ideal" in text
    assert "ten square-`-4`" in text
    assert "do not extend" in text
    assert "full nef cone" in text
