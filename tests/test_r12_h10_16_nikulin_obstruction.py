from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_r12_h10_16_nikulin_obstruction import build_certificate


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "r12_h10_16_nikulin_obstruction.json"
)
REPORT = ROOT / "results" / "r12_h10_16_nikulin_obstruction.md"


def _archived() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_certificate_replays_exactly() -> None:
    archived = _archived()
    assert build_certificate() == archived
    assert archived["all_checks_pass"]
    assert archived["case_id"] == "R12-D8A1x2-C2b-H10-16"
    assert archived["certificate_payload_sha256"] == (
        "31255fd53197baca745e0ec5b831929628869efa59222f3f8b6c1f883e5e9960"
    )


def test_archived_polytope_remains_positive_and_finite_volume() -> None:
    polytope = _archived()["archived_positive_polytope"]
    assert polytope == {
        "bad_edge_count": 0,
        "canonical_full_audit_sha256": (
            "900df6d000ca91669cea117474cda559e4188600fb2744bafa0f69fd1541c226"
        ),
        "compact": False,
        "dimension": 10,
        "edge_extension_histogram": {"2": 426},
        "elliptic_edge_count": 426,
        "endpoint_type_histogram": {"0,2": 10, "1,1": 132, "2,0": 284},
        "facet_count": 16,
        "finite_vertex_count": 70,
        "finite_volume": True,
        "ideal_vertex_count": 6,
        "parabolic_neighborhood_audit_count": 14,
        "simple_roots_complete": True,
    }


def test_discriminant_gate_passes_before_primitive_embedding_fails() -> None:
    archived = _archived()
    finite = archived["finite_discriminant_gate"]
    assert finite["status"] == "passes"
    assert finite["actions_equivariant"]
    assert finite["all_16_classes_checked"]
    assert len(finite["quadratic_form_table"]) == 16
    assert finite["anti_isometry_source_to_target_mod_two"] == [
        [1, 0, 1, 1],
        [0, 1, 1, 1],
        [0, 0, 1, 1],
        [0, 0, 1, 0],
    ]

    obstruction = archived["primitive_E7_minus_2_obstruction"]
    assert obstruction["polar_rank"] == 6
    assert obstruction["polar_radical_dimension"] == 4
    assert obstruction["qbar_is_nonzero_on_polar_radical"]
    assert obstruction["singular_radical_dimension"] == 3
    assert obstruction["maximal_totally_singular_dimension"] == 6
    assert obstruction["bound_is_sharp"]
    assert len(obstruction["sharp_six_dimensional_basis"]) == 6
    assert obstruction["E7_minus_2_mod_two_dimension"] == 7
    assert obstruction["E7_minus_2_qbar_is_identically_zero"]
    assert not obstruction["primitive_embedding_exists"]


def test_k3_scope_is_exact_and_higher_order_lane_stays_open() -> None:
    archived = _archived()
    action = archived["picard_action"]
    assert action["coinvariant_lattice"] == "<-4>"
    assert action["coinvariant_divisibility"] == 2
    assert action["fixed_lattice_smith_diagonal_absolute"] == [
        *[1] * 8,
        2,
        2,
        4,
    ]
    assert not action["fixed_lattice_is_two_elementary"]

    consequence = archived["k3_consequence"]
    assert not consequence["symplectic_Nikulin_lift_exists"]
    assert not consequence["nonsymplectic_involution_with_exact_fixed_NS_exists"]
    assert consequence["minimum_remaining_automorphism_order"] == 4
    assert "order-four or higher" in consequence["remaining_positive_lane"]
    assert "no global no-go" in archived["novelty_scope"]


def test_report_distinguishes_polytope_from_carrier_obstruction() -> None:
    report = " ".join(REPORT.read_text(encoding="utf-8").lower().split())
    assert "h10/16 polytope survives" in report
    assert "finite-volume noncompact" in report
    assert "cannot be the restriction of a symplectic nikulin" in report
    assert "order at least four" in report
    assert "invalidate the hyperbolic polytope" in report
    assert "not a global no-go for higher-order automorphisms" in report
