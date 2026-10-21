from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_r12_h10_16_h10_45_monodromy_carriers import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "data"
    / "enumerations"
    / "r12_h10_16_h10_45_monodromy_carriers.json"
)
REPORT = ROOT / "results" / "r12_h10_16_h10_45_monodromy_carriers.md"


def _archived() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_certificate_replays_exactly() -> None:
    archived = _archived()
    assert build_certificate() == archived
    assert archived["all_checks_pass"]
    assert archived["certificate_payload_sha256"] == (
        "ce457d3a0b1718e73391eb68ec7f7ca0aadb8f2e1f32c4b3eb12a470c3aef9c0"
    )


def test_scope_is_an_algebraic_local_system_not_a_fiber_automorphism() -> None:
    archived = _archived()
    scope = archived["scope"]
    assert scope["generic_fiber_picard_exact"]
    assert not scope["single_fiber_automorphism_claimed"]
    assert not scope["new_abstract_polytope_claimed"]
    assert "parallel transport" in scope["realization_type"]
    assert "equation-to-monodromy-to-Coxeter" in scope["novelty"]


def test_h10_16_explicit_transverse_component_braid() -> None:
    case = _archived()["H10_16"]
    family = case["family"]
    assert family["equation"] == "y^2=x^3+P(t^2)*x+t*Q(t^2)"
    assert family["generic_NS"] == "U+D8(-1)+A1(-1)^2"
    assert family["generic_Picard_rank"] == 12
    assert family["generic_moduli_dimension"] == 4
    assert family["Mordell_Weil_rank"] == family["Mordell_Weil_torsion"] == 0
    assert "trivial lattice is saturated" in family["primitive_saturation_proof"]

    meridian = case["transverse_meridian"]
    assert meridian["slice"] == "r=1, alpha=2, beta=1, gamma=-2/3+z"
    assert meridian["normal_discriminant"] == (
        "alpha^2-12*beta-12*gamma=-12*z"
    )
    assert meridian["third_blowup_polynomial"] == "3*q**2 + 2*q + z + 1/3"
    assert meridian["root_monodromy"] == (
        "q_+ <-> q_- after one positive meridian"
    )
    assert meridian["central_fiber_enhancement"] == "I4*/D8 -> I6*/D10 at z=0"
    assert case["component_transport"]["exchanged"] == [
        "D8 spinor component S_+",
        "D8 spinor component S_-",
    ]
    assert all(case["checks"].values())


def test_h10_16_full_integral_discriminant_and_vinberg_identification() -> None:
    case = _archived()["H10_16"]
    action = case["integral_action_identification"]
    assert action["basis_determinant"] == -1
    assert action["basis_carries_spinors_to_archived_walls"] == [10, 11]
    assert action["full_integral_action_pairs_conjugate"]
    assert action["fixed_rank_signature_determinant"] == [11, [1, 10], 16]
    assert action["fixed_smith"] == [*[1] * 8, 2, 2, 4]
    assert action["coinvariant"] == "<-4>, divisibility 2"
    finite = action["finite_discriminant_form"]
    assert finite["group"] == "(Z/2)^4"
    assert finite["actions_equivariant"]
    assert finite["all_q_classes_checked"] == 16
    assert len(finite["q_class_table"]) == 16

    vinberg = case["vinberg"]
    assert vinberg == {
        "bad_edges": 0,
        "canonical_full_audit_sha256": (
            "900df6d000ca91669cea117474cda559e4188600fb2744bafa0f69fd1541c226"
        ),
        "classification": "finite_volume_noncompact",
        "cusps": 6,
        "dimension": 10,
        "edge_extension_histogram": {"2": 426},
        "elliptic_edges": 426,
        "facets": 16,
        "finite_vertices": 70,
        "parabolic_neighborhood_audits": 14,
    }


def test_h10_45_explicit_transverse_component_braid() -> None:
    case = _archived()["H10_45"]
    family = case["family"]
    assert family["equation"] == "y^2=x^3+t^2*P(t^2)*x+t^3*Q(t^2)"
    assert family["generic_NS"] == "U+D4(-1)^2+A1(-1)^2"
    assert family["generic_Picard_rank"] == 12
    assert family["generic_moduli_dimension"] == 4
    assert family["Mordell_Weil_rank"] == family["Mordell_Weil_torsion"] == 0
    assert "trivial lattice is saturated" in family["primitive_saturation_proof"]

    meridian = case["transverse_meridian"]
    assert meridian["slice"] == "alpha=beta=gamma=1, delta=2+z"
    assert meridian["zero_Tate_cubic"] == "rho**3 - 3*rho - z + 2"
    assert meridian["normal_discriminant"] == "27*z*(4-z)"
    assert meridian["boundary_factorization"] == "(rho-1)^2*(rho+2)"
    assert meridian["central_fiber_enhancement"] == "I0*/D4 -> IV*/E6 at z=0"
    assert case["component_transport"]["exchanged"] == [
        "first D4 outer component E_+",
        "first D4 outer component E_-",
    ]
    assert all(case["checks"].values())


def test_h10_45_full_integral_discriminant_and_vinberg_identification() -> None:
    case = _archived()["H10_45"]
    action = case["integral_action_identification"]
    assert action["basis_determinant"] == -1
    assert action["basis_carries_outer_pair_to_archived_walls"] == [2, 3]
    assert action["full_integral_action_pairs_conjugate"]
    assert action["fixed_rank_signature_determinant"] == [11, [1, 10], 64]
    assert action["fixed_smith"] == [*[1] * 6, 2, 2, 2, 2, 4]
    assert action["coinvariant"] == "<-4>, divisibility 2"
    finite = action["finite_discriminant_form"]
    assert finite["group"] == "(Z/2)^6"
    assert finite["actions_equivariant"]
    assert finite["all_q_classes_checked"] == 64
    assert len(finite["q_class_table"]) == 64

    vinberg = case["vinberg"]
    assert vinberg["classification"] == "finite_volume_noncompact"
    assert (vinberg["dimension"], vinberg["facets"]) == (10, 45)
    assert (vinberg["finite_vertices"], vinberg["cusps"]) == (1048, 93)
    assert vinberg["elliptic_edges"] == 7332
    assert vinberg["endpoint_type_histogram"] == {
        "0,2": 578,
        "1,1": 3028,
        "2,0": 3726,
    }
    assert vinberg["parabolic_neighborhood_audits"] == 302
    assert vinberg["bad_edges"] == vinberg["bad_parabolic_neighborhoods"] == 0
    assert vinberg["inherited_gate_sha256"] == (
        "99ce1cc8684c5811fc3d684a013629a13820b7dc493c4abea414c4faeef9daea"
    )


def test_report_states_the_positive_result_and_its_exact_scope() -> None:
    report = " ".join(REPORT.read_text(encoding="utf-8").lower().split())
    assert "positive theorem" in report
    assert "two explicit algebraic carriers" in report
    assert "parallel transport" in report
    assert "not an automorphism of a single k3 surface" in report
    assert "finite-volume noncompact" in report
    assert "no new abstract-polytope claim" in report
