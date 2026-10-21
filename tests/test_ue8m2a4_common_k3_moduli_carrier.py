from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
import sympy as sp

from scripts.certify_ue8m2a4_common_k3_moduli_carrier import (
    OUTPUT,
    SOURCES,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_full_discriminant_form_gives_exact_alphabet(certificate: dict) -> None:
    lattice = certificate["lattice"]
    assert lattice["name"] == "S=U+E8(-2)+<-4>"
    assert lattice["signature"] == [1, 10, 0]
    assert lattice["determinant"] == 1024
    assert lattice["smith"] == [1, 1] + [2] * 8 + [4]
    assert lattice["discriminant_form"] == "q_plus(F2^8) + <-1/4>"
    assert lattice["primitive_negative_integral_reflection_alphabet"] == [
        {"square": -2, "divisibility": 1},
        {"square": -4, "divisibility": 2},
        {"square": -4, "divisibility": 4},
    ]

    audit = lattice["full_discriminant_form_audit"]
    assert audit["representative_count"] == 1024
    assert audit["distinct_classes_mod_S"] == 1024
    assert audit["all_representatives_lie_in_S_dual"] is True
    assert audit["exact_order_q_mod_2_histogram"] == {
        "order_1_q_0": 1,
        "order_2_q_0": 255,
        "order_2_q_1": 256,
        "order_4_q_3/4": 240,
        "order_4_q_7/4": 272,
    }


def test_naive_exponent_four_superset_is_strictly_filtered(
    certificate: dict,
) -> None:
    audit = certificate["lattice"]["full_discriminant_form_audit"]
    candidates = {
        (row["square"], row["divisibility"]): row
        for row in audit["candidate_filter"]
    }
    assert set(candidates) == {
        (-2, 1),
        (-2, 2),
        (-4, 2),
        (-2, 4),
        (-4, 4),
        (-8, 4),
    }
    assert {
        pair
        for pair, row in candidates.items()
        if row["survives_all_arithmetic_filters"]
    } == {(-2, 1), (-4, 2), (-4, 4)}
    assert {
        pair
        for pair, row in candidates.items()
        if row["self_pairing_divisibility_compatible"]
    } == {(-2, 1), (-2, 2), (-4, 2), (-4, 4), (-8, 4)}
    assert candidates[(-2, 2)]["required_discriminant_q_mod_2"] == "3/2"
    assert candidates[(-2, 4)]["required_discriminant_q_mod_2"] == "15/8"
    assert candidates[(-8, 4)]["required_discriminant_q_mod_2"] == "3/2"
    assert all(
        row["status"] == "excluded_by_full_discriminant_form"
        for pair, row in candidates.items()
        if pair in {(-2, 2), (-8, 4)}
    )
    assert candidates[(-2, 4)]["status"] == (
        "excluded_by_d_divides_r_square_before_discriminant_form"
    )
    assert candidates[(-2, 4)][
        "self_pairing_divisibility_compatible"
    ] is False
    assert audit["after_d_divides_self_pairing_filter"] == [
        {"square": -2, "divisibility": 1},
        {"square": -2, "divisibility": 2},
        {"square": -4, "divisibility": 2},
        {"square": -4, "divisibility": 4},
        {"square": -8, "divisibility": 4},
    ]


def test_all_surviving_types_have_literal_S_witnesses(
    certificate: dict,
) -> None:
    witnesses = certificate["lattice"]["explicit_negative_type_witnesses"]
    assert {
        (row["square"], row["divisibility"])
        for row in witnesses.values()
    } == {(-2, 1), (-4, 2), (-4, 4)}
    assert all(
        row["primitive"]
        and row["dual_class_exact_order"] == row["divisibility"]
        and row["reflection_integral"]
        and row["reflection_preserves_S"]
        and row["hyperbolic_mirror_normal"]
        for row in witnesses.values()
    )

    gram = sp.Matrix(
        json.loads(SOURCES["shared_k3_carrier"].read_text(encoding="utf-8"))
        ["shared_Picard_lattice"]["gram"]
    )
    for row in witnesses.values():
        root = sp.Matrix(row["coordinates"])
        reflection = (
            sp.eye(11)
            - sp.Rational(2, row["square"]) * root * (root.T * gram)
        )
        assert all(value.is_Integer for value in reflection)
        assert reflection.T * gram * reflection == gram


def test_positive_integral_reflection_is_explicitly_out_of_scope(
    certificate: dict,
) -> None:
    boundary = certificate["lattice"][
        "positive_norm_integral_reflection_boundary_witness"
    ]
    assert boundary["coordinates"] == [1, 1] + [0] * 9
    assert (boundary["square"], boundary["divisibility"]) == (2, 1)
    assert boundary["primitive"] is True
    assert boundary["reflection_integral"] is True
    assert boundary["reflection_preserves_S"] is True
    assert boundary["hyperbolic_mirror_normal"] is False


def test_carrier_is_period_level_not_equation_level(certificate: dict) -> None:
    carrier = certificate["period_carrier"]
    assert carrier["moduli_dimension"] == 9
    assert carrier["very_general_full_Picard"] == "NS(X)=S"
    assert carrier["orthogonal_complement"] == "U+E8(-2)+<4>"
    assert carrier["level"].startswith("period-level")
    assert "no explicit common" in carrier["level"]
    assert "Equations" in carrier["wall_scope"]
    assert certificate["scope"]["equation_family"] == "not asserted"
    assert "period-level" in certificate["theorem"][
        "common_period_moduli_criterion"
    ]


def test_full_group_indices_and_noncompactness_gate(certificate: dict) -> None:
    group = certificate["full_reflection_group"]
    assert group["fundamental_facets"] == 13
    assert group["finite_vertices"] == 29
    assert group["cusps"] == 3
    assert group["nonuniform"] is True
    assert group["rational_euler_characteristic"] == "-31/11147673600"
    proof = group["proof_of_subgroup_criterion"]
    assert "fundamental domain" in proof
    assert "disjoint union" in proof
    assert "finite iff the index is finite" in proof
    assert "powers of parabolic" in proof

    calibrations = certificate["calibrations"]
    assert calibrations["F3-V4-H10-123"]["index_in_W_S"] == 65_802_240
    assert calibrations["JK-H10-0007"]["index_in_W_S"] == 16_450_560
    assert calibrations["F3-V4-H10-123"][
        "rational_euler_characteristic"
    ] == "-527/2880"
    assert calibrations["JK-H10-0007"][
        "rational_euler_characteristic"
    ] == "-527/11520"
    assert certificate["theorem"]["compact_no_go_inside_S"].startswith(
        "No compact H10"
    )


def test_calibrations_and_geometric_promotion_are_bounded(
    certificate: dict,
) -> None:
    rows = certificate["calibrations"]
    assert rows["F3-V4-H10-123"]["wall_types"] == {
        "-2,1": 10,
        "-4,2": 31,
        "-4,4": 82,
    }
    assert rows["JK-H10-0007"]["wall_types"] == {
        "-2,1": 8,
        "-4,2": 33,
        "-4,4": 43,
    }
    assert [rows[key]["wall_span_index"] for key in rows] == [2, 2]
    assert [rows[key]["cusps"] for key in rows] == [1732, 862]
    boundary = certificate["promotion_boundary"]
    assert boundary["full_nef_criterion"].startswith(
        "After chamber-specific Weyl re-marking"
    )
    assert "(-4,div2) or (-4,div4)" in boundary["minus_four_obstruction"]
    assert "complete Vinberg" in boundary["volume_gate"]
    assert "finite-index" in boundary["volume_gate"]
    assert certificate["scope"]["surface_automorphisms"] == "not asserted"
    assert certificate["scope"]["new_abstract_polytopes"] == "none"


def test_all_sources_and_self_payload_are_frozen(certificate: dict) -> None:
    assert set(certificate["sources"]) == set(SOURCES)
    for key, path in SOURCES.items():
        source = certificate["sources"][key]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == source[
            "file_sha256"
        ]
        document = json.loads(path.read_text(encoding="utf-8"))
        canonical = hashlib.sha256(
            json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        assert canonical == source["payload_sha256"]

    unsealed = dict(certificate)
    expected = unsealed.pop("payload_sha256")
    observed = hashlib.sha256(
        json.dumps(unsealed, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert observed == expected
    assert certificate["all_checks_pass"] is True
    assert all(certificate["checks"].values())


def test_report_states_exact_inverse_carrier_boundary(certificate: dict) -> None:
    report = (
        ROOT / "results/ue8m2a4_common_k3_moduli_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "every genuine finite-sided coxeter polyhedron" in normalized
    assert "full discriminant form" in normalized
    assert "(-2,2)" in normalized and "(-8,4)" in normalized
    assert "period-level, not equation-level" in normalized
    assert "finite volume, finite index" in normalized
    assert "no compact `h10` coxeter polytope" in normalized
    assert "not a new abstract-polytope claim" in normalized
    assert certificate["payload_sha256"] in report
