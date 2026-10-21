from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_ss2_2a_h13_50_ud8d4_weierstrass_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    assert rebuilt["all_checks_pass"] is True
    return rebuilt


def test_explicit_six_parameter_weierstrass_family(certificate: dict) -> None:
    family = certificate["weierstrass_family"]
    assert family["parameters"] == ["a", "b", "c", "d", "e", "f"]
    assert family["affine_equation"] == "y^2=x^3+t^2*P(t)*x+t^3*Q(t)"
    assert family["fiber_configuration"] == "I0*+I4*+8I1"
    assert [row["kodaira"] for row in family["fibers"]] == ["I0*", "I4*", "I1"]
    assert family["nonempty_witness"]["parameters"] == [-2] * 6
    assert family["nonempty_witness"][
        "residual_discriminant_polynomial_discriminant"
    ] == -308186920225388788221712962876472295424
    assert family["nonempty_witness"]["resultant_with_P"] == 11489984737344
    assert all(family["checks"].values())


def test_generic_picard_rank_and_mordell_weil_saturation(certificate: dict) -> None:
    picard = certificate["very_general_picard_lattice"]
    torsion = picard["torsion_saturation"]
    assert picard["isometry_class"] == "U+D8(-1)+D4(-1)"
    assert picard["rank_signature_determinant"] == [14, [1, 13, 0], -16]
    assert picard["smith"] == [1] * 10 + [2] * 4
    assert torsion["D8_accessible_component_corrections"] == ["1", "2", "2"]
    assert torsion["D4_accessible_component_corrections"] == ["1", "1", "1"]
    assert torsion["maximum_total_correction"] == 3
    assert torsion["nonzero_section_height_lower_bound"] == 1
    assert "torsion is trivial" in torsion["conclusion"]
    assert all(picard["checks"].values())


def test_exact_equation_to_fifty_mirror_bridge(certificate: dict) -> None:
    bridge = certificate["explicit_lattice_bridge"]
    mirrors = certificate["transported_H13_50_mirrors"]
    assert bridge["determinant_T"] == -1
    assert bridge["N_to_standard_fiber_basis"]["determinant"] == 1
    assert bridge["identity"] == (
        "T*Gram(Roulleau)*T^t=Gram(U+D8(-1)+D4(-1))"
    )
    assert bridge["affine_D8_fiber"]["walls_zero_based"] == [
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        15,
        16,
    ]
    assert all(bridge["checks"].values())
    assert mirrors["count"] == 50
    assert mirrors["square_divisibility_histogram"] == {
        "-2,1": 40,
        "-4,2": 10,
    }
    assert mirrors["generate_full_NS"] is True
    assert mirrors["all_reflections_integral"] is True
    assert mirrors["strict_interior_square"] == 2012
    assert mirrors["minimum_wall_pairing"] == 6
    assert all(mirrors["checks"].values())
    assert all(all(row["checks"].values()) for row in mirrors["records"])


def test_primary_lineage_literal_vinberg_and_scope(certificate: dict) -> None:
    lineage = certificate["primary_lineage"]
    vinberg = certificate["literal_vinberg_handoff"]
    scope = certificate["scope"]
    assert lineage["arxiv_version"] == "2009.09635v6"
    assert lineage["doi"].endswith("10.1002/mana.202200197")
    assert "I4*+I0*+8I1" in lineage["exact_locations"]["Theorem_3_7"]
    assert vinberg["classification"] == "finite_volume_noncompact"
    assert (vinberg["dimension"], vinberg["facets"]) == (13, 50)
    assert (vinberg["finite_vertices"], vinberg["cusps"]) == (45954, 2598)
    assert vinberg["elliptic_edges"] == 391614
    assert vinberg["bad_edges"] == vinberg["bad_parabolic_neighborhoods"] == 0
    assert vinberg["independent_exact_cone_rays"]["hyperideal"] == 0
    assert all(vinberg["checks"].values())
    assert "published" in scope["known_algebraic_family"]
    assert "square-minus-four" in scope["not_full_nef"]
    assert "no surface automorphism" in scope["not_automorphism_fixed"]
    assert "not a new abstract polytope" in scope["novelty"]


def test_report_freezes_positive_gain_and_caveats() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "ss2_2a_h13_50_ud8d4_weierstrass_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "equation-level" in normalized
    assert "i_4^*+i_0^*+8i_1" in normalized
    assert "mordell--weil torsion is trivial" in normalized
    assert "14-by-14" in normalized
    assert "45,954" in normalized and "2,598" in normalized
    assert "not the full nef cone" in normalized
    assert "no surface automorphism" in normalized
