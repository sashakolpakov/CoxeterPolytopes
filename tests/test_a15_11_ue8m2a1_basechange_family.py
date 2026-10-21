import json
from pathlib import Path

from scripts.certify_a15_11_ue8m2a1_basechange_family import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/a15_11_ue8m2a1_basechange_family.json"
)


def payload() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == payload()


def test_normalized_equation_open_witness_and_fibers() -> None:
    family = payload()["family"]
    assert family["parameter_count"] == 9
    assert family["k3_weierstrass_model"] == (
        "X: y^2=x^3+A(s^2)*x+B(s^2)"
    )
    assert family["deck_involution"] == "iota(s,x,y)=(-s,x,y)"
    assert family["two_form_multiplier"] == "-1 (non-symplectic)"
    assert family["normalized_i1"] == {
        "D_order": 1,
        "D_prime": 108,
        "c4": 144,
        "place": "t=0",
    }
    assert family["fibers"] == {
        "R": {"I1": 12, "smooth_marked_at_infinity": 1},
        "X": {"I1": 22, "I2": 1, "smooth_at_infinity": 1},
        "euler_sum_X": 24,
    }
    witness = family["nonempty_open_witness"]
    assert witness["A"] == "t**4 + t - 3"
    assert witness["B"] == "t**6 + 2"
    assert witness["Q"] == (
        "31*t**11 + 12*t**8 - 36*t**7 + 120*t**5 - 72*t**4 + "
        "108*t**3 + 4*t**2 - 36*t + 108"
    )
    assert witness["Q_discriminant"] != 0
    assert witness["resultant_A_Q"] != 0
    assert witness["infinity_discriminant"] == 31
    assert "not asserted to be generic" in witness["claim_scope"]


def test_fixed_locus_and_literal_integral_lattice_identification() -> None:
    fixed = payload()["fixed_locus_and_invariant_lattice"]
    assert fixed["fixed_locus_g_k"] == [1, 1]
    assert fixed["nikulin_invariants"] == {"a": 9, "delta": 1, "r": 11}
    assert fixed["standard_lattice"] == "F=U+E8(-2)+<-2>"
    assert fixed["signature"] == [1, 10, 0]
    assert fixed["determinant"] == 512
    assert fixed["smith"] == [1, 1, *([2] * 9)]
    assert fixed["gram"][:2] == [
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ]
    assert fixed["gram"][-1][-1] == -2
    geometric = fixed["geometric_integral_identification"]
    assert "doubles every pairing" in geometric["E8(-2)"]
    assert "nonidentity component" in geometric["<-2>"]
    assert "index in H^2(X,Z)^iota is one" in geometric["saturation"]
    assert "not the exceptional component" in fixed[
        "ramified_i1_local_model"
    ]


def test_very_general_full_ns_is_not_bare_period_map_existence() -> None:
    generic = payload()["generic_full_picard"]
    assert generic["claim"] == (
        "NS(X)=U+E8(-2)+<-2> for a very general parameter"
    )
    assert generic["parameter_moduli_dimension"] == 9
    assert generic["two_elementary_pair_moduli_dimension"] == 9
    assert generic["parameter_recovery"]["jacobian_determinant"] == 1
    assert generic["parameter_recovery"]["symbolic_jacobian"] == [
        [1 if row == column else 0 for column in range(9)]
        for row in range(9)
    ]
    assert "recover R" in generic["quotient_recovery"]
    assert "zero-dimensional generic fibers" in generic[
        "quotient_recovery"
    ]
    assert "dominates" in generic["noether_lefschetz_step"]
    assert "very general" in generic["scope"]


def test_literal_a15_11_h10_21_handoff_and_finite_volume() -> None:
    handoff = payload()["coxeter_handoff"]
    assert handoff["source_polytope"] == "A15_11_1_C02 / H10/21"
    assert handoff["same_standard_gram_literal"]
    assert handoff["walls"] == 21
    assert handoff["all_reflections_integral_on_very_general_full_NS"]
    assert handoff["finite_volume"] == {
        "classification": "finite_volume_noncompact_Coxeter_polytope",
        "elliptic_edges": 4159,
        "finite_vertices": 607,
        "ideal_vertices": 50,
        "zero_bad_audits": True,
    }
    assert "full NS" in handoff["meaning"]


def test_lineage_and_claim_boundaries_are_explicit() -> None:
    data = payload()
    lineage = data["lineage"]
    assert "Garbagnati" not in json.dumps(lineage)  # keyed, not hand-waved
    assert lineage["garbagnati_salgado_double_cover"]["primary"] == (
        "https://arxiv.org/abs/1703.02783"
    )
    assert "Proposition 2.8" in lineage[
        "garbagnati_salgado_genus_one_classification"
    ]["use"]
    assert "no worldwide priority" in lineage["novelty_scope"]
    scope = data["scope"]
    assert scope["not_a_new_abstract_polytope"]
    assert scope["not_claimed_full_nef_cone"]
    assert scope["reflections_not_claimed_surface_automorphisms"]
    assert scope["square_minus_four_walls_not_claimed_curve_walls"]
    assert data["all_checks_pass"]
