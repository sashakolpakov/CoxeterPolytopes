import json
from pathlib import Path

from scripts.certify_ud4x2a1x2_order4_weierstrass_family import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/ud4x2a1x2_order4_weierstrass_family.json"
)


def payload() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == payload()


def test_explicit_equation_kodaira_fibers_and_open_witness() -> None:
    family = payload()["family"]
    assert family["affine_equation"] == (
        "y^2=x^3+t^2*P(t^2)*x+t^3*Q(t^2)"
    )
    assert family["automorphism"] == "sigma(x,y,t)=(-x,i*y,-t)"
    assert family["global_order"] == 4
    assert [row["kodaira"] for row in family["fibers"]] == [
        "I0*",
        "I0*",
        "I2",
        "I1",
    ]
    witness = family["nonempty_open_witness"]
    assert witness["parameters"] == {
        "alpha": 0,
        "beta": 1,
        "delta": 1,
        "gamma": 0,
    }
    assert witness["P_in_u_equals_s_minus_1"] == "u^2-3"
    assert witness["Q_in_u_equals_s_minus_1"] == "u^3+2"
    assert witness["R"] == "31*u**4 - 36*u**2 + 108*u + 108"
    assert witness["R_discriminant"] == -1900537572096
    assert witness["zero_i0star_coefficient"] == -5
    assert witness["infinity_i0star_coefficient"] == 31
    assert witness["resultant_R_P"] == 42849


def test_both_d4_diagram_actions_are_identity_not_triality() -> None:
    action = payload()["d4_component_actions"]
    at_zero = action["at_t_zero"]
    at_infinity = action["at_t_infinity"]
    assert at_zero["tate_root_coordinate"] == "xi=x/t"
    assert at_zero["sigma_action"] == "xi->(-x)/(-t)=xi"
    assert "individually fixed" in at_zero["conclusion"]
    assert "RHS=t^3" in at_zero["actual_blowup_chart"]
    assert "(t,xi=x/t)->(-t,xi)" in at_zero["actual_blowup_chart"]
    assert at_infinity["tate_root_coordinate"] == "xi=X/u"
    assert at_infinity["sigma_action"] == "xi->(-X)/(-u)=xi"
    assert "individually fixed" in at_infinity["conclusion"]
    assert "RHS=v^3" in at_infinity["actual_blowup_chart"]
    assert "(v,xi=X/v)->(-v,xi)" in at_infinity["actual_blowup_chart"]
    assert action["global_conclusion"] == (
        "both central components and all six outer arms are fixed; "
        "neither D4 triality group is used"
    )


def test_exact_generic_ns_and_mordell_weil_saturation() -> None:
    ns = payload()["neron_severi"]
    assert ns["rank_signature_determinant"] == [12, [1, 11, 0], -64]
    assert ns["smith"] == [*[1] * 6, *[2] * 6]
    generic = ns["generic_exactness"]
    assert generic["parameter_moduli_dimension"] == 4
    assert generic["parameter_recovery"]["symbolic_jacobian"] == [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
    assert generic["parameter_recovery"]["jacobian_determinant"] == 1
    assert generic["normalized_coordinate_stabilizer"]["finite"]
    assert generic["marked_to_unmarked_zero_dimensional"][
        "fiber_dimension"
    ] == 0
    assert "rho>=14" in generic["period_argument"]
    assert generic["mordell_weil_rank"] == 0
    assert "h(P)>=4-(1+1+1/2+1/2)=1" in generic[
        "torsion_height_argument"
    ]
    assert generic["conclusion"] == (
        "very general NS=U+D4(-1)^2+A1(-1)^2"
    )


def test_literal_picard_action_fixed_coinvariant_and_glue() -> None:
    action = payload()["picard_action"]
    assert action["a1_swap"] == [[0, 1], [1, 0]]
    assert action["characteristic_polynomial"] == (
        "(lambda - 1)**11*(lambda + 1)"
    )
    assert action["fixed_lattice"] == {
        "a1_diagonal_gram": [[-4]],
        "class": "U+D4(-1)^2+<-4>",
        "determinant": 64,
        "discriminant_q_histogram": {
            "0": 16,
            "1": 16,
            "3/4": 12,
            "7/4": 20,
        },
        "rank": 11,
        "signature": [1, 10, 0],
        "smith": [*[1] * 6, *[2] * 4, 4],
    }
    assert action["coinvariant"] == {
        "class": "<-4>",
        "divisibility_in_NS": 2,
        "generator": "R_+ - R_-",
        "gram": [[-4]],
        "rootless": True,
    }
    assert action["primitive_glue"]["index"] == 2


def test_bh_chamber_and_literal_vinberg_handoffs() -> None:
    data = payload()
    bh = data["brandhorst_hofmann_handoff"]
    assert bh["id"] == "0.4.4.6"
    assert bh["full_NS_action_pairs_integrally_conjugate"]
    assert bh["all_glue_pairs_carried"]
    assert data["chamber_handoff"]["complete_Roulleau_nef_chamber"]
    assert data["hyperbolic_handoff"] == {
        "abstract_polytope_novelty": (
            "bounded-corpus novelty candidate only; this theorem is a new "
            "explicit algebraic realization in the project"
        ),
        "bad_audits": 0,
        "audit_digest_schema": (
            "canonical_json_v1_finite_vertices_ideal_components_edges_"
            "parabolic_audits"
        ),
        "canonical_full_audit_sha256": (
            "79c800f246b68e6eff1908dd1922abdf7b31fa3ac436e43920fce8b22a47cf03"
        ),
        "classification": "finite_volume_noncompact",
        "cusps": 141,
        "dimension": 10,
        "elliptic_edges": 4518,
        "facets": 38,
        "finite_vertices": 526,
        "polytope_id": "R12-D4x2A1x2-p67-H10-38",
    }
    fixed = data["cohomology_and_fixed_loci"]
    assert fixed["fixed_euler_powers_1_2_3"] == {"1": 12, "2": 4, "3": 12}
    assert fixed["bh_generator_fingerprint"]["meaning"] == (
        "8 isolated points and 2 rational curves"
    )
    assert fixed["bh_square_fingerprint"]["meaning"] == (
        "3 rational curves and one genus-2 curve"
    )
