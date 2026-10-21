import json
from pathlib import Path

from scripts.certify_ud8a1x2_order4_weierstrass_family import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/ud8a1x2_order4_weierstrass_family.json"


def payload() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == payload()


def test_explicit_equation_and_kodaira_fibers() -> None:
    family = payload()["family"]
    assert family["affine_equation"] == "y^2=x^3+P(t^2)*x+t*Q(t^2)"
    assert family["automorphism"] == "sigma(x,y,t)=(-x,i*y,-t)"
    assert family["global_order"] == 4
    assert [row["kodaira"] for row in family["fibers"]] == [
        "I4*",
        "I2",
        "I1",
    ]
    witness = family["nonempty_open_witness"]
    assert witness["parameters"] == {"alpha": 0, "beta": 1, "gamma": 0, "r": 1}
    assert witness["P"] == "-3*s^3+s-1"
    assert witness["Q"] == "2*s^4-s+1"
    assert witness["R"] == "108*s**5 - 36*s**3 + 31*s - 4"
    assert witness["R_discriminant"] == 6477927509458944
    assert witness["R_at_s_zero"] == -4
    assert witness["resultant_R_P"] == -1440747


def test_d8_diagram_action_is_identity_not_spinor_swap() -> None:
    action = payload()["d8_component_action"]
    assert "xi=X/u is invariant" in action["repeated_tate_root_coordinate"]
    assert "all three branch germs" in action["repeated_tate_root_coordinate"]
    assert action["branchwise_fixed"] is True
    assert "spinor nodes are not exchanged" in action["resolution_conclusion"]
    assert action["standard_d8_equation"] == "Y^2+u*V^2+kappa*u^7=0"
    actual = action["actual_weierstrass_germ"]
    assert actual["translation"] == "V=X-u"
    assert "q0=V/u is fixed" in actual["first_blowup"]
    assert "q1=V/u^2 is negated" in actual["second_blowup"]
    assert "q2=V/u^3 is fixed" in actual["third_blowup"]
    assert actual["slope_discriminant"] == "alpha^2-12*beta-12*gamma"


def test_exact_generic_ns_and_mordell_weil_saturation() -> None:
    ns = payload()["neron_severi"]
    assert ns["rank_signature_determinant"] == [12, [1, 11, 0], -16]
    assert ns["smith"] == [*[1] * 8, *[2] * 4]
    generic = ns["generic_exactness"]
    assert generic["parameter_moduli_dimension"] == 4
    assert generic["parameter_recovery"]["rational_identity"] == (
        "r=-3*Q(1)/(2*P(1)) on r!=0"
    )
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
    assert "h(P)>=4-3=1" in generic["torsion_height_argument"]
    assert generic["conclusion"] == "very general NS=U+D8(-1)+A1(-1)^2"


def test_literal_picard_action_fixed_coinvariant_and_glue() -> None:
    action = payload()["picard_action"]
    assert action["a1_swap"] == [[0, 1], [1, 0]]
    assert action["characteristic_polynomial"] == (
        "(lambda - 1)**11*(lambda + 1)"
    )
    assert action["fixed_lattice"] == {
        "a1_diagonal_gram": [[-4]],
        "class": "U+D8(-1)+<-4>",
        "determinant": 16,
        "discriminant_q_histogram": {
            "0": 4,
            "1": 4,
            "3/4": 2,
            "7/4": 6,
        },
        "rank": 11,
        "signature": [1, 10, 0],
        "smith": [*[1] * 8, 2, 2, 4],
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
    assert bh["id"] == "0.4.4.2"
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
            "00977c1aa97e62df0ac6bbd22e635cd3b37dbf38e22b3cf52770ecd7a42c15fb"
        ),
        "classification": "finite_volume_noncompact",
        "cusps": 10,
        "dimension": 10,
        "elliptic_edges": 294,
        "facets": 15,
        "finite_vertices": 42,
        "polytope_id": "R12-D8A1x2-C2a-H10-15",
    }
    fixed = data["cohomology_and_fixed_loci"]
    assert fixed["fixed_euler_powers_1_2_3"] == {"1": 12, "2": 4, "3": 12}
    assert fixed["bh_square_fingerprint"]["meaning"] == (
        "4 rational curves and one genus-3 curve"
    )
