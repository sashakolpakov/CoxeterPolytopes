import json
from pathlib import Path

from scripts.certify_ue8a1x2_order4_weierstrass_family import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT / "data/enumerations/ue8a1x2_order4_weierstrass_family.json"
)


def payload() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == payload()


def test_explicit_equation_and_kodaira_fibers() -> None:
    data = payload()
    family = data["family"]
    assert family["affine_equation"] == "y^2=x^3+A(t^2)*x+t*C(t^2)"
    assert family["automorphism"] == "sigma(x,y,t)=(-x,i*y,-t)"
    assert family["global_order"] == 4
    assert family["two_form_multiplier"] == "-i (a primitive fourth root)"
    assert [row["kodaira"] for row in family["fibers"]] == [
        "II*",
        "I2",
        "I1",
    ]
    assert family["nonempty_open_witness"] == {
        "Q_over_27_at_s_zero": -4,
        "Q_over_27_discriminant": 283904,
        "Q_over_27_in_u_equals_s_minus_1": (
            "u**5 + u**4 - 2*u**3 + 2*u**2 + 5*u - 3"
        ),
        "parameters": {"alpha": 0, "beta": 0, "delta": 1, "gamma": 0},
    }


def test_exact_generic_ns_mw_and_saturation() -> None:
    ns = payload()["neron_severi"]
    assert ns["rank_signature_determinant"] == [12, [1, 11, 0], -4]
    assert ns["smith"] == [*[1] * 10, 2, 2]
    assert ns["generic_exactness"]["parameter_moduli_dimension"] == 4
    generic = ns["generic_exactness"]
    assert generic["parameter_recovery"] == {
        "jacobian_determinant": 1,
        "observables": ["A'(1)", "A''(1)/2", "C''(1)/2", "C'''(1)/6"],
        "recovered_parameters": ["alpha", "beta", "gamma", "delta"],
        "symbolic_jacobian": [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ],
    }
    assert generic["normalized_coordinate_stabilizer"]["finite"]
    assert generic["marked_to_unmarked_zero_dimensional"][
        "fiber_dimension"
    ] == 0
    assert "strict jump would have rho>=14" in ns["generic_exactness"][
        "argument"
    ]
    assert "signature (1,11)" in ns["primitive_at_rank_12"]
    assert ns["mordell_weil_very_general"] == {
        "rank": 0,
        "reason": (
            "Shioda--Tate after exact rho=12; the no-overlattice argument "
            "rules out a torsion saturation"
        ),
        "torsion": "trivial",
    }


def test_literal_a1_swap_fixed_coinvariant_and_glue() -> None:
    action = payload()["picard_action"]
    assert action["a1_swap"] == [[0, 1], [1, 0]]
    assert action["characteristic_polynomial"] == (
        "(lambda - 1)**11*(lambda + 1)"
    )
    assert action["fixed_lattice"] == {
        "a1_diagonal_gram": [[-4]],
        "class": "U+E8(-1)+<-4>",
        "determinant": 4,
        "discriminant_q_histogram": {"0": 1, "1": 1, "7/4": 2},
        "rank": 11,
        "signature": [1, 10, 0],
        "smith": [*[1] * 10, 4],
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
    assert data["brandhorst_hofmann_handoff"]["id"] == "0.4.4.1"
    assert data["brandhorst_hofmann_handoff"][
        "full_NS_action_pairs_integrally_conjugate"
    ]
    assert data["brandhorst_hofmann_handoff"]["all_glue_pairs_carried"]
    assert data["chamber_handoff"]["complete_Roulleau_nef_chamber"]
    assert data["hyperbolic_handoff"] == {
        "abstract_polytope_novelty": (
            "known; this is a new explicit algebraic realization in the project"
        ),
        "bad_audits": 0,
        "audit_digest_schema": (
            "canonical_json_v1_finite_vertices_ideal_components_edges_"
            "parabolic_audits"
        ),
        "canonical_full_audit_sha256": (
            "d20267d1d38090f54d650d20d88d71448aa0393686cf64705e9780f91c37bd5a"
        ),
        "classification": "finite_volume_noncompact",
        "cusps": 2,
        "dimension": 10,
        "elliptic_edges": 99,
        "facets": 12,
        "finite_vertices": 17,
        "polytope_id": "R12-E8A1x2-C2-H10-12",
    }
    fixed = data["cohomology_and_fixed_loci"]
    assert fixed["fixed_euler_powers_1_2_3"] == {"1": 12, "2": 4, "3": 12}
    assert fixed["bh_generator_fingerprint"]["meaning"] == (
        "8 isolated points and 2 rational curves"
    )
    assert fixed["bh_square_fingerprint"]["meaning"] == (
        "5 rational curves and one genus-4 curve"
    )


def test_lineage_is_positive_but_priority_bounded() -> None:
    lineage = payload()["artebani_sarti_lineage"]
    assert "Example 6.3" in lineage["generic_even_odd_model"]
    assert "Example 7.2" in lineage["node_specialization"]
    assert "II*/E8" in lineage["our_fiber_root_specialization"]
    assert "no worldwide priority" in lineage["novelty_scope"]
