import json
from pathlib import Path

from scripts.certify_ue8d6_order4_weierstrass_family import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/ue8d6_order4_weierstrass_family.json"


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_equation_parameter_open_and_fibers() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    family = payload["family"]
    assert family["affine_equation"] == (
        "y^2 = x^3 + t^2*(t^2-3)*x + t^3*(2+b*t^2+c*t^4)"
    )
    assert family["automorphism"] == "sigma(x,y,t)=(-x,i*y,-t)"
    assert family["global_order"] == 4
    assert family["two_form_multiplier"] == "-i"
    assert [row["kodaira"] for row in family["fibers"]] == [
        "I2*",
        "I1",
        "II*",
    ]
    assert family["open_nonempty_witness"] == {
        "Q": "27*u**3 + 4*u**2 + 72*u + 108",
        "Q_at_3": 1089,
        "b": 0,
        "c": 1,
        "disc_u_Q": -254721456,
    }
    assert "432*(3*b + 9*c + 2)**3" in family["disc_u_Q"]


def test_exact_generic_ns_and_outer_d6_action() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    ns = payload["neron_severi"]
    action = payload["picard_action"]
    assert ns["rank_signature_determinant"] == [16, [1, 15, 0], -4]
    assert ns["generic_exactness"]["parameter_moduli_dimension"] == 2
    assert "lambda^2=1" in ns["generic_exactness"]["normalization_argument"]
    assert "rho>=18" in ns["generic_exactness"]["period_argument"]
    assert action["characteristic_polynomial"] == (
        "(lambda - 1)**15*(lambda + 1)"
    )
    assert action["fixed_lattice"]["rank_signature_determinant"] == [
        15,
        [1, 14, 0],
        4,
    ]
    assert action["coinvariant"] == {
        "class": "<-4>",
        "gram": [[-4]],
        "rootless": True,
    }
    assert "Proposition 1" in action["description"]


def test_discriminant_and_primitive_glue_match() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    action = payload["picard_action"]
    assert action["discriminant_action"]["our_permutation"] == [0, 1, 3, 2]
    assert action["discriminant_action"]["archived_permutation"] == [0, 2, 1, 3]
    assert action["discriminant_action"][
        "our_to_archived_class_identification"
    ] == [0, 3, 1, 2]
    assert action["primitive_glue"] == {
        "adjoined_vector": "alpha_5=(s+d)/2",
        "coinvariant_order_two_class": "(alpha_5-alpha_6)/2",
        "fixed_order_two_class": "(alpha_5+alpha_6)/2",
        "index": 2,
        "source_smith_pair": [2, 2],
    }


def test_bh_fixed_locus_chamber_and_vinberg_handoffs() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert payload["brandhorst_hofmann_handoff"]["id"] == "0.4.2.2"
    assert payload["brandhorst_hofmann_handoff"]["full_group_row"] == "1.4.2.1"
    loci = payload["cohomology_and_fixed_loci"]
    assert loci["generator"]["isolated_points"] == 10
    assert loci["generator"]["fixed_rational_curves"] == 3
    assert loci["square"]["fixed_rational_curves"] == 7
    assert loci["square"]["fixed_positive_genus_curves"] == [
        {"count": 1, "genus": 2}
    ]
    assert payload["kahler_chamber_handoff"]["fixed_section"] == (
        "R16-UE8D6-C2-H14-17"
    )
    hyperbolic = payload["hyperbolic_handoff"]
    assert hyperbolic["classification"] == "finite_volume_noncompact"
    assert (hyperbolic["dimension"], hyperbolic["facets"]) == (14, 17)
    assert hyperbolic["bad_edges"] == hyperbolic["hyperideal_vertices"] == 0
    assert payload["artebani_sarti_lineage"]["locator"] == (
        "Example 6.3 and Proposition 1"
    )

