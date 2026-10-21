import json
from pathlib import Path

from scripts.certify_ue8d4_order6_weierstrass_family import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/ue8d4_order6_weierstrass_family.json"


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_explicit_weierstrass_family_and_fibers() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert payload["family"]["affine_equation"] == (
        "y^2 = x^3 + t^3*(t-1)*(t-a)*(t-b)*(t-c)"
    )
    assert [row["kodaira"] for row in payload["family"]["fibers"]] == [
        "I0*",
        "II",
        "II*",
    ]
    assert payload["dillies_lineage"]["noether_lefschetz_collision"] == (
        "II + IV -> I0*"
    )
    lineage = payload["dillies_lineage"]
    assert "zeta and zeta^2" in lineage["pre_collision_match_scope"]
    assert "does not determine the zeta^3" in lineage["pre_collision_match_scope"]
    assert lineage["zeta3_collision_change"]["before_collision_genus"] == 5
    assert lineage["zeta3_collision_change"]["after_collision_genus"] == 3


def test_exact_generic_ns_and_triality() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    ns = payload["neron_severi"]
    action = payload["picard_action"]
    assert ns["rank_signature_determinant"] == [14, [1, 13, 0], -4]
    assert ns["generic_exactness"]["parameter_moduli_dimension"] == 3
    generic_argument = ns["generic_exactness"]["argument"]
    assert "rho=22-rank(T) are even" in generic_argument
    assert "any strict jump would therefore have rho>=16" in generic_argument
    assert action["characteristic_polynomial"] == (
        "(lambda - 1)**12*(lambda**2 + lambda + 1)"
    )
    assert action["fixed_lattice"] == {
        "class": "U+E8(-1)+A2(-1)",
        "d4_fixed_gram": [[-2, 3], [3, -6]],
        "determinant": -3,
        "rank": 12,
    }
    assert action["coinvariant"] == {
        "class": "A2(-2)",
        "gram": [[-4, 2], [2, -4]],
        "rootless": True,
    }


def test_bh_and_vinberg_handoffs() -> None:
    payload = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert payload["brandhorst_hofmann_handoff"]["id"] == "0.6.3.1"
    assert payload["brandhorst_hofmann_handoff"]["same_d4_action_entrywise"]
    assert payload["cohomology_and_fixed_loci"]["fixed_euler_powers_1_2_3"] == {
        "1": 17,
        "2": 9,
        "3": 8,
    }
    assert payload["hyperbolic_handoff"] == {
        "abstract_polytope_novelty": (
            "known; this is a new explicit realization in the project"
        ),
        "bad_edges": 0,
        "classification": "finite_volume_noncompact",
        "cusps": 1,
        "dimension": 11,
        "elliptic_edges": 162,
        "facets": 13,
        "finite_vertices": 27,
    }
