import json
from pathlib import Path

import pytest

from scripts.certify_ue8d4_named_c12_member import (
    build_certificate,
    point_count_audit,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/ue8d4_named_c12_member.json"


def payload() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == payload()


def test_named_member_point_counts_and_exact_ns() -> None:
    certificate = payload()
    assert certificate["surface"]["affine_equation"] == (
        "y^2 = x^3 + t^3*(t^2-1)*(t^2-4)"
    )
    arithmetic = certificate["point_count_and_tate_certificate"]
    assert [row["resolved_surface_points"] for row in arithmetic["counts"]] == [
        127,
        2913,
        122452,
        5803953,
    ]
    assert arithmetic["quotient_power_traces_n1_to_n4"] == [0, -28, 0, 12740]
    assert arithmetic["frobenius_polynomial"] == (
        "T**8 + 14*T**6 - 3087*T**4 + 33614*T**2 + 5764801"
    )
    assert not any(arithmetic["cyclotomic_gcd_degrees"].values())
    assert certificate["picard_rank_conclusion"]["characteristic_zero_picard_rank"] == 14
    assert certificate["picard_rank_conclusion"]["exact_ns"] == (
        "U+E8(-1)+D4(-1)"
    )


def test_c12_relations_and_resolved_component_actions() -> None:
    certificate = payload()
    automorphisms = certificate["automorphisms"]
    assert automorphisms["relations"]["cyclic_group"] == (
        "<zeta,iota>=<h>=C12; zeta=h^10 and iota=h^3"
    )
    assert automorphisms["picard_actions"] == {
        "order12_h": "identity on U+E8(-1), outer triality on D4(-1)",
        "order4_iota": "identity on NS",
        "order6_zeta": "identity on U+E8(-1), outer triality on D4(-1)",
    }
    resolved = automorphisms["resolved_fiber_action"]
    assert "X->X" in resolved["iota"]
    assert "3-cycle" in resolved["zeta"]
    assert certificate["symmetric_c12_subfamily"]["moduli_dimension"] == 1


def test_fixed_coinvariant_glue_and_bh_matches() -> None:
    certificate = payload()
    triality = certificate["fixed_coinvariant_glue"]["triality_action"]
    assert triality["fixed_lattice"]["rank_signature_determinant"] == [
        12,
        [1, 11, 0],
        -3,
    ]
    assert triality["coinvariant_lattice"] == {
        "class": "A2(-2)",
        "gram": [[-4, 2], [2, -4]],
        "rootless": True,
    }
    assert triality["fixed_plus_coinvariant_index"] == 3
    assert triality["quotient_smith"] == [*[1] * 13, 3]

    bh = certificate["brandhorst_hofmann"]
    assert [row["id"] for row in bh["rows"]] == ["0.4.3.1", "0.6.3.1", "0.12.1.1"]
    assert [row["member_row_count"] for row in bh["rows"]] == [79, 150, 48]
    assert bh["order6_order12_exact_pair_equality"]["same_ns_gram_entrywise"]
    assert bh["order6_order12_exact_pair_equality"]["same_ns_action_entrywise"]


def test_hyperbolic_handoffs_are_positive_but_not_novelty_overclaims() -> None:
    handoff = payload()["coxeter_handoff"]
    assert handoff["order4_full_chamber"]["dimension"] == 13
    assert handoff["order4_full_chamber"]["facets"] == 15
    triality = handoff["orders6_and12_triality_chamber"]
    assert triality["classification"] == "finite_volume_noncompact"
    assert (triality["dimension"], triality["facets"]) == (11, 13)
    assert (triality["finite_vertices"], triality["cusps"]) == (27, 1)
    assert "abstract chambers are known" in handoff["novelty"]


def test_optional_python_flint_recount() -> None:
    pytest.importorskip("flint")
    recounted = point_count_audit(recount=True)
    assert [row["resolved_surface_points"] for row in recounted["counts"]] == [
        127,
        2913,
        122452,
        5803953,
    ]
