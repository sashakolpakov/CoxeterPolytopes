import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from scripts.certify_og10_golden_h9_infinite_stabilizer import build


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/candidates/og10_golden_h9_infinite_stabilizer.json"


def test_frozen_infinite_stabilizer_certificate_is_reproducible() -> None:
    result = build()
    assert result == json.loads(ARTIFACT.read_text())
    assert result["verified"] is True
    assert all(result["checks"].values())


def test_bugaenko_reversal_is_a_complete_chamber_symmetry() -> None:
    result = build()
    reversal = result["bugaenko_facet_reversal"]
    assert reversal["diagram_permutation"] == [
        8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 9
    ]
    assert "exactly one adjacent chamber on each side" in reversal[
        "complete_chamber_argument"
    ]
    assert result["checks"][
        "bugaenko_reversal_maps_all_11_oriented_facet_walls"
    ]
    assert result["checks"][
        "bugaenko_reversal_fixes_facet_normal_and_preserves_u_side"
    ]


def test_product_has_infinite_order_and_forces_infinite_volume() -> None:
    result = build()
    product = result["infinite_order_product"]
    assert product["trace"] == Golden(8, 8).to_json()
    assert Golden.from_json(product["trace_minus_dimension"]).sign() > 0
    assert result["group_theoretic_conclusion"][
        "complete_chamber_stabilizer"
    ] == "infinite"
    assert result["group_theoretic_conclusion"][
        "P9_hyperbolic_volume"
    ] == "infinite"
    assert result["ambient_consequence"][
        "golden_H10_chamber_volume"
    ] == "infinite"
    assert result["ambient_consequence"][
        "enhanced_OG10_birational_centralizer"
    ] == "infinite"
