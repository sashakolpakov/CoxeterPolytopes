import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden, PHI
from scripts.certify_og10_golden_h9_fixed_slice import (
    HISTORICAL_FULL_CENSUS,
    build,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/candidates/og10_golden_h9_fixed_slice.json"


def test_frozen_h9_fixed_slice_certificate_is_reproducible() -> None:
    frozen = json.loads(ARTIFACT.read_text())
    result = build()
    assert result == frozen
    assert result["verified"] is True
    assert all(result["checks"].values())


def test_complete_fixed_slice_and_global_ultraparallel_lemma() -> None:
    result = build()
    lattice = result["fixed_lattice"]
    assert lattice["decomposition"] == "<-2*phi> + E8 + <2>"
    assert lattice["identity_inertia"] == [9, 1, 0]
    assert lattice["conjugate_inertia"] == [10, 0, 0]

    theorem = result["global_fixed_slice_theorem"]
    assert theorem["nonfixed_simple_orbits"] == (
        "globally_ultraparallel_and_redundant_on_H9"
    )
    assert theorem["arithmetic_lemma"]["verified"] is True
    assert theorem["arithmetic_lemma"]["sqrt_comparison_squared_margin"] == {
        "a": "25/121",
        "b": "0",
    }


def test_current_and_historical_hyperideal_escape_certificates() -> None:
    result = build()
    current = result["known_61_wall_outer_cone"]
    assert current["facet_count"] == 61
    assert current["hyperideal_witness_active_facets"] == [
        0, 1, 2, 4, 5, 6, 12, 48, 49
    ]
    assert current["hyperideal_witness_ray"] == [
        value.to_json()
        for value in (
            Golden(1),
            2 * PHI,
            4 * PHI,
            6 * PHI,
            Golden(0, 64) / 13,
            Golden(0, 48) / 13,
            Golden(0, 32) / 13,
            Golden(0, 16) / 13,
            Golden(0, 40) / 13,
            PHI / 13,
        )
    ]
    assert current["hyperideal_witness_norm"] == (
        Golden(210, -128) / 169
    ).to_json()

    historical = result["historical_44_wall_outer_cone"]
    assert historical["facet_count"] == 44
    assert historical["hyperideal_witness_norm"] == Golden(178, -110).to_json()
    assert historical["full_exact_census"] == HISTORICAL_FULL_CENSUS


def test_bugaenko_polytope_is_a_complete_right_angled_facet() -> None:
    facet = build()["right_angled_bugaenko_facet"]
    assert facet["verified_diagram_gram"] is True
    assert facet["first_nine_root_coordinate_determinant"] == {
        "a": "-3",
        "b": "2",
    }
    assert facet["first_nine_root_determinant_norm"] == "-1"
    assert facet["verified_integral_lattice_identification"] is True
    assert facet["complete_facet_wall_count"] == 11
    assert facet["complete_facet_f_vector"] == [
        41, 164, 316, 374, 294, 156, 54, 11, 1
    ]
    assert "not finite volume" in facet["scope"]


def test_full_h2_lattices_and_kahler_regularization_gate() -> None:
    geometry = build()["algebro_geometric_involution"]
    eigenspaces = geometry["full_H2_eigensublattices"]
    assert eigenspaces["verified"] is True
    assert eigenspaces["fixed_lattice"] == {
        "rank": 20,
        "inertia": [1, 19, 0],
        "determinant": -16,
        "smith_diagonal": [1] * 16 + [2] * 4,
        "two_elementary_invariants": {"r": 20, "a": 4, "delta": 1},
        "classification_isometry": "<2> + <-2>^3 + E8(-1)^2",
    }
    assert eigenspaces["anti_invariant_lattice"] == {
        "rank": 4,
        "inertia": [2, 2, 0],
        "determinant": 48,
        "smith_diagonal": [2, 2, 2, 6],
        "classification_isometry": "<2>^2 + <-2> + <-6>",
    }
    regularization = geometry["kahler_wall_containment_audit"]
    assert regularization["verified"] is True
    assert regularization["rational_norms_of_ratios"] == [
        "1/11", "4/11", "5/11"
    ]
