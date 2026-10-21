import json

import sympy as sp

from scripts.certify_looijenga_h4_inverse_carrier import (
    BOUNDARY,
    BOUNDARY_ORDER,
    CANONICAL,
    OUTPUT,
    PERIOD_BASIS,
    WALLS,
    WALL_ORDER,
    _extreme_rays,
    _numerical_exceptional_classes,
    _period_roots,
    _reflection_matrix,
    build_certificate,
    pair,
    square,
    verify_archive,
    wall_gram,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_polytope_reconstructs_unimodular_picard_canonical_and_boundary() -> None:
    payload = build_certificate()
    assert payload["all_checks_pass"]
    assert wall_gram().rank() == 5
    assert payload["coxeter_input"]["wall_gram_inertia"] == [1, 4, 3]
    assert payload["coxeter_input"]["unimodular_wall_basis_determinant"] in {-1, 1}
    assert square(CANONICAL) == 5
    assert tuple(
        sum(BOUNDARY[name][index] for name in BOUNDARY_ORDER)
        for index in range(5)
    ) == tuple(-value for value in CANONICAL)
    assert payload["surface_reconstruction"]["boundary_gram"] == [
        [-1, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ]


def test_period_lattice_and_exact_root_stratum() -> None:
    assert pair(PERIOD_BASIS["A"], PERIOD_BASIS["A"]) == -2
    assert pair(PERIOD_BASIS["delta"], PERIOD_BASIS["delta"]) == -2
    assert pair(PERIOD_BASIS["A"], PERIOD_BASIS["delta"]) == 1
    assert _period_roots() == [
        (-1, -1),
        (-1, 0),
        (0, -1),
        (0, 1),
        (1, 0),
        (1, 1),
    ]
    period = build_certificate()["period_problem"]
    assert period["exact_carrier_stratum"] == "u=1 and t!=1"
    assert period["stratum_isomorphism"] == "C* minus {1}"
    assert period["period_trivial_roots_on_stratum"] == ["A", "-A"]
    assert period["chosen_nef_branch"].startswith("x.A>=0")


def test_complete_negative_curve_enumeration_and_decompositions() -> None:
    numerical = _numerical_exceptional_classes()
    assert len(numerical) == 10
    assert all(square(vector) == -1 for vector in numerical)
    assert all(pair(CANONICAL, vector) == -1 for vector in numerical)
    certificate = build_certificate()["negative_curve_completeness"]
    assert certificate["irreducible_negative_curves"] == list(WALL_ORDER)
    assert certificate["reducible_numerical_exceptional_decompositions"] == {
        "H-E1-E3": ["A", "E4"],
        "H-E1-E4": ["A", "E3"],
        "H-E3-E4": ["A", "E1"],
    }
    assert certificate["redundant_full_nef_inequalities"] == {
        "D2": ["A", "E1", "E4"],
        "D3": ["A", "E1", "E3"],
        "H-E1-E3": ["A", "E4"],
        "H-E1-E4": ["A", "E3"],
        "H-E3-E4": ["A", "E1"],
    }


def test_all_angles_and_reflections_are_exactly_coxeter_arithmetic() -> None:
    payload = build_certificate()["coxeter_input"]
    profile = {
        relation: sum(
            row["relation"] == relation for row in payload["coxeter_pairs"]
        )
        for relation in ("m2", "m4", "parallel")
    }
    assert profile == {"m2": 19, "m4": 3, "parallel": 6}
    ambient = sp.diag(1, -1, -1, -1, -1)
    for name in WALL_ORDER:
        reflection = _reflection_matrix(WALLS[name])
        assert all(value.q == 1 for value in reflection)
        assert reflection.det() == -1
        assert reflection.T * ambient * reflection == ambient
        assert reflection * reflection == sp.eye(5)


def test_all_facets_and_extreme_rays_prove_exact_finite_volume() -> None:
    payload = build_certificate()
    witnesses = payload["cone_certificate"]["facet_witnesses"]
    assert len(witnesses) == 8
    assert all(row["facet_pairing"] == 0 for row in witnesses)
    assert all(row["minimum_other_pairing"] > 0 for row in witnesses)
    rays = _extreme_rays()
    assert len(rays) == 9
    assert [ray["kind"] for ray in rays].count("finite") == 5
    assert [ray["kind"] for ray in rays].count("ideal") == 4
    assert all(ray["square"] >= 0 for ray in rays)
    assert payload["cone_certificate"]["decision"] == "finite_volume_noncompact"


def test_frozen_archive_replays_without_drift() -> None:
    archived = _archive()
    assert archived == build_certificate()
    assert archived["certificate_sha256"] == (
        "8bc33c023b1a2a6f95b5ade6d7c57dc760b54d6bd5017e262ec3cb47407324d7"
    )
    assert verify_archive() == {
        "certificate_sha256": archived["certificate_sha256"],
        "decision": "finite_volume_noncompact",
        "facet_count": 8,
        "finite_vertex_count": 5,
        "ideal_vertex_count": 4,
        "vertex_count": 9,
    }
