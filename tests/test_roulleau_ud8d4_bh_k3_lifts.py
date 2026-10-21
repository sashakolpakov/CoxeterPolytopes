from __future__ import annotations

import json

import sympy as sp

from scripts.certify_roulleau_ud8d4_bh_k3_lifts import (
    CASES,
    OUTPUT,
    build_certificate,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


CERTIFICATE = json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_two_genuine_k3_lifts_replay_offline() -> None:
    assert build_certificate() == CERTIFICATE
    assert CERTIFICATE["all_checks_pass"] is True
    assert CERTIFICATE["summary"] == {
        "compact": 0,
        "dimensions": [11, 10],
        "finite_volume_noncompact": 2,
        "genuine_k3_automorphism_lifts": 2,
        "new_abstract_diagram_claims": 0,
    }


def test_database_rows_and_picard_restriction_orders_are_exact() -> None:
    assert [case["brandhorst_hofmann_row"]["id"] for case in CERTIFICATE["cases"]] == [
        "0.8.1.3",
        "0.16.0.3",
    ]
    assert [
        case["kahler_torelli_handoff"]["global_order_vs_picard_order"]
        for case in CERTIFICATE["cases"]
    ] == [[8, 2], [16, 4]]
    for case in CERTIFICATE["cases"]:
        row = case["brandhorst_hofmann_row"]
        assert row["generic_transcendental_lattice"]["rank"] == 8
        assert row["generic_transcendental_lattice"]["signature"] == [2, 6, 0]
        assert row["generic_picard_lattice"]["nikulin_invariants"] == [14, 4, 0]
        assert case["kahler_torelli_handoff"][
            "database_row_is_a_genuine_k3_automorphism_class"
        ]


def test_fixed_discriminant_maps_and_primitive_glue_are_exact() -> None:
    expected = [([4, 4], [[1, 1], [1, 0]], 16), ([8], [[7]], 8)]
    for case, (invariants, matrix, count) in zip(CERTIFICATE["cases"], expected):
        conjugacy = case["exact_action_conjugacy"]
        discr = conjugacy["fixed_lattices"]["discriminant_isometry"]
        assert discr["source_invariant_factors"] == invariants
        assert discr["matrix_on_nontrivial_smith_coordinates"] == matrix
        assert discr["elements_checked"] == count
        assert discr["bijective"] and discr["quadratic_form_preserved_mod_2Z"]
        assert conjugacy["primitive_glue"]["index_of_fixed_plus_coinvariant"] == 4
        assert conjugacy["primitive_glue"]["all_four_pairs_carried_exactly"]
        assert len(conjugacy["primitive_glue"]["archived_pairs_in_smith_coordinates"]) == 4
        assert conjugacy["conclusion"].endswith(
            "conjugacy of the full Picard-lattice actions"
        )


def test_coinvariants_are_integrally_intertwined_and_rootless() -> None:
    for specification, case in zip(CASES, CERTIFICATE["cases"]):
        coin = case["exact_action_conjugacy"]["coinvariant_lattices"]
        archived = sp.Matrix(coin["archived_gram"])
        bh = sp.Matrix(coin["bh_gram"])
        intertwiner = sp.Matrix(
            coin["basis_intertwiner_P_with_PGarchivePt_equal_Gbh"]
        )
        assert abs(intertwiner.det()) == 1
        assert intertwiner * archived * intertwiner.T == bh
        assert exact_inertia(bh) == (0, bh.rows, 0)
        assert all((bh / 2)[i, i] % 2 == 0 for i in range(bh.rows))
        assert coin["action_intertwined"] and coin["rootless"]
        assert specification["expected"]["coinvariant_smith"] == coin["smith"]


def test_vinberg_handoffs_are_finite_volume_not_merely_lorentzian() -> None:
    assert [
        (
            case["hyperbolic_polytope"]["dimension"],
            case["hyperbolic_polytope"]["facets"],
            case["hyperbolic_polytope"]["finite_vertices"],
            case["hyperbolic_polytope"]["cusps"],
            case["hyperbolic_polytope"]["elliptic_edges"],
            case["hyperbolic_polytope"]["bad_edges"],
        )
        for case in CERTIFICATE["cases"]
    ] == [(11, 16, 112, 9, 750, 0), (10, 14, 44, 3, 254, 0)]
    for case in CERTIFICATE["cases"]:
        assert case["hyperbolic_polytope"]["classification"] == (
            "finite_volume_noncompact"
        )
        assert case["scope"]["new_positive_surface_realization_in_this_project"]
        assert not case["scope"]["new_abstract_coxeter_diagram_claimed"]
        assert not case["scope"]["worldwide_surface_realization_priority_claimed"]
