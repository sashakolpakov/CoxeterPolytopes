from __future__ import annotations

import json

from scripts.certify_shimada_aut3singk3_x2_rational_involution_h10_37 import (
    OUTPUT,
    SCHEMA,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_x2_h10_37_archive_rebuilds_exactly() -> None:
    assert build_certificate() == _archive()


def test_complete_high_subgroup_denominator_has_only_two_C2_classes() -> None:
    payload = _archive()
    denominator = payload["complete_high_subgroup_denominator"]
    assert payload["schema"] == SCHEMA
    assert denominator["fixed_rank_histogram"] == {
        "3": 324,
        "4": 144,
        "6": 90,
        "8": 80,
        "11": 36,
        "12": 45,
        "20": 1,
    }
    assert denominator["high_involution_classes"] == [
        {
            "surface_type": "rational",
            "source_cardinality": 36,
            "fixed_rank": 11,
            "hyperbolic_dimension": 10,
            "representative_zero_based": 596,
        },
        {
            "surface_type": "symplectic",
            "source_cardinality": 45,
            "fixed_rank": 12,
            "hyperbolic_dimension": 11,
            "representative_zero_based": 470,
        },
    ]
    assert denominator["v4_subgroup_count"] == 120
    assert denominator["v4_common_fixed_rank_histogram"] == {"7": 90, "8": 30}
    assert denominator["v4_element_rank_pattern_histogram"] == {
        "11,11,12": 90,
        "12,12,12": 30,
    }
    assert all(denominator["checks"].values())


def test_rational_class_is_strict_integral_finite_volume_H10_37() -> None:
    payload = _archive()
    fold = payload["rational_involution_fold"]
    lattice = payload["saturated_fixed_divisor_lattice"]
    vinberg = payload["literal_vinberg_volume_certificate"]
    assert fold["wall_orbit_count"] == 242
    assert fold["positive_ray_deduplicated_restriction_count"] == 207
    assert fold["restriction_square_sign_histogram"] == {
        "negative": 37,
        "zero": 45,
        "positive": 125,
    }
    assert fold["strict_facet_count"] == 37
    assert len(fold["strict_facet_witnesses"]) == 37
    assert fold["pair_gate"]["coxeter_edge_counts"] == {
        "2": 405,
        "3": 15,
        "4": 70,
        "6": 0,
        "infinity": 121,
        "ultraparallel": 55,
        "forbidden": 0,
    }
    assert lattice["nikulin_invariants"] == {"r": 11, "a": 9, "delta": 1}
    assert lattice["classification"] == (
        "By uniqueness of indefinite even 2-elementary lattices with these "
        "invariants: U + E8(-2) + <-2>."
    )
    assert lattice["mirror_type_histogram"] == {
        "-2,div1": 6,
        "-2,div2": 21,
        "-4,div2": 10,
    }
    assert lattice["facet_mirror_coordinate_smith_diagonal"] == [1] * 11
    stable = lattice["exhaustive_reflection_stable_overlattice_audit"]
    assert stable["stable_integral_subgroup_orders"] == [1]
    assert stable["stable_even_subgroup_orders"] == [1]
    assert vinberg["elliptic_edge_count"] == 24350
    assert vinberg["finite_vertex_count"] == 3236
    assert vinberg["ideal_vertex_count"] == 261
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0
    assert vinberg["finite_volume"] is True
    assert vinberg["compact"] is False


def test_symplectic_class_fails_and_full_group_lineage_is_exact() -> None:
    payload = _archive()
    negative = payload["symplectic_involution_negative_screen"]
    lineage = payload["arithmetic_and_literal_lineage"]
    assert negative["fixed_rank"] == 12
    assert negative["strict_facet_count"] == 90
    assert negative["pair_gate"]["coxeter_edge_counts"]["forbidden"] == 352
    assert negative["pair_gate"]["all_pairs_coxeter"] is False
    assert lineage["solomon"] == {
        "candidate": {
            "connected_elliptic_special_subdiagram_count": 472,
            "rational_euler_characteristic": "-527/138240",
        },
        "full_Wr_chamber": {
            "connected_elliptic_special_subdiagram_count": 79,
            "rational_euler_characteristic": "-527/33443020800",
        },
        "index_candidate_in_full_Wr": 241920,
        "index_formula": "chi(W_candidate)/chi(W_r)",
    }
    assert lineage["novelty_status"] == {
        **lineage["novelty_status"],
        "new_arithmetic_lattice_class": False,
        "new_full_reflection_group": False,
        "new_abstract_polytope_worldwide": "open",
        "new_algebro_geometric_realization_candidate": True,
    }
    assert payload["all_checks_pass"] is True

