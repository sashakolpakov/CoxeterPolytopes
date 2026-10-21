from __future__ import annotations

from functools import lru_cache
import json

from scripts.audit_alexeev_parked_h9_doubling_products_hostile_audit import (
    OUTPUT,
    audit,
)


@lru_cache(maxsize=1)
def _rebuilt():
    return audit()


def test_hostile_product_audit_archive_rebuilds_exactly() -> None:
    payload = _rebuilt()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert payload["independence"] == {
        "imports_product_certifier": False,
        "imports_product_scout": False,
        "calls_product_build": False,
        "all_integer_translation_exponents_solved_symbolically": True,
        "fresh_GMP_double_description_for_six_rows": True,
        "independent_ray_incidence_time_facet_and_Vinberg_replay": True,
    }


def test_all_parameter_formula_gates_are_closed() -> None:
    payload = _rebuilt()
    paths = payload["infinite_path_cube_families"]
    assert paths["A14_8_1_C02"]["component_formula"]["simplified"] == "4*N+13"
    assert paths["A14_8_1_C02"]["component_formula"]["all_N_at_least_one"]
    assert paths["A14_10_1_C02"]["component_formula"]["simplified"] == "3*N+14"
    assert paths["A14_10_1_C02"]["component_formula"]["all_N_at_least_one"]
    assert all(
        family["global_wall_orbit_audit"][
            "stabilizers_exactly_generated_by_orthogonal_selected_reflections"
        ]
        and family["global_wall_orbit_audit"][
            "cross_seed_wall_orbit_identification_count"
        ]
        == 0
        and family["global_wall_orbit_audit"][
            "opposite_orientation_identification_count"
        ]
        == 0
        for family in paths.values()
    )
    rectangle = payload["infinite_rectangle_family"]
    assert rectangle["direct_product"] == {
        "group": "D_infinity x D_infinity",
        "injective": True,
        "proof": (
            "the y factor fixes a,b; an x normal form fixing a,b is the identity; "
            "then the y normal form is the identity"
        ),
    }
    assert rectangle["all_parameter_component_count"]["facet_formula"] == "m*n+m+n+12"
    assert rectangle["corner_orbit_separation"][
        "all_integer_k_ell_covered_symbolically"
    ]
    assert len(rectangle["corner_orbits"]) == 4


def test_six_direct_replays_and_bounded_corpus_gate_are_green() -> None:
    payload = _rebuilt()
    assert len(payload["direct_exact_replays"]) == 6
    assert all(
        row["all_exact_checks_pass"]
        and row["independent_group_boundary_reconstruction"]
        and row["hyperideal_vertex_count"] == 0
        and row["independent_strict_timelike_facet_witness_count"]
        == row["facet_count"]
        for row in payload["direct_exact_replays"]
    )
    assert payload["coxiter_bounded_comparison"][
        "total_exact_encoded_labelled_graph_match_count"
    ] == 0
    theorem = payload["theorem_status_split"]
    assert theorem["six_explicit_descendants"]["count"] == 6
    assert theorem["finite_C2_cubed_census"]["marked_subset_count"] == 7
    assert theorem["finite_C2_cubed_census"]["labelled_type_count"] == 3
    assert theorem["three_all_parameter_families"]["status"] == "universal_closed_orbit_proofs"
