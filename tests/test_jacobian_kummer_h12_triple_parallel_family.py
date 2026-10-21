from __future__ import annotations

from functools import lru_cache
import json

from scripts.certify_jacobian_kummer_h12_triple_parallel_family import (
    OUTPUT,
    build,
    facet_formula,
)


@lru_cache(maxsize=1)
def _payload() -> dict:
    return build()


def test_archive_rebuilds_exactly() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["independence"]["exploratory_scout_imported_or_invoked"] is False
    assert payload["independence"][
        "production_path_cube_builder_imported_or_invoked"
    ] is False
    assert payload["independence"]["finite_parameter_interpolation_used"] is False


def test_exact_triple_parallel_structures_and_scoped_maxima() -> None:
    payload = _payload()
    graph = payload["scoped_maximum_product_theorem"]
    assert graph["doubling_wall_count"] == 124
    assert graph["doubling_relation_histogram"] == {
        "m2": 264,
        "parallel": 1138,
        "ultraparallel": 6224,
    }
    assert graph["parallel_pair_count"] == 1138
    assert graph["parallel_factor_compatibility_edge_count"] == 6
    assert graph["maximum_D_infinity_factor_count"] == 3
    assert graph["D_infinity_cubed_configurations_one_based"] == [
        [[48, 42], [51, 39], [56, 36]],
        [[58, 33], [59, 32], [60, 31]],
    ]
    assert graph["triple_common_orthogonal_doubling_walls"] == [[], []]
    assert graph["D_infinity_squared_common_cube_rank_histogram"] == {"1": 6}
    assert graph["D_infinity_common_cube_rank_histogram"] == {
        "0": 452,
        "1": 240,
        "2": 438,
        "3": 8,
    }
    assert graph["maximum_orthogonal_cube_rank"] == 4
    assert graph["maximum_C2_four_cube_count"] == 16
    assert all(graph["checks"].values())

    structures = payload["triple_parallel_structures"]["configurations"]
    assert len(structures) == 2
    expected_gram = [
        [-4, 4, 0, 0, 0, 0],
        [4, -4, 0, 0, 0, 0],
        [0, 0, -4, 4, 0, 0],
        [0, 0, 4, -4, 0, 0],
        [0, 0, 0, 0, -4, 4],
        [0, 0, 0, 0, 4, -4],
    ]
    assert all(row["selected_wall_Gram"] == expected_gram for row in structures)
    assert all(all(row["checks"].values()) for row in structures)


def test_all_Z3_collision_and_endpoint_proofs() -> None:
    payload = _payload()
    collision = payload["all_integer_Z3_collision_audit"]
    assert collision["bounded_parameter_scan_used"] is False
    assert collision["ordered_signed_equation_count"] == 440_896
    assert collision["moving_branch_source_count"] == 1232
    assert collision["quotient_affine_rank_histogram"] == {"4": 1232}
    assert collision["finite_solution_count"] == 462
    assert collision["finite_solution_branch_histogram"] == {
        "000": 154,
        "001": 60,
        "010": 60,
        "011": 36,
        "100": 60,
        "101": 36,
        "110": 36,
        "111": 20,
    }
    assert collision["all_Z3_identity_count"] == 96
    assert collision["negative_root_collision_count"] == 0
    assert collision["cross_facet_collision_count"] == 0
    assert all(collision["checks"].values())

    endpoints = payload["endpoint_orbit_audit"]
    assert endpoints["outer_endpoint_count_for_every_finite_box"] == 6
    assert endpoints["cross_factor_endpoint_quotient_collision_count"] == 0
    assert all(endpoints["checks"].values())


def test_phase_refined_formula_and_direct_reconstructions() -> None:
    payload = _payload()
    theorem = payload["theorem"]
    assert theorem["selected_subgroup"] == "D_infinity^3"
    assert theorem["seed_chamber_count_and_reflection_subgroup_index"] == "N1*N2*N3"
    assert theorem["volume"] == "Vol(P_N1_N2_N3)=N1*N2*N3*Vol(P_seed)"
    assert theorem["finite_volume"] and not theorem["compact"]

    assert facet_formula(1, 1, 1, 0) == 172
    assert facet_formula(2, 1, 1, 0) == 296
    assert facet_formula(2, 2, 1, 0) == 523
    assert facet_formula(2, 2, 2, 0) == 941
    assert facet_formula(2, 2, 2, 1) == 942
    assert facet_formula(3, 3, 3, 0) == facet_formula(3, 3, 3, 1) == 2860
    # Regression gate: the preliminary phase-blind formula is false.
    assert facet_formula(2, 2, 2, 0) != facet_formula(2, 2, 2, 1)

    formula = payload["formula_audit"]
    assert len(formula["symbolic_parity_identities"]) == 64
    assert all(formula["checks"].values())
    direct = payload["direct_grid_reconstructions"]
    assert len(direct["rows"]) == 30
    assert sum(len(row["phase_reconstructions"]) for row in direct["rows"]) == 240
    row_222 = next(row for row in direct["rows"] if row["parameters"] == [2, 2, 2])
    assert row_222["even_phase_expected_facet_count"] == 941
    assert row_222["odd_phase_expected_facet_count"] == 942
    assert {
        phase["facet_count"] for phase in row_222["phase_reconstructions"]
    } == {941, 942}
    assert all(direct["checks"].values())


def test_symmetry_full_carrier_volume_and_lineage_scope() -> None:
    payload = _payload()
    symmetry = payload["symmetry_orbit_audit"]
    assert symmetry["published_chamber_group_order"] == 23_040
    assert symmetry["involution_centralizer_order"] == 768
    assert symmetry["configuration_orbit_size"] == 2
    assert symmetry["maps_from_representative_to_each_configuration"] == {
        "1": 384,
        "2": 384,
    }
    assert symmetry["configuration_stabilizer_order"] == 384
    assert len(symmetry["labelled_factor_endpoint_actions"]) == 24
    assert all(symmetry["checks"].values())

    lattice = payload["full_lattice_core_audit"]
    assert len(lattice["phase_cores"]) == 8
    assert all(row["Smith_diagonal"] == [1] * 13 for row in lattice["phase_cores"])
    carrier = payload["theorem"]["carrier"]
    assert carrier["carrier_lattice"] == "U(2)+D5(-1)^2+<-4>"
    assert carrier["carrier_rank"] == 13
    assert carrier["carrier_moduli_dimension"] == 7
    assert all(carrier["checks"].values())

    lineage = payload["bounded_coxiter_corpus_replay"]
    assert lineage["ordinary_file_count"] == 1214
    assert lineage["index_two_file_count"] == 2959
    assert lineage["literal_seed_or_descendant_headers_absent"] is True
    assert lineage["worldwide_priority"] == "open"
    assert "not a publication-complete denominator" in lineage["scope_warning"]

    boundary = payload["claim_boundary"]
    assert boundary["new_project_family"] is True
    assert boundary["worldwide_priority"] == "open"
    assert boundary["weak_full_Picard_divisor_carrier"] is True
    assert boundary["full_nef_claimed"] is False
    assert boundary["surface_automorphism_claimed"] is False
    assert boundary["scoped_product_maximum_not_arbitrary_redoubling_classification"]
