import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/enumerations/a3_h3_plane_arrangement_coble_factories.json"


def _load():
    return json.loads(DATA.read_text(encoding="utf-8"))


def test_h3_intersection_lattice_is_exact():
    row = _load()["arrangements"]["H3"]
    assert row["field"] == "Q(sqrt(5))"
    assert row["line_count"] == 15
    assert row["point_count"] == 31
    assert row["intersection_multiplicity_histogram"] == {"2": 15, "3": 10, "5": 6}
    assert row["pair_identity"] == row["expected_pair_count"] == 105
    assert row["projective_arrangement_automorphism_order"] == 60


def test_a3_intersection_lattice_is_exact():
    row = _load()["arrangements"]["A3"]
    assert row["field"] == "Q"
    assert row["line_count"] == 6
    assert row["point_count"] == 7
    assert row["intersection_multiplicity_histogram"] == {"2": 3, "3": 4}
    assert row["pair_identity"] == row["expected_pair_count"] == 15
    assert row["projective_arrangement_automorphism_order"] == 24


def test_a3_exhaustive_one_boundary_no_go():
    census = _load()["A3_hilbert_census"]
    assert census["scope"] == "all parent-closed A3 proper/first-near clusters of sizes 10..10"
    assert sum(census["rank_census_by_size"]["10"].values()) == 48_822
    assert census["strong_coble_hilbert_survivor_count_by_size"] == {"10": 9_925}
    orbits = census["survivor_orbits"]["10"]
    assert orbits["orbit_count"] == 442
    assert orbits["orbit_count_by_sextic_rank"] == {"24": 4, "25": 46, "26": 162, "27": 230}
    assert orbits["marked_count_by_sextic_rank"] == {
        "24": 25,
        "25": 948,
        "26": 3_678,
        "27": 5_274,
    }
    audit = census["unique_boundary_orbit_audit"]
    assert audit["orbit_count"] == 230
    assert audit["factor_partition_orbit_histogram"] == {
        "1+1+1+1+1+1": 160,
        "1+1+1+1+2": 35,
        "1+1+1+3": 16,
        "1+1+2+2": 2,
        "1+1+4": 13,
        "1+2+3": 1,
        "6": 3,
    }
    assert audit["irreducible_orbit_count"] == 3
    assert audit["local_resolution_pass_orbit_count"] == 0
    assert audit["nondegenerate_one_boundary_strong_coble_orbit_count"] == 0


def test_a3_full_six_mirror_surface_and_natural_cone():
    row = _load()["A3_full_six_mirror_surface"]
    assert row["field"] == "Q, characteristic zero"
    assert row["proper_center_count"] == 7
    assert row["first_near_center_count"] == 12
    assert row["blowup_count"] == 19
    assert row["picard_rank"] == 20
    assert row["hyperbolic_dimension"] == 19
    assert row["boundary_component_count"] == 10
    assert row["boundary_norms"] == [-4] * 10
    assert row["boundary_components_pairwise_disjoint"]
    assert row["obvious_exceptional_count"] == 15
    assert row["obvious_exceptional_norms"] == [-1] * 15
    assert row["strong_coble"]
    cone = row["declared_wall_cone"]
    assert cone["wall_count"] == 25
    assert cone["wall_rank"] == 20
    assert cone["extreme_ray_count"] == 3_986
    assert cone["finite_vertex_count"] == 2_865
    assert cone["ideal_vertex_count"] == 331
    assert cone["hyperideal_vertex_count"] == 790
    assert cone["all_facets_irredundant"]
    witness = cone["sample_hyperideal_rays"][0]
    assert witness["square"] < 0
    assert witness["active_facet_rank"] == 19
    assert witness["all_wall_pairings_nonnegative"]


def test_winger_hilbert_and_singularity_gates():
    row = _load()["H3_Winger_surface"]
    assert row["field"] == "Q(sqrt(5)), characteristic zero"
    assert row["no_three_centers_collinear"]
    assert row["cubic_rank"] == 10
    assert row["h0_minus_K"] == 0
    assert row["sextic_rank"] == 27
    assert row["h0_minus_2K"] == 1
    assert row["all_ten_prescribed_singularities_ordinary_nodes"]
    singular = row["singular_exact_audit"]
    assert singular["projective_jacobian_degree"] == 10
    assert singular["projective_radical_degree"] == 10
    assert singular["minimal_associated_prime_count"] == 10
    assert singular["nonunit_factor_count_over_Qsqrt5"] == 1
    assert row["absolute_irreducibility_argument"]["absolutely_irreducible"]
    assert row["surface_conclusion"]["strong_coble"]


def test_winger_lineage_is_not_overclaimed():
    lineage = _load()["H3_Winger_surface"]["lineage"]
    assert lineage["classical_name"] == "Winger sextic / A5 Coble surface"
    assert lineage["new_surface_claim"] is False


def test_natural_declared_wall_cones_are_infinite_volume():
    carriers = _load()["H3_winger_ten_node_lane"]["carriers"]
    expected = {
        "petersen": (25, 294, 1, 10, 283),
        "petersen_plus_boundary": (26, 173, 1, 10, 162),
        "all_pair_lines_plus_boundary": (56, 2103, 121, 220, 1762),
    }
    for name, values in expected.items():
        row = carriers[name]
        assert (
            row["wall_count"],
            row["extreme_ray_count"],
            row["finite_vertex_count"],
            row["ideal_vertex_count"],
            row["hyperideal_vertex_count"],
        ) == values
        assert row["hyperideal_vertex_count"] > 0


def test_literal_hyperideal_witnesses_have_full_active_rank():
    lane = _load()["H3_winger_ten_node_lane"]
    triples = lane["boundary_compatible_independent_triple_hyperideal_witnesses"]
    assert len(triples) == 30
    assert all(row["square"] == -2 for row in triples)
    assert all(row["boundary_pairing"] == 0 for row in triples)
    assert all(row["all_wall_pairings_nonnegative"] for row in triples)
    assert all(row["active_wall_rank"] == 10 for row in triples)

    degree_two = lane["degree_two_truncation_seven_subset_hyperideal_witnesses"]
    assert len(degree_two) == 120
    assert all(row["square"] == -3 for row in degree_two)
    assert all(row["primitive_boundary_pairing"] == 1 for row in degree_two)
    assert all(row["all_wall_pairings_nonnegative"] for row in degree_two)
    assert all(row["active_wall_rank"] == 10 for row in degree_two)


def test_six_mirror_strong_coble_factory_census():
    summary = _load()["six_mirror_boundary_orbits"]["H3_summary"]
    assert summary["orbit_count"] == 103
    assert summary["subset_count"] == 5005
    assert summary["admissible_orbit_count"] == 76
    assert summary["admissible_subset_count"] == 3595
    assert summary["admissible_subsets_by_triple_count"] == {
        "0": 20,
        "1": 260,
        "2": 1890,
        "3": 1300,
        "4": 125,
    }
    rows = _load()["six_mirror_boundary_orbits"]["H3"]
    t2 = [row for row in rows if row["ordinary_reduced_boundary_admissible"] and row["triple_count"] == 2]
    assert len(t2) == 39
    assert sum(row["orbit_size"] for row in t2) == 1890
    assert all(row["boundary_component_count"] == 8 for row in t2)
    assert all(row["blowup_count"] == row["hyperbolic_dimension"] == 17 for row in t2)
