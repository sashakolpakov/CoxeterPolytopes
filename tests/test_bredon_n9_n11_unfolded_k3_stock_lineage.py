from scripts.certify_bredon_n9_n11_unfolded_k3_stock_lineage import analyze


def test_n9_exact_complete_stock_nonmatch() -> None:
    result = analyze()["n9"]
    assert result["bredon_row"] == "ADEG-N9-01-branched-path"
    assert result["unfolded_rank_root_count"] == [12, 13]
    assert [row["key"] for row in result["classified_stock_candidates_with_same_rank_and_wall_count"]] == [
        "R12_UA2E8"
    ]
    assert result["exact_weighted_gram_isomorphism_count"] == 0
    assert result["target_graph_invariants"] == {
        "vertices": 13,
        "edges": 12,
        "connected_components": 1,
        "cycle_rank": 0,
        "triangle_count": 0,
        "degree_histogram": {"1": 5, "2": 6, "3": 1, "4": 1},
        "weighted_graph_automorphism_order": 12,
    }
    assert result["stock_graph_invariants"] == {
        "vertices": 13,
        "edges": 13,
        "connected_components": 1,
        "cycle_rank": 1,
        "triangle_count": 1,
        "degree_histogram": {"1": 2, "2": 9, "3": 2},
        "weighted_graph_automorphism_order": 2,
    }
    assert result["target_root_lattice_determinant"] == -12
    assert result["stock_full_NS_determinant"] == -3
    assert not result["integral_lattice_isometry"]
    assert not result["triality_compatible_stock_automorphism"]


def test_n11_exact_graph_lattice_fold_and_surface_match() -> None:
    result = analyze()["n11"]
    assert result["bredon_row"] == "ADEG-N11-01-branched-path"
    assert result["unfolded_rank_root_count"] == [14, 15]
    assert [row["key"] for row in result["classified_stock_candidates_with_same_rank_and_wall_count"]] == [
        "R14_UE8D4"
    ]
    assert result["exact_weighted_gram_isomorphism_count"] == 6
    assert result["canonical_target_to_source_wall_permutation_one_based"] == [
        2,
        3,
        5,
        6,
        7,
        4,
        8,
        10,
        9,
        1,
        11,
        12,
        13,
        14,
        15,
    ]
    assert result["long_target_wall_unfolds_to_target_roots_one_based"] == [13, 14, 15]
    assert result["their_source_images_one_based"] == [13, 14, 15]
    assert result["source_c3_orbit_one_based"] == [13, 14, 15]
    assert result["triality_orbits_respected"]
    assert result["folded_integral_gram_entrywise_equal"]
    assert result["folded_normalized_relation_matrix_entrywise_equal"]
    assert abs(result["unfolded_basis_change_determinant"]) == 1
    assert result["unfolded_integral_lattice_isometry"]
    assert result["c3_integral_action_entrywise_intertwined"]
    assert abs(result["fixed_lattice_basis_change_determinant"]) == 1
    assert result["fixed_lattice_determinant"] == -3
    assert result["fixed_lattice"] == "U+E8(-1)+A2(-1)"
    assert result["equation_level_surface"]["global_order"] == 6
    assert result["equation_level_surface"]["order_on_NS"] == 3
    assert result["same_polytope_vinberg_certificate"] == {
        "abstract_polytope_novelty": "known; this is a new explicit realization in the project",
        "bad_edges": 0,
        "classification": "finite_volume_noncompact",
        "cusps": 1,
        "dimension": 11,
        "elliptic_edges": 162,
        "facets": 13,
        "finite_vertices": 27,
    }
    assert result["known_lineage"]["reflection_subgroup_index"] == 1


def test_aggregate_lineage_verdict() -> None:
    result = analyze()
    assert result["summary"] == {
        "exact_stock_matches": 1,
        "exact_stock_failures": 1,
        "new_abstract_polytope_count": 0,
        "new_named_same_P_polytope_to_surface_identifications": 1,
    }
    assert len(result["canonical_sha256"]) == 64
