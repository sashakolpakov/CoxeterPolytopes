from scripts.verify_k3_quotient_nef_slices import analyze


def test_rank18_quotient_nef_slice_is_a_finite_volume_h9_simplex() -> None:
    report = analyze()
    assert report["all_checks_pass"] is True
    case = next(case for case in report["cases"] if case["ambient_rho"] == 18)
    assert case["fixed_rank"] == 10
    assert case["hyperbolic_dimension"] == 9
    assert len(case["facet_orbits"]) == 10
    assert case["coxeter_edge_counts"] == {"3": 8, "4": 1, "infinity": 0}
    assert case["finite_vertex_count"] == 8
    assert case["ideal_vertex_count"] == 2
    assert case["finite_volume"] is True
    assert case["compact"] is False
    literal = case["literal_coxiter_prop_1_10"]
    assert literal["spherical_rank_n_minus_one_count"] == 45
    assert literal["spherical_rank_n_count"] == 8
    assert literal["parabolic_rank_n_minus_one_count"] == 2
    assert literal["extension_histogram"] == {"2": 45}
    assert literal["endpoint_type_histogram"] == {
        "2,0": 28,
        "1,1": 16,
        "0,2": 1,
    }
    assert literal["has_required_vertex"] is True
    assert literal["bad_edges"] == []
    assert literal["finite_covolume"] is True
    assert literal["matches_optimized_endpoint_sets"] is True
    assert case["vinberg_volume_criterion"] == {
        "simple_roots_complete": True,
        "decision": "finite_volume_noncompact",
        "elliptic_edge_count": 45,
        "finite_vertex_count": 8,
        "ideal_vertex_count": 2,
        "finite_volume_bad_edges": [],
        "ideal_vertices": case["vinberg_volume_criterion"]["ideal_vertices"],
        "parabolic_neighborhood_audit_count": 2,
    }


def test_rank19_quotient_nef_slice_is_a_finite_volume_h10_polytope() -> None:
    report = analyze()
    case = next(case for case in report["cases"] if case["ambient_rho"] == 19)
    assert case["fixed_rank"] == 11
    assert case["hyperbolic_dimension"] == 10
    assert len(case["facet_orbits"]) == 13
    assert case["redundant_timelike_orbits"] == [[23, 24]]
    assert case["coxeter_edge_counts"] == {"3": 9, "4": 2, "infinity": 2}
    assert case["finite_vertex_count"] == 34
    assert case["ideal_vertex_count"] == 4
    assert case["finite_volume"] is True
    assert case["compact"] is False
    criterion = case["vinberg_volume_criterion"]
    assert criterion["decision"] == "finite_volume_noncompact"
    assert criterion["elliptic_edge_count"] == 208
    assert criterion["finite_vertex_count"] == 34
    assert criterion["ideal_vertex_count"] == 4
    assert criterion["finite_volume_bad_edges"] == []
    # Three of the four cusps are non-simple.  Their full parabolic
    # neighbourhoods have two affine components and must each count as one
    # endpoint, rather than once for each component or newly added node.
    assert sorted(
        len(vertex["affine_components"])
        for vertex in criterion["ideal_vertices"]
    ) == [1, 2, 2, 2]
    assert criterion["parabolic_neighborhood_audit_count"] == 7
    literal = case["literal_coxiter_prop_1_10"]
    assert literal["spherical_rank_n_minus_one_count"] == 208
    assert literal["spherical_rank_n_count"] == 34
    assert literal["parabolic_rank_n_minus_one_count"] == 4
    assert literal["extension_histogram"] == {"2": 208}
    assert literal["endpoint_type_histogram"] == {
        "2,0": 137,
        "1,1": 66,
        "0,2": 5,
    }
    assert literal["has_required_vertex"] is True
    assert literal["bad_edges"] == []
    assert literal["finite_covolume"] is True
    assert literal["matches_optimized_endpoint_sets"] is True


def test_every_quotient_cusp_has_an_exact_parabolic_link() -> None:
    report = analyze()
    assert [
        len(case["ideal_link_checks"]) for case in report["cases"]
    ] == [2, 4]
    assert all(
        link["negative_semidefinite_of_rank_dimension_minus_one"]
        for case in report["cases"]
        for link in case["ideal_link_checks"]
    )
