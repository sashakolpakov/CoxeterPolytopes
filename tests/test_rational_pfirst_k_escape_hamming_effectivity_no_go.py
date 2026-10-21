from scripts.certify_rational_pfirst_k_escape_hamming_effectivity_no_go import (
    payload,
)


def test_exact_hamming_k_escape_effectivity_obstruction() -> None:
    data = payload()
    assert data["rank"] == 10
    assert data["canonical_square"] == 0
    assert data["wall_count"] == 32
    assert data["wall_kind_histogram"] == {"minus_one": 16, "minus_two": 16}
    assert data["coxeter_order_histogram"] == {
        "2": 352,
        "4": 128,
        "infinity": 16,
    }
    assert data["unimodular_wall_subset_determinant"] == -1
    assert data["canonical_escape"]["K_is_not_in_P"] is True
    assert data["polyhedron"]["extreme_ray_count"] == 2306
    assert data["polyhedron"]["primitive_ray_square_histogram"] == {
        "0": 226,
        "1": 1040,
        "2": 1040,
    }
    assert data["polyhedron"]["active_histogram"] == {
        "finite_9": 2080,
        "ideal_12": 224,
        "ideal_16": 2,
    }
    assert data["polyhedron"]["facet_active_ray_counts"] == [670] * 32
    assert data["polyhedron"]["facet_active_span_ranks"] == [9] * 32
    assert (
        data["polyhedron"]["sorted_primitive_ray_sha256"]
        == "33ad7a1bb8df2dda4c4fd82b7cea62577f9f64d1f9972222707923104e248b09"
    )
    assert data["simultaneous_effectivity_obstruction"][
        "minimum_euler_contribution_of_eight_fibers"
    ] == 16
    assert data["simultaneous_effectivity_obstruction"]["rational_surface_c2"] == 12
    assert data["explicit_missing_reflection"]["vector"] == [
        1,
        0,
        0,
        0,
        0,
        -1,
        0,
        -1,
        0,
        -1,
    ]
    assert data["explicit_missing_reflection"]["square"] == -2
    assert data["explicit_missing_reflection"]["strict_interior_pairing"] == 0
    assert data["scope"]["relative_reflection_maximality_of_this_packet"] == "DISPROVED"
    assert (
        data["scope"]["maximal_plus_K_escape_nonrealizability_counterexample"]
        == "OPEN"
    )
