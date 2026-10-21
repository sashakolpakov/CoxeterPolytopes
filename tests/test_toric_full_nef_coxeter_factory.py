import json

import sympy as sp

from scripts.certify_toric_full_nef_coxeter_factory import (
    ENDPOINT_CANONICAL,
    ENDPOINT_INTERIOR,
    ENDPOINT_WALL_ORDER,
    ENDPOINT_WALLS,
    J,
    OUTPUT,
    build_certificate,
    canonical_cycle,
    endpoint_anticanonical_lattice_points,
    endpoint_congruence_audit,
    endpoint_cusp_audits,
    endpoint_extreme_rays,
    endpoint_fan_self_intersections,
    endpoint_gram,
    endpoint_solomon_euler,
    endpoint_wall_matrix,
    enumerate_weak_fano_cycles,
    exact_inertia,
    pair,
    p6_clique_euler,
    reflection_matrix,
    square,
    verify_archive,
)


def test_exact_toric_weak_fano_fan_classification() -> None:
    cycles = enumerate_weak_fano_cycles()
    assert len(cycles) == 16
    assert {
        length: sum(len(cycle) == length for cycle in cycles)
        for length in range(3, 10)
    } == {3: 1, 4: 3, 5: 2, 6: 4, 7: 2, 8: 3, 9: 1}
    assert max(map(len, cycles)) == 9
    assert [cycle for cycle in cycles if len(cycle) == 9] == [
        canonical_cycle((-2, -1, -2) * 3)
    ]
    assert all(min(cycle) >= -2 for cycle in cycles)
    assert all(sum(cycle) == 12 - 3 * len(cycle) for cycle in cycles)


def test_h6_endpoint_is_full_primitive_lorentzian_wall_stock() -> None:
    walls = endpoint_wall_matrix()
    gram = endpoint_gram()
    assert walls.rank() == 7
    assert gram.rank() == 7
    assert exact_inertia(gram) == (1, 6, 2)
    assert sorted(square(ENDPOINT_WALLS[name]) for name in ENDPOINT_WALL_ORDER) == [
        -2,
        -2,
        -2,
        -2,
        -2,
        -2,
        -1,
        -1,
        -1,
    ]
    assert square(ENDPOINT_CANONICAL) == 3
    assert all(
        pair(ENDPOINT_CANONICAL, ENDPOINT_WALLS[name])
        == -2 - square(ENDPOINT_WALLS[name])
        for name in ENDPOINT_WALL_ORDER
    )
    assert square(ENDPOINT_INTERIOR) > 0
    assert all(
        pair(ENDPOINT_INTERIOR, ENDPOINT_WALLS[name]) > 0
        for name in ENDPOINT_WALL_ORDER
    )


def test_h6_endpoint_angles_and_reflections_are_integral_coxeter() -> None:
    payload = build_certificate()["h6_endpoint"]
    assert payload["coxeter_label_histogram"] == {
        "2": 27,
        "3": 3,
        "4": 6,
        "infinity": 0,
    }
    for name in ENDPOINT_WALL_ORDER:
        reflection = reflection_matrix(ENDPOINT_WALLS[name])
        assert all(value.q == 1 for value in reflection)
        assert reflection.det() == -1
        assert reflection * reflection == sp.eye(7)
        assert reflection.T * J * reflection == J


def test_h6_endpoint_complete_exact_ray_and_facet_audit() -> None:
    rays = endpoint_extreme_rays()
    assert len(rays) == 21
    assert sum(ray["kind"] == "finite" for ray in rays) == 18
    assert sum(ray["kind"] == "ideal" for ray in rays) == 3
    assert all(ray["square"] >= 0 for ray in rays)
    assert {
        tuple(ray["vector"])
        for ray in rays
        if ray["kind"] == "ideal"
    } == {
        (1, -1, 0, 0, 0, 0, 0),
        (1, 0, -1, 0, 0, 0, 0),
        (1, 0, 0, -1, 0, 0, 0),
    }
    endpoint = build_certificate()["h6_endpoint"]
    assert endpoint["active_subset_trials"] == 84
    assert endpoint["decision"] == "finite_volume_noncompact"
    assert all(endpoint["checks"].values())


def test_h6_endpoint_fan_cubic_model_and_literal_cusps() -> None:
    assert endpoint_fan_self_intersections() == (-2, -1, -2) * 3
    assert endpoint_anticanonical_lattice_points() == (
        (-1, -1),
        (0, -1),
        (0, 0),
        (1, 2),
    )
    cusps = endpoint_cusp_audits()
    assert len(cusps) == 3
    assert all(cusp["active_facet_count"] == 7 for cusp in cusps)
    assert all(cusp["active_gram_inertia"] == [0, 5, 2] for cusp in cusps)
    assert all(cusp["component_rank_sum"] == 5 for cusp in cusps)
    assert all(cusp["maximal_parabolic_h6_cusp"] for cusp in cusps)


def test_h6_endpoint_a2_cubed_congruence_and_p6_gluing() -> None:
    congruence = endpoint_congruence_audit()
    assert congruence["minus_two_diagram_is_three_disjoint_A2_edges"]
    assert congruence["finite_group_order"] == 216
    assert congruence["mod_two_image_order"] == 216
    assert congruence["minus_one_orbit_sizes"] == [9, 9, 9]
    assert congruence["minus_one_orbits_are_disjoint"]
    assert congruence["orbit_union_count"] == 27
    assert congruence[
        "orbit_union_is_standard_6_plus_15_plus_6_cubic_lines"
    ]
    assert endpoint_solomon_euler() == {
        "subset_count": 512,
        "spherical_subset_count": 322,
        "euler_characteristic": "-1/1728",
        "gauss_bonnet_factor_absolute": "8*pi^3/15",
        "orbifold_volume": "pi^3/3240",
    }
    assert p6_clique_euler() == {
        "clique_counts_by_size": {
            "0": 1,
            "1": 27,
            "2": 216,
            "3": 720,
            "4": 1080,
            "5": 648,
            "6": 72,
        },
        "euler_characteristic": "-1/8",
        "gauss_bonnet_factor_absolute": "8*pi^3/15",
        "orbifold_volume": "pi^3/15",
    }


def test_frozen_toric_archive_replays_without_drift() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archived == build_certificate()
    assert archived["certificate_sha256"] == (
        "bd3cd228d883510444e2aa7948b859f76984498048e4ddfdcb848c5ab445f027"
    )
    assert verify_archive() == {
        "certificate_sha256": archived["certificate_sha256"],
        "cycle_count": 16,
        "maximum_hyperbolic_dimension": 6,
        "facet_count": 9,
        "finite_vertex_count": 18,
        "ideal_vertex_count": 3,
        "decision": "finite_volume_noncompact",
    }
