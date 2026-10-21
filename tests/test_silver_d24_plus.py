from dataclasses import replace

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_d24_plus import (
    analyze_silver_d24_plus_cone,
    d24_cartan_gram,
    d24_height_zero_simple_roots,
    d24_plus_pair,
    run_silver_d24_plus_first_shell,
    run_silver_d24_plus_unit_prefix,
)


def test_height_zero_walls_are_actual_d24_simple_roots() -> None:
    roots = d24_height_zero_simple_roots()
    search = run_silver_d24_plus_first_shell(10)
    assert len(roots) == 24
    assert tuple(
        tuple(d24_plus_pair(search.seed, left, right) for right in roots)
        for left in roots
    ) == tuple(
        tuple(Silver(value) for value in row) for row in d24_cartan_gram()
    )
    assert all(d24_plus_pair(search.seed, root, root) == 2 for root in roots)


def test_corrected_first_packet_exact_cone_geometry() -> None:
    search = run_silver_d24_plus_first_shell(10)
    assert len(search.accepted_roots) == 2
    assert search.wall_count == len(search.walls) == 26
    assert all(
        d24_plus_pair(search.seed, root.vector, wall).sign(0) <= 0
        for root in search.accepted_roots
        for wall in search.walls[:24]
    )

    geometry = analyze_silver_d24_plus_cone(search)
    assert len(geometry.extreme_rays) == 48
    assert (
        geometry.finite_vertex_count,
        geometry.cusp_count,
        geometry.hyperideal_ray_count,
    ) == (5, 0, 43)
    assert geometry.irredundant
    assert not geometry.finite_volume and not geometry.compact
    assert len(geometry.insertion_profiles) == 1
    profile = geometry.insertion_profiles[0]
    assert (
        profile.negative_ray_count,
        profile.zero_ray_count,
        profile.positive_ray_count,
        profile.output_ray_count,
    ) == (2, 0, 23, 48)
    # The two spinor endpoints of the highest-root alcove survive and have
    # square 24-2*sqrt(2)>0.
    assert sum(
        ray.vector[0]
        and ray.norm / (ray.vector[0] * ray.vector[0]) == 24 - 2 * SQRT2
        for ray in geometry.extreme_rays
    ) == 2


def test_unit_coefficient_glue_complete_h30_prefix() -> None:
    search = run_silver_d24_plus_unit_prefix(30)
    assert len(search.control_levels) == 12
    assert len(search.packets) == 8
    assert [
        (packet.height, packet.candidate_count, len(packet.roots))
        for packet in search.packets
    ] == [
        (Silver(1, 1) / 2, 0, 0),
        (Silver(1, 1), 0, 0),
        (Silver(2, 2), 1, 1),
        (Silver(7, 5) / 2, 0, 0),
        (Silver(7, 5), 1, 1),
        (Silver(17, 13) / 2, 0, 0),
        (Silver(23, 17) / 2, 0, 0),
        (Silver(14, 10), 13, 11),
    ]
    assert len(search.accepted_roots) == 13
    assert search.wall_count == 37

    # The complete wall set is too large for the generic exact DD routine,
    # but the exact prefix through its first two positive walls is small.
    first_two = replace(search, walls=search.walls[:26], wall_count=26)
    geometry = analyze_silver_d24_plus_cone(first_two)
    assert len(geometry.extreme_rays) == 67
    assert (
        geometry.finite_vertex_count,
        geometry.cusp_count,
        geometry.hyperideal_ray_count,
    ) == (3, 0, 64)
    assert geometry.irredundant and not geometry.finite_volume
