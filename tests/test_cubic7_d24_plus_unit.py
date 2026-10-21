from coxeter_polytopes.exact.cubic7 import THETA7
from coxeter_polytopes.search.cubic7_d24_plus_unit import (
    analyze_unit_d24_plus_prefix,
    run_unit_d24_plus_prefix,
    unit_d24_plus_first_packet,
)


def test_cubic7_unit_d24_plus_affine_alcove_geometry() -> None:
    audit = unit_d24_plus_first_packet()
    assert audit.complete_control_height_bound == 100
    assert audit.control_level_count == 21
    assert audit.first_control == 2 + 5 * THETA7 + 2 * THETA7**2
    assert audit.first_height == 5 + 14 * THETA7 + 6 * THETA7**2
    assert audit.root_norm == 1
    assert audit.spatial_coset == "D"
    assert audit.glue_basis_integral and audit.reflection_integral
    assert audit.affine_marks_verified
    assert audit.cap == 1 + THETA7 - THETA7**2
    assert len(audit.rays) == 25
    assert audit.finite_vertex_count == 2
    assert audit.cusp_count == 0
    assert audit.hyperideal_ray_count == 23
    assert not audit.finite_volume and not audit.compact


def test_cubic7_unit_d24_plus_glue_complete_h100_prefix() -> None:
    search = run_unit_d24_plus_prefix(100)
    assert search.control_level_count == 21
    assert len(search.accepted_roots) == 4
    nonempty = [packet for packet in search.packets if packet.candidate_count]
    assert [(packet.candidate_count, len(packet.accepted_roots)) for packet in nonempty] == [
        (1, 1),
        (3, 3),
    ]
    assert all(root.coset == "D" for root in search.accepted_roots)
    assert [root.norm for root in search.accepted_roots] == [1, 2, 2, 2]

    geometry = analyze_unit_d24_plus_prefix(search)
    assert [
        (
            packet.wall_count,
            packet.extreme_ray_count,
            packet.finite_vertex_count,
            packet.cusp_count,
            packet.hyperideal_ray_count,
        )
        for packet in geometry.packet_geometries
    ] == [
        (25, 25, 2, 0, 23),
        (28, 134, 4, 0, 130),
    ]
    assert [
        (
            profile.negative_ray_count,
            profile.zero_ray_count,
            profile.positive_ray_count,
            profile.output_ray_count,
        )
        for profile in geometry.insertion_profiles
    ] == [(3, 0, 22, 69), (23, 0, 46, 111), (110, 0, 1, 134)]
    assert geometry.finite_vertex_count == 4
    assert geometry.cusp_count == 0
    assert geometry.hyperideal_ray_count == 130
    assert not geometry.finite_volume and not geometry.compact
