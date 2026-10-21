from coxeter_polytopes.exact.cubic7 import THETA7
from coxeter_polytopes.search.cubic7_a24_plus_unit import unit_a24_plus_first_packet


def test_cubic7_unit_a24_niemeier_affine_alcove_geometry() -> None:
    audit = unit_a24_plus_first_packet()
    assert audit.complete_control_height_bound == 100
    assert audit.control_level_count == 21
    assert audit.first_control == 2 + 5 * THETA7 + 2 * THETA7**2
    assert audit.first_height == 5 + 14 * THETA7 + 6 * THETA7**2
    assert audit.root_norm == 1
    assert audit.spatial_coset == "A24"
    assert audit.reflection_integral
    assert audit.cap == 1 + THETA7 - THETA7**2
    assert len(audit.rays) == 25
    assert audit.finite_vertex_count == 1
    assert audit.cusp_count == 0
    assert audit.hyperideal_ray_count == 24
    assert not audit.finite_volume and not audit.compact
