from coxeter_polytopes.exact.cubic7 import THETA7, Cubic7
from coxeter_polytopes.search.cubic7_a1_golay import (
    A_UNIT,
    UNIT_A1_GOLAY_SEED,
    height_zero_patterns,
    unit_prefix_through_h100,
)
from coxeter_polytopes.search.cubic7_timelike import unit_signature_rank


def test_cubic7_unit_a1_golay_complete_h100_prefix() -> None:
    assert [A_UNIT.sign(place) for place in range(3)] == [1, -1, -1]
    assert A_UNIT.norm() == 1
    assert unit_signature_rank() == 3
    assert UNIT_A1_GOLAY_SEED.canonical_mirror_norms == (Cubic7(1), Cubic7(2))
    assert len(height_zero_patterns(UNIT_A1_GOLAY_SEED)) == 1

    audit = unit_prefix_through_h100()
    assert audit.control_level_count == 21
    assert audit.raw_square_partition_count == 1329
    assert audit.nonempty_level_count == 2

    assert audit.first_level.norm == 1
    assert audit.first_level.control == 2 + 5 * THETA7 + 2 * THETA7**2
    assert audit.first_level.height == 5 + 14 * THETA7 + 6 * THETA7**2
    assert audit.first_pattern.labeled_count == 24
    assert audit.first_geometry.extreme_ray_count == 2**24
    assert audit.first_geometry.finite_vertex_count == 25
    assert audit.first_geometry.cusp_count == 0
    assert audit.first_geometry.hyperideal_ray_count == 2**24 - 25

    assert audit.second_level.norm == 2
    assert audit.second_level.control == 2 + 8 * THETA7 + 4 * THETA7**2
    assert audit.second_level.height == 6 + 18 * THETA7 + 8 * THETA7**2
    assert audit.second_pattern.labeled_count == 24 * 23
    assert all(pairing.sign(0) <= 0 for pairing in audit.first_to_second_pairings)
    assert all(pairing.sign(0) <= 0 for pairing in audit.second_same_height_pairings)

    geometry = audit.second_geometry
    assert geometry.total_wall_count == 24 + 24 + 24 * 23
    assert geometry.extreme_ray_count == 218_103_784
    assert geometry.finite_vertex_count == 301
    assert geometry.cusp_count == 0
    assert geometry.hyperideal_ray_count == 218_103_483
    assert geometry.uniform_coordinate == THETA7 - 1
    assert geometry.uniform_norm == 49 - 97 * THETA7 + 48 * THETA7**2
    assert geometry.active_wall_count == 552
    assert geometry.active_spatial_rank == 24
    assert not geometry.finite_volume and not geometry.compact
