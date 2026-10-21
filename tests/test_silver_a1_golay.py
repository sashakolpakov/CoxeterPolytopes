from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_a1_golay import (
    EPSILON,
    SQRT2_A1_GOLAY_SEED,
    TWO_SQRT2_A1_GOLAY_SEED,
    UNIT_A1_GOLAY_SEED,
    enumerate_control_levels,
    enumerate_level_patterns,
    first_packet_cube_geometry,
    height_zero_patterns,
    unit_dense_prefix_through_h20,
    unit_symmetric_slope_upper_bound,
)


def test_complete_height_zero_and_first_cube_packets() -> None:
    expected = {
        UNIT_A1_GOLAY_SEED: (2 + 2 * SQRT2, Silver(1), 47 - SQRT2),
        SQRT2_A1_GOLAY_SEED: (1 + SQRT2, SQRT2 / 2, 24 - SQRT2),
        TWO_SQRT2_A1_GOLAY_SEED: (1 + SQRT2, Silver(1), 48 - 2 * SQRT2),
    }
    for seed, (height, cap, uniform_norm) in expected.items():
        height_zero = height_zero_patterns(seed)
        assert len(height_zero) == 1
        assert height_zero[0].labeled_count == 24
        geometry = first_packet_cube_geometry(seed)
        assert geometry.height == height
        assert geometry.cap_coordinate == cap
        assert geometry.accepted_wall_count == 24
        assert geometry.total_wall_count == 48
        assert geometry.extreme_ray_count == 2**24
        assert (
            geometry.finite_vertex_count,
            geometry.cusp_count,
            geometry.hyperideal_ray_count,
        ) == (25, 0, 2**24 - 25)
        assert geometry.uniform_ray_norm == uniform_norm
        assert not geometry.finite_volume and not geometry.compact


def test_unit_first_two_complete_nonempty_packets() -> None:
    audit = unit_dense_prefix_through_h20()
    assert audit.complete_height_bound == 20
    assert audit.first_height == 2 + 2 * SQRT2
    assert audit.first_family_count == 24
    assert audit.second_height == 7 + 5 * SQRT2
    assert audit.second_control == 2 + SQRT2
    assert audit.second_pattern == (
        ((3 + SQRT2) / 2, 1),
        (EPSILON / 2, 7),
        (Silver(0), 16),
    )
    assert audit.second_family_count == 759 * 8 == 6072
    assert audit.second_family_rank == 24
    assert all(margin.sign(0) <= 0 for margin in audit.second_same_height_margins)
    assert sum(not margin for margin in audit.second_same_height_margins) == 1
    assert audit.symmetric_coordinate == Silver(2, 0) / 7 + SQRT2 / 14
    assert audit.symmetric_norm == Silver(167, 47) / 49
    assert audit.symmetric_extreme and audit.symmetric_hyperideal


def test_unit_h100_relaxed_slope_proves_prefix_hyperideal_witness() -> None:
    audit = unit_symmetric_slope_upper_bound(100)
    assert audit.relaxed_pattern_count == 210
    assert audit.maximum_slope == 7 - 3 * SQRT2 / 2
    assert audit.maximizing_height == 41 + 29 * SQRT2
    assert audit.maximizing_control == 4 + 3 * SQRT2
    assert audit.critical_slope_square == 12 * EPSILON
    assert audit.slope_square_below_critical
    assert audit.feasible_uniform_coordinate == Silver(10, 0) / 89 + 17 * SQRT2 / 178
    assert audit.feasible_uniform_norm == Silver(3815, 239) / 7921
    assert audit.hyperideal_prefix_witness


def test_two_sqrt2_nonreflective_glue_patterns_are_rejected() -> None:
    first_level = enumerate_control_levels(TWO_SQRT2_A1_GOLAY_SEED, 10)[0]
    # The raw square equation has pair/octet/single spatial patterns, but the
    # pi^3 norm class requires y in pi*N. Only the coordinate caps survive.
    patterns = enumerate_level_patterns(TWO_SQRT2_A1_GOLAY_SEED, first_level)
    assert len(patterns) == 1
    assert patterns[0].labeled_count == 24

