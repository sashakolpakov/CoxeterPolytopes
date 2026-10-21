from collections import Counter
from itertools import groupby

from coxeter_polytopes.exact.cubic7 import THETA7, Cubic7
from coxeter_polytopes.search.cubic7_bna1 import (
    DELTA7,
    L2_BNA1_SEED,
    L7_BNA1_SEED,
    analyze_bna1_cone_incremental,
    bna1_diagonal_pair,
    dominant_bna1_representative,
    enumerate_bna1_root_pool,
    enumerate_height_zero_patterns,
    height_zero_is_bn_minus_1_times_a1,
    run_bna1_vinberg_search,
)
from coxeter_polytopes.search.cubic7_nonunimodular import (
    L2_SEED,
    L7_SEED,
    enumerate_nonunimodular_root_pool,
)


def _stabilize_rank_six_vector(vector, target_dimension):
    # H5 has four ordinary coordinates.  New ordinary coordinates are zero
    # and are inserted immediately before the distinguished final slot.
    return (*vector[:-1], *(Cubic7(0) for _ in range(target_dimension - 5)), vector[-1])


def _zero_pad_vector(vector, source_dimension, target_dimension):
    return (
        *vector[:-1],
        *(Cubic7(0) for _ in range(target_dimension - source_dimension)),
        vector[-1],
    )


def test_dimension_parametric_local_norm_classes_are_complete_and_admissible():
    assert DELTA7 == 2 - THETA7
    assert DELTA7.norm() == 7
    assert all(DELTA7.sign(embedding) > 0 for embedding in range(3))
    assert L2_BNA1_SEED.canonical_mirror_norms == (Cubic7(1), Cubic7(2))
    assert L7_BNA1_SEED.canonical_mirror_norms == (
        Cubic7(1),
        Cubic7(2),
        DELTA7,
        2 * DELTA7,
    )
    for dimension in (5, 22, 29):
        for seed in (L2_BNA1_SEED, L7_BNA1_SEED):
            coefficients = seed.coefficients(dimension)
            assert len(coefficients) == dimension + 1
            assert coefficients[0] == -THETA7
            assert coefficients[-1] == seed.last_coefficient
            assert all(coefficient.sign(0) > 0 for coefficient in coefficients[1:])
            assert all(
                coefficient.sign(embedding) > 0
                for coefficient in coefficients
                for embedding in (1, 2)
            )


def test_height_zero_exhaustion_is_exact_b_n_minus_1_times_a1_through_h29():
    for dimension in (5, 22, 29):
        ordinary_count = dimension - 1
        expected_root_count = 2 * ordinary_count * ordinary_count + 2
        for seed in (L2_BNA1_SEED, L7_BNA1_SEED):
            patterns = enumerate_height_zero_patterns(seed, dimension)
            assert len(patterns) == 3
            assert sum(pattern.orbit_size for pattern in patterns) == expected_root_count
            assert height_zero_is_bn_minus_1_times_a1(seed, dimension)


def test_dimension_four_specialization_replays_the_independent_rank_five_backend():
    for new_seed, old_seed in (
        (L2_BNA1_SEED, L2_SEED),
        (L7_BNA1_SEED, L7_SEED),
    ):
        new = enumerate_bna1_root_pool(new_seed, dimension=4, height_bound=5)
        old = enumerate_nonunimodular_root_pool(old_seed, height_bound=5)
        assert new.control_counts_by_norm == old.control_counts_by_norm
        assert {
            (root.norm, root.height, root.vector) for root in new.roots
        } == {
            (root.norm, root.height, root.vector) for root in old.roots
        }


def test_rank_six_h5_prefixes_and_all_extreme_rays_are_exact():
    expected = {
        "L2-BnA1": {
            "pool": Counter({Cubic7(2): 3, Cubic7(1): 1}),
            "accepted": Counter({Cubic7(2): 1}),
            "walls": 6,
            "rays": (3, 0, 3),
        },
        "L7-BnA1": {
            "pool": Counter({Cubic7(2): 3, 2 * DELTA7: 1, Cubic7(1): 1}),
            "accepted": Counter({Cubic7(2): 1, 2 * DELTA7: 1}),
            "walls": 7,
            "rays": (5, 0, 5),
        },
    }
    for seed in (L2_BNA1_SEED, L7_BNA1_SEED):
        search = run_bna1_vinberg_search(seed, dimension=5, height_bound=5)
        datum = expected[seed.seed_id]
        assert Counter(root.norm for root in search.pool.roots) == datum["pool"]
        assert Counter(root.norm for root in search.accepted_roots) == datum["accepted"]
        assert len(search.walls) == datum["walls"]
        assert search.geometry is not None
        assert (
            search.geometry.finite_vertex_count,
            search.geometry.cusp_count,
            search.geometry.hyperideal_ray_count,
        ) == datum["rays"]
        assert search.geometry.irredundant
        assert search.geometry.consistent_time_sheet
        assert not search.geometry.finite_volume
        assert not search.geometry.compact
        assert all(
            bna1_diagonal_pair(seed, 5, left, right).sign(0) <= 0
            for index, left in enumerate(search.walls)
            for right in search.walls[index + 1 :]
        )


def test_target_dimension_h5_pools_are_certified_stabilizations_not_new_closures():
    for seed in (L2_BNA1_SEED, L7_BNA1_SEED):
        rank_six = enumerate_bna1_root_pool(seed, dimension=5, height_bound=5)
        target = enumerate_bna1_root_pool(seed, dimension=22, height_bound=5)
        assert target.control_counts_by_norm == rank_six.control_counts_by_norm
        assert {
            (root.norm, root.height, root.vector) for root in target.roots
        } == {
            (
                root.norm,
                root.height,
                _stabilize_rank_six_vector(root.vector, target_dimension=22),
            )
            for root in rank_six.roots
        }


def test_same_height_survivors_are_mutually_acute_and_order_independent():
    for dimension, height in ((5, 20), (22, 5)):
        for seed in (L2_BNA1_SEED, L7_BNA1_SEED):
            search = run_bna1_vinberg_search(
                seed, dimension=dimension, height_bound=height, analyze_geometry=False
            )
            walls = list(search.pool.cone_roots)
            accepted_forward = []
            accepted_reverse = []
            for _, level_iterator in groupby(
                search.pool.roots, key=lambda root: root.height
            ):
                level = tuple(level_iterator)
                forward = tuple(
                    root
                    for root in level
                    if all(
                        bna1_diagonal_pair(seed, dimension, root.vector, wall).sign(0)
                        <= 0
                        for wall in walls
                    )
                )
                reverse = tuple(
                    root
                    for root in reversed(level)
                    if all(
                        bna1_diagonal_pair(seed, dimension, root.vector, wall).sign(0)
                        <= 0
                        for wall in walls
                    )
                )
                assert set(forward) == set(reverse)
                assert all(
                    bna1_diagonal_pair(seed, dimension, left.vector, right.vector).sign(0)
                    <= 0
                    for index, left in enumerate(forward)
                    for right in forward[index + 1 :]
                )
                accepted_forward.extend(forward)
                accepted_reverse.extend(reverse)
                walls.extend(root.vector for root in forward)
            assert set(accepted_forward) == set(search.accepted_roots)
            assert set(accepted_reverse) == set(search.accepted_roots)


def test_ordered_square_partitions_are_the_full_height_zero_weyl_quotient():
    """The high-dimensional search is exhaustive modulo B_(d-1) x A1.

    This is stronger than observing that enumerated roots happen to lie in an
    acute sector: arbitrary signed permutations, the final A1 reflection,
    and reversal of root orientation all return the same exact representative.
    """

    for seed in (L2_BNA1_SEED, L7_BNA1_SEED):
        pool = enumerate_bna1_root_pool(seed, dimension=5, height_bound=20)
        assert all(
            dominant_bna1_representative(seed, 5, root.vector) == root.vector
            for root in pool.roots
        )
        vector = pool.roots[-1].vector
        transformed = (
            -vector[0],
            vector[4],
            -vector[2],
            vector[1],
            -vector[3],
            -vector[5],
        )
        assert dominant_bna1_representative(seed, 5, transformed) == vector
        assert bna1_diagonal_pair(seed, 5, transformed, transformed) == (
            bna1_diagonal_pair(seed, 5, vector, vector)
        )


def test_incremental_rays_match_independent_support_oracle_in_h5_and_h22():
    def projective_key(vector):
        pivot = next(entry for entry in vector if entry)
        return tuple(entry / pivot for entry in vector)

    for dimension, height in ((5, 20), (22, 5)):
        for seed in (L2_BNA1_SEED, L7_BNA1_SEED):
            search = run_bna1_vinberg_search(seed, dimension, height)
            assert search.geometry is not None
            incremental = analyze_bna1_cone_incremental(
                seed, dimension, search.walls
            )
            assert {
                projective_key(ray.vector) for ray in incremental.extreme_rays
            } == {
                projective_key(ray.vector) for ray in search.geometry.extreme_rays
            }
            assert (
                incremental.finite_vertex_count,
                incremental.cusp_count,
                incremental.hyperideal_ray_count,
                incremental.finite_volume,
                incremental.compact,
            ) == (
                search.geometry.finite_vertex_count,
                search.geometry.cusp_count,
                search.geometry.hyperideal_ray_count,
                search.geometry.finite_volume,
                search.geometry.compact,
            )


def test_h50_patterns_stabilize_exactly_from_h22_through_h29():
    """Endpoint equality certifies every intermediate dimension.

    A d-dimensional ordered partition injects into dimension d+1 by adding
    a zero.  Conversely, equality of the complete d=22 and d=29 pools shows
    that no solution using 22--28 nonzero ordinary coordinates first appears
    in an intermediate rank.  Acute selection is then identical as well.
    """

    for seed in (L2_BNA1_SEED, L7_BNA1_SEED):
        low = run_bna1_vinberg_search(seed, 22, 50, analyze_geometry=False)
        high = run_bna1_vinberg_search(seed, 29, 50, analyze_geometry=False)

        def padded(roots):
            return {
                (root.norm, root.height, _zero_pad_vector(root.vector, 22, 29))
                for root in roots
            }

        assert padded(low.pool.roots) == {
            (root.norm, root.height, root.vector) for root in high.pool.roots
        }
        assert padded(low.accepted_roots) == {
            (root.norm, root.height, root.vector) for root in high.accepted_roots
        }
        assert max(
            sum(bool(coordinate) for coordinate in root.vector[1:-1])
            for root in high.pool.roots
        ) == 5
