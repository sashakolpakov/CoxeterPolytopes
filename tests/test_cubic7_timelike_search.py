from collections import Counter
from fractions import Fraction

from coxeter_polytopes.exact.cubic7 import THETA7, Cubic7
from coxeter_polytopes.search.cubic7_seed_manifest import alvin_rc7_to_cubic7
from coxeter_polytopes.search.cubic7_timelike import (
    DELTA7,
    MINKOWSKI_BOUND_DENOMINATOR,
    MINKOWSKI_BOUND_NUMERATOR,
    T2_TIMELIKE_SEED,
    T7_TIMELIKE_SEED,
    UNIT_SIGNATURE_BITS,
    UNIT_SIGNATURE_GENERATORS,
    enumerate_height_zero_mirrors,
    enumerate_timelike_cone_rays_bruteforce,
    enumerate_timelike_cone_rays_incremental,
    height_zero_is_b4,
    run_timelike_vinberg_search,
    timelike_diagonal_pair,
    unit_signature_rank,
)


def test_complete_norm_class_number_theory_certificate():
    # Minkowski's bound for a totally real cubic field of discriminant 49.
    assert Fraction(MINKOWSKI_BOUND_NUMERATOR, MINKOWSKI_BOUND_DENOMINATOR) == Fraction(14, 9)
    assert Fraction(14, 9) < 2

    # These three units surject onto every real sign pattern modulo squares.
    assert UNIT_SIGNATURE_GENERATORS == (-Cubic7(1), THETA7, THETA7**2 - 2)
    assert UNIT_SIGNATURE_BITS == ((1, 1, 1), (0, 1, 1), (1, 1, 0))
    assert unit_signature_rank() == 3
    assert all(abs(unit.norm()) == 1 for unit in UNIT_SIGNATURE_GENERATORS)

    # x^3+x^2+1 has no F2 root, so (2) is inert.  Delta generates the
    # ramified prime over seven.  The local pairing-ideal proof then leaves
    # exactly these possible squarefree ideal classes.
    assert all((x**3 + x**2 + 1) % 2 for x in (0, 1))
    assert DELTA7 == 2 - THETA7
    assert DELTA7.norm() == 7
    assert all(DELTA7.sign(embedding) > 0 for embedding in range(3))
    assert T2_TIMELIKE_SEED.canonical_mirror_norms == (Cubic7(1), Cubic7(2))
    assert T7_TIMELIKE_SEED.canonical_mirror_norms == (
        Cubic7(1),
        Cubic7(2),
        DELTA7,
        2 * DELTA7,
    )
    assert "q(r) divides 2I(r)" in T2_TIMELIKE_SEED.complete_norm_class_proof
    assert "q(r) divides 2I(r)" in T7_TIMELIKE_SEED.complete_norm_class_proof


def test_height_zero_exhaustion_is_exactly_b4_for_both_seeds():
    for seed in (T2_TIMELIKE_SEED, T7_TIMELIKE_SEED):
        roots = enumerate_height_zero_mirrors(seed)
        assert len(roots) == 32
        assert height_zero_is_b4(seed)
        assert Counter(timelike_diagonal_pair(seed, root, root) for root in roots) == Counter(
            {Cubic7(2): 24, Cubic7(1): 8}
        )
        assert all(root[0] == 0 for root in roots)


def test_exact_height_100_prefix_counts_and_simultaneous_packets():
    expected = {
        "T2": (13, 7, 11, (11, 0, 15)),
        "T7": (16, 6, 10, (6, 0, 15)),
    }
    for seed in (T2_TIMELIKE_SEED, T7_TIMELIKE_SEED):
        search = run_timelike_vinberg_search(seed, 100)
        pool_count, accepted_count, wall_count, ray_counts = expected[seed.seed_id]
        assert len(search.pool.roots) == pool_count
        assert len(search.accepted_roots) == accepted_count
        assert len(search.walls) == wall_count
        assert all(
            timelike_diagonal_pair(seed, left, right).sign(0) <= 0
            for index, left in enumerate(search.walls)
            for right in search.walls[index + 1 :]
        )
        assert max(Counter(root.height for root in search.accepted_roots).values()) >= 2
        for root in search.pool.roots:
            assert timelike_diagonal_pair(seed, root.vector, root.vector) == root.norm
            assert (Cubic7(100) - root.height).sign(0) >= 0
            assert root.control.sign(0) > 0
            assert all(
                (2 * coefficient * coordinate / root.norm).is_integral
                for coefficient, coordinate in zip(seed.coefficients, root.vector)
            )
        geometry = enumerate_timelike_cone_rays_incremental(search)
        assert (
            geometry.finite_count,
            geometry.cusp_count,
            geometry.hyperideal_count,
        ) == ray_counts
        assert geometry.irredundant
        assert not geometry.finite_volume


def test_incremental_ray_engine_matches_independent_support_enumeration():
    for seed in (T2_TIMELIKE_SEED, T7_TIMELIKE_SEED):
        search = run_timelike_vinberg_search(seed, 50)
        incremental = enumerate_timelike_cone_rays_incremental(search)
        brute = enumerate_timelike_cone_rays_bruteforce(search)
        assert (
            incremental.finite_count,
            incremental.cusp_count,
            incremental.hyperideal_count,
        ) == (
            brute.finite_count,
            brute.cusp_count,
            brute.hyperideal_count,
        )
        assert {
            tuple(entry / next(value for value in ray.vector if value) for entry in ray.vector)
            for ray in incremental.rays
        } == {
            tuple(entry / next(value for value in ray.vector if value) for entry in ray.vector)
            for ray in brute.rays
        }


def _rc7_vector(*coordinates):
    return tuple(
        Cubic7(coordinate)
        if isinstance(coordinate, int)
        else alvin_rc7_to_cubic7(coordinate)
        for coordinate in coordinates
    )


def test_first_five_walls_match_live_alvin_packets_exactly():
    t2_alvin = {
        _rc7_vector((0, 0, -1), (0, -1, -1), (0, -1, -1), 0, 0),
        _rc7_vector((0, -1, -2), (0, -2, -3), (0, 0, -1), 0, 0),
        _rc7_vector((0, -1, -2), (0, -1, -2), (0, -1, -2), (0, -1, -1), (0, -1, -1)),
        _rc7_vector((-1, -2, -3), (0, -2, -4), 0, 0, 0),
        _rc7_vector((-1, -2, -3), (0, -1, -2), (0, -1, -2), (0, -1, -2), (0, -1, -2)),
    }
    t7_alvin = {
        _rc7_vector((0, 0, -1), (0, -1, -1), 0, 0, 0),
        _rc7_vector((-1, -3, -4), (-1, -2, -3), (-1, -2, -3), 1, 0),
        _rc7_vector((-1, -3, -4), (-1, -2, -3), (0, -1, -2), (0, -1, -1), (0, -1, -1)),
        _rc7_vector((-1, -3, -5), (-1, -2, -3), (-1, -2, -3), (-1, -2, -3), (0, 0, -1)),
        _rc7_vector((-1, -3, -5), (-1, -2, -4), (0, -1, -2), (0, -1, -2), 0),
    }
    assert {
        root.vector for root in run_timelike_vinberg_search(T2_TIMELIKE_SEED, 50).accepted_roots
    } == t2_alvin
    assert {
        root.vector for root in run_timelike_vinberg_search(T7_TIMELIKE_SEED, 50).accepted_roots
    } == t7_alvin
