from coxeter_polytopes.exact.cubic7 import Cubic7, THETA7
from coxeter_polytopes.search.cubic7_e8_cubes import (
    DELTA_THETA_E8_CUBE_SEED,
    SUPPORTED_E8_CUBE_SEEDS,
    THETA_E8_CUBE_SEED,
    TWO_THETA_E8_CUBE_SEED,
    e8_cube_pair,
    height_zero_is_e8_cubed,
    run_h24_e8_cube_vinberg,
)
from coxeter_polytopes.search.cubic7_root_blocks import (
    DELTA_THETA_ROOT_BLOCK_SEED,
    E6_RADICAL_MOD_3,
    E6_RADICAL_SQUARE,
    SUPPORTED_ROOT_BLOCK_SEEDS,
    THETA_ROOT_BLOCK_SEED,
    TWO_THETA_ROOT_BLOCK_SEED,
    _e6_cartan_gram,
    height_zero_is_e8_e8_e6,
    root_block_pair,
    run_h22_root_block_vinberg,
)


def test_e6_discriminant_radical_local_certificate() -> None:
    gram = _e6_cartan_gram()
    radical = E6_RADICAL_MOD_3
    covector = tuple(sum(gram[i][j] * radical[j] for j in range(6)) for i in range(6))
    square = sum(radical[i] * gram[i][j] * radical[j] for i in range(6) for j in range(6))
    assert all(entry % 3 == 0 for entry in covector)
    assert square == E6_RADICAL_SQUARE == 12
    assert square % 18 == 12


def test_complete_norm_representative_lists() -> None:
    delta = 2 - THETA7
    assert THETA_ROOT_BLOCK_SEED.canonical_mirror_norms == tuple(map(Cubic7, (1, 2, 3, 6)))
    assert TWO_THETA_ROOT_BLOCK_SEED.canonical_mirror_norms == tuple(map(Cubic7, (1, 2, 3, 6)))
    assert DELTA_THETA_ROOT_BLOCK_SEED.canonical_mirror_norms == (
        Cubic7(1), Cubic7(2), Cubic7(3), Cubic7(6), delta, 2 * delta, 3 * delta, 6 * delta
    )
    assert THETA_E8_CUBE_SEED.canonical_mirror_norms == (Cubic7(1), Cubic7(2))
    assert TWO_THETA_E8_CUBE_SEED.canonical_mirror_norms == (Cubic7(1), Cubic7(2))
    assert DELTA_THETA_E8_CUBE_SEED.canonical_mirror_norms == (
        Cubic7(1), Cubic7(2), delta, 2 * delta
    )


def test_height_zero_dominant_censuses_are_complete() -> None:
    assert all(height_zero_is_e8_e8_e6(seed) for seed in SUPPORTED_ROOT_BLOCK_SEEDS)
    assert all(height_zero_is_e8_cubed(seed) for seed in SUPPORTED_E8_CUBE_SEEDS)


def test_h22_two_theta_first_packet() -> None:
    search = run_h22_root_block_vinberg(
        TWO_THETA_ROOT_BLOCK_SEED, 10, analyze_geometry=True
    )
    assert len(search.control_levels) == 1
    assert len(search.accepted_roots) == 3
    assert len(search.walls) == 25
    assert {root.height for root in search.accepted_roots} == {2 * THETA7 + THETA7**2}
    assert all(root_block_pair(search.seed, root.vector, root.vector) == 2 for root in search.accepted_roots)
    geometry = search.geometry
    assert geometry is not None
    assert len(geometry.extreme_rays) == 9 * 9 * 7 == 567
    assert len({ray.vector for ray in geometry.extreme_rays}) == 567
    assert (geometry.finite_vertex_count, geometry.cusp_count, geometry.hyperideal_ray_count) == (4, 0, 563)
    assert all(len(ray.active_facets) == 22 for ray in geometry.extreme_rays)
    assert all(ray.active_rank_verified for ray in geometry.extreme_rays)


def test_h24_h10_prefixes_and_exact_product_of_simplices_cone() -> None:
    theta = run_h24_e8_cube_vinberg(THETA_E8_CUBE_SEED, 10)
    twice = run_h24_e8_cube_vinberg(TWO_THETA_E8_CUBE_SEED, 10)
    delta = run_h24_e8_cube_vinberg(DELTA_THETA_E8_CUBE_SEED, 10)
    assert (len(theta.control_levels), len(theta.accepted_roots)) == (5, 0)
    assert (len(twice.control_levels), len(twice.accepted_roots)) == (1, 3)
    assert (len(delta.control_levels), len(delta.accepted_roots)) == (7, 0)
    assert all(e8_cube_pair(twice.seed, root.vector, root.vector) == 2 for root in twice.accepted_roots)
    geometry = twice.geometry
    assert geometry is not None
    assert len(geometry.extreme_rays) == 9**3 == 729
    assert len({ray.vector for ray in geometry.extreme_rays}) == 729
    assert (geometry.finite_vertex_count, geometry.cusp_count, geometry.hyperideal_ray_count) == (4, 0, 725)
    assert geometry.irredundant
    assert not geometry.finite_volume
    assert all(len(ray.active_facets) == 24 for ray in geometry.extreme_rays)
    assert all(ray.active_rank_verified for ray in geometry.extreme_rays)
