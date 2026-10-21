from collections import Counter
from fractions import Fraction

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_e8_cubes import (
    INADMISSIBLE_TWO_PLUS_SQRT2,
    SQRT2_E8_CUBE_SEED,
    TWO_SQRT2_E8_CUBE_SEED,
    bounded_silver_integers,
    coefficient_is_arithmetically_admissible,
    first_packet_product_geometry,
    height_zero_dominant_roots,
    height_zero_is_e8_cubed,
    run_silver_h24_e8_cube_vinberg,
    silver_e8_cube_pair,
)


def test_coefficient_signature_gate_and_complete_norm_classes():
    assert coefficient_is_arithmetically_admissible(SQRT2)
    assert coefficient_is_arithmetically_admissible(2 * SQRT2)
    assert not coefficient_is_arithmetically_admissible(INADMISSIBLE_TWO_PLUS_SQRT2)
    assert SQRT2_E8_CUBE_SEED.canonical_mirror_norms == (
        Silver(2),
        Silver(2, 1),
    )
    assert TWO_SQRT2_E8_CUBE_SEED.canonical_mirror_norms == (
        Silver(2),
        Silver(4, 2),
    )
    for seed in (SQRT2_E8_CUBE_SEED, TWO_SQRT2_E8_CUBE_SEED):
        assert seed.admissible
        assert "every totally-positive unit is a square" in seed.complete_norm_class_proof
        assert all(norm.sign(place) > 0 for norm in seed.canonical_mirror_norms for place in (0, 1))


def test_bounded_silver_minkowski_rectangle_is_exact():
    values = bounded_silver_integers((Fraction(4), Fraction(1)))
    brute = {
        Silver(a, b)
        for a in range(-5, 6)
        for b in range(-5, 6)
        if all((Silver(cap) - Silver(a, b) ** 2).sign(place) >= 0 for place, cap in enumerate((4, 1)))
    }
    assert set(values) == brute


def test_height_zero_is_exactly_e8_cubed():
    for seed in (SQRT2_E8_CUBE_SEED, TWO_SQRT2_E8_CUBE_SEED):
        roots = height_zero_dominant_roots(seed)
        assert height_zero_is_e8_cubed(seed)
        assert len(roots) == 3
        assert all(root.norm == 2 and not root.control for root in roots)


def test_h10_screen_and_product_of_simplices_geometry_are_exact():
    expected = {
        SQRT2_E8_CUBE_SEED.seed_id: (5, 6, 3, Counter({Silver(2, 1): 3})),
        TWO_SQRT2_E8_CUBE_SEED.seed_id: (3, 6, 3, Counter({Silver(4, 2): 3})),
    }
    for seed in (SQRT2_E8_CUBE_SEED, TWO_SQRT2_E8_CUBE_SEED):
        search = run_silver_h24_e8_cube_vinberg(seed, 10)
        levels, candidates, accepted, norms = expected[seed.seed_id]
        assert len(search.control_levels) == levels
        assert sum(packet.candidate_count for packet in search.packets) == candidates
        assert len(search.accepted_roots) == accepted
        assert Counter(root.norm for root in search.accepted_roots) == norms
        assert len(search.walls) == 27
        assert all(
            silver_e8_cube_pair(seed, left, right).sign(0) <= 0
            for index, left in enumerate(search.walls)
            for right in search.walls[index + 1 :]
        )
        geometry = first_packet_product_geometry(search)
        assert geometry.extreme_ray_count == 729
        assert (
            geometry.finite_vertex_count,
            geometry.cusp_count,
            geometry.hyperideal_ray_count,
        ) == (7, 0, 722)
        assert not geometry.finite_volume
        assert not geometry.compact
