from pathlib import Path
import subprocess
import sys

from scripts.certify_rational_single_chain_standard_simplex_gate import (
    augmented_extreme_rays,
    build_certificate,
    canonical_class,
    doubled_stage_ten_wall,
    dual_rays,
    f2_augmented_extreme_rays,
    f2_canonical_class,
    f2_doubled_stage_nine_wall,
    f2_dual_rays,
    f2_pair,
    f2_standard_walls,
    f2_tail_ray,
    pair,
    standard_walls,
)


ROOT = Path(__file__).parents[1]
ARCHIVE = (
    ROOT / "data/enumerations/rational_single_chain_standard_simplex_gate.json"
)


def test_closed_form_duals_and_exact_visible_threshold() -> None:
    for m in range(3, 33):
        walls = standard_walls(m)
        rays = dual_rays(m)
        assert [
            [pair(ray, wall) for wall in walls]
            for ray in rays
        ] == [
            [int(i == j) for j in range(m + 1)]
            for i in range(m + 1)
        ]
        expected = [1] + [min(i, 3) ** 2 - i for i in range(1, m + 1)]
        assert [pair(ray, ray) for ray in rays] == expected
        assert all(value >= 0 for value in expected) is (m <= 9)


def test_escape_rays_have_genus_one_and_start_at_minus_one() -> None:
    for m in range(10, 33):
        canonical = canonical_class(m)
        for j in range(10, m + 1):
            ray = dual_rays(m)[j]
            assert pair(ray, ray) == 9 - j
            assert pair(canonical, ray) == j - 9
            assert 1 + (pair(ray, ray) + pair(canonical, ray)) // 2 == 1


def test_doubled_stage_ten_wall_cuts_every_old_tail_ray() -> None:
    for m in range(10, 33):
        wall = doubled_stage_ten_wall(m)
        half_root = tuple(value // 2 for value in wall)
        assert pair(wall, wall) == -4
        assert pair(half_root, half_root) == -1
        assert pair(canonical_class(m), wall) == 2
        assert [pair(wall, dual_rays(m)[j]) for j in range(10, m + 1)] == [
            -2
        ] * (m - 9)


def test_one_wall_repair_stops_at_m_eleven() -> None:
    status = {}
    for m in range(10, 20):
        squares = [pair(ray, ray) for ray in augmented_extreme_rays(m)]
        status[m] = all(value >= 0 for value in squares)
    assert [m for m, finite in status.items() if finite] == [10, 11]
    v_1 = dual_rays(12)[1]
    v_12 = dual_rays(12)[12]
    first_mixed = tuple(a + 2 * b for a, b in zip(v_1, v_12))
    assert pair(first_mixed, first_mixed) == -4


def test_f2_closed_form_duals_and_exact_visible_threshold() -> None:
    for m in range(2, 33):
        walls = f2_standard_walls(m)
        rays = f2_dual_rays(m)
        assert [
            [f2_pair(ray, wall) for wall in walls]
            for ray in rays
        ] == [
            [int(i == j) for j in range(m + 2)]
            for i in range(m + 2)
        ]
        expected = [0, 2] + [
            2 * min(i, 2) ** 2 - i for i in range(1, m + 1)
        ]
        assert [f2_pair(ray, ray) for ray in rays] == expected
        assert all(value >= 0 for value in expected) is (m <= 8)


def test_f2_escape_rays_are_genus_one_and_stage_nine_is_anticanonical() -> None:
    for m in range(9, 33):
        canonical = f2_canonical_class(m)
        for j in range(9, m + 1):
            ray = f2_tail_ray(m, j)
            assert f2_pair(ray, ray) == 8 - j
            assert f2_pair(canonical, ray) == j - 8
            assert 1 + (
                f2_pair(ray, ray) + f2_pair(canonical, ray)
            ) // 2 == 1
        if m == 9:
            assert tuple(-value for value in canonical) == f2_tail_ray(m, 9)


def test_f2_doubled_stage_nine_wall_cuts_every_old_tail_ray() -> None:
    for m in range(9, 33):
        wall = f2_doubled_stage_nine_wall(m)
        half_root = tuple(value // 2 for value in wall)
        assert f2_pair(wall, wall) == -4
        assert f2_pair(half_root, half_root) == -1
        assert f2_pair(f2_canonical_class(m), wall) == 2
        assert [
            f2_pair(wall, f2_tail_ray(m, j)) for j in range(9, m + 1)
        ] == [-2] * (m - 8)


def test_f2_one_wall_repair_stops_at_m_ten() -> None:
    status = {}
    for m in range(9, 20):
        squares = [f2_pair(ray, ray) for ray in f2_augmented_extreme_rays(m)]
        status[m] = all(value >= 0 for value in squares)
    assert [m for m, finite in status.items() if finite] == [9, 10]
    fibre = f2_dual_rays(11)[0]
    u_11 = f2_tail_ray(11, 11)
    first_mixed = tuple(a + 2 * b for a, b in zip(fibre, u_11))
    assert f2_pair(first_mixed, first_mixed) == -4


def test_scope_does_not_overclaim_a_geometric_counterfamily() -> None:
    payload = build_certificate(32)
    assert payload["status"] == "PASS_EXACT_SCOPED_NO_GO"
    scope = payload["scope"]
    assert scope["alpha_0_not_universal"] is True
    assert scope["geometric_effectivity_of_doubled_numerical_wall_claimed"] is False
    assert scope["uniform_chain_length_bound_claimed"] is False
    assert scope["classification_or_census_used"] is False
    assert payload["one_wall_augmented_cone"]["finite_volume_exactly_m"] == [10, 11]
    f2 = payload["F2_minimal_base_companion"]
    assert f2["visible_simplex"]["finite_volume_exactly_m_at_most"] == 8
    assert f2["one_wall_augmented_cone"]["finite_volume_exactly_m"] == [9, 10]
    assert f2["geometric_effectivity_of_doubled_numerical_wall_claimed"] is False
    assert f2["uniform_chain_length_bound_claimed"] is False


def test_archive_replays_exactly() -> None:
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/certify_rational_single_chain_standard_simplex_gate.py"),
            "--check",
            str(ARCHIVE),
        ],
        cwd=ROOT,
        check=True,
    )


def test_certifier_is_independent_of_existing_root_enumerators() -> None:
    source = (
        ROOT / "scripts/certify_rational_single_chain_standard_simplex_gate.py"
    ).read_text(encoding="utf-8")
    assert "certify_rational_i1n_inverse_design" not in source
    assert "enumerate_simple_roots" not in source
