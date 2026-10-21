"""Independent hostile replay of the direct full-NS rational n<=11 theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_direct_full_ns_n11_dimension_theorem.md"
CERTIFIER = (
    ROOT / "scripts/certify_rational_direct_full_ns_n11_dimension_theorem.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/rational_direct_full_ns_n11_dimension_theorem.json"
)
OWNER_TESTS = (
    ROOT / "tests/test_rational_direct_full_ns_n11_dimension_theorem.py"
)
PINS = {
    REPORT: "a5dbd6e2fbbd7cd55025f5af324de6b778a196cf4f72cf06903a1dce2e4ae27f",
    CERTIFIER: "23715fc4ed8ab4505cb9925b57b770eba388b48ef74f7f9bf261f5cdb5b576d2",
    ARCHIVE: "faea105582a7973363688f49710251e26269d861da5fb3bd13c629171b9ebbed",
    OWNER_TESTS: "db92c6f9078f19b7e48dc084c4c36d90ecd81c2de907b144fec54e17a11dd987",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plane_pair(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def f2_pair(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return (
        -2 * left[0] * right[0]
        + left[0] * right[1]
        + left[1] * right[0]
        - sum(a * b for a, b in zip(left[2:], right[2:]))
    )


def monotone_patterns(length: int, maximum: int):
    for row in combinations_with_replacement(range(maximum + 1), length):
        yield tuple(reversed(row))


def test_owner_artifacts_are_exactly_pinned() -> None:
    assert {path: sha256(path) for path in PINS} == PINS


def test_scope_allows_square_zero_tangent_rays_but_not_extra_negative_species() -> None:
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    scope = payload["scope"]
    assert scope["physical_negative_extremal_alphabet"] == (
        "smooth rational (-1),(-2),(-4)"
    )
    assert scope["square_zero_tangent_Mori_rays_allowed"] is True
    assert scope["host"] == "direct full Neron-Severi lattice only"
    assert scope["adapted_half_boundary_host_claimed"] is False
    conclusion = payload["conclusion"]
    assert conclusion["plane_blowups_at_most"] == 11
    assert conclusion["F2_blowups_at_most"] == 10
    assert conclusion["hyperbolic_dimension_at_most"] == 11
    assert conclusion["projective_nef_dimension"].startswith("rho-1")


def test_satellite_final_exceptional_prime_is_always_forbidden() -> None:
    alphabet = {-1, -2, -4}
    for children in range(2, 80):
        root = (1, *((-1,) * children))
        square = -sum(value * value for value in root)
        assert square == -(children + 1)
        assert root[0] == 1  # primitive
        if square == -4:
            # Reflection of the first coordinate basis vector has halves.
            coefficient = Fraction(-2 * (-1), square)
            image = tuple(
                Fraction(int(index == 0)) + coefficient * value
                for index, value in enumerate(root)
            )
            assert any(entry.denominator == 2 for entry in image)
        else:
            assert square not in alphabet


def test_degree_one_adjunction_gates_are_unique_and_half_minus_four_is_empty() -> None:
    for length in (3, 5, 12, 20):
        def solve(total: int, squares: int):
            return [
                row
                for row in monotone_patterns(length, min(total, squares))
                if sum(row) == total
                and sum(value * value for value in row) == squares
            ]

        assert solve(2, 2) == [(1, 1, *((0,) * (length - 2)))]
        assert solve(3, 3) == [(1, 1, 1, *((0,) * (length - 3)))]
        assert solve(4, 2) == []
        assert solve(1, 1) == [(1, *((0,) * (length - 1)))]
        assert solve(3, 1) == []


def test_f2_ample_degree_does_not_hide_the_negative_section() -> None:
    ample = (1, 3)
    section = (1, 0)
    fibre = (0, 1)
    assert f2_pair(ample, section) == 1
    assert f2_pair(ample, fibre) == 1
    delta = (0, 1, -1, -1)
    ample_pullback = (1, 3, 0, 0)
    assert f2_pair(ample_pullback, delta) == 1
    # Pushing f-e1-e2 to F2 gives f, never the section s.
    assert delta[:2] == fibre
    assert delta[:2] != section


def test_exceptional_remainders_have_a_negative_free_chain_coefficient() -> None:
    for m in range(3, 40):
        plane = tuple(-1 if index == 2 else 0 for index in range(m))
        ruled = tuple(-1 if index == 1 else 0 for index in range(m))

        def coefficients(target: tuple[int, ...]) -> tuple[int, ...]:
            running = 0
            answer = []
            for entry in target:
                running += entry
                answer.append(running)
            return tuple(answer)

        plane_coefficients = coefficients(plane)
        ruled_coefficients = coefficients(ruled)
        assert plane_coefficients[2] == -1
        assert ruled_coefficients[1] == -1
        assert not all(value >= 0 for value in plane_coefficients)
        assert not all(value >= 0 for value in ruled_coefficients)


def test_first_mixed_rays_reconstruct_independently() -> None:
    for m in range(12, 45):
        ray = (7, -3, *((-2,) * 11), *((0,) * (m - 12)))
        alpha = (1, -1, -1, -1, *((0,) * (m - 3)))
        walls = [alpha]
        for index in range(1, m):
            wall = [0] * (m + 1)
            wall[index], wall[index + 1] = 1, -1
            walls.append(tuple(wall))
        walls.append((0, *((0,) * (m - 1)), 1))
        pairings = [plane_pair(ray, wall) for wall in walls]
        assert plane_pair(ray, ray) == -4
        assert pairings == [
            1 if index == 1 else 2 if index == 12 else 0
            for index in range(m + 1)
        ]

    for m in range(11, 45):
        ray = (4, 9, *((-2,) * 11), *((0,) * (m - 11)))
        walls = [
            (1, 0, *((0,) * m)),
            (0, 1, -1, -1, *((0,) * (m - 2))),
        ]
        for index in range(1, m):
            wall = [0] * (m + 2)
            wall[index + 1], wall[index + 2] = 1, -1
            walls.append(tuple(wall))
        walls.append((0, 0, *((0,) * (m - 1)), 1))
        pairings = [f2_pair(ray, wall) for wall in walls]
        assert f2_pair(ray, ray) == -4
        assert pairings == [
            1 if index == 0 else 2 if index == 12 else 0
            for index in range(m + 2)
        ]


def test_low_pairing_minus_four_branches_have_no_solution() -> None:
    plane_solutions = []
    for R in range(4):
        for tail in monotone_patterns(11, R):
            if tail[0] + tail[1] > R:
                continue
            numerator = sum(tail) - R - 2
            if numerator < 0 or numerator % 2:
                continue
            d = numerator // 2
            if d - R < tail[0]:
                continue
            if sum(value * value for value in tail) == 2 * d * R - R * R + 4:
                plane_solutions.append((R, d, tail))
    assert plane_solutions == []

    ruled_solutions = []
    for a in range(4):
        for row in monotone_patterns(11, a):
            if row[0] + row[1] > a:
                continue
            numerator = sum(row) - 2
            if numerator < 0 or numerator % 2:
                continue
            b = numerator // 2
            if b < 2 * a:
                continue
            if sum(value * value for value in row) == -2 * a * a + 2 * a * b + 4:
                ruled_solutions.append((a, b, row))
    assert ruled_solutions == []


def test_report_closes_nonsimplicial_and_short_chain_gaps() -> None:
    report = REPORT.read_text(encoding="utf-8")
    required = (
        "same** Coxeter group",
        "symmetrizable generalized",
        "facet generators give the Coxeter",
        "choose a generic",
        "point of its mirror",
        "nonsimpliciality causes no gap",
        r"\(m<3\) over \(\mathbf P^2\)",
        r"\(m<2\) over \(\mathbb F_2\)",
        r"use the ample degree \((s+3f).x\)",
        "pushing (3.6) down to the base",
        "Square-zero tangent rays are allowed",
    )
    for phrase in required:
        assert phrase in report


def test_owner_archive_replays_without_importing_owner_code() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
    )
