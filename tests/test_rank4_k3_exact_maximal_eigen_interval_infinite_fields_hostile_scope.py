"""Standalone hostile replay of the maximal rank-four eigen-interval theorem."""

from __future__ import annotations

import hashlib
from math import isqrt, sqrt
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = (
    ROOT
    / "results/rank4_k3_exact_maximal_eigen_interval_infinite_fields.md"
)
CERTIFIER = (
    ROOT
    / "scripts/certify_rank4_k3_exact_maximal_eigen_interval_infinite_fields.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/rank4_k3_exact_maximal_eigen_interval_infinite_fields.json"
)
OWNER_TESTS = (
    ROOT
    / "tests/test_rank4_k3_exact_maximal_eigen_interval_infinite_fields.py"
)
AUDIT = (
    ROOT
    / "results/rank4_k3_exact_maximal_eigen_interval_infinite_fields_hostile_scope_audit.md"
)
PINS = {
    REPORT: "9351875964fbc73da81b55ec3ff68db54af3dbe1ec90272555cff1938ab13c95",
    CERTIFIER: "1b1e6c60ef6e7987adac7768dc9fdf87345482e6986d8966618c10c1a5eaa48d",
    ARCHIVE: "a238ea60701c86d9c2c688ee6fa889b73c9953632a5ef1f344e2ef5441b900cd",
    OWNER_TESTS: "9cf8b5c551764c9f456c2d204eee6524d2330c2640f7c28b13b6032913adad69",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def transpose(matrix):
    return tuple(
        tuple(row[column] for row in matrix)
        for column in range(len(matrix[0]))
    )


def multiply(left, right):
    columns = transpose(right)
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def gram(a: int):
    return (
        (-2, 0, 0, 0),
        (0, 2, 0, 0),
        (0, 0, 2 - 4 * a, -2),
        (0, 0, -2, -2),
    )


def companion(a: int):
    return ((0, -1), (1, 2 * a))


def lift_reflections(a: int):
    c = companion(a)
    identity = ((1, 0), (0, 1))

    def linear(c0, c1):
        return tuple(
            tuple(
                c0 * int(row == column) + c1 * c[row][column]
                for column in range(2)
            )
            for row in range(2)
        )

    def blocks(top_left, top_right, bottom_left, bottom_right):
        return tuple(
            tuple(
                (top_left if row < 2 and column < 2 else
                 top_right if row < 2 else
                 bottom_left if column < 2 else bottom_right)
                [row % 2][column % 2]
                for column in range(4)
            )
            for row in range(4)
        )

    s0 = blocks(identity, ((0, 0), (0, 0)), ((0, 0), (0, 0)), ((-1, 0), (0, -1)))
    s1 = blocks(
        linear(-1, 2),
        linear(2, -2),
        linear(0, 2),
        linear(1, -2),
    )
    return s0, s1


def q(left, right, d):
    return (
        -left[0] * right[0]
        + left[1] * right[1]
        + left[2] * right[2]
        + 2 * d * left[3] * right[3]
    )


def sector_roots(d: int, max_height: int):
    for height in range(1, max_height + 1):
        target = height * height + 1
        max_z = isqrt(target // (2 * d))
        for z in range(-max_z, max_z + 1):
            remaining = target - 2 * d * z * z
            if remaining < 0:
                continue
            for m in range(isqrt(remaining) + 1):
                u = isqrt(remaining - m * m)
                if u * u == remaining - m * m:
                    yield height, -m, u, z


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_integral_lifts_reconstruct_without_producer_imports():
    identity4 = tuple(
        tuple(int(row == column) for column in range(4))
        for row in range(4)
    )
    for a in range(2, 401):
        c = companion(a)
        tau = (
            (c[0][0], c[0][1], 0, 0),
            (c[1][0], c[1][1], 0, 0),
            (0, 0, c[0][0], c[0][1]),
            (0, 0, c[1][0], c[1][1]),
        )
        for reflection in lift_reflections(a):
            assert multiply(reflection, reflection) == identity4
            assert multiply(
                multiply(transpose(reflection), gram(a)), reflection
            ) == gram(a)
            assert multiply(reflection, tau) == multiply(tau, reflection)


def test_exact_root_free_interval_replays_on_a_large_shell():
    for a in range(2, 61):
        d = a - 1
        t = a + sqrt(a * a - 1)
        equalities = []
        for height, x, u, z in sector_roots(d, 120):
            p0 = -t * height - x
            p1 = -t * height - x + (t - 1) * u - 2 * d * z
            assert p0 < -1e-10
            assert p1 <= 1e-10
            if abs(p1) <= 1e-10:
                equalities.append((height, x, u, z))
        assert equalities == [(1, -1, 1, 0)]


def test_maximality_trace_window_and_square_obstructions_are_empty():
    for a in range(2, 401):
        d = a - 1
        t = a + sqrt(a * a - 1)
        other = 1 / t
        candidates = []
        for m in range(-50, 51):
            for n in range(-50, 51):
                if (
                    1 < m + n * t < 2 * t - 1
                    and abs(m + n * other) <= 1 + 1e-12
                ):
                    candidates.append((m, n))
        assert candidates == [(-1, 1), (0, 1)]

        for q0 in range(-80, 81):
            for q1 in range(-80, 81):
                square = (
                    q0 * q0 - q1 * q1,
                    2 * q1 * (q0 + a * q1),
                )
                assert square != (2 * d, 2 * d)
                assert square != (1, 2 * d - 1)


def test_ambient_infinite_facets_and_fixed_polarization_are_independent():
    for d in range(1, 81):
        polarization = (3, -1, 1, 0)
        assert -2 * q(polarization, polarization, d) == 14
        for root in sector_roots(d, 60):
            assert q(polarization, root, d) < 0
        for k in range(1, 31):
            height = d * k * k
            root = (height, 0, height - 1, k)
            support = (
                1 + height * height,
                0,
                height * (height - 1),
                height * k,
            )
            assert q(root, root, d) == 1
            assert q(support, root, d) == 0
            assert q(support, support, d) == -1 - height * height


def test_archive_scope_does_not_promote_dimension_or_absolute_maximality():
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert payload["hyperbolic_dimension"] == 1
    conclusions = payload["conclusions"]
    assert conclusions["finite_volume_eigen_interval"] is True
    assert conclusions["complete_restricted_root_saturation"] is True
    assert conclusions["reflection_maximal_in_declared_Z_t_host"] is True
    assert conclusions["fixed_ample_polarization_square"] == 14
    assert conclusions["bounded_polarized_correspondence_height"] is False
    assert conclusions["finite_field_set"] is False
    scope = payload["scope"]
    assert scope["n_at_least_two_claimed"] is False
    assert scope["maximal_across_overorders_claimed"] is False
    assert scope["maximal_across_commensurable_overlattices_claimed"] is False
    assert scope["maximal_in_full_O_L_claimed"] is False


def test_reports_state_the_scope_boundary_literally():
    theorem = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "hyperbolic-dimension-one counterfamily",
        "relative to the declared integral full-trace order and lattice",
        "No maximality across arbitrary overorders",
        "does not settle the required \\(n\\ge2\\) field theorem",
        "primitive ample polarization of fixed degree \\(14\\)",
    ):
        assert phrase in theorem
    for phrase in (
        "PASS for \\(n=1\\)",
        "No \\(n\\ge2\\) or overorder-maximal theorem",
        "literal host-relative maximality statement",
        "promoting the stated result to full-\\(O(L_a)\\) maximality would be false",
    ):
        assert phrase in audit


def test_archive_replay_passes_in_a_fresh_process():
    completed = subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "PASS archive replay" in completed.stdout
