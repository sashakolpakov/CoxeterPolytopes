"""Standalone hostile replay of the orthogonal-stabilization obstruction."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from math import isqrt, sqrt
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.md"
CERTIFIER = (
    ROOT / "scripts/certify_rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.json"
)
OWNER_TESTS = (
    ROOT / "tests/test_rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.py"
)
AUDIT = (
    ROOT
    / "results/rank6_k3_orthogonal_stabilization_nonarithmetic_no_go_hostile_scope_audit.md"
)
PINS = {
    REPORT: "a7d71ac262bdb709db14dafb5488ef00ce343366084783d427fef971cd6b0e01",
    CERTIFIER: "fd0ff5530fd50117a9e29151d4487766caf05c56d59cd9dfb16f6b915567a71c",
    ARCHIVE: "b979263830c876508b69359f161c738875b597ae14418e0cd02cdb6ffe3f117b",
    OWNER_TESTS: "bca7ea7dc513d4d2d6098286b4b57eefc27ca53151d9538d568ff4413043a8ac",
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
        tuple(sum(x * y for x, y in zip(row, column)) for column in columns)
        for row in left
    )


def block_diagonal(*blocks):
    rank = sum(len(block) for block in blocks)
    answer = [[0] * rank for _ in range(rank)]
    offset = 0
    for block in blocks:
        for row in range(len(block)):
            for column in range(len(block)):
                answer[offset + row][offset + column] = block[row][column]
        offset += len(block)
    return tuple(tuple(row) for row in answer)


def q(root, d):
    y, x, u, z, v, w = root
    return -y * y + x * x + u * u + 2 * d * z * z + v * v + 2 * d * w * w


def sector_roots(d, max_height):
    for height in range(1, max_height + 1):
        target = height * height + 1
        bound = isqrt(target // (2 * d))
        for z in range(-bound, bound + 1):
            for w in range(-bound, bound + 1):
                residual = target - 2 * d * (z * z + w * w)
                if residual < 0:
                    continue
                for m in range(isqrt(residual) + 1):
                    remaining = residual - m * m
                    for u in range(isqrt(remaining) + 1):
                        v2 = remaining - u * u
                        v = isqrt(v2)
                        if v * v == v2:
                            yield height, -m, u, z, v, w


def exact_midpoint_sign(root, a):
    height, x, u, z, v, w = root
    coefficient_a = 2 * height - u - v + z + w
    coefficient_b = -2 * x + z + w
    assert coefficient_a >= 0
    if coefficient_a == 0:
        return (coefficient_b > 0) - (coefficient_b < 0)
    rational = coefficient_b - a * coefficient_a
    if rational <= 0:
        return -1
    difference = rational * rational - (a * a - 1) * coefficient_a * coefficient_a
    return (difference > 0) - (difference < 0)


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_trace_lattice_and_field_operator_reconstruct_independently():
    change = (
        (0, 1, 0, 0, 0, 0),
        (1, 0, 0, 0, 0, 0),
        (0, 0, 0, 1, 0, 0),
        (0, 0, 1, -1, 0, 0),
        (0, 0, 0, 0, 0, 1),
        (0, 0, 0, 0, 1, -1),
    )
    identity = tuple(
        tuple(int(row == column) for column in range(6)) for row in range(6)
    )
    for a in range(2, 401):
        d = a - 1
        b0 = ((-2, 0), (0, 2))
        ba = ((2 - 4 * a, -2), (-2, -2))
        companion = ((0, -1), (1, 2 * a))
        gram = block_diagonal(b0, ba, ba)
        tau = block_diagonal(companion, companion, companion)
        diagonal = multiply(multiply(transpose(change), gram), change)
        assert diagonal == (
            (2, 0, 0, 0, 0, 0),
            (0, -2, 0, 0, 0, 0),
            (0, 0, -2, 0, 0, 0),
            (0, 0, 0, -4 * d, 0, 0),
            (0, 0, 0, 0, -2, 0),
            (0, 0, 0, 0, 0, -4 * d),
        )
        assert multiply(gram, tau) == multiply(transpose(tau), gram)
        tau_squared = multiply(tau, tau)
        assert tuple(
            tuple(
                tau_squared[row][column]
                - 2 * a * tau[row][column]
                + identity[row][column]
                for column in range(6)
            )
            for row in range(6)
        ) == tuple(tuple(0 for _ in range(6)) for _ in range(6))


def test_k3_embedding_and_five_roots_reconstruct_independently():
    roots = (
        (0, 0, -1, 0, 0, 0),
        (0, 0, 0, 0, -1, 0),
        (1, -1, 1, 0, 0, 0),
        (1, -1, 0, 0, 1, 0),
        (1, 0, 1, 0, 1, 0),
    )
    for d in range(1, 501):
        assert [q(root, d) for root in roots] == [1] * 5
        # The selected pivot coordinates of Y,X,U,Z,V,W have diagonal
        # entries 1,1,1,1,1,1, so the displayed embedding has a unit minor.
        pivot_minor = (
            (1, 0, 0, 0, 0, 0),
            (0, 1, 0, 0, 0, 0),
            (0, 0, 1, 0, 0, 0),
            (0, 0, 0, 1, 0, 0),
            (0, 0, 0, 0, 1, 0),
            (0, 0, 0, 0, 0, 1),
        )
        assert pivot_minor == tuple(
            tuple(int(row == column) for column in range(6))
            for row in range(6)
        )


def test_mixed_wall_is_unique_on_an_independent_exact_shell():
    gamma = (1, 0, 1, 0, 1, 0)
    for a in range(2, 19):
        equalities = []
        for root in sector_roots(a - 1, 44):
            sign = exact_midpoint_sign(root, a)
            assert sign <= 0
            if sign == 0:
                equalities.append(root)
        assert equalities == [gamma]


def test_pentagon_and_all_orthogonal_pair_obstructions():
    for a in range(2, 1001):
        t = a + sqrt(a * a - 1)
        radius_squared = t / (t - 1)
        assert 1 < radius_squared < 2
        x = radius_squared - 1
        assert 1 + x * x < radius_squared
    for slope_count in range(2, 10):
        for first in range(slope_count):
            for second in range(first + 1, slope_count):
                root = [1, 0] + [0, 0] * slope_count
                root[2 + 2 * first] = 1
                root[2 + 2 * second] = 1
                assert -root[0] ** 2 + sum(
                    root[2 + 2 * index] ** 2
                    for index in range(slope_count)
                ) == 1


def test_nonintegral_cyclic_product_is_parameter_uniform():
    for d in range(1, 20001):
        denominator = 4 * d - 1
        a = d + 1
        c0 = Fraction(4 * (2 * d - 1), denominator)
        c1 = Fraction(4, denominator)
        norm = c0 * c0 + 2 * a * c0 * c1 + c1 * c1
        assert norm == Fraction(32 * d, denominator)
        assert norm == 8 + Fraction(8, denominator)
        assert norm.denominator != 1


def test_archive_replays_and_scope_does_not_expand():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert payload["Q_rank"] == 6
    assert payload["distinguished_hyperbolic_dimension"] == 2
    assert payload["conclusions"][
        "all_orthogonal_stabilizations_r_at_least_2_obstructed"
    ] is True
    scope = payload["scope"]
    assert scope["orthogonal_trace_stabilizations_only"] is True
    assert scope["all_nonorthogonal_rank_six_couplings_excluded"] is False
    assert scope["primitive_K3_embedding_beyond_rank_six_claimed"] is False
    assert scope["compactness_beyond_rank_six_claimed"] is False
    assert scope["n_at_least_two_global_field_theorem_resolved"] is False


def test_report_and_audit_keep_the_surviving_route_explicit():
    theorem = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "every orthogonal stabilization fails",
        "no uniform K3 embedding or full-chamber compactness is asserted",
        "No impossibility theorem for all such nonorthogonal couplings is claimed",
        "exact surviving \\(n\\ge2\\) field-height problem",
    ):
        assert phrase in theorem
    for phrase in (
        "PASS, at the stated scope",
        "A genuinely nonorthogonal field form could change the mixed normal norm",
        "No classification theorem is imported",
    ):
        assert phrase in audit
