"""Independent hostile checks for the genus-three vertical-tail no-go."""

from __future__ import annotations

from collections import deque
import hashlib
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT
    / "results"
    / "rational_genus3_vertical_blowup_invariant_no_go_hostile_audit.md"
)
SOURCE_HASHES = {
    "results/rational_genus3_vertical_blowup_invariant_no_go.md":
        "a2ee22bea4b5e9ccfafcc2ad2e99f12da9d2d26162c5e2e9c10ef58d7907adbb",
    "scripts/certify_rational_genus3_vertical_blowup_invariant_no_go.py":
        "6a8d5bafb99ed1f279f48716ee16bacf6bc5d2d68ed7231547a6153760651a06",
    "tests/test_rational_genus3_vertical_blowup_invariant_no_go.py":
        "9f45165f365e24c08ffd44bb61cfe2f51c2fa169054a8fa8382f9454fa2ccc9b",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def negative_pair(left, right):
    return -sum(a * b for a, b in zip(left, right))


def tail_roots(m: int):
    roots = []
    for index in range(m - 1):
        root = [0] * m
        root[index], root[index + 1] = 1, -1
        roots.append(tuple(root))
    endpoint = [0] * m
    endpoint[-1] = 1
    roots.append(tuple(endpoint))
    return tuple(roots)


def apply_simple_generator(frame, generator: int):
    """Left-compose a signed-permutation frame by a B_m generator."""

    result = []
    m = len(frame)
    for image in frame:
        sign = 1 if image > 0 else -1
        coordinate = abs(image) - 1
        if generator < m - 1:
            if coordinate == generator:
                coordinate = generator + 1
            elif coordinate == generator + 1:
                coordinate = generator
        elif coordinate == m - 1:
            sign *= -1
        result.append(sign * (coordinate + 1))
    return tuple(result)


def signed_permutation_group(m: int):
    identity = tuple(range(1, m + 1))
    found = {identity}
    queue = deque([identity])
    while queue:
        frame = queue.popleft()
        for generator in range(m):
            image = apply_simple_generator(frame, generator)
            if image not in found:
                found.add(image)
                queue.append(image)
    return found


def test_producer_sources_are_the_repaired_frozen_bytes():
    assert {
        relative: sha256(ROOT / relative)
        for relative in SOURCE_HASHES
    } == SOURCE_HASHES


def test_relative_noether_translation_is_symbolically_exact():
    for rho in range(1, 81):
        k_f_square = 26 - rho
        e_f = rho + 10
        assert k_f_square + e_f == 36 == 12 * 3
        for blowups in (1, 2, 7, 29, 101):
            new_k = 26 - (rho + blowups)
            new_e = rho + blowups + 10
            assert new_k - k_f_square == -blowups
            assert new_e - e_f == blowups
            assert new_k + new_e == 36


def test_physical_tail_gram_and_canonical_colors_include_the_mixed_curve():
    for m in range(1, 31):
        # Restricted orthogonal coordinates (C0,e1,...,em), all square -1.
        canonical = (1, *([1] * m))
        mixed = (1, -1, *([0] * (m - 1)))
        pure = tail_roots(m)
        curves = (mixed, *[(0, *root) for root in pure])
        assert [negative_pair(curve, curve) for curve in curves] == [
            *([-2] * m),
            -1,
        ]
        assert [negative_pair(canonical, curve) for curve in curves] == [
            *([0] * m),
            -1,
        ]
        gram = [
            [negative_pair(left, right) for right in curves]
            for left in curves
        ]
        assert all(gram[index][index + 1] == 1 for index in range(m))
        assert all(
            gram[left][right] == 0
            for left in range(m + 1)
            for right in range(m + 1)
            if abs(left - right) > 1
        )


def test_tail_reflections_generate_the_full_signed_permutation_group():
    for m in range(1, 6):
        group = signed_permutation_group(m)
        assert len(group) == (2 ** m) * math.factorial(m)
        assert all(
            sorted(abs(value) for value in frame) == list(range(1, m + 1))
            for frame in group
        )


def test_mixed_projection_orbit_is_exactly_all_signed_coordinates():
    for m in range(1, 65):
        orbit = {-1}
        queue = deque([-1])
        while queue:
            value = queue.popleft()
            sign = 1 if value > 0 else -1
            coordinate = abs(value) - 1
            for generator in range(m):
                new_sign = sign
                new_coordinate = coordinate
                if generator < m - 1:
                    if coordinate == generator:
                        new_coordinate = generator + 1
                    elif coordinate == generator + 1:
                        new_coordinate = generator
                elif coordinate == m - 1:
                    new_sign *= -1
                image = new_sign * (new_coordinate + 1)
                if image not in orbit:
                    orbit.add(image)
                    queue.append(image)
        assert orbit == set(range(1, m + 1)) | set(range(-m, 0))
        # The corresponding coordinate vectors contain the standard basis.
        assert len({abs(value) for value in orbit}) == m


def test_stable_scope_and_counterfamily_nonclaims_are_explicit():
    producer = (ROOT / "results/rational_genus3_vertical_blowup_invariant_no_go.md").read_text(
        encoding="utf-8"
    )
    assert "a smooth free centre gives a tail" in producer
    assert "a nodal centre gives a bridge" in producer
    assert "the stabilized dual graph" in producer
    assert "does not prove that the full\nMori cone" in producer
    assert "reflection-maximal" in producer

    text = AUDIT.read_text(encoding="utf-8")
    required = [
        "**PASS after one scope-precision repair.**",
        "relative canonical algebras",
        "stabilized fibre data",
        r"\(2^m m!\) signed permutation matrices",
        "finite-volume attachment checks only",
        "not advertised as\na universal fibre configuration",
        "genuinely global",
    ]
    for phrase in required:
        assert phrase in text
    forbidden = set(range(0, 9)) | {11, 12} | set(range(14, 32)) | {127}
    assert not any(byte in forbidden for byte in AUDIT.read_bytes())
