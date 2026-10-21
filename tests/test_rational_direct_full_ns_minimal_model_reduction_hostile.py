"""Independent hostile checks for the direct full-NS MMP reduction."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT
    / "results"
    / "rational_direct_full_ns_minimal_model_reduction_hostile_audit.md"
)
SOURCE_HASHES = {
    "results/rational_direct_full_ns_minimal_model_reduction.md":
        "de94807a390d6d8ec7238b16b00deb304077cd942c757aa96dc1cba529b32abc",
    "scripts/certify_rational_direct_full_ns_minimal_model_reduction.py":
        "d70f4d7df8a08aa7ee219b3e71eb0382269896aad55eaf055d6faacf770f444e",
    "tests/test_rational_direct_full_ns_minimal_model_reduction.py":
        "305c8a9a9b0bdef6ac6de34e649a7048a97e2a3b3b142c4f2210cb34a4766d87",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ruled_pair(left, right, e: int):
    return (
        -e * left[0] * right[0]
        + left[0] * right[1]
        + left[1] * right[0]
        - sum(a * b for a, b in zip(left[2:], right[2:]))
    )


def reflect(vector, root, e: int):
    scalar = Fraction(-2 * ruled_pair(vector, root, e), ruled_pair(root, root, e))
    return tuple(
        Fraction(value) + scalar * root_value
        for value, root_value in zip(vector, root)
    )


def basis(rank: int):
    return tuple(
        tuple(int(row == column) for row in range(rank))
        for column in range(rank)
    )


def test_producer_sources_are_the_frozen_bytes():
    assert {
        relative: sha256(ROOT / relative)
        for relative in SOURCE_HASHES
    } == SOURCE_HASHES


def test_F0_swap_is_an_integral_isometry_in_every_exceptional_rank():
    for exceptional_rank in range(0, 41):
        rank = exceptional_rank + 2
        vectors = basis(rank)
        c, d = vectors[:2]
        root = tuple(a - b for a, b in zip(c, d))
        assert ruled_pair(root, root, 0) == -2
        images = tuple(reflect(vector, root, 0) for vector in vectors)
        assert images[0] == d
        assert images[1] == c
        assert images[2:] == vectors[2:]
        assert all(value.denominator == 1 for image in images for value in image)
        assert [
            [ruled_pair(left, right, 0) for right in images]
            for left in images
        ] == [
            [ruled_pair(left, right, 0) for right in vectors]
            for left in vectors
        ]

        # Every exceptional-supported divisor, hence any relative ample
        # choice, is orthogonal to the swap root.  So is c+d.
        exceptional_part = (0, 0, *[-(index + 1) for index in range(exceptional_rank)])
        ample_shape = tuple(
            97 * (a + b) + value
            for a, b, value in zip(c, d, exceptional_part)
        )
        assert ruled_pair(ample_shape, root, 0) == 0


def test_strict_section_is_primitive_and_integral_exactly_when_e_plus_t_divides_two():
    for e in range(2, 18):
        for hits in range(0, 14):
            extra = 3
            rank = 2 + hits + extra
            section = (1, 0, *[-1] * hits, *[0] * extra)
            fibre = (0, 1, *[0] * (hits + extra))
            n = e + hits
            assert math.gcd(*(abs(value) for value in section)) == 1
            assert ruled_pair(section, section, e) == -n
            assert ruled_pair(fibre, section, e) == 1
            image = reflect(fibre, section, e)
            assert image == tuple(
                Fraction(value) + Fraction(2, n) * root_value
                for value, root_value in zip(fibre, section)
            )
            assert all(value.denominator == 1 for value in image) is (n in (1, 2))

            images = tuple(reflect(vector, section, e) for vector in basis(rank))
            integral = all(
                value.denominator == 1 for row in images for value in row
            )
            assert integral is (n in (1, 2))


def test_exponent_two_section_bound_and_ruling_parity_are_exactly_scoped():
    # From (2/n)v in (1/2)L and primitive v, exactly n|4 is possible.
    for n in range(1, 65):
        coefficient = Fraction(2, n)
        belongs_to_half_lattice = (2 * coefficient).denominator == 1
        assert belongs_to_half_lattice is (4 % n == 0)

    # On F0, the shift of B/2 is (B.r/2)r.  Odd pairing would force
    # r/2 into an integral doubled host, where it pairs half-integrally with
    # B/2; even pairing gives an integral shift.
    root = (1, -1, 0, 0)
    boundaries = (
        (3, 1, 0, 0),   # B.r=-2
        (2, 1, 0, 0),   # B.r=-1
        (-1, 3, 2, -4), # B.r=4
    )
    for boundary in boundaries:
        pairing = ruled_pair(boundary, root, 0)
        shift = tuple(Fraction(pairing, 2) * value for value in root)
        assert all(value.denominator == 1 for value in shift) is (pairing % 2 == 0)
        if pairing % 2:
            doubled_host_pairing = Fraction(pairing, 2)
            assert doubled_host_pairing.denominator == 2


def test_hostile_report_records_the_gain_and_the_nonclaims():
    text = AUDIT.read_text(encoding="utf-8")
    required = [
        "**PASS.**",
        "no missing minimal base",
        "infinitely-near centres",
        r"exactly when \(n\mid2\)",
        r"\(n\mid4\)",
        "only necessary",
        "No depth bound",
        "adapted-host minimal-base classification",
    ]
    for phrase in required:
        assert phrase in text
    forbidden = set(range(0, 9)) | {11, 12} | set(range(14, 32)) | {127}
    assert not any(byte in forbidden for byte in AUDIT.read_bytes())
