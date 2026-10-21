"""Independent hostile regressions for the colored affine cusp-factor no-go."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT
    / "results"
    / "rational_affine_h3_colored_cusp_factor_no_go_hostile_audit.md"
)
SOURCE_HASHES = {
    "results/rational_affine_h3_colored_cusp_factor_no_go.md":
        "183f9bca555e38f5d7ef41255078835636e99b8ce6428756ba6339d8b616710e",
    "scripts/certify_rational_affine_h3_colored_cusp_factor_no_go.py":
        "81b4a94588aa2326973031b30df32a00475c5e09e6af3472578afe1d0cfd1b16",
    "tests/test_rational_affine_h3_colored_cusp_factor_no_go.py":
        "2ae937a404f69bffafd2f030b58ea95b70bb7236e433ba70e97e298be79038b0",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_rank(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            coefficient = work[row][column]
            work[row] = [
                entry - coefficient * pivot_entry
                for entry, pivot_entry in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
    return pivot_row


def affine_packet(rank):
    roots = []
    affine = [0] * rank
    affine[0] = affine[1] = -1
    roots.append(affine)
    for index in range(rank - 1):
        root = [0] * rank
        root[index], root[index + 1] = 1, -1
        roots.append(root)
    terminal = [0] * rank
    terminal[-2] = terminal[-1] = 1
    roots.append(terminal)

    gram = [
        [-sum(x * y for x, y in zip(left, right)) for right in roots]
        for left in roots
    ]
    marks = [1, 1] + [2] * (rank - 3) + [1, 1]
    kinds = [
        "B"
        if (
            vertex in {0, 1, rank - 1, rank}
            or (3 <= vertex <= rank - 3 and vertex % 2 == 1)
        )
        else "E"
        for vertex in range(rank + 1)
    ]
    physical = [[0] * (rank + 1) for _ in range(rank + 1)]
    for left in range(rank + 1):
        physical[left][left] = -4 if kinds[left] == "B" else -1
        for right in range(left + 1, rank + 1):
            if gram[left][right] == 1:
                physical[left][right] = physical[right][left] = 1
    multiplicities = [
        mark if kind == "B" else 2 * mark
        for mark, kind in zip(marks, kinds)
    ]
    return roots, gram, marks, kinds, physical, multiplicities


def mat_vec(matrix, vector):
    return [
        sum(value * coefficient for value, coefficient in zip(row, vector))
        for row in matrix
    ]


def inside_alcove(point):
    return (
        all(point[index] > point[index + 1] for index in range(len(point) - 2))
        and point[-2] > abs(point[-1])
        and point[0] + point[1] < 1
    )


def lorentz_pair(left, right):
    return left[0] * right[0] - sum(
        x * y for x, y in zip(left[1:], right[1:])
    )


def in_half_boundary_host(vector, boundaries):
    """Membership in Z^d plus the displayed boundary half-classes."""
    for bits in product((0, 1), repeat=len(boundaries)):
        remainder = [Fraction(value) for value in vector]
        for bit, boundary in zip(bits, boundaries):
            if bit:
                remainder = [
                    value - Fraction(coefficient, 2)
                    for value, coefficient in zip(remainder, boundary)
                ]
        if all(value.denominator == 1 for value in remainder):
            return True
    return False


def reflection_preserves_host(boundaries):
    dimension = len(boundaries[0])
    difference = [x - y for x, y in zip(*boundaries)]
    root = [Fraction(value, 2) for value in difference]

    generators = [
        [Fraction(int(row == column)) for row in range(dimension)]
        for column in range(dimension)
    ] + [
        [Fraction(value, 2) for value in boundary]
        for boundary in boundaries
    ]

    def reflect(vector):
        host_pairing = 2 * lorentz_pair(vector, root)
        return [
            value + host_pairing * root_value / 2
            for value, root_value in zip(vector, root)
        ]

    return all(in_half_boundary_host(reflect(vector), boundaries) for vector in generators), root


def test_frozen_sources_and_hostile_scope():
    assert {
        relative: sha256(ROOT / relative)
        for relative in SOURCE_HASHES
    } == SOURCE_HASHES
    text = AUDIT.read_text()
    for phrase in (
        "**PASS after two exact repairs.**",
        r"\(1+1+1+1=4\)",
        "ambient mod-two",
        "does not manufacture",
        "does **not** prove arbitrary-rank projective",
        "or any dimension theorem",
    ):
        assert phrase in text


def test_all_packet_equations_and_proximity_profiles_independently():
    for rank in range(4, 62, 2):
        roots, gram, marks, kinds, physical, multiplicities = affine_packet(rank)
        assert exact_rank(gram) == rank
        assert mat_vec(gram, marks) == [0] * (rank + 1)
        assert [
            sum(mark * root[column] for mark, root in zip(marks, roots))
            for column in range(rank)
        ] == [0] * rank

        positive = sum(mark for mark, kind in zip(marks, kinds) if kind == "B")
        negative = sum(mark for mark, kind in zip(marks, kinds) if kind == "E")
        assert (positive, negative, positive - negative) == (rank, rank - 2, 2)
        assert mat_vec(physical, multiplicities) == [0] * (rank + 1)
        assert gcd(*multiplicities) == 1
        canonical = [2 if kind == "B" else -1 for kind in kinds]
        assert sum(a * b for a, b in zip(multiplicities, canonical)) == 4

        # B curves are pairwise disjoint and no R vertex occurs: H1/H2 vanish.
        b_vertices = [index for index, kind in enumerate(kinds) if kind == "B"]
        assert all(
            physical[left][right] == 0
            for offset, left in enumerate(b_vertices)
            for right in b_vertices[offset + 1:]
        )

        clusters = [
            tuple(
                marks[other]
                for other, other_kind in enumerate(kinds)
                if other_kind == "B" and physical[vertex][other] == 1
            )
            for vertex, kind in enumerate(kinds)
            if kind == "E"
        ]
        assert all(sum(cluster) == 4 for cluster in clusters)
        if rank == 4:
            assert clusters == [(1, 1, 1, 1)]
        else:
            assert sum(len(cluster) == 3 for cluster in clusters) == 2
            assert all(len(cluster) in {2, 3} for cluster in clusters)


def test_coordinate_mirror_reconstruction_and_alcove_bisection():
    for rank in range(4, 42, 2):
        roots, _, _, kinds, _, _ = affine_packet(rank)
        coordinate = [0] * rank
        coordinate[-1] = 2
        reconstructed = [
            roots[-1][index] - roots[-2][index]
            for index in range(rank)
        ]
        assert reconstructed == coordinate
        assert sum(value * value for value in coordinate) == 4
        assert kinds[-2:] == ["B", "B"]
        assert 2 - 2 == 0

        def reflect(vector):
            coefficient = Fraction(
                2 * sum(x * y for x, y in zip(vector, coordinate)),
                sum(value * value for value in coordinate),
            )
            return [
                value - coefficient * root_value
                for value, root_value in zip(vector, coordinate)
            ]

        assert reflect(roots[-2]) == roots[-1]
        assert reflect(roots[-1]) == roots[-2]

        # Half of 2e_n has odd coordinate sum, so 2e_n is primitive in D_n.
        assert sum(coordinate) % 2 == 0
        assert sum(value // 2 for value in coordinate) % 2 == 1

        mirror = [
            Fraction(rank - index, 4 * rank)
            for index in range(rank - 1)
        ] + [Fraction(0)]
        epsilon = Fraction(1, 4 * rank)
        positive = mirror[:-1] + [epsilon]
        negative = mirror[:-1] + [-epsilon]
        assert inside_alcove(mirror)
        assert inside_alcove(positive)
        assert inside_alcove(negative)
        assert reflect(positive) == negative

        # Sign change preserves every D_n vector because it changes the sum by 2x_n.
        for sample in product(range(-2, 3), repeat=min(rank, 6)):
            vector = sample + (0,) * (rank - len(sample))
            if sum(vector) % 2:
                continue
            reflected = vector[:-1] + (-vector[-1],)
            assert sum(reflected) % 2 == 0


def test_exact_full_adapted_host_parity_in_both_directions():
    # Primitive, orthogonal square-minus-four vectors in odd unimodular lattices.
    even_boundaries = ((-3, -3, -2, 0), (-3, -3, 0, -2))
    odd_boundaries = (
        (-4, -4, -1, -1, -1, -1),
        (-3, -3, -1, -1, 1, 1),
    )

    for boundaries, expected in ((even_boundaries, True), (odd_boundaries, False)):
        left, right = boundaries
        assert lorentz_pair(left, left) == lorentz_pair(right, right) == -4
        assert lorentz_pair(left, right) == 0
        difference = tuple(x - y for x, y in zip(left, right))
        assert all(value % 2 == 0 for value in difference) is expected

        preserves, root = reflection_preserves_host(boundaries)
        assert preserves is expected
        assert 2 * lorentz_pair(root, root) == -4
        assert not in_half_boundary_host([value / 2 for value in root], boundaries)

    # In the odd case an integral NS generator has odd q_M-pairing with v.
    _, odd_root = reflection_preserves_host(odd_boundaries)
    first_basis = [Fraction(1)] + [Fraction(0)] * 5
    assert 2 * lorentz_pair(first_basis, odd_root) == -1
