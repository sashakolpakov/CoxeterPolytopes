#!/usr/bin/env python3
"""No-project-import replay of the silver infinite-stabilizer theorem."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT / "data/candidates/silver_k3_norm2_vinberg_extension_6_plus_3sqrt2.json"
)
PRODUCTION = ROOT / "data/candidates/silver_k3_norm2_infinite_stabilizer.json"
REPLAY_ARTIFACT = (
    ROOT
    / "data/candidates/silver_k3_norm2_infinite_stabilizer_independent_replay.json"
)
EXPECTED_INPUT_FILE_SHA256 = (
    "2fff5cc20e1e97504071afb226cafa4c8732e39650b1611fa5bd02ab45537208"
)
EXPECTED_INPUT_CERTIFICATE_SHA256 = (
    "01c8c722b0db00f757d67789d82674f2e9724ec2e81e4b7e720870e1700a67f5"
)
EXPECTED_SUMMARY_SHA256 = (
    "edbd1429942e9baa6d0f7e76b345cc8ee086bb7c6baaa3b210c0c31995d6b497"
)


@dataclass(frozen=True)
class Q2:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __init__(self, a: int | Fraction = 0, b: int | Fraction = 0) -> None:
        object.__setattr__(self, "a", Fraction(a))
        object.__setattr__(self, "b", Fraction(b))

    @staticmethod
    def coerce(value: int | Fraction | "Q2") -> "Q2":
        return value if isinstance(value, Q2) else Q2(value)

    def __add__(self, other: int | Fraction | "Q2") -> "Q2":
        other = Q2.coerce(other)
        return Q2(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> "Q2":
        return Q2(-self.a, -self.b)

    def __sub__(self, other: int | Fraction | "Q2") -> "Q2":
        return self + (-Q2.coerce(other))

    def __rsub__(self, other: int | Fraction | "Q2") -> "Q2":
        return Q2.coerce(other) - self

    def __mul__(self, other: int | Fraction | "Q2") -> "Q2":
        other = Q2.coerce(other)
        return Q2(
            self.a * other.a + 2 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    __rmul__ = __mul__

    def __truediv__(self, other: int | Fraction | "Q2") -> "Q2":
        other = Q2.coerce(other)
        norm = other.a * other.a - 2 * other.b * other.b
        if not norm:
            raise ZeroDivisionError
        return Q2(
            (self.a * other.a - 2 * self.b * other.b) / norm,
            (self.b * other.a - self.a * other.b) / norm,
        )

    def sign(self, conjugate: bool = False) -> int:
        a, b = (self.a, -self.b) if conjugate else (self.a, self.b)
        if not a:
            return (b > 0) - (b < 0)
        if not b or (a > 0 and b > 0) or (a < 0 and b < 0):
            return (a > 0) - (a < 0)
        comparison = a * a - 2 * b * b
        if a > 0:  # a - |b| sqrt(2)
            return (comparison > 0) - (comparison < 0)
        # -|a| + b sqrt(2)
        return (comparison < 0) - (comparison > 0)

    @property
    def integral(self) -> bool:
        return self.a.denominator == self.b.denominator == 1

    def to_json(self) -> dict[str, str]:
        return {"a": str(self.a), "b": str(self.b)}

    @staticmethod
    def from_json(value: object) -> "Q2":
        if not isinstance(value, dict):
            raise TypeError("malformed quadratic value")
        return Q2(Fraction(value["a"]), Fraction(value["b"]))


ZERO = Q2()
ONE = Q2(1)
S = Q2(0, 1)
Vector = tuple[Q2, ...]
Matrix = tuple[Vector, ...]


def canonical_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def identity(size: int = 10) -> Matrix:
    return tuple(
        tuple(ONE if row == column else ZERO for column in range(size))
        for row in range(size)
    )


def mat_vec(matrix: Matrix, vector: Vector) -> Vector:
    size = len(vector)
    return tuple(
        sum((matrix[row][column] * vector[column] for column in range(size)), ZERO)
        for row in range(size)
    )


def dot(left: Vector, right: Vector) -> Q2:
    if len(left) != len(right):
        raise ValueError("hostile replay dot-product dimension mismatch")
    return sum((left[index] * right[index] for index in range(len(left))), ZERO)


def compose(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum((left[row][middle] * right[middle][column] for middle in range(10)), ZERO)
            for column in range(10)
        )
        for row in range(10)
    )


def transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[row][column] for row in range(10)) for column in range(10))


def inverse(matrix: Matrix) -> Matrix:
    work = [list(matrix[row]) + list(identity()[row]) for row in range(10)]
    for column in range(10):
        pivot = next(row for row in range(column, 10) if work[row][column] != ZERO)
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [entry / scale for entry in work[column]]
        for row in range(10):
            if row == column or work[row][column] == ZERO:
                continue
            scale = work[row][column]
            work[row] = [
                work[row][index] - scale * work[column][index] for index in range(20)
            ]
    return tuple(tuple(row[10:]) for row in work)


def quadratic(vector: Vector, gram: Matrix) -> Q2:
    return dot(vector, mat_vec(gram, vector))


def matrix_digest(matrix: Matrix) -> str:
    return canonical_digest([[entry.to_json() for entry in row] for row in matrix])


def ldlt_pivots(matrix: Matrix) -> tuple[Q2, ...]:
    """Exact field-form pivots, independent of project inertia routines."""

    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise AssertionError("hostile replay field form is not square")
    if any(matrix[i][j] != matrix[j][i] for i in range(size) for j in range(i)):
        raise AssertionError("hostile replay field form is not symmetric")
    work = [list(row) for row in matrix]
    output: list[Q2] = []
    for index in range(size):
        pivot = work[index][index]
        if pivot == ZERO:
            raise AssertionError("hostile replay exact LDL decomposition hit zero")
        output.append(pivot)
        updated = {
            (row, column): (
                work[row][column]
                - work[row][index] * work[index][column] / pivot
            )
            for row in range(index + 1, size)
            for column in range(index + 1, size)
        }
        for (row, column), value in updated.items():
            work[row][column] = value
    return tuple(output)


def signature(pivots: tuple[Q2, ...], *, conjugate: bool) -> tuple[int, int, int]:
    signs = tuple(pivot.sign(conjugate=conjugate) for pivot in pivots)
    return signs.count(1), signs.count(-1), signs.count(0)


def fraction_inverse(matrix: tuple[tuple[Fraction, ...], ...]) -> tuple[tuple[Fraction, ...], ...]:
    size = len(matrix)
    work = [
        [*matrix[row], *(Fraction(row == column) for column in range(size))]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [entry / scale for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                work[row][index] - scale * work[column][index]
                for index in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in work)


def fraction_quadratic(
    vector: tuple[int, ...], matrix: tuple[tuple[Fraction, ...], ...]
) -> Fraction:
    return sum(
        Fraction(vector[i]) * matrix[i][j] * Fraction(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def fraction_prefix_schur(
    matrix: tuple[tuple[Fraction, ...], ...]
) -> tuple[tuple[tuple[Fraction, ...], ...], ...]:
    size = len(matrix)
    output = []
    for prefix in range(1, size + 1):
        leading = tuple(tuple(matrix[i][j] for j in range(prefix)) for i in range(prefix))
        if prefix == size:
            output.append(leading)
            continue
        cross = tuple(
            tuple(matrix[i][j] for j in range(prefix, size))
            for i in range(prefix)
        )
        tail = tuple(
            tuple(matrix[i][j] for j in range(prefix, size))
            for i in range(prefix, size)
        )
        tail_inverse = fraction_inverse(tail)
        output.append(
            tuple(
                tuple(
                    leading[i][j]
                    - sum(
                        cross[i][a] * tail_inverse[a][b] * cross[j][b]
                        for a in range(size - prefix)
                        for b in range(size - prefix)
                    )
                    for j in range(prefix)
                )
                for i in range(prefix)
            )
        )
    return tuple(output)


def floor_sqrt_fraction(value: Fraction) -> int:
    candidate = isqrt(value.numerator // value.denominator)
    while Fraction((candidate + 1) ** 2) <= value:
        candidate += 1
    while Fraction(candidate**2) > value:
        candidate -= 1
    return candidate


def integer_short_vectors(
    gram: tuple[tuple[int, ...], ...], target: int = 2
) -> tuple[tuple[int, ...], ...]:
    """Exhaust a positive integral form using exact Schur-complement pruning."""

    rational = tuple(tuple(Fraction(entry) for entry in row) for row in gram)
    inverse_gram = fraction_inverse(rational)
    bounds = tuple(
        floor_sqrt_fraction(Fraction(target) * inverse_gram[index][index])
        for index in range(len(gram))
    )
    prefixes = fraction_prefix_schur(rational)
    vectors: list[tuple[int, ...]] = []

    def extend(prefix: tuple[int, ...]) -> None:
        depth = len(prefix)
        if depth and fraction_quadratic(prefix, prefixes[depth - 1]) > target:
            return
        if depth == len(gram):
            if fraction_quadratic(prefix, prefixes[-1]) == target:
                vectors.append(prefix)
            return
        for entry in range(-bounds[depth], bounds[depth] + 1):
            extend((*prefix, entry))

    extend(())
    return tuple(vectors)


def positive_trace_gram(block: Matrix) -> tuple[tuple[int, ...], ...]:
    """Gram of ``-Tr((-1/2-sqrt(2)/4)h)`` in the interleaved ring basis."""

    scalar = Q2(Fraction(-1, 2), Fraction(-1, 4))
    basis = (ONE, S)
    size = len(block)
    rows: list[tuple[int, ...]] = []
    for i in range(size):
        for basis_i in basis:
            row = []
            for j in range(size):
                for basis_j in basis:
                    value = -(scalar * basis_i * block[i][j] * basis_j).a * 2
                    if value.denominator != 1:
                        raise AssertionError("hostile replay trace form is nonintegral")
                    row.append(value.numerator)
            rows.append(tuple(row))
    return tuple(rows)


def height_zero_exhaustion(
    gram: Matrix, archived_roots: tuple[Vector, ...]
) -> tuple[tuple[int, ...], tuple[int, ...], set[Vector]]:
    """Independently exhaust the finite norm-two subsystem at ``x=0``."""

    intervals = ((1, 5), (5, 9), (9, 10))
    trace_counts: list[int] = []
    field_counts: list[int] = []
    exhaustive: set[Vector] = set()
    for start, stop in intervals:
        block = tuple(tuple(gram[i][j] for j in range(start, stop)) for i in range(start, stop))
        trace_gram = positive_trace_gram(block)
        if any(trace_gram[index][index] % 2 for index in range(len(trace_gram))):
            raise AssertionError("hostile replay height-zero trace block is not even")
        short = integer_short_vectors(trace_gram)
        trace_counts.append(len(short))
        field_roots = tuple(
            tuple(Q2(vector[2 * i], vector[2 * i + 1]) for i in range(stop - start))
            for vector in short
        )
        field_roots = tuple(root for root in field_roots if quadratic(root, block) == Q2(2))
        field_counts.append(len(field_roots))
        for local in field_roots:
            vector = [ZERO for _ in range(10)]
            vector[start:stop] = local
            exhaustive.add(tuple(vector))

    simple = archived_roots[:9]
    closure = set(simple)
    closure.update(tuple(-entry for entry in root) for root in simple)
    frontier = list(closure)
    while frontier:
        vector = frontier.pop()
        for root in simple:
            pairing = dot(vector, mat_vec(gram, root))
            reflected = tuple(
                vector[index] - pairing * root[index]
                for index in range(10)
            )
            if reflected not in closure:
                closure.add(reflected)
                frontier.append(reflected)
    if closure != exhaustive:
        raise AssertionError("hostile replay height-zero closure is not exhaustive")
    return tuple(trace_counts), tuple(field_counts), closure


def field_gram() -> Matrix:
    z = ZERO
    return (
        (Q2(2), Q2(-1), z, z, -S, z, z, z, z, z),
        (Q2(-1), Q2(2), -S, z, Q2(-1, 1), z, z, z, z, z),
        (z, -S, Q2(2), Q2(-1), Q2(-1, 1), z, z, z, z, z),
        (z, z, Q2(-1), Q2(2), -S, z, z, z, z, z),
        (-S, Q2(-1, 1), Q2(-1, 1), -S, Q2(3, -1), z, z, z, z, z),
        (z, z, z, z, z, Q2(2), Q2(-1), z, z, z),
        (z, z, z, z, z, Q2(-1), Q2(2), -S, z, z),
        (z, z, z, z, z, z, -S, Q2(2), Q2(-1), z),
        (z, z, z, z, z, z, z, Q2(-1), Q2(2), z),
        (z, z, z, z, z, z, z, z, z, Q2(2)),
    )


def control_vector() -> Vector:
    return (Q2(-1, 1), Q2(-1, 1), Q2(2, -1), ONE, ONE, ZERO, ZERO, ZERO, ZERO, ZERO)


def columns(vectors: tuple[Vector, ...]) -> Matrix:
    return tuple(tuple(vectors[column][row] for column in range(10)) for row in range(10))


def seed_reversal(roots: tuple[Vector, ...]) -> Matrix:
    permutation = (1, 0, 9, 3, 4, 5, 6, 7, 8, 2)
    return compose(
        columns(tuple(roots[permutation[index]] for index in range(10))),
        inverse(columns(roots[:10])),
    )


def wall_swap(roots: tuple[Vector, ...], gram: Matrix) -> Matrix:
    difference = tuple(roots[4][index] - roots[18][index] for index in range(10))
    image = mat_vec(gram, difference)
    coefficient = Q2(2) / dot(difference, image)
    unit = identity()
    return tuple(
        tuple(
            unit[row][column] - coefficient * difference[row] * image[column]
            for column in range(10)
        )
        for row in range(10)
    )


def seed_witness() -> Vector:
    return (
        ONE, Q2(Fraction(147, 71), Fraction(-30, 71)),
        Q2(Fraction(5, 71), Fraction(112, 71)), Q2(1, 2),
        Q2(Fraction(141, 71), Fraction(77, 71)),
        Q2(Fraction(54, 71), Fraction(31, 71)),
        Q2(Fraction(96, 71), Fraction(63, 71)),
        Q2(Fraction(96, 71), Fraction(63, 71)),
        Q2(Fraction(54, 71), Fraction(31, 71)),
        Q2(Fraction(6, 71), Fraction(-1, 142)),
    )


def swap_boundary() -> Vector:
    return (
        ONE, ONE, S, Q2(1, 1), Q2(1, 1),
        Q2(Fraction(2, 3), Fraction(2, 3)), Q2(1, 1),
        Q2(Fraction(4, 3), Fraction(2, 3)),
        Q2(Fraction(2, 3), Fraction(1, 3)),
        Q2(Fraction(1, 3), Fraction(1, 6)),
    )


def swap_strict() -> Vector:
    return (
        ONE, Q2(Fraction(856, 713), Fraction(45, 713)),
        Q2(Fraction(166, 713), Fraction(850, 713)),
        Q2(Fraction(41, 31), Fraction(40, 31)),
        Q2(Fraction(955, 713), Fraction(844, 713)),
        Q2(Fraction(486, 713), Fraction(487, 713)),
        Q2(Fraction(882, 713), Fraction(831, 713)),
        Q2(Fraction(1140, 713), Fraction(623, 713)),
        Q2(Fraction(586, 713), Fraction(329, 713)),
        Q2(Fraction(226, 713), Fraction(227, 1426)),
    )


def separator_bound(control: Vector, witness: Vector, gram: Matrix) -> Q2:
    p2 = quadratic(control, gram)
    v2 = quadratic(witness, gram)
    pv = dot(control, mat_vec(gram, witness))
    return Q2(2) * (pv * pv + (-p2) * v2) / (-v2)


def charpoly(matrix: Matrix) -> tuple[Q2, ...]:
    work = identity()
    output = [ONE]
    for degree in range(1, 11):
        product = compose(matrix, work)
        coefficient = -sum((product[i][i] for i in range(10)), ZERO) / degree
        output.append(coefficient)
        work = tuple(
            tuple(
                product[row][column] + (coefficient if row == column else ZERO)
                for column in range(10)
            )
            for row in range(10)
        )
    return tuple(output)


def admissible_below_four() -> tuple[Q2, ...]:
    values = set()
    upper = Q2(4)
    conjugate_square_bound = Q2(2, 2)
    for a in range(-8, 9):
        for b in range(-8, 9):
            value = Q2(a, b)
            if (
                value.sign() > 0
                and (upper - value).sign() >= 0
                and (conjugate_square_bound - Q2(a, -b) * Q2(a, -b)).sign() >= 0
            ):
                values.add(value)
    def compare(left: Q2, right: Q2) -> int:
        return (left - right).sign()
    return tuple(sorted(values, key=cmp_to_key(compare)))


def build_summary() -> dict[str, object]:
    if sha256(INPUT.read_bytes()).hexdigest() != EXPECTED_INPUT_FILE_SHA256:
        raise AssertionError("hostile replay input file changed")
    document = json.loads(INPUT.read_text(encoding="utf-8"))
    if document.get("certificate_sha256") != EXPECTED_INPUT_CERTIFICATE_SHA256:
        raise AssertionError("hostile replay input certificate changed")
    roots = tuple(
        tuple(Q2.from_json(entry) for entry in row)
        for row in document["ordered_prefix_roots"]
    )
    if len(roots) != 1291:
        raise AssertionError("hostile replay root count changed")
    gram = field_gram()
    pivots = ldlt_pivots(gram)
    identity_signature = signature(pivots, conjugate=False)
    conjugate_signature = signature(pivots, conjugate=True)
    determinant = ONE
    for pivot in pivots:
        determinant *= pivot
    field_anisotropic = conjugate_signature == (10, 0, 0)
    if not (
        len(pivots) == 10
        and identity_signature == (9, 1, 0)
        and conjugate_signature == (10, 0, 0)
        and determinant == Q2(-2, -2)
        and field_anisotropic
    ):
        raise AssertionError("hostile replay rejected the rank-ten field form")
    root_images = tuple(mat_vec(gram, root) for root in roots)
    if any(quadratic(root, gram) != Q2(2) for root in roots):
        raise AssertionError("hostile replay found a non-root")
    trace_counts, field_counts, zero_roots = height_zero_exhaustion(gram, roots)
    if not (
        trace_counts == (240, 240, 4)
        and field_counts == (32, 48, 2)
        and len(zero_roots) == 82
        and all(root[0] == ZERO for root in roots[:9])
        and all(root[0].sign() > 0 for root in roots[9:])
        and all(
            dot(roots[i], root_images[j]).sign() <= 0
            for i in range(9)
            for j in range(i)
        )
    ):
        raise AssertionError("hostile replay rejected height-zero completeness")

    g = seed_reversal(roots)
    h = wall_swap(roots, gram)
    unit = identity()
    if not (
        all(entry.integral for row in g for entry in row)
        and all(entry.integral for row in h for entry in row)
        and compose(transpose(g), compose(gram, g)) == gram
        and compose(transpose(h), compose(gram, h)) == gram
        and compose(g, g) == unit
        and compose(h, h) == unit
    ):
        raise AssertionError("hostile replay rejected the two involutions")

    control = control_vector()
    vg = seed_witness()
    vh = tuple((left + right) / 2 for left, right in zip(swap_boundary(), swap_strict(), strict=True))
    pairings_g = tuple(dot(vg, image) for image in root_images)
    pairings_h = tuple(dot(vh, image) for image in root_images)
    bound_g = separator_bound(control, vg, gram)
    bound_h = separator_bound(control, vh, gram)
    future_g = (
        mat_vec(g, vg) == vg
        and quadratic(vg, gram).sign() < 0
        and dot(control, mat_vec(gram, vg)).sign() < 0
    )
    future_h = (
        mat_vec(h, vh) == vh
        and quadratic(vh, gram).sign() < 0
        and dot(control, mat_vec(gram, vh)).sign() < 0
    )
    boundary = swap_boundary()
    strict = swap_strict()
    boundary_pairings = tuple(dot(boundary, image) for image in root_images)
    strict_pairings = tuple(dot(strict, image) for image in root_images)
    boundary_zero_set = tuple(
        index for index, pairing in enumerate(boundary_pairings) if pairing == ZERO
    )
    if not (
        future_g
        and future_h
        and all(value.sign() < 0 for value in pairings_g)
        and all(value.sign() < 0 for value in pairings_h)
        and mat_vec(h, boundary) == boundary
        and quadratic(boundary, gram) == Q2(0, Fraction(-1, 3))
        and all(value.sign() <= 0 for value in boundary_pairings)
        and boundary_zero_set == (0, 1, 2, 3, 5, 6, 7)
        and mat_vec(h, strict) == strict
        and all(value.sign() < 0 for value in strict_pairings)
        and (Q2(16) - bound_g).sign() > 0
        and (Q2(9) - bound_h).sign() > 0
    ):
        raise AssertionError("hostile replay rejected a complete-chamber witness")

    short_coordinates = admissible_below_four()
    swap_coordinates = tuple(value for value in short_coordinates if (Q2(3) - value).sign() >= 0)
    archived_coordinates = tuple(Q2.from_json(value) for value in document["candidate_coordinates"])
    if not (
        short_coordinates == archived_coordinates[: len(short_coordinates)]
        and swap_coordinates == archived_coordinates[: len(swap_coordinates)]
    ):
        raise AssertionError("hostile replay found an unprocessed short coordinate")

    product = compose(g, h)
    trace = sum((product[index][index] for index in range(10)), ZERO)
    coefficients = charpoly(product)
    expected = (
        Q2(1), Q2(-10, -4), Q2(39, 24), Q2(-76, -56), Q2(70, 56),
        Q2(0), Q2(-70, -56), Q2(76, 56), Q2(-39, -24), Q2(10, 4), Q2(-1),
    )
    if trace != Q2(10, 4) or coefficients != expected:
        raise AssertionError("hostile replay rejected the infinite-order product")

    summary: dict[str, object] = {
        "kind": "silver_k3_norm2_infinite_stabilizer_no_project_import_replay",
        "input_wall_count": len(roots),
        "all_roots_norm_two": True,
        "field_form_ldlt_pivots": [value.to_json() for value in pivots],
        "field_form_determinant": determinant.to_json(),
        "field_form_identity_signature": list(identity_signature),
        "field_form_conjugate_signature": list(conjugate_signature),
        "field_anisotropic_over_Qsqrt2_derived": field_anisotropic,
        "height_zero_trace_short_vector_counts": list(trace_counts),
        "height_zero_field_norm_two_counts": list(field_counts),
        "height_zero_root_count": len(zero_roots),
        "height_zero_closure_exhaustive": True,
        "short_admissible_coordinates": [value.to_json() for value in short_coordinates],
        "seed_matrix_sha256": matrix_digest(g),
        "seed_witness_norm": quadratic(vg, gram).to_json(),
        "seed_largest_wall_pairing": max(pairings_g, key=cmp_to_key(lambda a, b: (a - b).sign())).to_json(),
        "seed_separator_bound_squared": bound_g.to_json(),
        "seed_preserves_future_component": future_g,
        "swap_matrix_sha256": matrix_digest(h),
        "swap_witness_norm": quadratic(vh, gram).to_json(),
        "swap_largest_wall_pairing": max(pairings_h, key=cmp_to_key(lambda a, b: (a - b).sign())).to_json(),
        "swap_separator_bound_squared": bound_h.to_json(),
        "swap_preserves_future_component": future_h,
        "swap_boundary_zero_wall_indices": list(boundary_zero_set),
        "product_matrix_sha256": matrix_digest(product),
        "product_trace": trace.to_json(),
        "product_charpoly_coefficients": [value.to_json() for value in coefficients],
        "separator_hypothesis_includes_equality": True,
        "complete_chamber_symmetries": True,
        "infinite_order_product": True,
    }
    production = json.loads(PRODUCTION.read_text(encoding="utf-8")) if PRODUCTION.exists() else None
    if production is not None and not (
        production["seed_C5_reversal"]["matrix_sha256"] == summary["seed_matrix_sha256"]
        and production["wall_4_18_swap"]["matrix_sha256"] == summary["swap_matrix_sha256"]
        and production["infinite_order_product"]["matrix_sha256"] == summary["product_matrix_sha256"]
        and production["infinite_order_product"]["trace"] == summary["product_trace"]
    ):
        raise AssertionError("hostile and production theorem cores disagree")
    summary["summary_sha256"] = canonical_digest(summary)
    if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN" and summary["summary_sha256"] != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("hostile replay summary changed")
    return summary


def verify_frozen() -> dict[str, object]:
    actual = build_summary()
    expected = json.loads(REPLAY_ARTIFACT.read_text(encoding="utf-8"))
    if actual != expected:
        raise AssertionError("hostile replay artifact changed")
    return actual


if __name__ == "__main__":
    print(json.dumps(verify_frozen(), indent=2, sort_keys=True))
