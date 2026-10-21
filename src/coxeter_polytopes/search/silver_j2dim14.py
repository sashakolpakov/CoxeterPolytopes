"""Exact arithmetic substrate for silver ``J2dim14`` hosts in ``H^28``.

``J2dim14`` is the extremal 5-modular rank-28 lattice of minimum eight in
the Nebe--Sloane catalogue.  This module keeps the catalogue lower triangle
literally and supplies the local data needed for Vinberg searches over
``Z[sqrt(2)]``.

The search routines deliberately return *necessary* binary Gram types.
Actual shell realization and same-height Vinberg acuteness are separate
checks; callers must not turn a nonempty necessary packet into an accepted
wall claim without those checks.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import isqrt
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse


IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]

J2DIM14_SOURCE_URL = (
    "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/J2dim14.html"
)
J2DIM14_RANK = 28
J2DIM14_DETERMINANT = 5**14
J2DIM14_SMITH_INVARIANTS = (1,) * 14 + (5,) * 14
J2DIM14_MINIMUM = 8
J2DIM14_KISSING_NUMBER = 23_310
J2DIM14_NORM_TEN_SHELL_COUNT = 374_976
J2DIM14_GROUP_NAME = "C2xJ2"

# Lines 44--71 of the primary catalogue page, stored as its lower triangle.
J2DIM14_LOWER_TRIANGLE: tuple[tuple[int, ...], ...] = (
    (8,),
    (-4, 8),
    (-4, 0, 8),
    (3, 1, -4, 8),
    (-1, 3, 2, 2, 8),
    (3, -1, -4, 0, -4, 8),
    (2, 0, -4, 4, 2, 2, 8),
    (2, 2, -4, 2, 2, 0, 4, 8),
    (-3, -1, 0, 1, -3, -1, -2, -4, 10),
    (-4, 2, 1, -1, -2, -2, -2, -2, 5, 8),
    (2, -4, 1, -1, 1, -1, -1, -1, -2, -4, 8),
    (-2, 4, 0, 2, 0, -2, -1, 0, 1, 4, -2, 8),
    (-4, 2, 1, 0, -1, -3, -1, 0, 2, 2, 0, 4, 8),
    (-4, 2, 0, 1, -1, -1, 0, 0, 3, 4, -2, 2, 4, 8),
    (2, -2, 1, 2, -1, -1, -1, 0, -1, 0, 0, 0, -1, 0, 8),
    (-4, 4, 2, -3, 3, -1, 0, 2, -2, 1, -2, 0, 1, 2, -2, 8),
    (-4, 2, 2, -2, 2, -4, 0, 2, -1, 2, -1, 2, 2, 2, -1, 2, 8),
    (-1, 3, 2, 2, 2, -2, -2, 0, -1, -1, 0, 3, 1, -1, 3, -1, 0, 8),
    (4, -2, -2, 4, -2, 0, -1, 1, 0, -2, 0, 0, 1, 1, 5, -3, -2, 1, 10),
    (2, -4, 0, -3, -4, 4, -2, -1, 1, -2, 1, -3, -2, -3, 1, -2, -3, 0, 0, 10),
    (-4, 4, 1, 1, 1, -2, -1, 0, 3, 4, -4, 4, 2, 2, 1, 2, 2, 2, -1, -2, 8),
    (3, -1, -2, 4, 2, -1, 1, 2, 1, -2, 1, -2, -2, -1, 2, -1, -2, 0, 4, 0, 0, 8),
    (-4, 2, 2, 1, 3, -1, 2, -2, 3, 0, -2, 0, 2, 1, -2, 3, 0, 0, -2, -2, 3, 0, 10),
    (-2, -2, 2, -4, -4, 2, -4, -4, 2, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, 4, 0, -2, 0, 8),
    (0, -1, 0, -1, -3, 3, -1, -3, 3, 0, -2, -3, -2, 0, 0, 0, -3, -1, 1, 3, -2, 1, 1, 3, 8),
    (-1, -1, 0, 0, -2, 2, 2, -2, 0, 1, 1, 1, 1, 2, 1, -1, -1, -1, -2, 1, 0, -3, 0, 0, 0, 8),
    (-2, 4, 0, -1, 3, 1, 0, 0, -1, 0, -2, 0, -2, 0, -4, 4, 0, 0, -4, -2, 1, -2, 3, 0, 0, -1, 8),
    (2, -4, -1, -1, -2, 2, 1, -2, 0, -2, 4, -1, -1, -2, -2, -2, -2, -2, -1, 2, -4, -1, -1, 2, 1, 2, 0, 8),
)


def _bareiss_determinant(matrix: Sequence[Sequence[int]]) -> int:
    work = [list(map(int, row)) for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("a determinant matrix must be square")
    previous = 1
    sign = 1
    for column in range(size - 1):
        pivot = next((i for i in range(column, size) if work[i][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        value = work[column][column]
        for row in range(column + 1, size):
            for index in range(column + 1, size):
                numerator = (
                    work[row][index] * value
                    - work[row][column] * work[column][index]
                )
                if numerator % previous:
                    raise AssertionError("Bareiss division ceased to be exact")
                work[row][index] = numerator // previous
        previous = value
        for row in range(column + 1, size):
            work[row][column] = 0
    return sign * work[-1][-1]


def _rank_mod_prime(matrix: IntegerMatrix, prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(row, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, prime)
        work[row] = [(inverse * value) % prime for value in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column]:
                continue
            factor = work[index][column]
            work[index] = [
                (work[index][j] - factor * work[row][j]) % prime
                for j in range(len(work[0]))
            ]
        row += 1
    return row


@lru_cache(maxsize=1)
def j2dim14_gram() -> IntegerMatrix:
    if tuple(map(len, J2DIM14_LOWER_TRIANGLE)) != tuple(range(1, 29)):
        raise AssertionError("the J2dim14 lower triangle has wrong row lengths")
    return tuple(
        tuple(J2DIM14_LOWER_TRIANGLE[max(i, j)][min(i, j)] for j in range(28))
        for i in range(28)
    )


def _kernel_mod_prime(matrix: IntegerMatrix, prime: int) -> tuple[tuple[int, ...], ...]:
    """Return a systematic basis for the right kernel over ``F_prime``."""

    work = [[entry % prime for entry in row] for row in matrix]
    pivot_columns: list[int] = []
    row = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(row, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, prime)
        work[row] = [(inverse * value) % prime for value in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column]:
                continue
            factor = work[index][column]
            work[index] = [
                (work[index][j] - factor * work[row][j]) % prime
                for j in range(len(work[0]))
            ]
        pivot_columns.append(column)
        row += 1
    free_columns = [i for i in range(len(work[0])) if i not in pivot_columns]
    output = []
    for free in free_columns:
        vector = [0] * len(work[0])
        vector[free] = 1
        for pivot_row, pivot_column in enumerate(pivot_columns):
            vector[pivot_column] = (-work[pivot_row][free]) % prime
        output.append(tuple(vector))
    return tuple(output)


def _pivot_columns_mod_prime(matrix: IntegerMatrix, prime: int) -> tuple[int, ...]:
    work = [[entry % prime for entry in row] for row in matrix]
    pivots: list[int] = []
    row = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(row, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], -1, prime)
        work[row] = [(inverse * value) % prime for value in work[row]]
        for index in range(row + 1, len(work)):
            if not work[index][column]:
                continue
            factor = work[index][column]
            work[index] = [
                (work[index][j] - factor * work[row][j]) % prime
                for j in range(len(work[0]))
            ]
        pivots.append(column)
        row += 1
    return tuple(pivots)


@lru_cache(maxsize=1)
def j2dim14_radical_five_basis() -> IntegerMatrix:
    """Column basis of ``{v in J : (v,J) subset 5 Z}``.

    In systematic kernel coordinates the free columns lift the mod-five
    radical and five times the pivot coordinate vectors complete a basis of
    its full integral preimage.
    """

    gram = j2dim14_gram()
    kernel = _kernel_mod_prime(gram, 5)
    pivot_columns = _pivot_columns_mod_prime(gram, 5)
    free_columns = [i for i in range(28) if i not in pivot_columns]
    # Recover the pivot complement directly; the kernel vectors are the first
    # columns, and 5e_i on the complementary coordinates are the remaining
    # columns.  The determinant audit below protects this construction.
    complement = list(pivot_columns)
    columns = list(kernel) + [tuple(5 if row == i else 0 for row in range(28)) for i in complement]
    return tuple(tuple(columns[column][row] for column in range(28)) for row in range(28))


def _transpose(matrix: Sequence[Sequence[int]]) -> IntegerMatrix:
    return tuple(tuple(matrix[j][i] for j in range(len(matrix))) for i in range(len(matrix[0])))


def _matmul(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> IntegerMatrix:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0])))
        for i in range(len(left))
    )


@lru_cache(maxsize=1)
def j2dim14_radical_five_gram() -> IntegerMatrix:
    basis = j2dim14_radical_five_basis()
    return _matmul(_matmul(_transpose(basis), j2dim14_gram()), basis)


def j2dim14_catalogue_audit() -> bool:
    gram = j2dim14_gram()
    radical_basis = j2dim14_radical_five_basis()
    radical_gram = j2dim14_radical_five_gram()
    return (
        all(gram[i][j] == gram[j][i] for i in range(28) for j in range(28))
        and all(gram[i][i] % 2 == 0 for i in range(28))
        and _bareiss_determinant(gram) == J2DIM14_DETERMINANT
        and _rank_mod_prime(gram, 5) == 14
        and abs(_bareiss_determinant(radical_basis)) == 5**14
        and all(value % 5 == 0 for row in radical_gram for value in row)
    )


@lru_cache(maxsize=None)
def dual_coordinate_direction(index: int = 13) -> tuple[Fraction, ...]:
    if not 0 <= index < 28:
        raise ValueError("the dual coordinate index must lie in [0,28)")
    inverse = _fraction_inverse(
        tuple(tuple(Fraction(value) for value in row) for row in j2dim14_gram())
    )
    return tuple(inverse[row][index] for row in range(28))


def j2_pair(left: Sequence[int | Fraction], right: Sequence[int | Fraction]) -> Fraction:
    if len(left) != 28 or len(right) != 28:
        raise ValueError("J2dim14 coordinate vectors have length 28")
    gram = j2dim14_gram()
    return sum(
        Fraction(left[i]) * gram[i][j] * Fraction(right[j])
        for i in range(28)
        for j in range(28)
        if gram[i][j]
    )


EPSILON = SILVER_ONE + SQRT2
RHO = 2 * SQRT2 - SILVER_ONE
P7_GENERATOR = 3 + SQRT2
DYADIC_ODD_GENERATOR = 2 + SQRT2
DYADIC_CUBIC_GENERATOR = 4 + 2 * SQRT2


@dataclass(frozen=True, slots=True)
class J2Seed:
    seed_id: str
    coefficient: Silver
    dyadic_case: str
    coefficient_p7: bool = False
    coefficient_five: bool = False

    @property
    def admissible(self) -> bool:
        return self.coefficient.sign(0) > 0 and self.coefficient.sign(1) < 0


UNIT_SEED = J2Seed("silver-H28-J2dim14-a-1+sqrt2", EPSILON, "unit")
SQRT2_SEED = J2Seed("silver-H28-J2dim14-a-sqrt2", SQRT2, "pi")
TWO_SQRT2_SEED = J2Seed("silver-H28-J2dim14-a-2sqrt2", 2 * SQRT2, "pi3")
P7_SEED = J2Seed("silver-H28-J2dim14-a-2sqrt2-minus1", RHO, "unit", True)
FIVE_UNIT_SEED = J2Seed("silver-H28-J2dim14-a-5-5sqrt2", 5 * EPSILON, "unit", False, True)
FIVE_SQRT2_SEED = J2Seed("silver-H28-J2dim14-a-5sqrt2", 5 * SQRT2, "pi", False, True)
FIVE_TWO_SQRT2_SEED = J2Seed("silver-H28-J2dim14-a-10sqrt2", 10 * SQRT2, "pi3", False, True)
FIVE_P7_SEED = J2Seed("silver-H28-J2dim14-a-10sqrt2-minus5", 5 * RHO, "unit", True, True)
SUPPORTED_SEEDS = (
    UNIT_SEED,
    SQRT2_SEED,
    TWO_SQRT2_SEED,
    P7_SEED,
    FIVE_UNIT_SEED,
    FIVE_SQRT2_SEED,
    FIVE_TWO_SQRT2_SEED,
    FIVE_P7_SEED,
)


@dataclass(frozen=True, slots=True)
class J2NormClass:
    norm: Silver
    five: bool
    p7: bool
    dyadic_spatial_gate: bool


def norm_classes(seed: J2Seed) -> tuple[J2NormClass, ...]:
    if seed.dyadic_case == "unit":
        dyadic = ((SILVER_ONE, False), (Silver(2), False))
    elif seed.dyadic_case == "pi":
        dyadic = ((DYADIC_ODD_GENERATOR, False), (Silver(2), False))
    elif seed.dyadic_case == "pi3":
        dyadic = ((Silver(2), False), (DYADIC_CUBIC_GENERATOR, True))
    else:
        raise ValueError("unsupported dyadic coefficient case")
    return tuple(
        J2NormClass(
            base * (5 if five else 1) * (P7_GENERATOR if p7 else 1),
            five,
            p7,
            spatial_gate,
        )
        for (base, spatial_gate), five, p7 in product(
            dyadic, (False, True), (False, True) if seed.coefficient_p7 else (False,)
        )
    )


@dataclass(frozen=True, slots=True)
class J2GramTriple:
    u_square: int
    v_square: int
    u_dot_v: int

    def combination_square(self, a: int, b: int) -> int:
        return a * a * self.u_square + b * b * self.v_square + 2 * a * b * self.u_dot_v


def candidate_gram_triples(target: Silver) -> tuple[J2GramTriple, ...]:
    """Necessary binary Gram types for two vectors in a min-eight even lattice."""

    if not target.is_integral or (target.b / 2).denominator != 1:
        return ()
    rational = int(target.a)
    dot = int(target.b / 2)
    output = []
    for v_square in range(0, rational // 2 + 1, 2):
        u_square = rational - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if u_square not in (0,) and u_square < J2DIM14_MINIMUM:
            continue
        if v_square not in (0,) and v_square < J2DIM14_MINIMUM:
            continue
        if dot * dot > u_square * v_square:
            continue
        triple = J2GramTriple(u_square, v_square, dot)
        if any(
            0 < triple.combination_square(a, b) < J2DIM14_MINIMUM
            for a in range(-3, 4)
            for b in range(-3, 4)
            if a or b
        ):
            continue
        output.append(triple)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class J2Control:
    seed: J2Seed
    norm_class: J2NormClass
    height: Silver
    control: Silver
    target: Silver


def _coefficients_divisible(value: Silver, prime: int) -> bool:
    return value.is_integral and int(value.a) % prime == 0 and int(value.b) % prime == 0


@lru_cache(maxsize=None)
def enumerate_controls(seed: J2Seed, height_bound: Silver) -> tuple[J2Control, ...]:
    if seed not in SUPPORTED_SEEDS or not seed.admissible:
        raise ValueError("unsupported or inadmissible J2dim14 seed")
    output = set()
    for norm_class in norm_classes(seed):
        q = norm_class.norm
        for x in _bounded_for_caps((q * height_bound / seed.coefficient, q / (-seed.coefficient))):
            if x.sign(0) <= 0 or not (2 * seed.coefficient * x / q).is_integral:
                continue
            if norm_class.five and not seed.coefficient_five and not _coefficients_divisible(x, 5):
                continue
            target = q + seed.coefficient * x * x
            height = seed.coefficient * x * x / q
            if target.sign(1) < 0 or height > height_bound:
                continue
            output.add(J2Control(seed, norm_class, height, x, target))
    return tuple(sorted(output, key=lambda row: (row.height, row.norm_class.norm, row.control)))


def reduced_target(control: J2Control) -> Silver:
    divisor = SILVER_ONE
    if control.norm_class.p7:
        divisor *= RHO * RHO
    if control.norm_class.five:
        divisor *= 5
    if control.norm_class.dyadic_spatial_gate:
        divisor *= 2
    return control.target / divisor


def gated_gram_triples(control: J2Control) -> tuple[J2GramTriple, ...]:
    return candidate_gram_triples(reduced_target(control))


def height_zero_gated_gram_triples(norm_class: J2NormClass) -> tuple[J2GramTriple, ...]:
    divisor = SILVER_ONE
    if norm_class.p7:
        divisor *= RHO * RHO
    if norm_class.five:
        divisor *= 5
    if norm_class.dyadic_spatial_gate:
        divisor *= 2
    return candidate_gram_triples(norm_class.norm / divisor)


UNIT_FIRST_HEIGHT = 14 + 10 * SQRT2
UNIT_FIRST_SLAB_BOUND = 2 + 2 * SQRT2
UNIT_WITNESS_INDEX = 13
UNIT_WITNESS_SCALE = Fraction(3, 2)


@dataclass(frozen=True, slots=True)
class J2UnitFirstPacketAudit:
    control_count: int
    distinct_height_count: int
    earlier_control_count: int
    first_control: J2Control
    first_triple: J2GramTriple
    norm_ten_shell_count: int
    wall_count: int
    dual_index: int
    dual_square: Fraction
    maximum_shell_pairing: int
    witness_scale: Fraction
    witness_spatial_square: Fraction
    witness_lorentz_norm: Silver
    wall_margin: Silver
    strict_hyperideal: bool


def unit_first_packet_audit(*, maximum_shell_pairing: int = 2) -> J2UnitFirstPacketAudit:
    controls = enumerate_controls(UNIT_SEED, UNIT_FIRST_HEIGHT)
    live = tuple((row, gated_gram_triples(row)) for row in controls if gated_gram_triples(row))
    expected = J2GramTriple(10, 10, 10)
    if len(live) != 1 or live[0][1] != (expected,):
        raise AssertionError("the J2dim14 unit first-packet census drifted")
    row = live[0][0]
    if row.height != UNIT_FIRST_HEIGHT or row.norm_class.norm != 2 or row.control != 2 + 2 * SQRT2:
        raise AssertionError("the J2dim14 unit first control drifted")
    h = dual_coordinate_direction(UNIT_WITNESS_INDEX)
    h_square = j2_pair(h, h)
    witness_square = UNIT_WITNESS_SCALE * UNIT_WITNESS_SCALE * h_square
    witness_norm = Silver(witness_square) - EPSILON
    margin = UNIT_FIRST_SLAB_BOUND - UNIT_WITNESS_SCALE * maximum_shell_pairing
    return J2UnitFirstPacketAudit(
        control_count=len(controls),
        distinct_height_count=len({item.height for item in controls}),
        earlier_control_count=sum(item.height < UNIT_FIRST_HEIGHT for item in controls),
        first_control=row,
        first_triple=expected,
        norm_ten_shell_count=J2DIM14_NORM_TEN_SHELL_COUNT,
        wall_count=J2DIM14_NORM_TEN_SHELL_COUNT,
        dual_index=UNIT_WITNESS_INDEX,
        dual_square=h_square,
        maximum_shell_pairing=maximum_shell_pairing,
        witness_scale=UNIT_WITNESS_SCALE,
        witness_spatial_square=witness_square,
        witness_lorentz_norm=witness_norm,
        wall_margin=margin,
        strict_hyperideal=witness_norm.sign(0) > 0 and margin.sign(0) > 0,
    )


# The first continuation is controlled by the dual vector h=G^-1 e_13.
# PARI proves the coordinate maxima two on the complete norm-eight and
# norm-ten shells.  The norm-twelve and norm-fourteen bounds below are
# independently recovered by the exact affine-slice enumerator.
SHORT_SHELL_COORDINATE_BOUNDS = {8: 2, 10: 2, 12: 4, 14: 3}
UNIT_CONTINUATION_HEIGHT = 137 + 97 * SQRT2
UNIT_CONTINUATION_CONTROL = 8 + 5 * SQRT2
UNIT_REPLACEMENT_SCALE = Fraction(5, 4)

UNIT_CUTTER_U: tuple[int, ...] = (
    2, 2, -11, -46, 31, 5, 2, 3, 25, -33, -23, 23, -2, 14,
    18, -14, -12, -5, -6, -13, -11, 6, -4, 8, -7, 3, -7, 2,
)
UNIT_CUTTER_V: tuple[int, ...] = (
    1, 2, -8, -33, 22, 3, 2, 2, 18, -24, -16, 17, -2, 10,
    13, -10, -9, -4, -4, -9, -8, 4, -3, 6, -5, 2, -5, 1,
)
UNIT_CUTTER_DIFFERENCE: tuple[int, ...] = tuple(
    left - right for left, right in zip(UNIT_CUTTER_U, UNIT_CUTTER_V)
)
UNIT_FIRST_PACKET_MAXIMIZER: tuple[int, ...] = (
    0, 1, -2, -7, 5, 0, 0, 1, 4, -5, -3, 4, -2, 2,
    2, -2, -2, -1, 0, -2, -2, 0, 0, 2, -1, 1, -2, 0,
)


@dataclass(frozen=True, slots=True)
class J2CoordinateSlice:
    norm: int
    coordinate: int
    node_count: int
    vectors: tuple[tuple[int, ...], ...]


@lru_cache(maxsize=1)
def _coordinate_slice_ldlt() -> tuple[
    tuple[tuple[Fraction, ...], ...],
    tuple[Fraction, ...],
    tuple[Fraction, ...],
]:
    """Exact LDL data for fixing coordinate 13 in the catalogue basis."""

    fixed = UNIT_WITNESS_INDEX
    indices = tuple(index for index in range(28) if index != fixed)
    gram = j2dim14_gram()
    matrix = tuple(
        tuple(Fraction(gram[i][j]) for j in indices) for i in indices
    )
    column = tuple(Fraction(gram[i][fixed]) for i in indices)
    inverse = _fraction_inverse(matrix)
    centre_unit = tuple(
        -sum(inverse[i][j] * column[j] for j in range(27))
        for i in range(27)
    )
    lower = [[Fraction(0) for _ in range(27)] for _ in range(27)]
    diagonal: list[Fraction] = []
    for j in range(27):
        lower[j][j] = Fraction(1)
        pivot = matrix[j][j] - sum(
            lower[j][k] * lower[j][k] * diagonal[k] for k in range(j)
        )
        if pivot <= 0:
            raise AssertionError("the J2dim14 coordinate complement is not positive")
        diagonal.append(pivot)
        for i in range(j + 1, 27):
            lower[i][j] = (
                matrix[i][j]
                - sum(
                    lower[i][k] * lower[j][k] * diagonal[k]
                    for k in range(j)
                )
            ) / pivot
    return tuple(map(tuple, lower)), tuple(diagonal), centre_unit


def _ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


@lru_cache(maxsize=None)
def coordinate_slice_vectors(norm: int, coordinate: int) -> J2CoordinateSlice:
    """Enumerate exactly ``d^2=norm, d[13]=coordinate``.

    This is an exact inhomogeneous Fincke--Pohst enumeration.  The interval
    endpoints are obtained with integer square roots of rational bounds; no
    floating-point decision enters completeness.
    """

    if norm < 0:
        raise ValueError("a spatial norm must be nonnegative")
    lower, diagonal, centre_unit = _coordinate_slice_ldlt()
    centre = tuple(coordinate * value for value in centre_unit)
    # The Schur complement is coordinate^2 / h^2 and h^2=8/5.
    residual = Fraction(norm) - Fraction(5 * coordinate * coordinate, 8)
    if residual < 0:
        return J2CoordinateSlice(norm, coordinate, 0, ())
    work = [0] * 27
    vectors: list[tuple[int, ...]] = []
    node_count = 0

    def descend(index: int, remaining: Fraction) -> None:
        nonlocal node_count
        node_count += 1
        offset = -centre[index] + sum(
            lower[j][index] * (work[j] - centre[j])
            for j in range(index + 1, 27)
        )
        radius_squared = remaining / diagonal[index]
        if radius_squared < 0:
            return
        offset_numerator = offset.numerator
        offset_denominator = offset.denominator
        radius_numerator = radius_squared.numerator
        radius_denominator = radius_squared.denominator
        bound = isqrt(
            (radius_numerator * offset_denominator * offset_denominator)
            // radius_denominator
        )
        low = _ceil_div(
            -bound - offset_numerator, offset_denominator
        )
        high = (bound - offset_numerator) // offset_denominator
        for value in range(low, high + 1):
            work[index] = value
            term = diagonal[index] * (Fraction(value) + offset) ** 2
            new_remaining = remaining - term
            if new_remaining < 0:
                continue
            if index:
                descend(index - 1, new_remaining)
            elif new_remaining == 0:
                iterator = iter(work)
                vector = tuple(
                    coordinate if i == UNIT_WITNESS_INDEX else next(iterator)
                    for i in range(28)
                )
                if j2_pair(vector, vector) != norm:
                    raise AssertionError("the exact coordinate slice produced a false vector")
                vectors.append(vector)

    descend(26, residual)
    return J2CoordinateSlice(norm, coordinate, node_count, tuple(vectors))


@dataclass(frozen=True, slots=True)
class J2DualPairingBound:
    maximum: Silver
    maximizing_u_coordinate: int
    maximizing_v_coordinate: int
    retained_count: int


@lru_cache(maxsize=None)
def dual_pairing_bound(triple: J2GramTriple) -> J2DualPairingBound:
    """Necessary exact upper bound for ``h.(u+sqrt(2)v)``.

    Besides exact positive semidefiniteness with ``h^2=8/5``, every short
    combination in ``[-3,3]^2`` uses the certified coordinate bound for its
    shell.  This retains a superset of actual pairs, so its maximum is safe.
    """

    s, t, dot = triple.u_square, triple.v_square, triple.u_dot_v
    h_square = Fraction(8, 5)
    p_bound = isqrt((8 * s) // 5) + 1
    q_bound = isqrt((8 * t) // 5) + 1
    retained: list[tuple[Silver, int, int]] = []
    for p in range(-p_bound, p_bound + 1):
        if p * p > h_square * s:
            continue
        for q in range(-q_bound, q_bound + 1):
            if q * q > h_square * t:
                continue
            determinant = (
                h_square * (s * t - dot * dot)
                - t * p * p
                - s * q * q
                + 2 * dot * p * q
            )
            if determinant < 0:
                continue
            admissible = True
            for a in range(-3, 4):
                for b in range(-3, 4):
                    if not (a or b):
                        continue
                    combination_norm = triple.combination_square(a, b)
                    combination_coordinate = a * p + b * q
                    if combination_norm == 0 and combination_coordinate:
                        admissible = False
                        break
                    bound = SHORT_SHELL_COORDINATE_BOUNDS.get(combination_norm)
                    if bound is not None and abs(combination_coordinate) > bound:
                        admissible = False
                        break
                if not admissible:
                    break
            if admissible:
                retained.append((Silver(p, q), p, q))
    if not retained:
        raise AssertionError("the augmented J2dim14 coordinate gate is empty")
    maximum, p, q = max(retained)
    return J2DualPairingBound(maximum, p, q, len(retained))


def unit_cutter_root_audit() -> bool:
    u = UNIT_CUTTER_U
    v = UNIT_CUTTER_V
    difference = UNIT_CUTTER_DIFFERENCE
    maximizer = UNIT_FIRST_PACKET_MAXIMIZER
    x = UNIT_CONTINUATION_CONTROL
    target = 276 + 194 * SQRT2
    spatial_norm = Silver(
        j2_pair(u, u) + 2 * j2_pair(v, v),
        2 * j2_pair(u, v),
    )
    return (
        len(u) == len(v) == 28
        and j2_pair(u, u) == 136
        and j2_pair(v, v) == 70
        and j2_pair(u, v) == 97
        and j2_pair(difference, difference) == 12
        and difference[UNIT_WITNESS_INDEX] == 4
        and u[UNIT_WITNESS_INDEX] == 14
        and v[UNIT_WITNESS_INDEX] == 10
        and spatial_norm == target
        and spatial_norm - EPSILON * x * x == 2
        and any(value % 2 for value in u)
        and j2_pair(maximizer, maximizer) == 10
        and j2_pair(maximizer, u) == 27
        and j2_pair(maximizer, v) == 20
    )


@dataclass(frozen=True, slots=True)
class J2ContinuationRow:
    control: J2Control
    triple: J2GramTriple
    pairing_bound: J2DualPairingBound
    old_witness_margin: Silver
    replacement_margin: Silver


@dataclass(frozen=True, slots=True)
class J2UnitContinuationAudit:
    control_count: int
    distinct_height_count: int
    nonempty_control_count: int
    gated_type_count: int
    rows: tuple[J2ContinuationRow, ...]
    first_old_witness_threat_height: Silver
    cutter_first_packet_maximum: Silver
    cutter_first_packet_acute_margin: Silver
    cutter_old_witness_pairing: Silver
    replacement_scale: Fraction
    replacement_lorentz_norm: Silver
    replacement_first_packet_margin: Silver
    replacement_all_margins_strict: bool


def unit_h137_continuation_audit(
    *, cutter_first_packet_maximum: Silver = Silver(27, 20)
) -> J2UnitContinuationAudit:
    if not unit_cutter_root_audit():
        raise AssertionError("the explicit J2dim14 continuation root failed")
    if len(coordinate_slice_vectors(12, 4).vectors) != 2:
        raise AssertionError("the norm-twelve coordinate-four slice drifted")
    if coordinate_slice_vectors(14, 4).vectors:
        raise AssertionError("the norm-fourteen coordinate-four slice is nonempty")
    controls = enumerate_controls(UNIT_SEED, UNIT_CONTINUATION_HEIGHT)
    rows = []
    for control in controls:
        for triple in gated_gram_triples(control):
            bound = dual_pairing_bound(triple)
            old_margin = (
                EPSILON * control.control
                - Silver(UNIT_WITNESS_SCALE) * bound.maximum
            )
            replacement_margin = (
                EPSILON * control.control
                - Silver(UNIT_REPLACEMENT_SCALE) * bound.maximum
            )
            rows.append(
                J2ContinuationRow(
                    control, triple, bound, old_margin, replacement_margin
                )
            )
    threats = tuple(
        row for row in rows
        if row.control.height > UNIT_FIRST_HEIGHT
        and row.old_witness_margin.sign(0) <= 0
    )
    if not threats:
        raise AssertionError("the J2dim14 continuation has no cutter threat")
    first_threat = min(row.control.height for row in threats)
    first_packet_threshold = 2 * EPSILON * UNIT_CONTINUATION_CONTROL
    acute_margin = first_packet_threshold - cutter_first_packet_maximum
    old_pairing = (
        Silver(UNIT_WITNESS_SCALE) * Silver(14, 10)
        - EPSILON * UNIT_CONTINUATION_CONTROL
    )
    replacement_norm = (
        Silver(UNIT_REPLACEMENT_SCALE * UNIT_REPLACEMENT_SCALE * Fraction(8, 5))
        - EPSILON
    )
    replacement_first_margin = (
        2 * EPSILON - Silver(2 * UNIT_REPLACEMENT_SCALE)
    )
    return J2UnitContinuationAudit(
        control_count=len(controls),
        distinct_height_count=len({control.height for control in controls}),
        nonempty_control_count=len({row.control for row in rows}),
        gated_type_count=len(rows),
        rows=tuple(rows),
        first_old_witness_threat_height=first_threat,
        cutter_first_packet_maximum=cutter_first_packet_maximum,
        cutter_first_packet_acute_margin=acute_margin,
        cutter_old_witness_pairing=old_pairing,
        replacement_scale=UNIT_REPLACEMENT_SCALE,
        replacement_lorentz_norm=replacement_norm,
        replacement_first_packet_margin=replacement_first_margin,
        replacement_all_margins_strict=(
            replacement_norm.sign(0) > 0
            and replacement_first_margin.sign(0) > 0
            and all(row.replacement_margin.sign(0) > 0 for row in rows)
        ),
    )


__all__ = [name for name in globals() if name.startswith("J2") or name.endswith("_SEED")] + [
    "SUPPORTED_SEEDS",
    "candidate_gram_triples",
    "coordinate_slice_vectors",
    "dual_coordinate_direction",
    "dual_pairing_bound",
    "enumerate_controls",
    "gated_gram_triples",
    "height_zero_gated_gram_triples",
    "j2dim14_catalogue_audit",
    "j2dim14_gram",
    "j2dim14_radical_five_basis",
    "j2dim14_radical_five_gram",
    "j2_pair",
    "norm_classes",
    "reduced_target",
    "unit_first_packet_audit",
    "unit_cutter_root_audit",
    "unit_h137_continuation_audit",
]
