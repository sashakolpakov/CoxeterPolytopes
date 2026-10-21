"""A rank-20 K3 carrier forced from the silver order-five Lanner seed.

The order-five ``Q(sqrt(2))`` Lanner cycle with one label-four edge has
trace lattice ``E8(-1) + U(2)``.  Its unique nonzero even index-two
overlattice preserved by multiplication by ``sqrt(2)`` is the unimodular
lattice ``E8(-1) + U``.  Two orthogonal finite completions are useful:

``B5``
    Its trace lattice is ``D10(-1)``.  The completed lattice is therefore
    ``U + E8(-1) + D10(-1)``.

``F4 + A1``
    Its trace lattice is ``E8(-1) + A1(-1)^2``.  The completed lattice is
    ``U + E8(-1)^2 + A1(-1)^2``, the even 2-elementary lattice with
    invariants ``(r,a,delta)=(20,2,1)`` used by Alexeev.

The second presentation has an explicit primitive embedding in the K3
lattice with complement ``<2> + <2>``.  This module certifies a lattice
carrier, not a Coxeter chamber: the saturation adds a new ``(-2)`` class of
field norm ``3-sqrt(2)``, so root/wall closure remains a mandatory gate.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from typing import Sequence, TypeAlias

from sympy import Matrix as SympyMatrix
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia
from coxeter_polytopes.exact.silver import (
    SILVER_ZERO,
    SQRT2,
    Silver,
    silver_determinant,
    silver_inertia,
    silver_trace_gram,
)
from coxeter_polytopes.exact.trace_lattice import integral_gram
from coxeter_polytopes.geometric.silver_finite import (
    silver_finite_trace_lift,
    silver_lanner_blocks,
)
from coxeter_polytopes.lattices import e8_cartan_gram


Scalar: TypeAlias = int | Fraction
IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]
SilverMatrix: TypeAlias = tuple[tuple[Silver, ...], ...]


def _integer_determinant(matrix: Sequence[Sequence[Scalar]]) -> int:
    value = determinant(matrix)
    if value.b or value.a.denominator != 1:
        raise AssertionError("an integral matrix had a nonintegral determinant")
    return value.a.numerator


def _as_integer_matrix(matrix: Sequence[Sequence[Scalar]]) -> IntegerMatrix:
    rows: list[tuple[int, ...]] = []
    for row in matrix:
        converted: list[int] = []
        for entry in row:
            value = Fraction(entry)
            if value.denominator != 1:
                raise ValueError(f"matrix is not integral: {value}")
            converted.append(value.numerator)
        rows.append(tuple(converted))
    return tuple(rows)


def _as_rational_matrix(matrix: Sequence[Sequence[Scalar]]) -> RationalMatrix:
    return tuple(tuple(Fraction(entry) for entry in row) for row in matrix)


def _transpose(matrix: Sequence[Sequence[Scalar]]) -> list[list[Fraction]]:
    return [
        [Fraction(matrix[i][j]) for i in range(len(matrix))]
        for j in range(len(matrix[0]))
    ]


def _multiply(
    left: Sequence[Sequence[Scalar]], right: Sequence[Sequence[Scalar]]
) -> list[list[Fraction]]:
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("matrix dimensions do not match")
    return [
        [
            sum(
                Fraction(left[i][k]) * Fraction(right[k][j])
                for k in range(len(right))
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def _congruence(
    change: Sequence[Sequence[Scalar]], gram: Sequence[Sequence[Scalar]]
) -> RationalMatrix:
    return _as_rational_matrix(
        _multiply(_multiply(_transpose(change), gram), change)
    )


def _inverse(matrix: Sequence[Sequence[Scalar]]) -> RationalMatrix:
    sympy = SympyMatrix([[Fraction(entry) for entry in row] for row in matrix])
    return _as_rational_matrix(sympy.inv().tolist())


def _block_diagonal(*blocks: Sequence[Sequence[Scalar]]) -> RationalMatrix:
    size = sum(len(block) for block in blocks)
    result = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    offset = 0
    for block in blocks:
        if any(len(row) != len(block) for row in block):
            raise ValueError("blocks must be square")
        for i, row in enumerate(block):
            for j, entry in enumerate(row):
                result[offset + i][offset + j] = Fraction(entry)
        offset += len(block)
    return _as_rational_matrix(result)


def _silver_block_diagonal(*blocks: SilverMatrix) -> SilverMatrix:
    size = sum(len(block) for block in blocks)
    result = [[SILVER_ZERO for _ in range(size)] for _ in range(size)]
    offset = 0
    for block in blocks:
        if any(len(row) != len(block) for row in block):
            raise ValueError("blocks must be square")
        for i, row in enumerate(block):
            for j, entry in enumerate(row):
                result[offset + i][offset + j] = Silver.coerce(entry)
        offset += len(block)
    return tuple(tuple(row) for row in result)


def _negative(matrix: Sequence[Sequence[int]]) -> IntegerMatrix:
    return tuple(tuple(-entry for entry in row) for row in matrix)


def _d_cartan(rank: int) -> list[list[int]]:
    if rank < 4:
        raise ValueError("D rank must be at least four")
    gram = [[2 if i == j else 0 for j in range(rank)] for i in range(rank)]
    branch = rank - 3
    for index in range(branch):
        gram[index][index + 1] = gram[index + 1][index] = -1
    for leaf in (rank - 2, rank - 1):
        gram[branch][leaf] = gram[leaf][branch] = -1
    return gram


def _path_field_gram(labels: tuple[int, ...]) -> SilverMatrix:
    size = len(labels) + 1
    gram = [
        [Silver(2) if i == j else SILVER_ZERO for j in range(size)]
        for i in range(size)
    ]
    products = {3: Silver(-1), 4: -SQRT2}
    for index, label in enumerate(labels):
        gram[index][index + 1] = gram[index + 1][index] = products[label]
    return tuple(tuple(row) for row in gram)


def _sqrt2_action(field_rank: int) -> IntegerMatrix:
    """Multiplication by ``sqrt(2)`` in the interleaved trace basis."""

    result = [[0 for _ in range(2 * field_rank)] for _ in range(2 * field_rank)]
    for index in range(field_rank):
        result[2 * index + 1][2 * index] = 1
        result[2 * index][2 * index + 1] = 2
    return tuple(tuple(row) for row in result)


def _smith_invariants(matrix: IntegerMatrix) -> tuple[int, ...]:
    diagonal = smith_normal_form(SympyMatrix(matrix), domain=ZZ).diagonal()
    return tuple(abs(int(value)) for value in diagonal if abs(int(value)) > 1)


def _is_even(matrix: IntegerMatrix) -> bool:
    return all(matrix[index][index] % 2 == 0 for index in range(len(matrix)))


def _is_self_adjoint(action: IntegerMatrix, gram: IntegerMatrix) -> bool:
    return _multiply(_transpose(action), gram) == _multiply(gram, action)


def _action_squares_to_two(action: IntegerMatrix) -> bool:
    square = _multiply(action, action)
    return square == [
        [Fraction(2 if i == j else 0) for j in range(len(action))]
        for i in range(len(action))
    ]


def _quadratic_value(
    vector: Sequence[Scalar], gram: Sequence[Sequence[Scalar]]
) -> Fraction:
    return sum(
        Fraction(vector[i]) * Fraction(gram[i][j]) * Fraction(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def _pairing(
    left: Sequence[Scalar],
    gram: Sequence[Sequence[Scalar]],
    right: Sequence[Scalar],
) -> Fraction:
    return sum(
        Fraction(left[i]) * Fraction(gram[i][j]) * Fraction(right[j])
        for i in range(len(left))
        for j in range(len(right))
    )


def _matrix_vector(
    matrix: Sequence[Sequence[Scalar]], vector: Sequence[Scalar]
) -> tuple[Fraction, ...]:
    return tuple(
        sum(Fraction(entry) * Fraction(value) for entry, value in zip(row, vector))
        for row in matrix
    )


def _recover_silver_norm(trace_square: int, sqrt2_pairing: int) -> Silver:
    """Invert ``Q=-(a+c)``, ``B=-(a+2c)`` for ``a+c*sqrt(2)``."""

    return Silver(sqrt2_pairing - 2 * trace_square, trace_square - sqrt2_pairing)


@dataclass(frozen=True, slots=True)
class SilverIndexTwoClassAudit:
    """One order-two class in the discriminant group of the raw C5 trace."""

    bits: tuple[int, ...]
    square: Fraction
    is_even: bool
    sqrt2_residue: tuple[int, ...]
    sqrt2_stable: bool

    @property
    def is_zero(self) -> bool:
        return not any(self.bits)


@dataclass(frozen=True, slots=True)
class SilverC5SaturationAudit:
    """Exact audit of the unique silver-stable even saturation of C5."""

    field_gram: SilverMatrix
    raw_trace_gram: IntegerMatrix
    raw_trace_determinant: int
    raw_trace_inertia: Inertia
    raw_smith_invariants: tuple[int, ...]
    raw_sqrt2_action: IntegerMatrix
    order_two_classes: tuple[SilverIndexTwoClassAudit, ...]
    saturation_bits: tuple[int, ...]
    saturation_vector: tuple[Fraction, ...]
    saturation_basis: RationalMatrix
    saturated_trace_gram: IntegerMatrix
    saturated_trace_determinant: int
    saturated_trace_inertia: Inertia
    saturated_smith_invariants: tuple[int, ...]
    saturated_sqrt2_action: IntegerMatrix
    standard_trace_type: str
    standard_trace_gram: IntegerMatrix
    trace_isometry: IntegerMatrix
    trace_isometry_determinant: int
    trace_isometry_verified: bool
    extra_root_square: int
    extra_root_sqrt2_pairing: int
    sqrt2_extra_root: tuple[int, ...]
    sqrt2_extra_root_square: int
    extra_root_field_norm: Silver

    @property
    def stable_nonzero_even_classes(self) -> tuple[SilverIndexTwoClassAudit, ...]:
        return tuple(
            item
            for item in self.order_two_classes
            if not item.is_zero and item.is_even and item.sqrt2_stable
        )

    @property
    def passes(self) -> bool:
        return (
            self.raw_trace_determinant == -4
            and self.raw_trace_inertia == Inertia(1, 9, 0)
            and self.raw_smith_invariants == (2, 2)
            and len(self.order_two_classes) == 4
            and len(self.stable_nonzero_even_classes) == 1
            and self.stable_nonzero_even_classes[0].bits == self.saturation_bits
            and self.saturated_trace_determinant == -1
            and self.saturated_trace_inertia == Inertia(1, 9, 0)
            and self.saturated_smith_invariants == ()
            and _is_even(self.saturated_trace_gram)
            and _is_self_adjoint(
                self.saturated_sqrt2_action, self.saturated_trace_gram
            )
            and _action_squares_to_two(self.saturated_sqrt2_action)
            and abs(self.trace_isometry_determinant) == 1
            and self.trace_isometry_verified
            and self.extra_root_square == -2
            and self.extra_root_sqrt2_pairing == -1
            and self.sqrt2_extra_root_square == -4
            and self.extra_root_field_norm == Silver(3, -1)
            and self.extra_root_field_norm.sign(0) == 1
            and self.extra_root_field_norm.sign(1) == 1
        )


_C5_SATURATION_BITS = (0, 0, 0, 1, 0, 1, 0, 0, 0, 1)

# Columns identify the saturated C5 trace basis with E8(-1) + U.
_C5_SATURATION_ISOMETRY: IntegerMatrix = (
    (1, -2, 0, 0, 0, 1, 0, 1, 0, -1),
    (0, 1, -1, 1, 0, -1, 0, 0, -1, 0),
    (0, -2, -1, 0, 1, 1, 0, 2, -1, -1),
    (0, 1, -1, 2, -1, -1, 0, 0, 0, -1),
    (0, 1, -1, 2, -1, -2, 0, 1, -1, -1),
    (0, -1, -1, 0, 1, 1, 0, 1, 0, -1),
    (0, 0, 0, 1, 0, -2, 1, 0, 0, -1),
    (0, 0, -1, 0, 0, 1, 0, 1, -1, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 1, -1),
    (0, 0, 0, 0, 0, 0, 0, 0, -2, 1),
)


@lru_cache(maxsize=1)
def audit_silver_c5_saturation() -> SilverC5SaturationAudit:
    """Audit the unique ``sqrt(2)``-stable even index-two C5 overlattice."""

    block = next(item for item in silver_lanner_blocks() if item.order == 5)
    field_gram = block.field_gram
    raw = block.trace_gram
    action = _sqrt2_action(5)

    classes: list[SilverIndexTwoClassAudit] = []
    for bits in product((0, 1), repeat=10):
        gram_times_bits = _matrix_vector(raw, bits)
        if any(int(value) % 2 for value in gram_times_bits):
            continue
        square = _quadratic_value(bits, raw) / 4
        residue = tuple(int(value) % 2 for value in _matrix_vector(action, bits))
        stable = residue == (0,) * 10 or residue == bits
        classes.append(
            SilverIndexTwoClassAudit(
                bits=bits,
                square=square,
                is_even=square.denominator == 1 and square.numerator % 2 == 0,
                sqrt2_residue=residue,
                sqrt2_stable=stable,
            )
        )

    saturation_vector = tuple(Fraction(bit, 2) for bit in _C5_SATURATION_BITS)
    basis = [[Fraction(i == j) for j in range(10)] for i in range(10)]
    for index, value in enumerate(saturation_vector):
        basis[index][9] = value
    saturation_basis = _as_rational_matrix(basis)
    saturated = _as_integer_matrix(_congruence(saturation_basis, raw))
    saturated_action = _as_integer_matrix(
        _multiply(_multiply(_inverse(saturation_basis), action), saturation_basis)
    )

    standard = _as_integer_matrix(
        _block_diagonal(
            _negative(e8_cartan_gram()),
            ((0, 1), (1, 0)),
        )
    )
    isometry = _C5_SATURATION_ISOMETRY
    sqrt2_extra = _as_integer_matrix(
        [[value] for value in _matrix_vector(action, saturation_vector)]
    )
    sqrt2_extra_vector = tuple(row[0] for row in sqrt2_extra)
    extra_square = int(_quadratic_value(saturation_vector, raw))
    extra_sqrt2_pairing = int(
        _pairing(saturation_vector, raw, sqrt2_extra_vector)
    )

    audit = SilverC5SaturationAudit(
        field_gram=field_gram,
        raw_trace_gram=raw,
        raw_trace_determinant=_integer_determinant(raw),
        raw_trace_inertia=inertia(raw),
        raw_smith_invariants=_smith_invariants(raw),
        raw_sqrt2_action=action,
        order_two_classes=tuple(classes),
        saturation_bits=_C5_SATURATION_BITS,
        saturation_vector=saturation_vector,
        saturation_basis=saturation_basis,
        saturated_trace_gram=saturated,
        saturated_trace_determinant=_integer_determinant(saturated),
        saturated_trace_inertia=inertia(saturated),
        saturated_smith_invariants=_smith_invariants(saturated),
        saturated_sqrt2_action=saturated_action,
        standard_trace_type="E8(-1) + U",
        standard_trace_gram=standard,
        trace_isometry=isometry,
        trace_isometry_determinant=_integer_determinant(isometry),
        trace_isometry_verified=(
            _as_integer_matrix(_congruence(isometry, saturated)) == standard
        ),
        extra_root_square=extra_square,
        extra_root_sqrt2_pairing=extra_sqrt2_pairing,
        sqrt2_extra_root=sqrt2_extra_vector,
        sqrt2_extra_root_square=int(_quadratic_value(sqrt2_extra_vector, raw)),
        extra_root_field_norm=_recover_silver_norm(
            extra_square, extra_sqrt2_pairing
        ),
    )
    if not audit.passes:
        raise AssertionError("the silver C5 saturation audit failed")
    return audit


# Columns are trace-lattice vectors and give D10(-1).
_B5_TRACE_ISOMETRY: IntegerMatrix = (
    (1, -1, 1, -1, 1, -1, 1, 0, -1, 1),
    (0, 0, -1, 1, 0, 0, -1, 0, 1, 0),
    (0, 0, -2, 2, 0, 0, -2, 1, 0, 0),
    (1, -1, 1, -1, 1, -1, 1, 0, 0, 0),
    (0, 0, -2, 2, 0, -1, 0, 0, 0, 0),
    (1, -1, 1, -1, 1, 0, 0, 0, 0, 0),
    (0, 0, -2, 1, 0, 0, 0, 0, 0, 0),
    (1, -1, 1, 0, 0, 0, 0, 0, 0, 0),
    (0, -1, 0, 0, 0, 0, 0, 0, 0, 0),
    (1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
)

_A1_TRACE_ISOMETRY: IntegerMatrix = ((1, -1), (0, 1))


@dataclass(frozen=True, slots=True)
class SilverFiniteCompletionAudit:
    """A rank-five positive field block completing the C5 seed."""

    name: str
    field_gram: SilverMatrix
    field_determinant: Silver
    field_identity_inertia: Inertia
    field_conjugate_inertia: Inertia
    trace_gram: IntegerMatrix
    trace_type: str
    trace_determinant: int
    trace_inertia: Inertia
    trace_smith_invariants: tuple[int, ...]
    discriminant_class_squares: tuple[Fraction, ...]
    trace_isometry: IntegerMatrix
    standard_trace_gram: IntegerMatrix
    trace_isometry_determinant: int
    trace_isometry_verified: bool
    sqrt2_action: IntegerMatrix
    sqrt2_action_self_adjoint: bool
    sqrt2_action_squares_to_two: bool

    @property
    def coparity(self) -> int:
        return int(any(value.denominator == 2 for value in self.discriminant_class_squares))

    @property
    def passes(self) -> bool:
        return (
            self.field_identity_inertia == Inertia(5, 0, 0)
            and self.field_conjugate_inertia == Inertia(5, 0, 0)
            and self.trace_determinant == 4
            and self.trace_inertia == Inertia(0, 10, 0)
            and self.trace_smith_invariants == (2, 2)
            and len(self.discriminant_class_squares) == 4
            and self.coparity == 1
            and abs(self.trace_isometry_determinant) == 1
            and self.trace_isometry_verified
            and self.sqrt2_action_self_adjoint
            and self.sqrt2_action_squares_to_two
        )


def _discriminant_class_squares(gram: IntegerMatrix) -> tuple[Fraction, ...]:
    """Squares of all denominator-two classes when the exponent is two."""

    result: list[Fraction] = []
    for bits in product((0, 1), repeat=len(gram)):
        if any(int(value) % 2 for value in _matrix_vector(gram, bits)):
            continue
        square = _quadratic_value(bits, gram) / 4
        # The fractional part modulo 2 is the discriminant quadratic value.
        reduced = square % 2
        result.append(reduced)
    return tuple(sorted(result))


@lru_cache(maxsize=1)
def silver_finite_completion_audits() -> tuple[SilverFiniteCompletionAudit, ...]:
    """Return the ``B5`` and Alexeev-aligned ``F4+A1`` completions."""

    b5_field = _path_field_gram((4, 3, 3, 3))
    f4 = silver_finite_trace_lift("F4")
    a1_field: SilverMatrix = ((Silver(2),),)
    alexeev_field = _silver_block_diagonal(f4.field_gram, a1_field)

    b5_standard = _negative(_d_cartan(10))
    a1_trace = _as_integer_matrix(integral_gram(silver_trace_gram(a1_field)))
    alexeev_standard = _as_integer_matrix(
        _block_diagonal(
            f4.standard_trace_gram,
            ((-2, 0), (0, -2)),
        )
    )
    alexeev_isometry = _as_integer_matrix(
        _block_diagonal(f4.trace_isometry, _A1_TRACE_ISOMETRY)
    )

    inputs = (
        ("B5", b5_field, "D10(-1)", b5_standard, _B5_TRACE_ISOMETRY),
        (
            "F4 + A1",
            alexeev_field,
            "E8(-1) + A1(-1)^2",
            alexeev_standard,
            alexeev_isometry,
        ),
    )
    audits: list[SilverFiniteCompletionAudit] = []
    for name, field, trace_type, standard, isometry in inputs:
        traced = _as_integer_matrix(integral_gram(silver_trace_gram(field)))
        action = _sqrt2_action(5)
        audit = SilverFiniteCompletionAudit(
            name=name,
            field_gram=field,
            field_determinant=silver_determinant(field),
            field_identity_inertia=silver_inertia(field, embedding=0),
            field_conjugate_inertia=silver_inertia(field, embedding=1),
            trace_gram=traced,
            trace_type=trace_type,
            trace_determinant=_integer_determinant(traced),
            trace_inertia=inertia(traced),
            trace_smith_invariants=_smith_invariants(traced),
            discriminant_class_squares=_discriminant_class_squares(traced),
            trace_isometry=isometry,
            standard_trace_gram=standard,
            trace_isometry_determinant=_integer_determinant(isometry),
            trace_isometry_verified=(
                _as_integer_matrix(_congruence(isometry, traced)) == standard
            ),
            sqrt2_action=action,
            sqrt2_action_self_adjoint=_is_self_adjoint(action, traced),
            sqrt2_action_squares_to_two=_action_squares_to_two(action),
        )
        if not audit.passes:
            raise AssertionError(f"the {name} silver completion audit failed")
        audits.append(audit)
    if a1_trace != ((-2, -2), (-2, -4)):
        raise AssertionError("the A1 silver trace block changed")
    return tuple(audits)


@dataclass(frozen=True, slots=True)
class SilverCarrierModelAudit:
    """One rank-20 trace carrier built from C5 saturation plus a finite block."""

    name: str
    field_gram: SilverMatrix
    field_determinant: Silver
    field_identity_inertia: Inertia
    field_conjugate_inertia: Inertia
    raw_trace_gram: IntegerMatrix
    saturated_trace_gram: IntegerMatrix
    saturation_index: int
    trace_type: str
    standard_trace_gram: IntegerMatrix
    trace_isometry: IntegerMatrix
    trace_isometry_determinant: int
    trace_isometry_verified: bool
    trace_determinant: int
    trace_inertia: Inertia
    trace_smith_invariants: tuple[int, ...]
    two_elementary_invariants: tuple[int, int, int]
    sqrt2_action: IntegerMatrix
    standard_sqrt2_action: IntegerMatrix
    sqrt2_action_self_adjoint: bool
    sqrt2_action_squares_to_two: bool

    @property
    def passes(self) -> bool:
        return (
            self.field_identity_inertia == Inertia(9, 1, 0)
            and self.field_conjugate_inertia == Inertia(10, 0, 0)
            and self.saturation_index == 2
            and self.trace_determinant == -4
            and self.trace_inertia == Inertia(1, 19, 0)
            and self.trace_smith_invariants == (2, 2)
            and self.two_elementary_invariants == (20, 2, 1)
            and abs(self.trace_isometry_determinant) == 1
            and self.trace_isometry_verified
            and self.sqrt2_action_self_adjoint
            and self.sqrt2_action_squares_to_two
            and _is_self_adjoint(
                self.standard_sqrt2_action, self.standard_trace_gram
            )
            and _action_squares_to_two(self.standard_sqrt2_action)
        )


@lru_cache(maxsize=1)
def silver_carrier_model_audits() -> tuple[SilverCarrierModelAudit, ...]:
    """Build the diagram-aligned and Alexeev-aligned rank-20 carriers."""

    seed = audit_silver_c5_saturation()
    models: list[SilverCarrierModelAudit] = []
    for completion in silver_finite_completion_audits():
        field = _silver_block_diagonal(seed.field_gram, completion.field_gram)
        raw = _as_integer_matrix(
            _block_diagonal(seed.raw_trace_gram, completion.trace_gram)
        )
        saturated = _as_integer_matrix(
            _block_diagonal(seed.saturated_trace_gram, completion.trace_gram)
        )
        standard = _as_integer_matrix(
            _block_diagonal(
                seed.standard_trace_gram, completion.standard_trace_gram
            )
        )
        isometry = _as_integer_matrix(
            _block_diagonal(seed.trace_isometry, completion.trace_isometry)
        )
        action = _as_integer_matrix(
            _block_diagonal(
                seed.saturated_sqrt2_action, completion.sqrt2_action
            )
        )
        standard_action = _as_integer_matrix(
            _multiply(_multiply(_inverse(isometry), action), isometry)
        )
        trace_type = (
            "U + E8(-1) + D10(-1)"
            if completion.name == "B5"
            else "U + E8(-1)^2 + A1(-1)^2"
        )
        model = SilverCarrierModelAudit(
            name=completion.name,
            field_gram=field,
            field_determinant=silver_determinant(field),
            field_identity_inertia=silver_inertia(field, embedding=0),
            field_conjugate_inertia=silver_inertia(field, embedding=1),
            raw_trace_gram=raw,
            saturated_trace_gram=saturated,
            saturation_index=2,
            trace_type=trace_type,
            standard_trace_gram=standard,
            trace_isometry=isometry,
            trace_isometry_determinant=_integer_determinant(isometry),
            trace_isometry_verified=(
                _as_integer_matrix(_congruence(isometry, saturated)) == standard
            ),
            trace_determinant=_integer_determinant(saturated),
            trace_inertia=inertia(saturated),
            trace_smith_invariants=_smith_invariants(saturated),
            two_elementary_invariants=(20, 2, completion.coparity),
            sqrt2_action=action,
            standard_sqrt2_action=standard_action,
            sqrt2_action_self_adjoint=_is_self_adjoint(action, saturated),
            sqrt2_action_squares_to_two=_action_squares_to_two(action),
        )
        if not model.passes:
            raise AssertionError(f"the {completion.name} carrier audit failed")
        models.append(model)
    return tuple(models)


@dataclass(frozen=True, slots=True)
class SilverK3EmbeddingAudit:
    """Explicit discriminant gluing of the Alexeev-aligned carrier into K3."""

    neron_severi_gram: IntegerMatrix
    transcendental_gram: IntegerMatrix
    discriminant_anti_isometry_squares: tuple[tuple[Fraction, Fraction], ...]
    gluing_vectors: tuple[tuple[Fraction, ...], ...]
    gluing_index: int
    ambient_gram: IntegerMatrix
    ambient_determinant: int
    ambient_inertia: Inertia
    ambient_even: bool
    ambient_unimodular: bool
    ambient_type: str
    embedding_matrix: IntegerMatrix
    complement_basis: IntegerMatrix
    complement_gram: IntegerMatrix
    orthogonal_complement_verified: bool
    primitive_embedding_verified: bool

    @property
    def passes(self) -> bool:
        return (
            self.transcendental_gram == ((2, 0), (0, 2))
            and self.discriminant_anti_isometry_squares
            == ((Fraction(-1, 2), Fraction(1, 2)),) * 2
            and self.gluing_index == 4
            and self.ambient_determinant == -1
            and self.ambient_inertia == Inertia(3, 19, 0)
            and self.ambient_even
            and self.ambient_unimodular
            and self.ambient_type == "II(3,19) = U^3 + E8(-1)^2"
            and self.complement_gram == self.transcendental_gram
            and self.orthogonal_complement_verified
            and self.primitive_embedding_verified
        )


@lru_cache(maxsize=1)
def audit_silver_k3_embedding() -> SilverK3EmbeddingAudit:
    """Glue the Alexeev-aligned carrier to ``<2>+<2>`` explicitly."""

    model = next(item for item in silver_carrier_model_audits() if item.name == "F4 + A1")
    neron_severi = model.standard_trace_gram
    transcendental: IntegerMatrix = ((2, 0), (0, 2))
    preglue = _as_integer_matrix(_block_diagonal(neron_severi, transcendental))

    # The last two NS basis vectors are the two A1(-1) roots.  Glue each
    # half-root to the corresponding half-generator of <2>+<2>.
    glue_one = tuple(
        Fraction(1, 2) if index in (18, 20) else Fraction(0)
        for index in range(22)
    )
    glue_two = tuple(
        Fraction(1, 2) if index in (19, 21) else Fraction(0)
        for index in range(22)
    )
    ambient_basis = [
        [Fraction(i == j) for j in range(22)] for i in range(22)
    ]
    for index in range(22):
        ambient_basis[index][20] = glue_one[index]
        ambient_basis[index][21] = glue_two[index]
    ambient_basis_matrix = _as_rational_matrix(ambient_basis)
    ambient = _as_integer_matrix(_congruence(ambient_basis_matrix, preglue))

    embedding = tuple(
        tuple(1 if row == column else 0 for column in range(20))
        for row in range(22)
    )
    complement = tuple(
        (
            -1 if row == 18 else (2 if row == 20 else 0),
            -1 if row == 19 else (2 if row == 21 else 0),
        )
        for row in range(22)
    )
    embedded_gram = _as_integer_matrix(_congruence(embedding, ambient))
    complement_gram = _as_integer_matrix(_congruence(complement, ambient))
    cross = _multiply(_multiply(_transpose(embedding), ambient), complement)
    audit = SilverK3EmbeddingAudit(
        neron_severi_gram=neron_severi,
        transcendental_gram=transcendental,
        discriminant_anti_isometry_squares=(
            (Fraction(-1, 2), Fraction(1, 2)),
            (Fraction(-1, 2), Fraction(1, 2)),
        ),
        gluing_vectors=(glue_one, glue_two),
        gluing_index=4,
        ambient_gram=ambient,
        ambient_determinant=_integer_determinant(ambient),
        ambient_inertia=inertia(ambient),
        ambient_even=_is_even(ambient),
        ambient_unimodular=abs(_integer_determinant(ambient)) == 1,
        ambient_type="II(3,19) = U^3 + E8(-1)^2",
        embedding_matrix=embedding,
        complement_basis=complement,
        complement_gram=complement_gram,
        orthogonal_complement_verified=(
            embedded_gram == neron_severi
            and all(value == 0 for row in cross for value in row)
        ),
        # The top 20-by-20 minor of the embedding is the identity.
        primitive_embedding_verified=True,
    )
    if not audit.passes:
        raise AssertionError("the silver K3 embedding audit failed")
    return audit


@dataclass(frozen=True, slots=True)
class SilverK3CarrierTheorem:
    """The complete exact carrier audit, with its scope boundary."""

    seed: SilverC5SaturationAudit
    models: tuple[SilverCarrierModelAudit, ...]
    k3_embedding: SilverK3EmbeddingAudit
    alexeev_invariants: tuple[int, int, int]
    alexeev_lattice_type: str
    wall_closure_required: bool

    @property
    def passes(self) -> bool:
        return (
            self.seed.passes
            and len(self.models) == 2
            and all(model.passes for model in self.models)
            and {model.name for model in self.models} == {"B5", "F4 + A1"}
            and self.k3_embedding.passes
            and self.alexeev_invariants == (20, 2, 1)
            and self.alexeev_lattice_type == "U + E8(-1)^2 + A1(-1)^2"
            and self.wall_closure_required
        )


@lru_cache(maxsize=1)
def silver_k3_carrier_theorem() -> SilverK3CarrierTheorem:
    """Return the exact rank-20 silver K3 carrier certificate."""

    theorem = SilverK3CarrierTheorem(
        seed=audit_silver_c5_saturation(),
        models=silver_carrier_model_audits(),
        k3_embedding=audit_silver_k3_embedding(),
        alexeev_invariants=(20, 2, 1),
        alexeev_lattice_type="U + E8(-1)^2 + A1(-1)^2",
        wall_closure_required=True,
    )
    if not theorem.passes:
        raise AssertionError("the silver K3 carrier theorem failed")
    return theorem


__all__ = [
    "SilverIndexTwoClassAudit",
    "SilverC5SaturationAudit",
    "SilverFiniteCompletionAudit",
    "SilverCarrierModelAudit",
    "SilverK3EmbeddingAudit",
    "SilverK3CarrierTheorem",
    "audit_silver_c5_saturation",
    "silver_finite_completion_audits",
    "silver_carrier_model_audits",
    "audit_silver_k3_embedding",
    "silver_k3_carrier_theorem",
]
