"""Crystallographic trace lifts of the finite golden Coxeter systems.

For the K3 scalar ``alpha=-phi/sqrt(5)``, restriction of scalars turns the
noncrystallographic finite root lattices over ``Z[phi]`` into familiar even
negative-definite integral lattices::

    H2  ->  A4(-1)
    H3  ->  D6(-1)
    H4  ->  E8(-1)

The change-of-basis matrices below are exact integral certificates, rather
than identifications inferred from floating-point spectra or determinants.
The field roots have norm two and the 5-edge product is ``-phi``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia
from coxeter_polytopes.exact.trace_lattice import integral_gram, trace_gram
from coxeter_polytopes.lattices import e8_cartan_gram

from .k3 import K3_TRACE_SCALAR


IntegerMatrix = tuple[tuple[int, ...], ...]
GoldenMatrix = tuple[tuple[Golden, ...], ...]


def _path_gram(labels: tuple[int, ...]) -> GoldenMatrix:
    size = len(labels) + 1
    gram = [
        [Golden(2) if i == j else Golden(0) for j in range(size)]
        for i in range(size)
    ]
    products = {3: Golden(-1), 5: -PHI}
    for index, label in enumerate(labels):
        try:
            product = products[label]
        except KeyError as error:
            raise ValueError("a golden finite path uses only labels 3 and 5") from error
        gram[index][index + 1] = gram[index + 1][index] = product
    return tuple(tuple(row) for row in gram)


def _a_cartan(rank: int) -> list[list[int]]:
    gram = [[2 if i == j else 0 for j in range(rank)] for i in range(rank)]
    for index in range(rank - 1):
        gram[index][index + 1] = gram[index + 1][index] = -1
    return gram


def _d6_cartan() -> list[list[int]]:
    gram = [[2 if i == j else 0 for j in range(6)] for i in range(6)]
    for left, right in ((0, 1), (1, 2), (2, 3), (3, 4), (3, 5)):
        gram[left][right] = gram[right][left] = -1
    return gram


def _negative(matrix: list[list[int]]) -> IntegerMatrix:
    return tuple(tuple(-entry for entry in row) for row in matrix)


_FIELD_GRAMS: dict[str, GoldenMatrix] = {
    "H2": _path_gram((5,)),
    "H3": _path_gram((5, 3)),
    "H4": _path_gram((5, 3, 3)),
}

_TRACE_TYPES = {"H2": "A4(-1)", "H3": "D6(-1)", "H4": "E8(-1)"}

_STANDARD_TRACE_GRAMS: dict[str, IntegerMatrix] = {
    "H2": _negative(_a_cartan(4)),
    "H3": _negative(_d6_cartan()),
    "H4": _negative(e8_cartan_gram()),
}

# For each type, Q_trace = S^t Q_standard S in the interleaved trace basis
# e_1, phi*e_1, ..., e_r, phi*e_r.  Every S has determinant +/-1.
_TRACE_ISOMETRIES: dict[str, IntegerMatrix] = {
    "H2": (
        (0, 0, 0, -1),
        (0, 1, 0, -1),
        (0, 0, 1, 0),
        (1, 1, 0, -1),
    ),
    "H3": (
        (0, -1, 0, 0, 1, 1),
        (0, -1, 0, 0, 1, 0),
        (0, -1, 0, -1, 1, 1),
        (0, -2, 1, 0, 0, 0),
        (0, -1, 0, 0, 0, 0),
        (1, 0, 0, -1, 0, 0),
    ),
    "H4": (
        (0, 0, 0, 0, 0, 0, -1, -1),
        (0, 0, 0, 0, -1, -1, 0, 0),
        (0, 0, 1, 0, -2, -2, 0, 0),
        (0, 1, 0, -1, -1, -1, 0, 0),
        (0, 0, 0, -1, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, -1, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, -1),
        (1, 1, 0, -1, -1, -1, 0, 0),
    ),
}


def _integer_determinant(matrix: IntegerMatrix) -> int:
    value = determinant(matrix)
    if value.b or value.a.denominator != 1:
        raise AssertionError("an integer matrix had a nonintegral determinant")
    return value.a.numerator


def _congruence(
    change: IntegerMatrix, gram: IntegerMatrix
) -> IntegerMatrix:
    size = len(change)
    return tuple(
        tuple(
            sum(
                change[k][i] * gram[k][ell] * change[ell][j]
                for k in range(size)
                for ell in range(size)
            )
            for j in range(size)
        )
        for i in range(size)
    )


@dataclass(frozen=True, slots=True)
class GoldenFiniteTraceLift:
    """An exact restriction-of-scalars certificate for ``H2``, ``H3`` or ``H4``."""

    field_type: str
    trace_type: str
    field_gram: GoldenMatrix
    trace_gram: IntegerMatrix
    standard_trace_gram: IntegerMatrix
    trace_isometry: IntegerMatrix
    field_determinant: Golden
    trace_determinant: int
    field_identity_inertia: Inertia
    field_conjugate_inertia: Inertia
    trace_inertia: Inertia
    isometry_determinant: int
    isometry_verified: bool

    @property
    def field_rank(self) -> int:
        return len(self.field_gram)

    @property
    def trace_rank(self) -> int:
        return len(self.trace_gram)

    @property
    def orthogonally_splits_every_integral_overlattice(self) -> bool:
        """Whether the trace core is unimodular, hence splits integrally."""

        return abs(self.trace_determinant) == 1

    @property
    def passes(self) -> bool:
        rank = self.field_rank
        return (
            self.field_identity_inertia == Inertia(rank, 0, 0)
            and self.field_conjugate_inertia == Inertia(rank, 0, 0)
            and self.trace_inertia == Inertia(0, 2 * rank, 0)
            and abs(self.isometry_determinant) == 1
            and self.isometry_verified
        )


def golden_finite_trace_lift(field_type: str) -> GoldenFiniteTraceLift:
    """Return the exact crystallographic trace lift of a golden finite type."""

    try:
        field_gram = _FIELD_GRAMS[field_type]
        standard = _STANDARD_TRACE_GRAMS[field_type]
        change = _TRACE_ISOMETRIES[field_type]
    except KeyError as error:
        raise ValueError("field_type must be one of H2, H3, H4") from error

    traced = tuple(
        tuple(row)
        for row in integral_gram(trace_gram(field_gram, K3_TRACE_SCALAR))
    )
    return GoldenFiniteTraceLift(
        field_type=field_type,
        trace_type=_TRACE_TYPES[field_type],
        field_gram=field_gram,
        trace_gram=traced,
        standard_trace_gram=standard,
        trace_isometry=change,
        field_determinant=determinant(field_gram),
        trace_determinant=_integer_determinant(traced),
        field_identity_inertia=inertia(field_gram, embedding=0),
        field_conjugate_inertia=inertia(field_gram, embedding=1),
        trace_inertia=inertia(traced),
        isometry_determinant=_integer_determinant(change),
        isometry_verified=traced == _congruence(change, standard),
    )


@dataclass(frozen=True, slots=True)
class ForcedGoldenCore:
    """Finite trace core forced inside a golden Lanner missing face."""

    missing_face_size: int
    field_type: str
    trace_type: str
    trace_rank: int
    trace_determinant: int
    orthogonally_splits: bool


def forced_golden_core(missing_face_size: int) -> ForcedGoldenCore:
    """Return the core forced by an order-four or order-five missing face.

    A Lanner diagram of order four contains ``H3`` or ``B3``; one of order
    five contains ``H4`` or ``F4``.  Equal-norm forms over ``Q(sqrt(5))`` do
    not admit a label-four edge, so only the ``H`` alternatives survive.
    """

    field_type_by_size = {4: "H3", 5: "H4"}
    try:
        lift = golden_finite_trace_lift(field_type_by_size[missing_face_size])
    except KeyError as error:
        raise ValueError(
            "the forced-core theorem applies to missing faces of size 4 or 5"
        ) from error
    return ForcedGoldenCore(
        missing_face_size=missing_face_size,
        field_type=lift.field_type,
        trace_type=lift.trace_type,
        trace_rank=lift.trace_rank,
        trace_determinant=lift.trace_determinant,
        orthogonally_splits=lift.orthogonally_splits_every_integral_overlattice,
    )


@dataclass(frozen=True, slots=True)
class H4LannerExtension:
    """A norm-two order-five Lanner extension of the standard ``H4`` path."""

    attachment_labels: tuple[int, int, int, int]
    field_determinant: Golden
    complement_norm: Golden
    trace_determinant: int
    trace_inertia: Inertia
    visible_curve_type: str


def h4_lanner_extensions() -> tuple[H4LannerExtension, ...]:
    """Enumerate all golden norm-two Lanner extensions of ``[5,3,3]``.

    The fifth root meets each ``H4`` root with label 2, 3 or 5.  All 81
    attachments are checked exactly: the full identity-place Gram matrix must
    have inertia ``(4,1)``, every deletion must be elliptic, and the conjugate
    Gram matrix must be positive definite.  Since ``H4`` is unimodular over
    ``Z[phi]``, the Schur-complement norm is represented integrally by the
    rank-six orthogonal complement in any ambient rank-ten candidate.
    """

    h4 = _FIELD_GRAMS["H4"]
    h4_determinant = determinant(h4)
    products = {2: Golden(0), 3: Golden(-1), 5: -PHI}
    extensions: list[H4LannerExtension] = []
    for labels in product((2, 3, 5), repeat=4):
        attachment = [products[label] for label in labels]
        gram = [
            [*h4[index], attachment[index]] for index in range(len(h4))
        ]
        gram.append([*attachment, Golden(2)])
        if inertia(gram, embedding=0) != Inertia(4, 1, 0):
            continue
        if inertia(gram, embedding=1) != Inertia(5, 0, 0):
            continue
        if any(
            inertia(
                [
                    [gram[i][j] for j in range(5) if j != omitted]
                    for i in range(5)
                    if i != omitted
                ],
                embedding=0,
            )
            != Inertia(4, 0, 0)
            for omitted in range(5)
        ):
            continue
        field_determinant = determinant(gram)
        norm_determinant = field_determinant.norm()
        if norm_determinant.denominator != 1:
            raise AssertionError("an integral trace determinant was fractional")
        traced = tuple(
            tuple(row)
            for row in integral_gram(trace_gram(gram, K3_TRACE_SCALAR))
        )
        trace_determinant = _integer_determinant(traced)
        if trace_determinant != norm_determinant.numerator:
            raise AssertionError("field norm and trace determinant disagree")
        trace_signature = inertia(traced)
        if trace_signature != Inertia(1, 9, 0):
            raise AssertionError("a Lanner trace lift had the wrong signature")
        if labels[:3] == (2, 2, 2) and labels[3] in (3, 5):
            visible_type = "A5(-1)"
        elif labels == (2, 2, 3, 2):
            visible_type = "D5(-1)"
        else:
            raise AssertionError("an unclassified golden H4 Lanner extension arose")
        extensions.append(
            H4LannerExtension(
                attachment_labels=labels,
                field_determinant=field_determinant,
                complement_norm=field_determinant / h4_determinant,
                trace_determinant=trace_determinant,
                trace_inertia=trace_signature,
                visible_curve_type=visible_type,
            )
        )
    return tuple(extensions)
