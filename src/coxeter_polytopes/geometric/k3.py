"""K3 wall tests for the Q(sqrt(5)) eigen-nef construction.

The trace lattice uses ``Tr(alpha*h)`` with
``alpha = -phi/sqrt(5)``.  If an integral field vector represents a K3
``(-2)``-class and its wall meets the indefinite identity eigenspace, its
field norm must be positive at both real places.  These requirements leave
exactly the five elements in :data:`K3_WALL_NORMS`.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

from coxeter_polytopes.exact.golden import Golden, PHI, RationalLike
from coxeter_polytopes.exact.trace_lattice import INV_SQRT5


K3_TRACE_SCALAR = -PHI * INV_SQRT5
K3_WALL_NORMS = tuple(Golden(2 - b, b) for b in range(-3, 2))
"""The five possible field norms of relevant trace-lattice (-2)-walls."""


def k3_trace_square(field_norm: Golden | RationalLike) -> Fraction:
    """Return ``Tr(alpha * field_norm)`` for the fixed K3 trace scalar."""

    return (K3_TRACE_SCALAR * Golden.coerce(field_norm)).trace()


def is_k3_wall_norm(field_norm: Golden | RationalLike) -> bool:
    """Test the necessary norm conditions for a wall of the eigen-nef cone.

    The test is exact: the trace-lattice square must be ``-2`` and the field
    norm must be strictly positive at both real embeddings.  For an element
    of ``Z[phi]`` this is equivalent to membership in ``K3_WALL_NORMS``.
    """

    value = Golden.coerce(field_norm)
    return (
        value.is_integral
        and k3_trace_square(value) == -2
        and value.sign(0) > 0
        and value.sign(1) > 0
    )


def _coerce_field_gram(
    field_gram: Sequence[Sequence[Golden | RationalLike]],
) -> list[list[Golden]]:
    gram = [[Golden.coerce(entry) for entry in row] for row in field_gram]
    size = len(gram)
    if any(len(row) != size for row in gram):
        raise ValueError("field Gram matrix must be square")
    if any(gram[i][j] != gram[j][i] for i in range(size) for j in range(i)):
        raise ValueError("field Gram matrix must be symmetric")
    return gram


def field_inner_product(
    field_gram: Sequence[Sequence[Golden | RationalLike]],
    left: Sequence[Golden | RationalLike],
    right: Sequence[Golden | RationalLike],
) -> Golden:
    """Evaluate an ``F``-bilinear Gram pairing exactly."""

    gram = _coerce_field_gram(field_gram)
    if len(left) != len(gram) or len(right) != len(gram):
        raise ValueError("vector lengths must equal the Gram-matrix rank")
    x = [Golden.coerce(entry) for entry in left]
    y = [Golden.coerce(entry) for entry in right]
    return sum(
        x[i] * gram[i][j] * y[j]
        for i in range(len(gram))
        for j in range(len(gram))
    )


def reflection_coefficients(
    field_gram: Sequence[Sequence[Golden | RationalLike]],
    root: Sequence[Golden | RationalLike],
) -> tuple[Golden, ...]:
    """Return the coefficients ``2*h(e_i, root)/h(root, root)``.

    They are the coordinates of the root covector controlling whether the
    reflection ``x -> x - 2*h(x,r)/h(r,r)*r`` preserves the free
    ``Z[phi]``-module.
    """

    gram = _coerce_field_gram(field_gram)
    if len(root) != len(gram):
        raise ValueError("root length must equal the Gram-matrix rank")
    vector = [Golden.coerce(entry) for entry in root]
    norm = field_inner_product(gram, vector, vector)
    if not norm:
        raise ValueError("a reflection root must have nonzero norm")
    return tuple(
        2 * sum(gram[i][j] * vector[j] for j in range(len(gram))) / norm
        for i in range(len(gram))
    )


def is_integral_reflection_root(
    field_gram: Sequence[Sequence[Golden | RationalLike]],
    root: Sequence[Golden | RationalLike],
) -> bool:
    """Whether ``root`` defines a reflection preserving the field lattice."""

    vector = [Golden.coerce(entry) for entry in root]
    if not all(entry.is_integral for entry in vector):
        return False
    try:
        coefficients = reflection_coefficients(field_gram, vector)
    except ValueError:
        return False
    return all(entry.is_integral for entry in coefficients)
