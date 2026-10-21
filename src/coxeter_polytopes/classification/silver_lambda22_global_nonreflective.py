"""Global nonreflectivity certificate for the silver ``Lambda22`` host.

Let ``K=Q(sqrt(2))``, ``O=Z[sqrt(2)]``, ``a=sqrt(2)-1``, and

``M = <-a> + Lambda22_O``.

The bounded Vinberg calculation through real height 150 is imported only as
a finite lemma.  An explicit infinite-order integral isometry carries the
basepoint to a point of the same reflection chamber.  Consequently that
chamber has an infinite stabilizer and hence infinite hyperbolic volume.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from coxeter_polytopes.exact.silver import SILVER_ONE, SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    LAMBDA22_MIRROR_NORMS,
    lambda22_gram,
    lambda22_pair_silver,
)


ROOT = Path(__file__).resolve().parents[3]
FRONTIER_ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h150_frontier.json"

DIRECTION = (1,) + (0,) * 21
DIRECTION_NORM = 4

# On the plane with coordinates ``(x,c)``, where ``y=c*d``, the isometry is
# ``[[A,16+12sqrt(2)],[2+sqrt(2),A]]``, with ``A=5+4sqrt(2)``.
# In integral lattice coordinates this is the formula in ``isometry`` below.
ISOMETRY_A = Silver(5, 4)
ISOMETRY_B = Silver(2, 1)
ISOMETRY_DOT = Silver(4, 3)
ISOMETRY_SHIFT = Silver(1, 1)

BASEPOINT = (SILVER_ONE,) + (SILVER_ZERO,) * 22
ORBIT_ONE = (ISOMETRY_A,) + tuple(ISOMETRY_B * entry for entry in DIRECTION)

SEPARATING_CONTROL_BOUND = Silver(8, 6)
SEPARATING_HEIGHT_BOUND = Silver(56, 40)
FRONTIER_HEIGHT_BOUND = Silver(150)


def lorentz_pair(
    left: Sequence[Silver | int], right: Sequence[Silver | int]
) -> Silver:
    """The exact field-valued form on ``<-a> + Lambda22``."""

    if len(left) != 23 or len(right) != 23:
        raise ValueError("silver Lambda22 vectors have length 23")
    return (
        -LAMBDA22_COEFFICIENT * Silver.coerce(left[0]) * Silver.coerce(right[0])
        + lambda22_pair_silver(left[1:], right[1:])
    )


def direction_pair(vector: Sequence[Silver | int]) -> Silver:
    """Return ``(d,vector)`` in the positive Lambda22 block."""

    if len(vector) != 22:
        raise ValueError("Lambda22 vectors have length 22")
    return lambda22_pair_silver(DIRECTION, vector)


def isometry(vector: Sequence[Silver | int]) -> tuple[Silver, ...]:
    """Apply the explicit infinite-order integral Lorentz isometry ``g``.

    ``g(x,y)=(A*x+(4+3sqrt(2))*(d,y),``
    ``y+(B*x+(1+sqrt(2))*(d,y))*d)``.
    """

    if len(vector) != 23:
        raise ValueError("silver Lambda22 vectors have length 23")
    x = Silver.coerce(vector[0])
    y = tuple(Silver.coerce(value) for value in vector[1:])
    dot = direction_pair(y)
    new_x = ISOMETRY_A * x + ISOMETRY_DOT * dot
    shift = ISOMETRY_B * x + ISOMETRY_SHIFT * dot
    new_y = tuple(y[index] + shift * DIRECTION[index] for index in range(22))
    return (new_x, *new_y)


def inverse_isometry(vector: Sequence[Silver | int]) -> tuple[Silver, ...]:
    """Apply ``g^{-1}``; the signs of the two off-diagonal terms reverse."""

    if len(vector) != 23:
        raise ValueError("silver Lambda22 vectors have length 23")
    x = Silver.coerce(vector[0])
    y = tuple(Silver.coerce(value) for value in vector[1:])
    dot = direction_pair(y)
    new_x = ISOMETRY_A * x - ISOMETRY_DOT * dot
    shift = -ISOMETRY_B * x + ISOMETRY_SHIFT * dot
    new_y = tuple(y[index] + shift * DIRECTION[index] for index in range(22))
    return (new_x, *new_y)


def is_integral_vector(vector: Sequence[Silver]) -> bool:
    return all(value.is_integral for value in vector)


def isometry_audit() -> bool:
    """Audit integrality, inverse, and ``g^TBg=B`` on the full basis."""

    basis = tuple(
        tuple(SILVER_ONE if index == column else SILVER_ZERO for index in range(23))
        for column in range(23)
    )
    images = tuple(isometry(vector) for vector in basis)
    return (
        lambda22_pair_silver(DIRECTION, DIRECTION) == DIRECTION_NORM
        and isometry(BASEPOINT) == ORBIT_ONE
        and all(is_integral_vector(vector) for vector in images)
        and all(inverse_isometry(images[index]) == basis[index] for index in range(23))
        and all(
            lorentz_pair(images[i], images[j]) == lorentz_pair(basis[i], basis[j])
            for i in range(23)
            for j in range(23)
        )
    )


def orbit_point(exponent: int) -> tuple[Silver, ...]:
    if exponent < 0:
        point = BASEPOINT
        for _ in range(-exponent):
            point = inverse_isometry(point)
        return point
    point = BASEPOINT
    for _ in range(exponent):
        point = isometry(point)
    return point


def orbit_recurrence(exponent: int) -> tuple[Silver, Silver]:
    """Return ``(X_n,C_n)`` for ``g^n v0=(X_n,C_n*d)``."""

    point = orbit_point(exponent)
    if any(point[index] for index in range(2, 23)):
        raise AssertionError("the orbit left the Lorentz plane")
    return point[0], point[1]


def infinite_order_audit() -> bool:
    """Certify an expanding real eigenvalue, hence infinite order.

    On ``span(v0,d/2)`` the matrix has trace ``2A`` and determinant one.
    Since ``A>1`` at the identity embedding, its real eigenvalues are
    ``A+-sqrt(A^2-1)`` and the larger is strictly greater than one.
    """

    determinant = (
        ISOMETRY_A**2
        - Silver(4) * ISOMETRY_B**2 / LAMBDA22_COEFFICIENT
    )
    return (
        determinant == SILVER_ONE
        and (ISOMETRY_A - SILVER_ONE).sign(0) > 0
        and orbit_recurrence(0) == (Silver(1), Silver(0))
        and orbit_recurrence(1) == (Silver(5, 4), Silver(2, 1))
        and orbit_recurrence(2) == (Silver(113, 80), Silver(36, 26))
    )


# For x=0, writing y=u+sqrt(2)v gives
# Tr(y^2)/2=u^2+2v^2=Tr(q)/2.  These are the exact lower bounds forced by
# the ambient minimum and the dyadic/ternary pairing-radical gates.
TIME_ZERO_NORM_CASES = (
    # q, rational target, minimum for nonzero u, minimum for nonzero v
    (Silver(1), 1, 4, 4),
    (Silver(2), 2, 4, 4),
    (Silver(4, 2), 4, 12, 4),
    (Silver(4), 4, 12, 12),
    (Silver(3), 3, 24, 24),
    (Silver(6), 6, 24, 24),
    (Silver(12, 6), 12, 96, 24),
    (Silver(12), 12, 96, 96),
)


def no_time_zero_mirror_audit() -> bool:
    """Prove that the basepoint lies on no reflective hyperplane.

    The eight norm classes are complete.  In every case a nonzero ``u``
    already exceeds ``Tr(q)/2``, or a nonzero ``v`` contributes at least
    twice its listed minimum.  Thus ``u=v=0``, which cannot have square q.
    """

    if tuple(row[0] for row in TIME_ZERO_NORM_CASES) != LAMBDA22_MIRROR_NORMS:
        return False
    return all(
        q.a == target and u_minimum > target and 2 * v_minimum > target
        for q, target, u_minimum, v_minimum in TIME_ZERO_NORM_CASES
    )


@dataclass(frozen=True, slots=True)
class SeparatingMirrorBound:
    control_bound: Silver
    height_bound: Silver
    frontier_bound: Silver


def separating_mirror_bound() -> SeparatingMirrorBound:
    """Exact endpoint/Cauchy bound for a mirror separating ``v0`` and ``gv0``.

    Write a root as ``r=(x,y)`` with square ``q>0``.  Orient it so
    ``(v0,r)=-a*x<0``.  If its mirror separates the endpoints then

    ``(gv0,r)=-(3+sqrt(2))*x+(2+sqrt(2))*(d,y) >= 0``.

    Thus ``(d,y)>=((4-sqrt(2))/2)*x``.  Cauchy and
    ``y^2=q+a*x^2`` give
    ``((9/2)-2sqrt(2))*x^2 <= 4(q+a*x^2)``, hence
    ``a*x^2/q <= 56+40*sqrt(2)``.
    Equivalently ``x^2/q <= 136+96*sqrt(2)=(8+6*sqrt(2))^2``.
    """

    return SeparatingMirrorBound(
        SEPARATING_CONTROL_BOUND,
        SEPARATING_HEIGHT_BOUND,
        FRONTIER_HEIGHT_BOUND,
    )


def separator_implies_positive_frontier_pairing() -> bool:
    """Audit the exact coefficient implication for ``z=(1,d/3)``.

    Separation gives ``(d,y)>=((4-sqrt(2))/2)x``.  Hence
    ``(z,r)>=(-a+(4-sqrt(2))/6)x=(10-7sqrt(2))x/6>0``.
    The frozen H150 frontier excludes every reflective root with this
    property, since the separating height bound is below 150.
    """

    lower_coefficient = -LAMBDA22_COEFFICIENT + Silver(4, -1) / 6
    return (
        lower_coefficient == Silver(10, -7) / 6
        and lower_coefficient.sign(0) > 0
        and SEPARATING_HEIGHT_BOUND < FRONTIER_HEIGHT_BOUND
    )


GLOBAL_NONREFLECTIVITY_PROOF = (
    "The integral isometry g preserves the chosen time sheet and maps v0 to "
    "gv0. If a reflective mirror separated these points, its Vinberg height "
    "would be at most 56+40*sqrt(2) and it would pair positively with "
    "z=(1,d/3). The complete H<=150 frontier certificate excludes such a "
    "root. The complete norm alphabet also excludes mirrors with time "
    "coordinate zero, so v0 is regular. Therefore v0 and gv0 lie in the same "
    "reflection chamber C. Every integral isometry permutes the complete set "
    "of primitive reflective roots, so g normalizes the reflection group and "
    "maps chambers to chambers; hence gC=C. The arithmetic orthogonal group is "
    "discrete and g has infinite order. Since the plane trace of g is >2, g is "
    "a hyperbolic translation. Proper discontinuity gives a small ball U whose "
    "distinct g-powers are pairwise disjoint. As gC=C, choosing U inside C puts "
    "every equal-volume translate g^nU inside C. Thus C has infinite volume. "
    "Equivalently, the standard reflection-chamber action identifies the "
    "chamber stabilizer with the quotient of the arithmetic group by its normal "
    "reflection subgroup; the infinite stabilizer gives infinite index."
)


__all__ = [
    "BASEPOINT",
    "DIRECTION",
    "FRONTIER_ARTIFACT",
    "FRONTIER_HEIGHT_BOUND",
    "GLOBAL_NONREFLECTIVITY_PROOF",
    "ISOMETRY_A",
    "ISOMETRY_B",
    "ISOMETRY_DOT",
    "ISOMETRY_SHIFT",
    "ORBIT_ONE",
    "SEPARATING_CONTROL_BOUND",
    "SEPARATING_HEIGHT_BOUND",
    "SeparatingMirrorBound",
    "direction_pair",
    "infinite_order_audit",
    "inverse_isometry",
    "isometry",
    "isometry_audit",
    "lorentz_pair",
    "no_time_zero_mirror_audit",
    "orbit_point",
    "orbit_recurrence",
    "separator_implies_positive_frontier_pairing",
    "separating_mirror_bound",
    "TIME_ZERO_NORM_CASES",
]
