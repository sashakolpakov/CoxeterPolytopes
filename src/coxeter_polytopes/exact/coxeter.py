"""Exact classification of pairs of Coxeter walls."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .golden import Golden, PHI, RationalLike


RelationKind = Literal["finite", "finite_unrecognized", "parallel", "ultraparallel"]


@dataclass(frozen=True, slots=True)
class CoxeterRelation:
    kind: RelationKind
    q: Golden
    label: int | None = None


_FINITE_Q_LABELS = {
    Golden(0): 2,
    Golden(1): 3,
    Golden(2): 4,
    PHI + 1: 5,
    Golden(3): 6,
    PHI + 2: 10,
}


def classify_root_pair(
    norm_i: Golden | RationalLike,
    inner_product: Golden | RationalLike,
    norm_j: Golden | RationalLike,
) -> CoxeterRelation:
    """Classify two outward-oriented spacelike root normals.

    The invariant is ``q = 4 (r_i,r_j)^2 / ((r_i,r_i)(r_j,r_j))``.
    No numerical tolerance or nearest Coxeter label is used.
    """

    norm_i = Golden.coerce(norm_i)
    inner_product = Golden.coerce(inner_product)
    norm_j = Golden.coerce(norm_j)
    if norm_i <= 0 or norm_j <= 0:
        raise ValueError("root norms must be positive in the identity embedding")
    if inner_product > 0:
        raise ValueError("simple roots must be oriented with nonpositive inner product")

    q = 4 * inner_product * inner_product / (norm_i * norm_j)
    label = _FINITE_Q_LABELS.get(q)
    if label is not None:
        return CoxeterRelation("finite", q, label)
    if q < 4:
        return CoxeterRelation("finite_unrecognized", q)
    if q == Golden(4):
        return CoxeterRelation("parallel", q)
    return CoxeterRelation("ultraparallel", q)

