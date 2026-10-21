"""Exact two-dotted four-node tree packets over ``Q(sqrt(5))``.

The full two-dotted four-node box has tens of billions of assignments, but
its tree tranche has no cycle equations and admits a closed exact count.
Modulo graph automorphisms there are three cases: a star, a path whose two
dotted edges are disjoint, and a path whose two dotted edges are adjacent.

In the star and adjacent-dotted cases selected-place Lorentzian inertia is
automatic, and conjugate positive definiteness reduces to one exact
inequality.  A disjoint-dotted path can also be Lorentzian, but only in a very
thin region: its determinant must be negative at the selected place and
positive at the conjugate place.  The counting routines below exploit exact
one- and two-dimensional dominance joins on the 6,869 dotted values.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cmp_to_key, lru_cache
from hashlib import sha256

from coxeter_polytopes.classification.golden_agcat_manifest import (
    golden_coxeter_frame_pair_alphabet,
    golden_finite_pair_values,
)
from coxeter_polytopes.exact.golden import Golden


def _conjugate_compare(left: Golden, right: Golden) -> int:
    return (left - right).sign(1)


@lru_cache(maxsize=1)
def conjugate_sorted_dotted_values() -> tuple[Golden, ...]:
    values = [
        value
        for value in golden_coxeter_frame_pair_alphabet()
        if (value - 4).sign(0) > 0
    ]
    values.sort(key=cmp_to_key(_conjugate_compare))
    if len(values) != 6_869:
        raise AssertionError("the golden dotted alphabet changed")
    if any((value - 4).sign(1) >= 0 or value.sign(1) <= 0 for value in values):
        raise AssertionError("a dotted conjugate left the open Cauchy interval")
    return tuple(values)


def count_star_pairs(values: tuple[Golden, ...], finite: Golden) -> int:
    """Count unordered ``x<=y`` with ``x'+y'+finite'<4``."""

    right = len(values) - 1
    total = 0
    for left, first in enumerate(values):
        while right >= left and (first + values[right] + finite - 4).sign(1) >= 0:
            right -= 1
        if right < left:
            break
        total += right - left + 1
    return total


def count_adjacent_path_pairs(values: tuple[Golden, ...], finite: Golden) -> int:
    """Count ordered remote/adjacent dotted pairs with conjugate PD Gram.

    For path edge products ``x,y,z`` in that order, where ``x,y`` are dotted
    and ``z`` is finite, the full determinant is

    ``(4-x)(4-z)-4y``.

    Its strict positivity already implies positivity of the preceding
    three-node determinant.
    """

    right = len(values) - 1
    total = 0
    for remote in values:
        right_hand_side = (4 - remote) * (4 - finite)
        while right >= 0 and (4 * values[right] - right_hand_side).sign(1) >= 0:
            right -= 1
        total += right + 1
    return total


class _Fenwick:
    def __init__(self, size: int) -> None:
        self.values = [0] * (size + 1)

    def add(self, index: int) -> None:
        index += 1
        while index < len(self.values):
            self.values[index] += 1
            index += index & -index

    def prefix(self, length: int) -> int:
        total = 0
        while length:
            total += self.values[length]
            length -= length & -length
        return total


def _conjugate_lower_bound(values: tuple[Golden, ...], bound: Golden) -> int:
    left, right = 0, len(values)
    while left < right:
        middle = (left + right) // 2
        if (values[middle] - bound).sign(1) < 0:
            left = middle + 1
        else:
            right = middle
    return left


def count_disjoint_path_pairs(values: tuple[Golden, ...], finite: Golden) -> int:
    """Count unordered endpoint-dotted paths passing both real places.

    For endpoint products ``x,z`` and middle finite product ``y``, the full
    determinant is ``(4-x)(4-z)-4y``.  Selected Lorentzian inertia is
    equivalent to its identity value being negative; conjugate positive
    definiteness is equivalent to its conjugate being positive.  For fixed
    ``x`` these are simultaneous strict upper bounds on the identity and
    conjugate embeddings of ``z``.  An offline Fenwick dominance join counts
    ordered pairs; symmetry then converts to unordered pairs.
    """

    identity_sorted = tuple(
        sorted(values, key=cmp_to_key(lambda a, b: (a - b).sign(0)))
    )
    conjugate_sorted = tuple(
        sorted(values, key=cmp_to_key(lambda a, b: (a - b).sign(1)))
    )
    conjugate_rank = {value: index for index, value in enumerate(conjugate_sorted)}
    queries = []
    for first in values:
        identity_bound = 4 + 4 * finite / (first - 4)
        conjugate_bound = 4 - 4 * finite / (4 - first)
        queries.append((identity_bound, conjugate_bound))
    queries.sort(key=cmp_to_key(lambda a, b: (a[0] - b[0]).sign(0)))

    tree = _Fenwick(len(values))
    inserted = 0
    ordered = 0
    for identity_bound, conjugate_bound in queries:
        while (
            inserted < len(identity_sorted)
            and (identity_sorted[inserted] - identity_bound).sign(0) < 0
        ):
            tree.add(conjugate_rank[identity_sorted[inserted]])
            inserted += 1
        ordered += tree.prefix(_conjugate_lower_bound(conjugate_sorted, conjugate_bound))

    diagonal = sum(
        ((value - 4) * (value - 4) - 4 * finite).sign(0) < 0
        and ((4 - value) * (4 - value) - 4 * finite).sign(1) > 0
        for value in values
    )
    if (ordered + diagonal) % 2:
        raise AssertionError("the symmetric disjoint-path count lost parity")
    return (ordered + diagonal) // 2


@dataclass(frozen=True, slots=True)
class GoldenTwoDottedTreeRow:
    finite_label: str
    finite_product: Golden
    star_canonical_before_gram: int
    star_retained: int
    disjoint_path_canonical_before_gram: int
    disjoint_path_retained: int
    adjacent_path_canonical_before_gram: int
    adjacent_path_retained: int


@dataclass(frozen=True, slots=True)
class GoldenTwoDottedTreeCensus:
    dotted_value_count: int
    finite_nonzero_value_count: int
    dotted_conjugate_order_sha256: str
    canonical_before_gram: int
    retained_count: int
    star_retained_count: int
    disjoint_path_retained_count: int
    adjacent_path_retained_count: int
    rows: tuple[GoldenTwoDottedTreeRow, ...]


@lru_cache(maxsize=1)
def golden_two_dotted_tree_census() -> GoldenTwoDottedTreeCensus:
    values = conjugate_sorted_dotted_values()
    finite_rows = tuple(row for row in golden_finite_pair_values() if row.product)
    unordered = len(values) * (len(values) + 1) // 2
    ordered = len(values) * len(values)
    rows = tuple(
        GoldenTwoDottedTreeRow(
            finite_label=row.label,
            finite_product=row.product,
            star_canonical_before_gram=unordered,
            star_retained=count_star_pairs(values, row.product),
            disjoint_path_canonical_before_gram=unordered,
            disjoint_path_retained=count_disjoint_path_pairs(values, row.product),
            adjacent_path_canonical_before_gram=ordered,
            adjacent_path_retained=count_adjacent_path_pairs(values, row.product),
        )
        for row in finite_rows
    )
    digest = sha256()
    for value in values:
        digest.update(f"{int(value.a)},{int(value.b)}\n".encode())
    star = sum(row.star_retained for row in rows)
    adjacent = sum(row.adjacent_path_retained for row in rows)
    census = GoldenTwoDottedTreeCensus(
        dotted_value_count=len(values),
        finite_nonzero_value_count=len(finite_rows),
        dotted_conjugate_order_sha256=digest.hexdigest(),
        canonical_before_gram=sum(
            row.star_canonical_before_gram
            + row.disjoint_path_canonical_before_gram
            + row.adjacent_path_canonical_before_gram
            for row in rows
        ),
        retained_count=star + sum(row.disjoint_path_retained for row in rows) + adjacent,
        star_retained_count=star,
        disjoint_path_retained_count=sum(row.disjoint_path_retained for row in rows),
        adjacent_path_retained_count=adjacent,
        rows=rows,
    )
    expected_rows = (
        ("3", 6_634_798, 7, 17_689_960),
        ("4", 2_948_447, 5, 11_792_207),
        ("5", 9_650_285, 44, 21_334_861),
        ("6", 737_201, 2, 5_895_415),
        ("10", 5_052_550, 23, 15_436_911),
    )
    observed_rows = tuple(
        (
            row.finite_label,
            row.star_retained,
            row.disjoint_path_retained,
            row.adjacent_path_retained,
        )
        for row in rows
    )
    if observed_rows != expected_rows:
        raise AssertionError(f"the two-dotted tree rows changed: {observed_rows!r}")
    if (
        census.canonical_before_gram,
        census.retained_count,
        census.star_retained_count,
        census.adjacent_path_retained_count,
    ) != (471_865_955, 97_172_716, 25_023_281, 72_149_354):
        raise AssertionError("the two-dotted tree census changed")
    return census


__all__ = [
    "GoldenTwoDottedTreeCensus",
    "GoldenTwoDottedTreeRow",
    "conjugate_sorted_dotted_values",
    "count_adjacent_path_pairs",
    "count_disjoint_path_pairs",
    "count_star_pairs",
    "golden_two_dotted_tree_census",
]
