r"""Small exact homology audits for finite simplicial complexes.

The H9/14 nerve search represents a complex by a Boolean table indexed by
vertex-subset masks.  This module deliberately uses no computer-algebra or
topology package: boundary ranks over :math:`\mathbf F_2` are computed by
integer bit elimination.  A SAT witness can therefore be checked independently
of the search backend.

Passing the tests here is necessary for a triangulated sphere, but not a PL
sphere certificate.  In particular, mod-two homology does not detect odd
torsion or certify that a homology sphere is simply connected.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log2
from typing import Iterable, Sequence


def _vertex_count(face_table: Sequence[bool]) -> int:
    """Infer and validate the ground-set size of a subset-indexed table."""

    length = len(face_table)
    if length == 0 or length & (length - 1):
        raise ValueError("face-table length must be a positive power of two")
    count = int(log2(length))
    if 1 << count != length:  # defensive against floating conversion
        raise ValueError("face-table length is not an exact power of two")
    if not face_table[0]:
        raise ValueError("the empty set must be a face")
    return count


def gf2_rank(columns: Iterable[int]) -> int:
    """Return the rank of binary column vectors encoded as Python integers."""

    pivots: dict[int, int] = {}
    for raw_column in columns:
        if raw_column < 0:
            raise ValueError("GF(2) vectors must be nonnegative integers")
        column = raw_column
        while column:
            pivot = column.bit_length() - 1
            previous = pivots.get(pivot)
            if previous is None:
                pivots[pivot] = column
                break
            column ^= previous
    return len(pivots)


def link_face_masks(
    face_table: Sequence[bool], face_mask: int
) -> tuple[tuple[int, ...], ...]:
    """Return the nonempty faces of ``link(face_mask)``, grouped by dimension.

    The masks retain the original ground-set coordinates.  The outer tuple is
    empty exactly for the link of a facet (the conventional ``S^-1`` case).
    """

    vertex_count = _vertex_count(face_table)
    full_mask = (1 << vertex_count) - 1
    if face_mask < 0 or face_mask & ~full_mask or not face_table[face_mask]:
        raise ValueError("the requested link base is not a face")

    remaining = full_mask ^ face_mask
    grouped: list[list[int]] = []
    subset = remaining
    while subset:
        if face_table[face_mask | subset]:
            dimension = subset.bit_count() - 1
            while len(grouped) <= dimension:
                grouped.append([])
            grouped[dimension].append(subset)
        subset = (subset - 1) & remaining
    return tuple(tuple(sorted(group)) for group in grouped)


def reduced_betti_numbers_mod2(
    faces_by_dimension: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    r"""Compute reduced Betti numbers over :math:`\mathbf F_2` exactly.

    ``faces_by_dimension[k]`` contains masks for the ``k``-faces.  The input
    must be downward closed.  The returned tuple runs from degree zero through
    the largest represented dimension.  The void-link/S^-1 convention is
    represented by an empty tuple and is handled by the calling sphere audit.
    """

    if not faces_by_dimension:
        return ()
    normalized = tuple(tuple(sorted(set(group))) for group in faces_by_dimension)
    if any(mask <= 0 or mask.bit_count() != dimension + 1
           for dimension, group in enumerate(normalized)
           for mask in group):
        raise ValueError("face masks do not match their declared dimensions")

    boundary_ranks = [0] * (len(normalized) + 1)
    # boundary_ranks[k] is rank(d_k : C_k -> C_{k-1}); d_0 has rank zero.
    for dimension in range(1, len(normalized)):
        lower_index = {
            face: index for index, face in enumerate(normalized[dimension - 1])
        }
        columns: list[int] = []
        for face in normalized[dimension]:
            column = 0
            bits = face
            while bits:
                bit = bits & -bits
                boundary_face = face ^ bit
                try:
                    row = lower_index[boundary_face]
                except KeyError as exc:
                    raise ValueError("face collection is not downward closed") from exc
                column |= 1 << row
                bits ^= bit
            columns.append(column)
        boundary_ranks[dimension] = gf2_rank(columns)

    ordinary = [
        len(normalized[dimension])
        - boundary_ranks[dimension]
        - boundary_ranks[dimension + 1]
        for dimension in range(len(normalized))
    ]
    if ordinary[0] <= 0:
        raise ValueError("a nonempty simplicial complex must have positive H_0")
    ordinary[0] -= 1
    return tuple(ordinary)


@dataclass(frozen=True)
class LinkHomologyFailure:
    """The first face whose link is not a mod-two homology sphere."""

    face_mask: int
    expected_sphere_dimension: int
    observed_reduced_betti: tuple[int, ...]


def first_mod2_homology_sphere_link_failure(
    face_table: Sequence[bool], complex_dimension: int
) -> LinkHomologyFailure | None:
    """Audit the complex and every proper face link over ``GF(2)``.

    A ``d``-dimensional homology manifold without boundary has link ``S^(d-k)``
    at every face with ``k`` vertices.  Facet links are the conventional
    ``S^-1`` and are checked directly to contain no nonempty faces.
    """

    vertex_count = _vertex_count(face_table)
    if complex_dimension < 0 or complex_dimension >= vertex_count:
        raise ValueError("invalid complex dimension")

    for face_mask, is_face in enumerate(face_table):
        if not is_face or face_mask.bit_count() > complex_dimension + 1:
            continue
        expected_dimension = complex_dimension - face_mask.bit_count()
        link_faces = link_face_masks(face_table, face_mask)
        if expected_dimension == -1:
            if link_faces:
                return LinkHomologyFailure(face_mask, -1, ())
            continue
        observed = reduced_betti_numbers_mod2(link_faces)
        expected_length = expected_dimension + 1
        padded = observed + (0,) * max(0, expected_length - len(observed))
        expected = tuple(
            1 if degree == expected_dimension else 0
            for degree in range(expected_length)
        )
        if padded[:expected_length] != expected or any(
            padded[expected_length:]
        ):
            return LinkHomologyFailure(face_mask, expected_dimension, observed)
    return None
