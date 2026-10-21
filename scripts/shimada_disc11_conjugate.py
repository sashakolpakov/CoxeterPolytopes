"""Exact ample-section certificate for a conjugate H11 marking.

The fixed H11 marking misses the nef cone.  This module archives a lattice
isometry ``T`` found by an exact walk through Shimada's induced chambers and
certifies a different statement: the identity eigenspace of the conjugated
field action meets the interior of Shimada chamber 68.

Shimada uses row vectors.  Thus ``T`` is checked by ``T G T^t = G``, while
the conjugated column action is ``T^t M T^-t``.  A wall row is recorded in
dual-basis coordinates, so its inequality on a row ``x`` is simply
``x v^t >= 0``.

This is a nef/ample nonemptiness certificate.  It makes no assertion that the
resulting section is a Coxeter polytope, or that it has finite volume or is
compact.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

from sympy import Matrix

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import rank
from coxeter_polytopes.geometric.k3 import (
    K3_TRACE_SCALAR,
    field_inner_product,
    is_integral_reflection_root,
    is_k3_wall_norm,
)
from scripts.analyze_shimada_disc11 import (
    expected_shimada_gram,
    field_vector_from_shimada,
    read_shimada_data,
    shimada_basis_change,
)
from scripts.coupled_discriminant11_candidate import field_gram, multiplication_by_phi


DISC11_SOURCE_SHA256 = "fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68"
D68_CHAMBER_INDEX = 68
D68_MINUS_TWO_ORBITS = tuple(range(9))


# Row isometry taking the original identity eigenspace to one meeting D_68.
CONJUGATING_ROW_ISOMETRY: tuple[tuple[int, ...], ...] = (
    (605, 296, -180, 60, -620, -418, -836, -1240, -1024, -808, -592, -296, -620, -418, -836, -1240, -1024, -808, -592, -296),
    (-309, -151, 92, -31, 316, 213, 426, 632, 522, 412, 302, 151, 316, 213, 426, 632, 522, 412, 302, 151),
    (760, 372, -226, 75, -780, -526, -1052, -1560, -1288, -1016, -744, -372, -780, -526, -1052, -1560, -1288, -1016, -744, -372),
    (-115, -56, 35, -11, 116, 78, 156, 232, 192, 152, 112, 56, 116, 78, 156, 232, 192, 152, 112, 56),
    (122, 60, -36, 12, -126, -85, -170, -252, -208, -164, -120, -60, -125, -84, -168, -250, -206, -163, -120, -60),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0),
    (-81, -40, 24, -8, 83, 56, 112, 166, 137, 108, 79, 40, 83, 56, 112, 167, 138, 109, 80, 40),
    (-1, 0, 0, 0, 3, 2, 4, 6, 5, 4, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    (21, 10, -6, 2, -23, -16, -31, -46, -38, -30, -22, -11, -21, -14, -28, -42, -35, -28, -20, -10),
    (-20, -10, 6, -2, 20, 14, 27, 40, 33, 26, 19, 9, 21, 14, 28, 42, 35, 28, 20, 10),
    (81, 40, -24, 8, -83, -56, -112, -166, -137, -108, -79, -39, -83, -56, -112, -166, -138, -109, -80, -40),
    (122, 60, -36, 12, -125, -84, -168, -250, -206, -163, -120, -60, -126, -85, -170, -252, -208, -164, -120, -60),
    (0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (-81, -40, 24, -8, 83, 56, 112, 167, 138, 109, 80, 40, 83, 56, 112, 166, 137, 108, 79, 40),
    (-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 4, 6, 5, 4, 3, 1),
    (21, 10, -6, 2, -21, -14, -28, -42, -35, -28, -20, -10, -23, -16, -31, -46, -38, -30, -22, -11),
    (-20, -10, 6, -2, 21, 14, 28, 42, 35, 28, 20, 10, 20, 14, 27, 40, 33, 26, 19, 9),
    (81, 40, -24, 8, -83, -56, -112, -166, -138, -109, -80, -40, -83, -56, -112, -166, -137, -108, -79, -39),
)


AMPLE_FIELD_COEFFICIENTS: tuple[Golden, ...] = (
    Golden(-54, 26),
    Golden(-1),
    Golden(3, -3),
    Golden(3, -3),
    Golden(5, -5),
    Golden(7, -7),
    Golden(11, -11),
    Golden(12, -12),
    Golden(15, -15),
    Golden(16, -16),
)


# Exact points on the eight root facets of the restricted D_68 section.
# At the witness for orbit i, both ambient members of orbit i vanish and all
# other archived D_68 wall inequalities are strictly positive.
D68_ROOT_FACET_WITNESSES: tuple[tuple[int, tuple[Golden, ...]], ...] = (
    (1, (
        Golden(-46, 22), Golden(-1), Golden(3, -3), Golden(3, -3),
        Golden(5, -5), Golden(7, -7), Golden(10, -10), Golden(11, -11),
        Golden(13, -13), Golden(14, -14),
    )),
    (2, (
        Golden(-54, 26), Golden(-1), Golden(1, -1), Golden(0),
        Golden(0), Golden(3, -3), Golden(8, -8), Golden(10, -10),
        Golden(14, -14), Golden(14, -14),
    )),
    (3, (
        Golden(-50, 24), Golden(-1), Golden(2, -2), Golden(2, -2),
        Golden(4, -4), Golden(6, -6), Golden(10, -10), Golden(11, -11),
        Golden(14, -14), Golden(15, -15),
    )),
    (4, (
        Golden(-48, 23), Golden(-1), Golden(4, -4), Golden(5, -5),
        Golden(8, -8), Golden(9, -9), Golden(12, -12), Golden(12, -12),
        Golden(14, -14), Golden(16, -16),
    )),
    (5, (
        Golden(-44, 21), Golden(-1), Golden(2, -2), Golden(1, -1),
        Golden(2, -2), Golden(4, -4), Golden(8, -8), Golden(9, -9),
        Golden(12, -12), Golden(12, -12),
    )),
    (6, (
        Golden(-52, 25), Golden(-1), Golden(1, -1), Golden(0),
        Golden(0), Golden(3, -3), Golden(8, -8), Golden(10, -10),
        Golden(14, -14), Golden(13, -13),
    )),
    (7, (
        Golden(-96, 47), Golden(-1), Golden(-6, 6), Golden(-15, 15),
        Golden(-22, 22), Golden(-12, 12), Golden(0), Golden(10, -10),
        Golden(21, -21), Golden(13, -13),
    )),
    (8, (
        Golden(-110, 54), Golden(-1), Golden(-7, 7), Golden(-17, 17),
        Golden(-25, 25), Golden(-13, 13), Golden(0), Golden(11, -11),
        Golden(24, -24), Golden(15, -15),
    )),
)


# Exact positive points on the three non-root facets of the restricted D_68
# section.  The first facet is the simultaneous restriction of two distinct
# ambient Borcherds walls, (9, 0) and (16, 0); their restricted functionals
# differ by the positive scalar phi.
D68_INTERNAL_FACET_WITNESSES: tuple[
    tuple[tuple[tuple[int, int], ...], tuple[Golden, ...]], ...
] = (
    (
        ((9, 0), (16, 0)),
        (
            Golden(-109, 53) / 2,
            Golden(0),
            Golden(3, -3),
            Golden(3, -3),
            Golden(5, -5),
            Golden(7, -7),
            Golden(11, -11),
            Golden(12, -12),
            Golden(15, -15),
            Golden(16, -16),
        ),
    ),
    (
        ((11, 1),),
        (
            Golden(-325, 156) / 6,
            Golden(-8, 1) / 6,
            Golden(13, -16) / 6,
            Golden(13, -16) / 6,
            Golden(25, -28) / 6,
            Golden(37, -40) / 6,
            Golden(61, -64) / 6,
            Golden(67, -70) / 6,
            Golden(85, -88) / 6,
            Golden(91, -94) / 6,
        ),
    ),
    (
        ((14, 0),),
        (
            Golden(-80, 24) / 3,
            Golden(52, -40) / 3,
            Golden(3, -3),
            Golden(3, -3),
            Golden(5, -5),
            Golden(7, -7),
            Golden(11, -11),
            Golden(12, -12),
            Golden(15, -15),
            Golden(16, -16),
        ),
    ),
)


D68_RESTRICTED_FACET_REPRESENTATIVES: tuple[tuple[int, int], ...] = (
    *((orbit_index, 0) for orbit_index in range(1, 9)),
    (9, 0),
    (11, 1),
    (14, 0),
)


# Positive exact dependences for every distinct restricted half-space that is
# not one of the eleven facets above.  Entries express the target functional
# as a positive O_F-linear combination of facet representatives.  Positive
# proportional copies inside facet classes are checked separately.
D68_RESTRICTED_REDUNDANCY_DEPENDENCES: tuple[
    tuple[
        tuple[int, int],
        tuple[tuple[Golden, tuple[int, int]], ...],
    ],
    ...,
] = (
    ((0, 0), (
        (PHI, (1, 0)), (2 + 2 * PHI, (2, 0)),
        (2 * PHI, (3, 0)), (1 + PHI, (4, 0)),
        (2 + 2 * PHI, (5, 0)), (2 + 2 * PHI, (6, 0)),
        (2 * PHI, (7, 0)), (1 + 2 * PHI, (9, 0)),
        (2 * PHI, (11, 1)),
    )),
    ((10, 0), ((PHI, (9, 0)), (PHI - 1, (14, 0)))),
    ((11, 0), (
        (Golden(8, -4), (1, 0)), (Golden(-3, 3), (2, 0)),
        (Golden(4, -2), (3, 0)), (Golden(-3, 3), (4, 0)),
        (Golden(-5, 5), (5, 0)), (Golden(-4, 4), (6, 0)),
        (Golden(12, -6), (7, 0)), (Golden(4, -2), (8, 0)),
        (Golden(1), (11, 1)),
    )),
    ((12, 0), (
        (Golden(4, 2), (1, 0)), (Golden(4, 6), (2, 0)),
        (Golden(1, 4), (3, 0)), (Golden(2, 4), (4, 0)),
        (Golden(4, 8), (5, 0)), (Golden(4, 7), (6, 0)),
        (Golden(5, 4), (7, 0)), (Golden(2, 1), (8, 0)),
        (4 * PHI, (11, 1)),
    )),
    ((12, 1), (
        (Golden(-4, 6), (1, 0)), (Golden(6, 4), (2, 0)),
        (Golden(-1, 5), (3, 0)), (Golden(4, 2), (4, 0)),
        (Golden(8, 4), (5, 0)), (Golden(7, 4), (6, 0)),
        (Golden(-5, 9), (7, 0)), (Golden(-2, 3), (8, 0)),
        (4 * PHI, (11, 1)),
    )),
    ((13, 0), (
        (Golden(-2, 4), (1, 0)), (Golden(3, 3), (2, 0)),
        (3 * PHI, (3, 0)), (Golden(3, 1), (4, 0)),
        (Golden(5, 3), (5, 0)), (Golden(4, 3), (6, 0)),
        (Golden(-3, 6), (7, 0)), (Golden(-1, 2), (8, 0)),
        (PHI, (9, 0)), (3 * PHI, (11, 1)),
    )),
    ((13, 1), (
        (Golden(2, 2), (1, 0)), (Golden(3, 3), (2, 0)),
        (3 * PHI, (3, 0)), (Golden(1, 3), (4, 0)),
        (Golden(3, 5), (5, 0)), (Golden(3, 4), (6, 0)),
        (Golden(3, 3), (7, 0)), (Golden(1, 1), (8, 0)),
        (PHI, (9, 0)), (3 * PHI, (11, 1)),
    )),
    ((15, 0), (
        (Golden(8, -2), (1, 0)), (6 * PHI, (2, 0)),
        (Golden(4, 1), (3, 0)), (Golden(-2, 4), (4, 0)),
        (Golden(-2, 8), (5, 0)), (Golden(-1, 7), (6, 0)),
        (Golden(12, -3), (7, 0)), (Golden(4, -1), (8, 0)),
        (Golden(1, 1), (9, 0)), (Golden(2, 2), (11, 1)),
    )),
)


# Exact positive dependence proving that the orbit-0 root inequality is
# redundant after restriction to the transported eigenspace.
D68_ORBIT0_RESTRICTED_DEPENDENCE: tuple[
    tuple[Golden, tuple[int, int]], ...
] = (
    (PHI, (1, 0)),
    (2 + 2 * PHI, (2, 0)),
    (2 * PHI, (3, 0)),
    (1 + PHI, (4, 0)),
    (2 + 2 * PHI, (5, 0)),
    (2 + 2 * PHI, (6, 0)),
    (2 * PHI, (7, 0)),
    (1 + 2 * PHI, (9, 0)),
    (2 * PHI, (11, 1)),
)


# Each defining wall is a tuple of its nonzero (zero-based index, value)
# pairs.  The outer keys and inner ordering are exactly those of Shimada's
# chamberdata[68][wallorbit].
D68_WALL_ORBITS_SPARSE: tuple[
    tuple[int, tuple[tuple[tuple[int, int], ...], ...]], ...
] = (
    (0, (((0, 1), (1, -2)),)),
    (1, (((10, 1), (11, -2)), ((18, 1), (19, -2)))),
    (2, (((4, 1), (6, 1), (7, -2), (8, 1)), ((12, 1), (14, 1), (15, -2), (16, 1)))),
    (3, (((12, -2), (15, 1)), ((4, -2), (7, 1)))),
    (4, (((13, -2), (14, 1)), ((5, -2), (6, 1)))),
    (5, (((8, 1), (9, -2), (10, 1)), ((16, 1), (17, -2), (18, 1)))),
    (6, (((7, 1), (8, -2), (9, 1)), ((15, 1), (16, -2), (17, 1)))),
    (7, (
        ((0, 10), (2, 10), (3, -6), (5, 1), (10, 2), (11, 1), (17, 1), (18, 2)),
        ((0, 10), (2, 10), (3, -6), (9, 1), (10, 2), (13, 1), (18, 2), (19, 1)),
    )),
    (8, (((1, 1), (11, 1)), ((1, 1), (19, 1)))),
    (9, (((0, -13), (2, -12), (3, 8), (6, -1), (10, -3), (14, -1), (18, -3)),)),
    (10, (((2, -1), (3, -3)),)),
    (11, (
        ((0, 3), (2, 3), (3, -1), (4, 1), (18, 1)),
        ((0, 3), (2, 3), (3, -1), (10, 1), (12, 1)),
    )),
    (12, (
        ((0, 2), (2, 2), (3, -2), (5, 1), (19, 1)),
        ((0, 2), (2, 2), (3, -2), (11, 1), (13, 1)),
    )),
    (13, (
        ((0, 13), (2, 12), (3, -8), (5, 1), (6, -1), (7, 1), (10, 3), (14, 1), (18, 3)),
        ((0, 13), (2, 12), (3, -8), (6, 1), (10, 3), (13, 1), (14, -1), (15, 1), (18, 3)),
    )),
    (14, (((0, 13), (1, 1), (2, 14), (3, -7), (6, 1), (10, 3), (14, 1), (18, 3)),)),
    (15, (((0, 8), (2, 8), (3, -3), (5, 1), (10, 2), (13, 1), (18, 2)),)),
    (16, (((0, 19), (2, 18), (3, -11), (6, 1), (10, 5), (14, 1), (18, 5)),)),
)


ARCHIVED_AMPLE_VECTOR: tuple[Golden, ...] = tuple(
    Golden(a, b)
    for a, b in (
        (1005, 1470), (496, 727), (-303, -426), (83, 151),
        (-1062, -1498), (-718, -1007), (-1436, -2013), (-2123, -2996),
        (-1747, -2481), (-1370, -1966), (-992, -1451), (-496, -726),
        (-1030, -1530), (-694, -1031), (-1387, -2062), (-2060, -3059),
        (-1703, -2525), (-1346, -1990), (-989, -1454), (-495, -727),
    )
)


def _integer_matrix(matrix: Matrix) -> list[list[int]]:
    if any(not entry.is_Integer for entry in matrix):
        raise ValueError("expected an integer matrix")
    return [[int(entry) for entry in row] for row in matrix.tolist()]


@lru_cache(maxsize=1)
def conjugating_row_isometry_inverse() -> tuple[tuple[int, ...], ...]:
    """Return and cache the integral inverse of the archived transport."""

    inverse = _integer_matrix(Matrix(CONJUGATING_ROW_ISOMETRY).inv())
    return tuple(tuple(row) for row in inverse)


def _row_product(row, matrix):
    return [
        sum((row[i] * matrix[i][j] for i in range(len(row))), Golden(0))
        for j in range(len(matrix[0]))
    ]


def _pair(left, gram, right):
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(len(left))
            for j in range(len(right))
        ),
        Golden(0),
    )


def d68_wall_orbits() -> tuple[tuple[int, tuple[tuple[int, ...], ...]], ...]:
    """Expand the archived sparse D_68 defining wall functionals."""

    result = []
    for orbit_index, sparse_orbit in D68_WALL_ORBITS_SPARSE:
        orbit = []
        for sparse_wall in sparse_orbit:
            wall = [0] * 20
            for coordinate, value in sparse_wall:
                wall[coordinate] = value
            orbit.append(tuple(wall))
        result.append((orbit_index, tuple(orbit)))
    return tuple(result)


def original_phi_action_shimada() -> list[list[int]]:
    """Return the original column action of phi in Shimada coordinates."""

    change = Matrix(shimada_basis_change())
    action = change.inv() * Matrix(multiplication_by_phi(10)) * change
    return _integer_matrix(action)


def conjugated_phi_action_shimada() -> list[list[int]]:
    """Return ``T^t M T^-t``, still as a column action."""

    transport = Matrix(CONJUGATING_ROW_ISOMETRY)
    action = transport.T * Matrix(original_phi_action_shimada()) * transport.T.inv()
    return _integer_matrix(action)


def transported_identity_eigenbasis() -> list[list[Golden]]:
    """Return a 10-by-20 row basis for the transported identity eigenspace."""

    inverse_change = Matrix(shimada_basis_change()).inv()
    original_basis = [
        [
            Golden(int(inverse_change[r, 2 * k]), int(inverse_change[r, 2 * k + 1]))
            for r in range(20)
        ]
        for k in range(10)
    ]
    transport = CONJUGATING_ROW_ISOMETRY
    return [
        [
            sum(
                (original_basis[k][r] * transport[r][j] for r in range(20)),
                Golden(0),
            )
            for j in range(20)
        ]
        for k in range(10)
    ]


def ample_witness() -> list[Golden]:
    """Generate the exact row vector lying strictly inside D_68."""

    basis = transported_identity_eigenbasis()
    return [
        sum(
            (AMPLE_FIELD_COEFFICIENTS[k] * basis[k][j] for k in range(10)),
            Golden(0),
        )
        for j in range(20)
    ]


def pullback_shimada_row(row: tuple[int, ...] | list[int]) -> list[int]:
    """Pull a physical divisor back through ``T`` to the old field marking."""

    if len(row) != 20:
        raise ValueError("a Shimada row must have 20 coordinates")
    inverse = conjugating_row_isometry_inverse()
    return [
        sum(row[i] * inverse[i][j] for i in range(20)) for j in range(20)
    ]


def conjugated_field_vector_from_shimada(
    row: tuple[int, ...] | list[int],
) -> list[Golden]:
    """Read a physical divisor in the conjugated ``O_F`` structure."""

    return field_vector_from_shimada(pullback_shimada_row(row))


def d68_minus_two_roots() -> tuple[dict[str, object], ...]:
    """Recover the 17 actual ``(-2)`` boundary roots of D_68.

    In Shimada's data the adjacent tag is ``_minustwo`` precisely for D_68
    wall orbits 0 through 8.  A stored wall row is a dual-basis functional;
    multiplication by ``G^-1`` recovers its primal root row.
    """

    gram = Matrix(expected_shimada_gram())
    inverse_gram = gram.inv()
    result = []
    for orbit_index, orbit in d68_wall_orbits():
        if orbit_index not in D68_MINUS_TWO_ORBITS:
            continue
        for member_index, wall in enumerate(orbit):
            root = _integer_matrix(Matrix(1, 20, wall) * inverse_gram)[0]
            square = (Matrix(1, 20, root) * gram * Matrix(20, 1, root))[0]
            if square != -2:
                raise ValueError("a D_68 _minustwo wall did not recover a (-2)-root")
            result.append(
                {
                    "orbit_index": orbit_index,
                    "member_index": member_index,
                    "dual_wall": list(wall),
                    "primal_root": root,
                }
            )
    return tuple(result)


def analyze_conjugated_roots(
    roots: list[list[int]] | tuple[tuple[int, ...], ...],
) -> dict[str, object]:
    """Apply the five-shell and reflection filters to physical ``(-2)`` rows."""

    h = field_gram()
    gram = expected_shimada_gram()
    witness = ample_witness()
    records = []
    for index, root in enumerate(roots):
        vector = conjugated_field_vector_from_shimada(root)
        norm = field_inner_product(h, vector, vector)
        ample_pairing = _pair(witness, gram, root)
        eigen_pairing = K3_TRACE_SCALAR * field_inner_product(
            h, AMPLE_FIELD_COEFFICIENTS, vector
        )
        records.append(
            {
                "index": index,
                "root": list(root),
                "field_vector": vector,
                "norm": norm,
                "trace_square": (K3_TRACE_SCALAR * norm).trace(),
                "meets_eigen_positive_cone": is_k3_wall_norm(norm),
                "integral_field_reflection": is_integral_reflection_root(h, vector),
                "ample_pairing": ample_pairing,
                "eigen_pairing": eigen_pairing,
            }
        )

    relevant = [record for record in records if record["meets_eigen_positive_cone"]]
    integral = [
        record for record in relevant if record["integral_field_reflection"]
    ]
    pair_signs: Counter[str] = Counter()
    for i, left in enumerate(integral):
        for right in integral[i + 1 :]:
            sign = field_inner_product(
                h, left["field_vector"], right["field_vector"]
            ).sign()
            pair_signs[{1: "positive", 0: "zero", -1: "negative"}[sign]] += 1

    ample_pairings = [record["ample_pairing"] for record in records]
    eigen_pairings = [record["eigen_pairing"] for record in records]
    return {
        "input_root_count": len(records),
        "trace_minus_two_verified": sum(
            record["trace_square"] == -2 for record in records
        ),
        "positive_at_both_places": len(relevant),
        "positive_shell_distribution": dict(
            sorted(Counter(str(record["norm"]) for record in relevant).items())
        ),
        "field_integral_reflection_roots": len(integral),
        "field_integral_shell_distribution": dict(
            sorted(Counter(str(record["norm"]) for record in integral).items())
        ),
        "field_span_rank": rank([record["field_vector"] for record in integral]),
        "pair_sign_distribution": dict(sorted(pair_signs.items())),
        "surviving_zero_based_indices": [record["index"] for record in integral],
        "all_ample_pairings_positive": all(
            pairing.sign() > 0 for pairing in ample_pairings
        ),
        "minimum_ample_pairing": (
            min(ample_pairings).to_json() if ample_pairings else None
        ),
        "all_eigen_pairings_positive": all(
            pairing.sign() > 0 for pairing in eigen_pairings
        ),
        "integral_root_data": [
            {
                "index": record["index"],
                "shimada_coordinates": record["root"],
                "field_coordinates": [
                    entry.to_json() for entry in record["field_vector"]
                ],
                "field_norm": record["norm"].to_json(),
                "ample_pairing": record["ample_pairing"].to_json(),
                "eigen_pairing": record["eigen_pairing"].to_json(),
            }
            for record in integral
        ],
    }


def analyze_d68_minus_two_roots() -> dict[str, object]:
    """Filter the 17 archived geometric boundary roots of D_68."""

    root_data = d68_minus_two_roots()
    result = analyze_conjugated_roots(
        tuple(tuple(item["primal_root"]) for item in root_data)
    )
    result["wall_labels"] = [
        [item["orbit_index"], item["member_index"]] for item in root_data
    ]
    return result


def restricted_d68_wall_coefficients() -> dict[
    tuple[int, int], tuple[Golden, ...]
]:
    """Restrict every archived D_68 dual wall to the transported eigenspace."""

    basis = transported_identity_eigenbasis()
    return {
        (orbit_index, member_index): tuple(
            sum((basis[k][i] * wall[i] for i in range(20)), Golden(0))
            for k in range(10)
        )
        for orbit_index, orbit in d68_wall_orbits()
        for member_index, wall in enumerate(orbit)
    }


def verify_d68_root_facet_certificate() -> dict[str, object]:
    """Verify the exact local root-facet structure of ``D_68 cap E``.

    The 16 roots in two-member orbits 1 through 8 restrict in positive
    proportional pairs.  Eight archived positive-norm face points prove that
    the resulting eight hyperplanes are genuine facets.  The remaining root
    orbit 0 is redundant by an archived positive linear dependence.
    """

    coefficients = restricted_d68_wall_coefficients()
    h = field_gram()
    facet_data = []
    all_facets_verified = True
    for orbit_index, witness in D68_ROOT_FACET_WITNESSES:
        first = coefficients[(orbit_index, 0)]
        second = coefficients[(orbit_index, 1)]
        pivot = next(i for i, entry in enumerate(first) if entry)
        ratio = second[pivot] / first[pivot]
        proportional = second == tuple(ratio * entry for entry in first)
        pairings = {
            label: sum(
                (coefficient[k] * witness[k] for k in range(10)), Golden(0)
            )
            for label, coefficient in coefficients.items()
        }
        target_zero = pairings[(orbit_index, 0)] == pairings[(orbit_index, 1)] == 0
        other_positive = all(
            pairing.sign() > 0
            for label, pairing in pairings.items()
            if label[0] != orbit_index
        )
        q_square = K3_TRACE_SCALAR * field_inner_product(h, witness, witness)
        verified = (
            proportional
            and ratio.sign() > 0
            and target_zero
            and other_positive
            and q_square.sign() > 0
        )
        all_facets_verified &= verified
        facet_data.append(
            {
                "orbit_index": orbit_index,
                "ambient_root_walls": [[orbit_index, 0], [orbit_index, 1]],
                "positive_proportionality_ratio": ratio.to_json(),
                "face_witness": [entry.to_json() for entry in witness],
                "face_witness_q_square": q_square.to_json(),
                "target_pairings_zero": target_zero,
                "all_other_D68_pairings_positive": other_positive,
            }
        )

    dependence_sum = tuple(
        sum(
            (
                scalar * coefficients[label][k]
                for scalar, label in D68_ORBIT0_RESTRICTED_DEPENDENCE
            ),
            Golden(0),
        )
        for k in range(10)
    )
    orbit_zero_redundant = (
        all(scalar.sign() > 0 for scalar, _ in D68_ORBIT0_RESTRICTED_DEPENDENCE)
        and dependence_sum == coefficients[(0, 0)]
    )
    return {
        "verified": all_facets_verified and orbit_zero_redundant,
        "ambient_minus_two_wall_count": 17,
        "ambient_roots_on_restricted_facets": 16,
        "distinct_restricted_root_facets": 8,
        "redundant_restricted_root_walls": [[0, 0]],
        "orbit_zero_positive_dependence_verified": orbit_zero_redundant,
        "facets": facet_data,
        "scope": "local_facets_of_D68_section_not_global_nef_wall_completeness",
    }


def verify_d68_restricted_facet_certificate() -> dict[str, object]:
    """Certify the complete irredundant half-space list of ``D_68 cap E``.

    This is stronger than the root-only audit above.  It partitions all 28
    ambient walls into 19 positive-proportional restricted half-spaces, gives
    positive face witnesses for exactly eleven of them, and supplies positive
    exact dependences for the remaining eight.  Hence the local eigensection
    has eight geometric ``(-2)`` facets and only three internal Borcherds
    facets.
    """

    coefficients = restricted_d68_wall_coefficients()
    h = field_gram()

    proportional_classes: list[list[tuple[int, int]]] = []
    for label, row in coefficients.items():
        for wall_class in proportional_classes:
            base = coefficients[wall_class[0]]
            pivot = next(index for index, entry in enumerate(base) if entry)
            ratio = row[pivot] / base[pivot]
            if ratio.sign() > 0 and row == tuple(ratio * entry for entry in base):
                wall_class.append(label)
                break
        else:
            proportional_classes.append([label])

    facet_classes = [
        ((orbit_index, 0), (orbit_index, 1))
        for orbit_index in range(1, 9)
    ] + [
        ((9, 0), (16, 0)),
        ((11, 1),),
        ((14, 0),),
    ]
    facet_class_set = {frozenset(wall_class) for wall_class in facet_classes}
    computed_class_set = {
        frozenset(wall_class) for wall_class in proportional_classes
    }

    witness_records = []
    witnesses = [
        (
            ((orbit_index, 0), (orbit_index, 1)),
            witness,
        )
        for orbit_index, witness in D68_ROOT_FACET_WITNESSES
    ] + list(D68_INTERNAL_FACET_WITNESSES)
    all_witnesses_verified = True
    for target_class, witness in witnesses:
        pairings = {
            label: sum(
                (row[index] * witness[index] for index in range(10)),
                Golden(0),
            )
            for label, row in coefficients.items()
        }
        target_set = set(target_class)
        target_zero = all(pairings[label] == 0 for label in target_class)
        others_positive = all(
            pairing.sign() > 0
            for label, pairing in pairings.items()
            if label not in target_set
        )
        q_square = K3_TRACE_SCALAR * field_inner_product(h, witness, witness)
        verified = target_zero and others_positive and q_square.sign() > 0
        all_witnesses_verified &= verified
        witness_records.append(
            {
                "restricted_wall_class": [list(label) for label in target_class],
                "face_witness": [entry.to_json() for entry in witness],
                "face_witness_q_square": q_square.to_json(),
                "target_pairings_zero": target_zero,
                "all_other_D68_pairings_positive": others_positive,
            }
        )

    dependence_records = []
    all_dependences_verified = True
    for target, terms in D68_RESTRICTED_REDUNDANCY_DEPENDENCES:
        dependence_sum = tuple(
            sum(
                (scalar * coefficients[label][index] for scalar, label in terms),
                Golden(0),
            )
            for index in range(10)
        )
        positive = all(scalar.sign() > 0 for scalar, _ in terms)
        verified = positive and dependence_sum == coefficients[target]
        all_dependences_verified &= verified
        dependence_records.append(
            {
                "redundant_wall": list(target),
                "terms": [
                    {
                        "positive_coefficient": scalar.to_json(),
                        "facet_wall": list(label),
                    }
                    for scalar, label in terms
                ],
                "verified": verified,
            }
        )

    redundant_classes = computed_class_set - facet_class_set
    expected_redundant_classes = {
        frozenset((target,))
        for target, _ in D68_RESTRICTED_REDUNDANCY_DEPENDENCES
    }
    classes_verified = (
        len(proportional_classes) == 19
        and len(facet_classes) == 11
        and facet_class_set <= computed_class_set
        and redundant_classes == expected_redundant_classes
    )
    return {
        "verified": (
            classes_verified
            and all_witnesses_verified
            and all_dependences_verified
        ),
        "ambient_wall_count": len(coefficients),
        "distinct_restricted_halfspaces": len(proportional_classes),
        "irredundant_restricted_facets": len(facet_classes),
        "root_facets": 8,
        "internal_Borcherds_facets": 3,
        "internal_facet_wall_classes": [
            [[9, 0], [16, 0]],
            [[11, 1]],
            [[14, 0]],
        ],
        "simultaneous_internal_crossing": [[9, 0], [16, 0]],
        "positive_face_witnesses": witness_records,
        "positive_redundancy_dependences": dependence_records,
        "scope": "complete_local_facets_of_D68_section_not_global_eigen_nef_completeness",
    }


def analyze_conjugated_disc11_file(path: Path) -> dict[str, object]:
    """Filter all 347 published curve representatives from external data."""

    gram, roots, _ = read_shimada_data(path)
    if gram != expected_shimada_gram():
        raise ValueError("unexpected disc11 Gram matrix")
    result = analyze_conjugated_roots(roots)
    result["source_file"] = str(path)
    result["source_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def read_disc11_d68_wall_orbits(
    path: Path,
) -> tuple[tuple[int, tuple[tuple[int, ...], ...]], ...]:
    """Extract D_68 wallorbit assignments from Shimada's external data."""

    prefix = "chamberdata[68][wallorbit]["
    result = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.replace("\\", "").strip()
            if not line.startswith(prefix):
                continue
            name, value = line.split(":=", 1)
            orbit_index = int(name[len(prefix) : -1])
            if not value.endswith(";"):
                raise ValueError("unexpected wrapped D_68 wallorbit assignment")
            orbit = ast.literal_eval(value[:-1])
            result.append((orbit_index, tuple(tuple(wall) for wall in orbit)))
    return tuple(sorted(result))


def verify_conjugate_marking_certificate() -> dict[str, object]:
    """Verify the archived nonempty ample-section certificate exactly."""

    gram = expected_shimada_gram()
    transport = Matrix(CONJUGATING_ROW_ISOMETRY)
    row_isometry = transport * Matrix(gram) * transport.T == Matrix(gram)
    # det(G)=-11.  The third column of G^-1 generates the cyclic
    # discriminant group; the column action T^t fixes it modulo the lattice.
    discriminant_generator = Matrix(gram).inv()[:, 2]
    discriminant_difference = transport.T * discriminant_generator - discriminant_generator
    discriminant_action_trivial = all(
        entry.is_Integer for entry in discriminant_difference
    )
    witness = ample_witness()
    action = conjugated_phi_action_shimada()
    action_matrix = Matrix(action)
    phi_polynomial = action_matrix * action_matrix == action_matrix + Matrix.eye(20)
    self_adjoint = action_matrix.T * Matrix(gram) == Matrix(gram) * action_matrix
    action_transpose = [list(row) for row in zip(*action)]
    eigenvector = _row_product(witness, action_transpose) == [
        PHI * entry for entry in witness
    ]
    square = _pair(witness, gram, witness)
    indexed_pairings = [
        (orbit_index, member_index, sum((witness[i] * wall[i] for i in range(20)), Golden(0)))
        for orbit_index, orbit in d68_wall_orbits()
        for member_index, wall in enumerate(orbit)
    ]
    minimum = min(pairing for _, _, pairing in indexed_pairings)
    minimum_witnesses = [
        [orbit_index, member_index]
        for orbit_index, member_index, pairing in indexed_pairings
        if pairing == minimum
    ]
    return {
        "verified": (
            transport.det() == 1
            and row_isometry
            and discriminant_action_trivial
            and witness == list(ARCHIVED_AMPLE_VECTOR)
            and phi_polynomial
            and self_adjoint
            and eigenvector
            and square == Golden(-2474, 2688)
            and square.sign() > 0
            and len(indexed_pairings) == 28
            and all(pairing.sign() > 0 for _, _, pairing in indexed_pairings)
            and minimum == Golden(1)
        ),
        "conclusion": "conjugated_identity_eigenspace_meets_D68_interior",
        "scope": "nef_ample_nonemptiness_only_not_a_Coxeter_or_compactness_claim",
        "row_isometry_determinant": str(transport.det()),
        "row_isometry_verified": row_isometry,
        "discriminant_action_trivial_verified": discriminant_action_trivial,
        "phi_polynomial_verified": phi_polynomial,
        "phi_self_adjoint_verified": self_adjoint,
        "phi_eigenvector_verified": eigenvector,
        "witness_square": square.to_json(),
        "defining_wall_count": len(indexed_pairings),
        "all_wall_pairings_strictly_positive": all(
            pairing.sign() > 0 for _, _, pairing in indexed_pairings
        ),
        "minimum_wall_pairing": minimum.to_json(),
        "minimum_pairing_walls": minimum_witnesses,
    }


def external_archive_check(path: Path) -> dict[str, object]:
    """Compare the archive with a downloaded copy of ``disc11.txt``."""

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "sha256": digest,
        "expected_sha256": DISC11_SOURCE_SHA256,
        "sha256_matches": digest == DISC11_SOURCE_SHA256,
        "d68_wall_orbits_match": read_disc11_d68_wall_orbits(path)
        == d68_wall_orbits(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "disc11",
        nargs="?",
        type=Path,
        help="optional path to Shimada's external disc11.txt",
    )
    parser.add_argument(
        "--full-roots",
        action="store_true",
        help="retain explicit integral-root records in JSON output",
    )
    args = parser.parse_args()

    local_filter = analyze_d68_minus_two_roots()
    if not args.full_roots:
        local_filter.pop("integral_root_data")
    result: dict[str, object] = {
        "conjugate_marking_certificate": verify_conjugate_marking_certificate(),
        "local_D68_minus_two_filter": local_filter,
        "local_D68_root_facets": verify_d68_root_facet_certificate(),
        "local_D68_all_facets": verify_d68_restricted_facet_certificate(),
    }
    if args.disc11 is not None:
        published_filter = analyze_conjugated_disc11_file(args.disc11)
        if not args.full_roots:
            published_filter.pop("integral_root_data")
        result["external_archive_check"] = external_archive_check(args.disc11)
        result["conjugated_published_B_filter"] = published_filter
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
