"""Exact first-shell feasibility audit for ``<-2*theta> + Lambda_24``.

Here ``theta = 2*cos(2*pi/7)`` and ``Lambda_24`` is the positive-definite
Leech lattice with minimum four.  The spatial lattice is even and unimodular.
Consequently every primitive reflection root has canonical norm two; unlike a
root-block Vinberg search, however, there are no height-zero roots from which
to choose a Weyl chamber.

This module deliberately stops at the first positive control shell.  It uses
the standard extended-Golay description of the 196,560 minimal Leech vectors
and an exact trace decomposition to prove that none has the required
field-valued square.  Continuing requires simultaneous Co0 orbit enumeration
of ordered Leech triples, not merely the ordinary one-vector theta series.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from typing import Iterator

from coxeter_polytopes.exact.cubic7 import CUBIC7_ZERO, THETA7, Cubic7
from coxeter_polytopes.search.cubic7_root_blocks import _bounded_for_caps


LEECH_MINIMUM = 4
LEECH_MINIMAL_VECTOR_COUNT = 196_560
CO0_ORDER = 8_315_553_613_086_720_000
CO2_ORDER = 42_305_421_312_000

TIMELIKE_COEFFICIENT = 2 * THETA7
CANONICAL_MIRROR_NORMS = (Cubic7(2),)

COMPLETE_MIRROR_NORM_PROOF = (
    "The lattice is even, so every root norm is divisible by the inert prime "
    "(2).  Away from (2), spatial unimodularity and primitivity make the "
    "pairing ideal a unit.  At (2), if the spatial component is not divisible "
    "by 2 then reflection integrality forces valuation at most one.  If it is "
    "divisible by 2, the control is a unit and q/2 is congruent to "
    "-theta*x^2 modulo 2, again giving valuation exactly one.  Thus the norm "
    "ideal is exactly (2).  Class number one and the full unit-signature map "
    "make every totally positive unit a square, leaving canonical norm 2."
)


# Tr(theta^(i+j)) in the integral power basis (1, theta, theta^2).
TRACE_GRAM = (
    (3, -1, 5),
    (-1, 5, -4),
    (5, -4, 13),
)

# If (lambda_0,lambda_1,lambda_2)^t = U (A,B,C)^t, then
# y=lambda_0+theta*lambda_1+theta^2*lambda_2
#   = -A + (-1+theta+theta^2)B - theta*C.
TRACE_CHANGE_OF_BASIS = (
    (-1, -1, 0),
    (0, 1, -1),
    (0, 1, 0),
)

REDUCED_TRACE_GRAM = (
    (3, -1, -1),
    (-1, 5, -2),
    (-1, -2, 5),
)

TRACE_DECOMPOSITION = (
    # coefficient, integral linear form in (A,B,C)
    (1, (1, -1, 0)),  # ||A-B||^2
    (1, (1, 0, -1)),  # ||A-C||^2
    (2, (0, 1, -1)),  # 2||B-C||^2
    (1, (1, 0, 0)),   # ||A||^2
    (2, (0, 1, 0)),   # 2||B||^2
    (2, (0, 0, 1)),   # 2||C||^2
)


def _transpose(matrix: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(matrix[row][column] for row in range(len(matrix))) for column in range(len(matrix[0])))


def _multiply(
    left: tuple[tuple[int, ...], ...], right: tuple[tuple[int, ...], ...]
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def reduced_trace_gram_from_change_of_basis() -> tuple[tuple[int, ...], ...]:
    """Return ``U^t T U`` for the archived integral basis change."""

    return _multiply(
        _multiply(_transpose(TRACE_CHANGE_OF_BASIS), TRACE_GRAM),
        TRACE_CHANGE_OF_BASIS,
    )


def trace_gram_from_norm_decomposition() -> tuple[tuple[int, ...], ...]:
    """Reconstruct the reduced trace Gram from six Leech norm terms."""

    gram = [[0, 0, 0] for _ in range(3)]
    for coefficient, linear_form in TRACE_DECOMPOSITION:
        for i in range(3):
            for j in range(3):
                gram[i][j] += coefficient * linear_form[i] * linear_form[j]
    return tuple(tuple(row) for row in gram)


def trace_minimum_and_equality_case() -> dict[str, object]:
    """Certify the trace minimum using only ``min(Lambda_24)=4``.

    The decomposition is

    ``||A-B||^2 + ||A-C||^2 + 2||B-C||^2``
    ``+ ||A||^2 + 2||B||^2 + 2||C||^2``.

    Exhausting which of ``A,B,C`` are equal or zero shows that every nonzero
    triple has at least three units of nonzero norm contribution.  Equality
    occurs only for ``B=C=0`` and a minimal ``A``.
    """

    minimum_weight = None
    equality_patterns: set[tuple[bool, bool, bool, bool, bool, bool]] = set()
    # Equality/non-equality of three abstract vectors can be represented by
    # small integer labels; this exhaustive finite model covers every zero and
    # coincidence pattern relevant to vanishing of the six linear forms.
    for a, b, c in product(range(4), repeat=3):
        if a == b == c == 0:
            continue
        values = (a, b, c)
        weight = sum(
            coefficient
            for coefficient, form in TRACE_DECOMPOSITION
            if sum(form[i] * values[i] for i in range(3)) != 0
        )
        if minimum_weight is None or weight < minimum_weight:
            minimum_weight = weight
            equality_patterns.clear()
        if weight == minimum_weight:
            equality_patterns.add(
                (a == b, a == c, b == c, a == 0, b == 0, c == 0)
            )
    expected_pattern = (False, False, True, False, True, True)
    if minimum_weight != 3 or equality_patterns != {expected_pattern}:
        raise AssertionError("the trace-minimum incidence audit drifted")
    return {
        "minimum_nonzero_term_weight": minimum_weight,
        "leech_minimum": LEECH_MINIMUM,
        "trace_minimum": minimum_weight * LEECH_MINIMUM,
        "equality_case": "B=C=0 and A is a minimal Leech vector",
    }


# One generator polynomial for the binary cyclic [23,12,7] Golay code.
_GOLAY23_GENERATOR = sum(1 << exponent for exponent in (11, 9, 7, 6, 5, 1, 0))


@lru_cache(maxsize=1)
def extended_binary_golay_codewords() -> tuple[int, ...]:
    """Return the standard extended binary Golay code as 24-bit words."""

    words: list[int] = []
    for message in range(1 << 12):
        word = 0
        for bit in range(12):
            if message & (1 << bit):
                word ^= _GOLAY23_GENERATOR << bit
        # Extend by the parity bit, producing the doubly-even [24,12,8] code.
        word |= (word.bit_count() & 1) << 23
        words.append(word)
    if len(set(words)) != 4096:
        raise AssertionError("the Golay generator lost dimension twelve")
    return tuple(words)


def golay_weight_distribution() -> dict[int, int]:
    return dict(sorted(Counter(word.bit_count() for word in extended_binary_golay_codewords()).items()))


def iter_minimal_leech_vectors_scaled() -> Iterator[tuple[str, tuple[int, ...]]]:
    """Yield all minimal Leech vectors in coordinates scaled by ``sqrt(8)``.

    The three standard coordinate types are ``(4^2,0^22)``,
    ``(2^8,0^16)`` on Golay octads with even sign parity, and
    ``(3,1^23)`` with Golay sign changes.  Every yielded integer vector has
    squared Euclidean length 32, hence Leech norm four after division by 8.
    """

    # Type 1: choose two positions and their independent signs.
    for first, second in combinations(range(24), 2):
        for sign_first, sign_second in product((-1, 1), repeat=2):
            vector = [0] * 24
            vector[first] = 4 * sign_first
            vector[second] = 4 * sign_second
            yield "4^2", tuple(vector)

    # Type 2: an octad support and an even number of negative coordinates.
    octads = tuple(word for word in extended_binary_golay_codewords() if word.bit_count() == 8)
    for octad in octads:
        support = tuple(index for index in range(24) if octad & (1 << index))
        for negative_mask in range(1 << 8):
            if negative_mask.bit_count() % 2:
                continue
            vector = [0] * 24
            for local_index, coordinate in enumerate(support):
                vector[coordinate] = -2 if negative_mask & (1 << local_index) else 2
            yield "2^8", tuple(vector)

    # Type 3: choose the coordinate of absolute value three and apply every
    # extended-Golay sign change.  The all-ones word includes global negation.
    for distinguished in range(24):
        for codeword in extended_binary_golay_codewords():
            vector = [-3 if index == distinguished else 1 for index in range(24)]
            for index in range(24):
                if codeword & (1 << index):
                    vector[index] = -vector[index]
            yield "3,1^23", tuple(vector)


@lru_cache(maxsize=1)
def minimal_leech_shell_audit() -> dict[str, object]:
    """Stream and hash the complete standard minimal shell."""

    counts: Counter[str] = Counter()
    digest = sha256()
    total = 0
    for kind, vector in iter_minimal_leech_vectors_scaled():
        if sum(coordinate * coordinate for coordinate in vector) != 32:
            raise AssertionError("a standard Leech minimal vector has wrong norm")
        counts[kind] += 1
        total += 1
        digest.update(bytes(coordinate + 4 for coordinate in vector))
    if total != LEECH_MINIMAL_VECTOR_COUNT:
        raise AssertionError("the standard Leech minimal shell count drifted")
    if CO0_ORDER // CO2_ORDER != total or CO0_ORDER % CO2_ORDER:
        raise AssertionError("the Co0/Co2 orbit-stabilizer count drifted")
    return {
        "coordinate_scale": "sqrt(8)",
        "scaled_squared_length": 32,
        "leech_norm": 4,
        "type_counts": dict(sorted(counts.items())),
        "total_count": total,
        "sha256_signed_coordinate_stream": digest.hexdigest(),
        "co0_orbit_count": 1,
        "co0_order": CO0_ORDER,
        "point_stabilizer": "Co2",
        "point_stabilizer_order": CO2_ORDER,
    }


@dataclass(frozen=True, slots=True)
class LeechControlLevel:
    norm: Cubic7
    control: Cubic7
    height: Cubic7
    target: Cubic7


@lru_cache(maxsize=None)
def leech_control_levels(height_bound: int) -> tuple[LeechControlLevel, ...]:
    """Exhaust all norm-two positive controls through exact height ``H``."""

    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")
    norm = Cubic7(2)
    a = TIMELIKE_COEFFICIENT
    caps = (norm * height_bound / a, norm / (-a), norm / (-a))
    levels: set[LeechControlLevel] = set()
    for control in _bounded_for_caps(caps):
        if control.sign(0) <= 0:
            continue
        if not (2 * a * control / norm).is_integral:
            continue
        target = norm + a * control * control
        if any(target.sign(place) < 0 for place in (1, 2)):
            continue
        height = a * control * control / norm
        if (Cubic7(height_bound) - height).sign(0) < 0:
            continue
        levels.add(LeechControlLevel(norm, control, height, target))
    return tuple(sorted(levels, key=lambda level: (level.height, level.control)))


def first_shell_audit() -> dict[str, object]:
    """Certify height zero and the first positive control packet."""

    trace_audit = trace_minimum_and_equality_case()
    if trace_audit["trace_minimum"] != 12:
        raise AssertionError("the Leech trace minimum drifted")
    levels = leech_control_levels(10)
    if len(levels) != 1:
        raise AssertionError("the first exact control window changed")
    level = levels[0]
    expected_control = Cubic7(-1, 1, 1)
    expected_height = Cubic7(0, 2, 1)
    expected_target = Cubic7(2, 4, 2)
    if (level.control, level.height, level.target) != (
        expected_control,
        expected_height,
        expected_target,
    ):
        raise AssertionError("the first Leech control packet changed")
    if level.target.trace() != 12:
        raise AssertionError("the first target no longer lies on the trace minimum")

    # Every trace-minimal spatial vector is constant and has field square 4,
    # whereas this packet requires 2+4*theta+2*theta^2.
    if level.target == Cubic7(4):
        raise AssertionError("the empty first packet unexpectedly acquired minimal vectors")
    return {
        "canonical_mirror_norms": [norm.to_json() for norm in CANONICAL_MIRROR_NORMS],
        "complete_mirror_norm_proof": COMPLETE_MIRROR_NORM_PROOF,
        "height_zero": {
            "target_norm": Cubic7(2).to_json(),
            "target_trace": str(Cubic7(2).trace()),
            "minimum_nonzero_spatial_trace": trace_audit["trace_minimum"],
            "root_count": 0,
        },
        "first_positive_control_window": {
            "exhaustive_height_bound": 10,
            "level_count": 1,
            "control": level.control.to_json(),
            "height": level.height.to_json(),
            "target": level.target.to_json(),
            "target_trace": str(level.target.trace()),
            "trace_minimal_spatial_vectors_examined": LEECH_MINIMAL_VECTOR_COUNT,
            "trace_minimal_spatial_field_square": Cubic7(4).to_json(),
            "candidate_root_count": 0,
        },
        "trace_lattice": {
            "power_basis_gram": TRACE_GRAM,
            "integral_change_of_basis": TRACE_CHANGE_OF_BASIS,
            "reduced_gram": REDUCED_TRACE_GRAM,
            "norm_decomposition": TRACE_DECOMPOSITION,
            **trace_audit,
        },
        "minimal_leech_shell": minimal_leech_shell_audit(),
        "continuation_blocker": (
            "The next targets require Co0 orbits of ordered triples (A,B,C) "
            "with prescribed six Leech norms and pairings. Individual Leech "
            "theta-series coefficients or one-vector Co0 shell orbits do not "
            "determine these simultaneous triple orbits."
        ),
    }


__all__ = [
    "CANONICAL_MIRROR_NORMS",
    "COMPLETE_MIRROR_NORM_PROOF",
    "LeechControlLevel",
    "extended_binary_golay_codewords",
    "first_shell_audit",
    "golay_weight_distribution",
    "iter_minimal_leech_vectors_scaled",
    "leech_control_levels",
    "minimal_leech_shell_audit",
    "reduced_trace_gram_from_change_of_basis",
    "trace_gram_from_norm_decomposition",
    "trace_minimum_and_equality_case",
]
