"""Exact first Vinberg packet for a rootless silver ``S26`` host.

Let ``S26`` be Conway--Borcherds' odd unimodular rank-26 lattice of
minimum three.  For ``rho=2*sqrt(2)-1`` this module audits

``<-rho> + S26``

over ``Z[sqrt(2)]``.  It has no height-zero roots.  The first nonempty
positive packet is a scalar copy of the norm-four shell of ``S26``.  A
fixed minimal vector supplies a strict hyperideal point for the whole packet.

The conclusion is a bounded-prefix obstruction only: later Vinberg walls may
cut the displayed point.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd, isqrt
from typing import TypeAlias

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps


IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]

S26_SOURCE_URL = "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/neb3.html"
S26_MINIMUM = 3
S26_KISSING_NUMBER = 3_120
S26_GROUP_ORDER = 18_720_000
S26_SHELL_COUNTS = {3: 3_120, 4: 102_180, 5: 1_482_624, 6: 13_191_360}

# The catalogue prints the lower triangle.  This literal tuple is expanded
# symmetrically by ``s26_gram``.
S26_LOWER_TRIANGLE: tuple[tuple[int, ...], ...] = (
    (3,),
    (1,3),
    (1,1,3),
    (1,1,1,3),
    (1,1,1,0,3),
    (1,1,1,1,1,3),
    (1,1,1,1,1,1,3),
    (1,1,1,1,1,1,1,3),
    (1,1,1,1,1,1,0,1,3),
    (1,1,1,1,1,1,1,1,1,3),
    (1,1,1,1,0,0,0,0,1,1,3),
    (1,1,1,0,1,1,0,0,1,1,1,3),
    (1,0,1,0,1,1,1,0,0,0,1,1,3),
    (1,1,1,0,1,1,1,1,0,1,0,1,1,3),
    (1,1,1,0,1,1,1,0,1,1,1,1,1,0,3),
    (1,1,1,1,0,0,1,0,0,0,0,0,0,1,0,3),
    (1,1,1,1,1,1,0,1,1,1,1,0,0,1,0,0,3),
    (1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,1,0,3),
    (1,1,1,1,1,1,1,0,1,1,0,1,0,1,1,1,1,0,3),
    (1,1,1,1,0,1,1,1,1,0,1,1,1,0,1,0,0,0,0,3),
    (1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,1,1,1,1,0,3),
    (0,1,1,1,1,1,1,1,1,1,0,1,0,1,0,0,1,0,1,1,0,3),
    (1,1,1,0,1,0,1,1,0,1,1,1,1,1,1,0,0,1,0,1,0,0,3),
    (1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,1,1,1,3),
    (1,0,1,1,0,0,1,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0,1,0,3),
    (1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,3),
)


def _bareiss_determinant(matrix: IntegerMatrix) -> int:
    work = [list(row) for row in matrix]
    size = len(work)
    previous = 1
    sign = 1
    for column in range(size - 1):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
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


@lru_cache(maxsize=1)
def s26_gram() -> IntegerMatrix:
    if tuple(map(len, S26_LOWER_TRIANGLE)) != tuple(range(1, 27)):
        raise AssertionError("the S26 lower triangle has the wrong row lengths")
    return tuple(
        tuple(
            S26_LOWER_TRIANGLE[max(i, j)][min(i, j)]
            for j in range(26)
        )
        for i in range(26)
    )


def s26_catalogue_audit() -> bool:
    gram = s26_gram()
    return (
        all(gram[i][j] == gram[j][i] for i in range(26) for j in range(26))
        and all(gram[i][i] == 3 for i in range(26))
        and _bareiss_determinant(gram) == 1
    )


RHO = 2 * SQRT2 - SILVER_ONE
EPSILON = SILVER_ONE + SQRT2
DYADIC_GENERATOR = 2 + SQRT2
P7_GENERATOR = 3 + SQRT2


@dataclass(frozen=True, slots=True)
class S26NormClass:
    norm: Silver
    dyadic_exponent: int
    p7: bool


NORM_CLASSES = tuple(
    S26NormClass(
        (SILVER_ONE, DYADIC_GENERATOR, Silver(2))[dyadic_exponent]
        * (P7_GENERATOR if p7 else SILVER_ONE),
        dyadic_exponent,
        p7,
    )
    for dyadic_exponent, p7 in product(range(3), (False, True))
)

COMPLETE_NORM_CLASS_PROOF = (
    "Put pi=sqrt(2), rho=2*pi-1 and delta=3+pi=epsilon*rho. "
    "For a primitive r=(x,y), reflection integrality is q(r)|2I(r), "
    "I(r)=(rho*x,(y,S26)). Since S26 is unimodular, every prime away "
    "from (pi) and (rho) is excluded by primitivity. At (pi), rho is a "
    "unit and I(r) is the unit ideal, so v_pi(q) is 0, 1, or 2. At "
    "p=(rho), the exponent is zero unless y lies in p*S26; in that case "
    "primitive x is a p-unit and q=-rho*x^2 mod p^2, so the exponent is "
    "exactly one. Class number one and U^+=U^2 give the six canonical "
    "totally-positive representatives (1,2+pi,2) times (1 or 3+pi)."
)


@dataclass(frozen=True, slots=True)
class S26GramTriple:
    u_square: int
    v_square: int
    u_dot_v: int

    def combination_square(self, a: int, b: int) -> int:
        return (
            a * a * self.u_square
            + b * b * self.v_square
            + 2 * a * b * self.u_dot_v
        )


def candidate_s26_gram_triples(target: Silver) -> tuple[S26GramTriple, ...]:
    """Necessary exact Gram data using integrality and minimum three.

    Positive semidefiniteness and the absence of norms one and two among
    combinations in ``[-3,3]^2`` define a rigorous finite superset of actual
    S26 pairs.  At the first packet this gate is exact.  At later heights it
    is deliberately only necessary, which is sufficient for upper-bounding
    every possible wall pairing.
    """

    if not target.is_integral or (target.b / 2).denominator != 1:
        return ()
    rational = int(target.a)
    dot = int(target.b / 2)
    output = []
    for v_square in range(rational // 2 + 1):
        u_square = rational - 2 * v_square
        if u_square < 0:
            continue
        if u_square not in (0,) and u_square < S26_MINIMUM:
            continue
        if v_square not in (0,) and v_square < S26_MINIMUM:
            continue
        if dot * dot > u_square * v_square:
            continue
        triple = S26GramTriple(u_square, v_square, dot)
        if any(
            triple.combination_square(a, b) in (1, 2)
            for a in range(-3, 4)
            for b in range(-3, 4)
            if (a, b) != (0, 0)
        ):
            continue
        output.append(triple)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class S26Control:
    norm_class: S26NormClass
    height: Silver
    control: Silver
    target: Silver


FIRST_PACKET_HEIGHT = 5 + 4 * SQRT2


@lru_cache(maxsize=None)
def enumerate_controls(height_bound: Silver = FIRST_PACKET_HEIGHT) -> tuple[S26Control, ...]:
    output = set()
    for norm_class in NORM_CLASSES:
        q = norm_class.norm
        for x in _bounded_for_caps((q * height_bound / RHO, q / (-RHO))):
            if x.sign(0) <= 0 or not (2 * RHO * x / q).is_integral:
                continue
            target = q + RHO * x * x
            height = RHO * x * x / q
            if target.sign(1) < 0 or height > height_bound:
                continue
            output.add(S26Control(norm_class, height, x, target))
    return tuple(sorted(output, key=lambda row: (row.height, row.norm_class.norm, row.control)))


def reduced_target(control: S26Control) -> Silver:
    return control.target / (RHO * RHO) if control.norm_class.p7 else control.target


def gated_gram_triples(control: S26Control) -> tuple[S26GramTriple, ...]:
    return candidate_s26_gram_triples(reduced_target(control))


def height_zero_gated_gram_triples(norm_class: S26NormClass) -> tuple[S26GramTriple, ...]:
    target = norm_class.norm / (RHO * RHO) if norm_class.p7 else norm_class.norm
    return candidate_s26_gram_triples(target)


@dataclass(frozen=True, slots=True)
class S26FirstPacketAudit:
    control_count: int
    distinct_height_count: int
    earlier_control_count: int
    earlier_surviving_control_count: int
    first_norm: Silver
    first_control: Silver
    first_target: Silver
    surviving_triple: S26GramTriple
    scalar_shell_norm: int
    wall_count: int
    slab_bound: Silver
    witness_scale: Silver
    witness_square: Silver
    witness_lorentz_norm: Silver
    pairing_bound: int
    wall_margin: Silver
    strict_hyperideal: bool
    finite_volume: bool
    compact: bool


def first_packet_audit(*, shell_four_count: int) -> S26FirstPacketAudit:
    controls = enumerate_controls()
    earlier = tuple(row for row in controls if row.height < FIRST_PACKET_HEIGHT)
    first_rows = tuple(row for row in controls if row.height == FIRST_PACKET_HEIGHT)
    if any(gated_gram_triples(row) for row in earlier):
        raise AssertionError("an earlier S26 necessary packet survived")
    surviving = tuple((row, gated_gram_triples(row)) for row in first_rows if gated_gram_triples(row))
    expected_triple = S26GramTriple(4, 4, 4)
    if len(surviving) != 1 or surviving[0][1] != (expected_triple,):
        raise AssertionError("the S26 first packet Gram census drifted")
    row = surviving[0][0]
    if row.norm_class.norm != 2 or row.control != 2 + SQRT2:
        raise AssertionError("the S26 first packet control drifted")

    # y=epsilon*d, d^2=4.  Dividing the wall inequality by epsilon gives
    # z.d <= 4-pi.  For fixed e^2=3, integrality gives e.d in Z.  The
    # minimum-three condition applied to e+-d improves Cauchy to |e.d|<=2:
    # a value 3 would make one of e+-d have norm one.
    slab_bound = 4 - SQRT2
    witness_scale = Silver(1)
    witness_square = 3 * witness_scale * witness_scale
    witness_lorentz_norm = witness_square - RHO
    wall_margin = slab_bound - 2 * witness_scale
    return S26FirstPacketAudit(
        control_count=len(controls),
        distinct_height_count=len({row.height for row in controls}),
        earlier_control_count=len(earlier),
        earlier_surviving_control_count=0,
        first_norm=row.norm_class.norm,
        first_control=row.control,
        first_target=row.target,
        surviving_triple=expected_triple,
        scalar_shell_norm=4,
        wall_count=shell_four_count,
        slab_bound=slab_bound,
        witness_scale=witness_scale,
        witness_square=witness_square,
        witness_lorentz_norm=witness_lorentz_norm,
        pairing_bound=2,
        wall_margin=wall_margin,
        strict_hyperideal=(
            witness_lorentz_norm.sign(0) > 0 and wall_margin.sign(0) > 0
        ),
        finite_volume=False,
        compact=False,
    )


def first_packet_root_is_reflective(d: tuple[int, ...]) -> bool:
    if len(d) != 26:
        raise ValueError("an S26 coordinate vector must have length 26")
    gram = s26_gram()
    square = sum(d[i] * gram[i][j] * d[j] for i in range(26) for j in range(26))
    if square != 4:
        return False
    q = Silver(2)
    x = 2 + SQRT2
    y_scale = EPSILON
    spatial_pairings = tuple(
        y_scale * sum(gram[i][j] * d[j] for j in range(26))
        for i in range(26)
    )
    coordinates = tuple([x, *(y_scale * value for value in d)])
    primitive = gcd(*(int(value.a) for value in coordinates), *(int(value.b) for value in coordinates)) == 1
    return (
        -RHO * x * x + y_scale * y_scale * square == q
        and (2 * RHO * x / q).is_integral
        and all((2 * value / q).is_integral for value in spatial_pairings)
        and primitive
    )


H1000_BOUND = Silver(1000)
H1000_WITNESS_SCALE = Silver(Fraction(43, 55))


@dataclass(frozen=True, slots=True)
class S26AugmentedPairingAudit:
    maximum_pairing: Silver
    maximizing_u_pairing: int
    maximizing_v_pairing: int
    retained_extension_count: int


def _augmented_combination_square(
    triple: S26GramTriple,
    p: int,
    q: int,
    i: int,
    j: int,
    k: int,
) -> int:
    return (
        3 * i * i
        + triple.u_square * j * j
        + triple.v_square * k * k
        + 2 * p * i * j
        + 2 * q * i * k
        + 2 * triple.u_dot_v * j * k
    )


@lru_cache(maxsize=None)
def augmented_pairing_audit(triple: S26GramTriple) -> S26AugmentedPairingAudit:
    """Maximise a fixed minimal vector's pairing with ``u+sqrt(2)v``.

    The enumeration is over a rigorous superset of actual S26 triples.  We
    impose exact positive semidefiniteness of the augmented 3-by-3 Gram and
    exclude norms one and two among all coefficient combinations in
    ``[-3,3]^3``.  This bounded short-combination test is only a necessary
    local gate, not a proof that a surviving abstract ternary Gram embeds in
    S26.  Maximising over the superset is therefore safe for the witness.

    The scalar rank-one case is handled separately, since the determinant is
    identically zero and actual pairings necessarily satisfy ``p=q``.
    """

    u_square = triple.u_square
    v_square = triple.v_square
    dot = triple.u_dot_v
    if u_square == v_square == dot:
        best_p: int | None = None
        count = 0
        for p in range(-isqrt(3 * u_square), isqrt(3 * u_square) + 1):
            if p * p > 3 * u_square:
                continue
            if any(
                3 * i * i + u_square * j * j + 2 * p * i * j in (1, 2)
                for i in range(-5, 6)
                for j in range(-5, 6)
                if (i, j) != (0, 0)
            ):
                continue
            count += 1
            if best_p is None or p > best_p:
                best_p = p
        if best_p is None:
            raise AssertionError("no scalar augmented pairing survived")
        return S26AugmentedPairingAudit(
            Silver(best_p, best_p), best_p, best_p, count
        )

    best: Silver | None = None
    best_pair: tuple[int, int] | None = None
    count = 0
    p_bound = isqrt(3 * u_square)
    q_bound = isqrt(3 * v_square)
    for p in range(-p_bound, p_bound + 1):
        if p * p > 3 * u_square:
            continue
        for q in range(-q_bound, q_bound + 1):
            if q * q > 3 * v_square:
                continue
            determinant = (
                3 * (u_square * v_square - dot * dot)
                - v_square * p * p
                - u_square * q * q
                + 2 * dot * p * q
            )
            if determinant < 0:
                continue
            if any(
                _augmented_combination_square(triple, p, q, i, j, k) in (1, 2)
                for i in range(-3, 4)
                for j in range(-3, 4)
                for k in range(-3, 4)
                if (i, j, k) != (0, 0, 0)
            ):
                continue
            count += 1
            value = Silver(p, q)
            if best is None or value > best:
                best = value
                best_pair = (p, q)
    if best is None or best_pair is None:
        raise AssertionError("no augmented S26 pairing survived")
    return S26AugmentedPairingAudit(best, *best_pair, count)


@dataclass(frozen=True, slots=True)
class S26H1000Row:
    control: S26Control
    triple: S26GramTriple
    pairing: S26AugmentedPairingAudit
    normalized_wall_bound: Silver
    witness_margin: Silver


@dataclass(frozen=True, slots=True)
class S26H1000Audit:
    complete_height_bound: Silver
    control_count: int
    distinct_height_count: int
    nonempty_control_count: int
    gated_binary_type_count: int
    nonempty_height_count: int
    rows: tuple[S26H1000Row, ...]
    worst_row: S26H1000Row
    witness_scale: Silver
    witness_square: Silver
    witness_lorentz_norm: Silver
    all_margins_strict: bool
    strict_hyperideal: bool
    finite_volume: bool
    compact: bool


@lru_cache(maxsize=1)
def h1000_minimal_direction_audit() -> S26H1000Audit:
    """Certify a fixed minimal direction through every possible H<=1000 wall."""

    return minimal_direction_audit_through(H1000_BOUND)


@lru_cache(maxsize=None)
def minimal_direction_audit_through(height_bound: Silver) -> S26H1000Audit:
    """Certify the fixed norm-three direction through an arbitrary exact bound.

    The return type retains its historical name for API compatibility.  All
    counts and margins are recomputed from the supplied algebraic height.
    """

    controls = enumerate_controls(height_bound)
    rows: list[S26H1000Row] = []
    for control in controls:
        triples = gated_gram_triples(control)
        for triple in triples:
            pairing = augmented_pairing_audit(triple)
            divisor = RHO if control.norm_class.p7 else SILVER_ONE
            bound = RHO * control.control / divisor
            margin = bound - H1000_WITNESS_SCALE * pairing.maximum_pairing
            rows.append(S26H1000Row(control, triple, pairing, bound, margin))
    if not rows:
        raise AssertionError("the S26 H1000 necessary census is empty")
    worst = min(rows, key=lambda row: row.witness_margin)
    witness_square = 3 * H1000_WITNESS_SCALE * H1000_WITNESS_SCALE
    witness_lorentz_norm = witness_square - RHO
    return S26H1000Audit(
        complete_height_bound=height_bound,
        control_count=len(controls),
        distinct_height_count=len({row.height for row in controls}),
        nonempty_control_count=sum(bool(gated_gram_triples(row)) for row in controls),
        gated_binary_type_count=len(rows),
        nonempty_height_count=len({row.control.height for row in rows}),
        rows=tuple(rows),
        worst_row=worst,
        witness_scale=H1000_WITNESS_SCALE,
        witness_square=witness_square,
        witness_lorentz_norm=witness_lorentz_norm,
        all_margins_strict=all(row.witness_margin.sign(0) > 0 for row in rows),
        strict_hyperideal=(
            witness_lorentz_norm.sign(0) > 0
            and all(row.witness_margin.sign(0) > 0 for row in rows)
        ),
        finite_volume=False,
        compact=False,
    )


__all__ = [
    "COMPLETE_NORM_CLASS_PROOF",
    "FIRST_PACKET_HEIGHT",
    "NORM_CLASSES",
    "RHO",
    "S26_GROUP_ORDER",
    "S26_KISSING_NUMBER",
    "S26_LOWER_TRIANGLE",
    "S26_MINIMUM",
    "S26_SHELL_COUNTS",
    "S26_SOURCE_URL",
    "candidate_s26_gram_triples",
    "enumerate_controls",
    "first_packet_audit",
    "first_packet_root_is_reflective",
    "gated_gram_triples",
    "h1000_minimal_direction_audit",
    "minimal_direction_audit_through",
    "height_zero_gated_gram_triples",
    "s26_catalogue_audit",
    "s26_gram",
]
