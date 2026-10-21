"""Exact finite front end for the golden AGCat rank-23--30 decision queue.

The global AGCat theorem proves that the queue is finite.  It does *not*
enumerate the queue or prove that it is empty.  This module freezes the first
universal finite set needed by an actual decision computation over
``K=Q(phi)``: the complete length-two Cartan-invariant box for a connected
Vinberg frame with bound ``A=62``.

If ``b_ij=a_ij*a_ji``, then at the distinguished place

    0 <= b_ij < 62^2,

and at the other real place positive definiteness and Cauchy--Schwarz give

    0 <= b_ij' < 4.

The strict conjugate inequality uses independence of two frame roots:
Cauchy equality in the positive-definite conjugate form would make them
proportional.  (The weak closed Cauchy box has one extra value, namely
``4``.)  The enumeration below is exact in ``Z[phi]``.  It then removes the
two subparallel values which are not Coxeter angle values in
``Q(sqrt(5))``.  Values greater than four at the distinguished place are the
dotted-edge possibilities and remain in the alphabet.

This is deliberately not presented as a chamber census: connected graph
enumeration, longer-cycle compatibility, ideal gauges, reflection orders,
conductor lattices, and complete Vinberg transcripts remain later finite
stages of the queue.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations_with_replacement
from math import isqrt

from coxeter_polytopes.exact.golden import Golden, PHI


FRAME_BOUND = 62
MIN_HYPERBOLIC_DIMENSION = 22
MAX_HYPERBOLIC_DIMENSION = 29


def _integral_coefficients(value: Golden) -> tuple[int, int]:
    if not value.is_integral:
        raise ValueError("the exact floor helper requires a golden integer")
    return int(value.a), int(value.b)


def floor_integral_golden(value: Golden, embedding: int = 0) -> int:
    """Return ``floor(sigma_embedding(value))`` without approximation."""

    if embedding not in (0, 1):
        raise ValueError("embedding must be zero or one")
    if embedding == 1:
        value = value.conjugate()
    a, b = _integral_coefficients(value)
    if b == 0:
        return a
    magnitude = abs(b)
    # sqrt(5)*magnitude is irrational for magnitude != 0, so this is its
    # exact floor.  Since phi=(1+sqrt(5))/2, parity causes no ambiguity.
    floor_positive = (magnitude + isqrt(5 * magnitude * magnitude)) // 2
    if b > 0:
        return a + floor_positive
    return a - floor_positive - 1


def ceil_integral_golden(value: Golden, embedding: int = 0) -> int:
    """Return ``ceil(sigma_embedding(value))`` without approximation."""

    floor = floor_integral_golden(value, embedding)
    conjugated = value if embedding == 0 else value.conjugate()
    _a, b = _integral_coefficients(conjugated)
    return floor if b == 0 else floor + 1


@lru_cache(maxsize=2)
def _golden_nonnegative_edge_box(*, strict_conjugate: bool) -> tuple[Golden, ...]:
    """Enumerate a nonnegative two-place length-two box."""

    selected_upper = FRAME_BOUND**2
    conjugate_upper = 2**2

    # If x=a+b*phi is nonnegative at both places, then
    # b*sqrt(5)=x-x'.  For b>0, b*sqrt(5)<selected_upper; for b<0,
    # (-b)*sqrt(5)<=conjugate_upper.  sqrt(5)>2 gives this exhaustive
    # integer interval.
    b_min = -((conjugate_upper - 1) // 2)
    b_max = (selected_upper - 1) // 2
    output: list[Golden] = []
    for b in range(b_min, b_max + 1):
        b_phi = b * PHI
        lower = max(
            ceil_integral_golden(-b_phi, 0),
            ceil_integral_golden(-b_phi, 1),
        )
        upper = min(
            # a+b*phi < selected_upper
            ceil_integral_golden(Golden(selected_upper) - b_phi, 0) - 1,
            # a+b*phi' < (or <=) conjugate_upper
            (
                ceil_integral_golden(Golden(conjugate_upper) - b_phi, 1) - 1
                if strict_conjugate
                else floor_integral_golden(Golden(conjugate_upper) - b_phi, 1)
            ),
        )
        for a in range(lower, upper + 1):
            value = Golden(a, b)
            if not (
                value.sign(0) >= 0
                and (value - selected_upper).sign(0) < 0
                and value.sign(1) >= 0
                and (value - conjugate_upper).sign(1) < (0 if strict_conjugate else 1)
            ):
                raise AssertionError("coefficient interval emitted an out-of-box value")
            output.append(value)
    if len(set(output)) != len(output):
        raise AssertionError("the golden edge box contains a duplicate")
    return tuple(output)


@lru_cache(maxsize=1)
def golden_closed_cauchy_edge_box() -> tuple[Golden, ...]:
    """The weak Cauchy box before frame-root independence is imposed."""

    return _golden_nonnegative_edge_box(strict_conjugate=False)


@lru_cache(maxsize=1)
def golden_frame_edge_box() -> tuple[Golden, ...]:
    """The exact two-place box for a pair of independent frame roots.

    At the positive-definite nonselected place, equality in Cauchy--Schwarz
    would make the two embedded roots proportional.  A k-independent frame
    pair stays independent after embedding, so its product is strictly less
    than four there.  In the golden closed box, equality occurs only for the
    rational element four, which is therefore absent here.  Compactness also
    independently excludes the corresponding parallel selected pair.
    """

    output = _golden_nonnegative_edge_box(strict_conjugate=True)
    closed = golden_closed_cauchy_edge_box()
    if set(closed) - set(output) != {Golden(4)}:
        raise AssertionError("the Cauchy equality boundary changed")
    return output


@dataclass(frozen=True, slots=True)
class GoldenCoxeterPairValue:
    label: str
    product: Golden


@lru_cache(maxsize=1)
def golden_finite_pair_values() -> tuple[GoldenCoxeterPairValue, ...]:
    """All finite Coxeter angle products lying in ``Q(sqrt(5))``.

    For finite label m, ``b_m=4*cos(pi/m)^2`` and
    ``Q(b_m)=Q(cos(2*pi/m))``.  Degree at most two forces
    ``m in {2,3,4,5,6,8,10,12}``; labels 8 and 12 generate the other
    quadratic fields ``Q(sqrt(2))`` and ``Q(sqrt(3))``.
    """

    return (
        GoldenCoxeterPairValue("2", Golden(0)),
        GoldenCoxeterPairValue("3", Golden(1)),
        GoldenCoxeterPairValue("4", Golden(2)),
        GoldenCoxeterPairValue("5", 1 + PHI),
        GoldenCoxeterPairValue("6", Golden(3)),
        GoldenCoxeterPairValue("10", 2 + PHI),
    )


@lru_cache(maxsize=1)
def golden_coxeter_frame_pair_alphabet() -> tuple[Golden, ...]:
    """Exact possible pair products after the Coxeter-angle filter."""

    finite = {row.product for row in golden_finite_pair_values()}
    return tuple(
        value
        for value in golden_frame_edge_box()
        if value in finite or (value - 4).sign(0) > 0
    )


def golden_integral_square_roots(value: Golden) -> tuple[Golden, ...]:
    """Enumerate every square root of a totally nonnegative golden integer."""

    if not value.is_integral or any(value.sign(place) < 0 for place in (0, 1)):
        return ()
    embedding_bounds = []
    for place in (0, 1):
        upper = 1
        while (Golden(upper * upper) - value).sign(place) <= 0:
            upper *= 2
        embedding_bounds.append(upper)
    # For z=a+b*phi, b*sqrt(5)=z-z'.  sqrt(5)>2 and phi<2.
    b_bound = sum(embedding_bounds) // 2 + 1
    a_bound = embedding_bounds[0] + 2 * b_bound
    roots = tuple(
        sorted(
            Golden(a, b)
            for a in range(-a_bound, a_bound + 1)
            for b in range(-b_bound, b_bound + 1)
            if Golden(a, b) ** 2 == value
        )
    )
    if len(roots) not in (0, 1, 2) or (len(roots) == 2 and roots[0] != -roots[-1]):
        raise AssertionError("a quadratic field element acquired unexpected square roots")
    return roots


@dataclass(frozen=True, slots=True)
class GoldenFiniteTrianglePacket:
    labels: tuple[int, int, int]
    edge_products: tuple[Golden, Golden, Golden]
    directed_three_cycle: Golden
    gram_determinant: Golden
    selected_determinant_sign: int
    conjugate_determinant_sign: int
    admissible_lorentzian: bool


@lru_cache(maxsize=1)
def golden_finite_angle_triangle_packets() -> tuple[GoldenFiniteTrianglePacket, ...]:
    """Exhaust finite-angle complete triangle cycle packets over Q(sqrt(5)).

    The input consists of the five nonorthogonal finite Coxeter products.
    Modulo the three vertex permutations there are 35 multisets.  A packet
    exists over the golden field precisely when the edge-product product has
    an integral square root; the acute three-cycle is the negative root at
    the distinguished place.  The determinant identity is

        det(G)=8-2(c12+c23+c31)+2*b123.

    Since frame normals are independent, admissibility requires a negative
    selected determinant and a positive conjugate determinant.
    """

    finite = {
        int(row.label): row.product
        for row in golden_finite_pair_values()
        if row.product
    }
    output = []
    for labels in combinations_with_replacement(sorted(finite), 3):
        edges = tuple(finite[label] for label in labels)
        product = edges[0] * edges[1] * edges[2]
        roots = golden_integral_square_roots(product)
        positive = next((root for root in roots if root.sign(0) > 0), None)
        if positive is None:
            continue
        cycle = -positive
        determinant = Golden(8) - 2 * sum(edges, start=Golden(0)) + 2 * cycle
        selected = determinant.sign(0)
        conjugate = determinant.sign(1)
        output.append(
            GoldenFiniteTrianglePacket(
                labels,
                edges,
                cycle,
                determinant,
                selected,
                conjugate,
                selected < 0 and conjugate > 0,
            )
        )
    return tuple(output)


def cycle_box_record(length: int) -> dict[str, object]:
    """A finite coefficient rectangle containing the AGCat cycle box."""

    if length < 2:
        raise ValueError("cycle length must be at least two")
    selected = FRAME_BOUND**length
    conjugate = 2**length
    # For x=a+b*phi with |x|<selected and |x'|<=conjugate,
    # |b|*sqrt(5)=|x-x'|<selected+conjugate.  Since sqrt(5)>2
    # and the sum is even, the displayed integral coefficient bounds follow.
    b_bound = (selected + conjugate) // 2 - 1
    a_bound = conjugate + b_bound
    return {
        "cycle_length": length,
        "set_definition": (
            "x in Z[phi], |sigma0(x)|<62^r, |sigma1(x)|<=2^r"
        ),
        "selected_absolute_bound_strict": selected,
        "conjugate_absolute_bound_closed": conjugate,
        "coefficient_rectangle": {
            "absolute_b_at_most": b_bound,
            "absolute_a_at_most": a_bound,
        },
        "rectangle_is_a_proved_superset": True,
        "box_materialized": length == 2,
        "length_two_materialization_kind": (
            "stronger nonnegative independent-edge specialization, with strict conjugate Cauchy bound"
            if length == 2
            else None
        ),
    }


def rank_manifest_rows() -> tuple[dict[str, object], ...]:
    pair_count = len(golden_coxeter_frame_pair_alphabet())
    rows = []
    for dimension in range(MIN_HYPERBOLIC_DIMENSION, MAX_HYPERBOLIC_DIMENSION + 1):
        rank = dimension + 1
        pair_slots = rank * (rank - 1) // 2
        rows.append(
            {
                "hyperbolic_dimension": dimension,
                "native_rank": rank,
                "connected_frame_wall_count": rank,
                "unordered_pair_slots": pair_slots,
                "exact_pair_alphabet_count": pair_count,
                "raw_labelled_pair_assignment_upper_bound": {
                    "base": pair_count,
                    "exponent": pair_slots,
                },
                "required_cycle_lengths": list(range(2, rank + 1)),
                "required_cycle_box_count": rank - 1,
                "pipeline_status": {
                    "bounded_connected_frame": "theorem_proved_A_equals_62",
                    "length_two_edge_box": "enumerated_exactly",
                    "connected_graph_and_complete_cycle_packets": "not_enumerated",
                    "golden_ideal_gauge": "not_implemented",
                    "reflection_order_and_conductor_lattices": "not_enumerated",
                    "complete_vinberg_host_transcripts": "not_run",
                    "empty": "not_proved",
                },
            }
        )
    return tuple(rows)


def carrier_logic_records() -> tuple[dict[str, object], ...]:
    """Freeze the two AG-carrier implications which must not be conflated."""

    return (
        {
            "carrier": "real_framed_projective_quadric_hyperplane_arrangement",
            "input_scope": "every finite_volume arithmetic Coxeter chamber in the Vinberg model",
            "exact_return": True,
            "essentially_surjective_onto_its_declared_image": True,
            "construction_is_tautological_arrangement_encoding": True,
            "independent_AG_dimension_bound": False,
            "what_it_proves_here": (
                "after the external n<=995 and A=62 geometry plus cycle/order/conductor "
                "arithmetic, only finitely many reflection-maximal hosts exist"
            ),
            "what_it_does_not_prove": (
                "that the Q(sqrt(5)), dimensions 22--29 finite set is empty"
            ),
        },
        {
            "carrier": "coherent_full_H2_complex_K3_golden_carrier",
            "input_scope": "only Coxeter data for which a coherent same-P K3 package exists",
            "exact_return": "proved only on scoped/constructed branches",
            "essential_surjectivity_for_Qsqrt5_dimensions_22_29": False,
            "false_means": "not_proved; no repository theorem supplies this universal arrow",
            "conditional_quadratic_dimension_bound": 8,
            "bound_reason": "2*(n+3)<=b2(K3)=22",
            "can_be_used_to_exclude_raw_dimensions_22_29": False,
            "why_not": (
                "the K3 carrier is impossible in those ranks, but there is no theorem that "
                "every raw golden chamber must admit that carrier"
            ),
        },
    )


__all__ = [
    "FRAME_BOUND",
    "MAX_HYPERBOLIC_DIMENSION",
    "MIN_HYPERBOLIC_DIMENSION",
    "GoldenCoxeterPairValue",
    "GoldenFiniteTrianglePacket",
    "ceil_integral_golden",
    "carrier_logic_records",
    "cycle_box_record",
    "floor_integral_golden",
    "golden_coxeter_frame_pair_alphabet",
    "golden_closed_cauchy_edge_box",
    "golden_finite_pair_values",
    "golden_finite_angle_triangle_packets",
    "golden_frame_edge_box",
    "golden_integral_square_roots",
    "rank_manifest_rows",
]
