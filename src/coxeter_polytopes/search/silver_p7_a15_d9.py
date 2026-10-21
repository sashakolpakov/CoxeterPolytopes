"""Prime-seven silver Vinberg prefix for ``N(A15 D9)``.

The timelike coefficient is ``a=2*sqrt(2)-1``.  With
``epsilon=1+sqrt(2)`` and ``delta=3+sqrt(2)``, one has
``epsilon*a=delta`` and ``N(delta)=7``.  Even unimodularity gives the
complete primitive norm classes ``1,2,delta,2*delta``.  A delta-supported
mirror has spatial part in ``delta*N(A15 D9)``; enumeration is performed
after dividing by this forced factor and then reconstructed in the full
cyclic index-eight glue lattice.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import groupby

from coxeter_polytopes.exact.silver import SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.silver_a15_d9 import (
    D9_COMPONENT,
    DominantWeight,
    SilverA15D9Candidate,
    SilverA15D9Witness,
    a15_d9_glue_membership,
    a15_residue,
    d9_residue,
    enumerate_dominant_a15_weights,
    enumerate_dominant_d9_weights,
)
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.niemeier_host_triage import cartan_gram
from coxeter_polytopes.search.silver_p7_d24_plus import (
    P7_DELTA,
    P7_TIMELIKE,
    SILVER_P7_D24_SEED,
    SilverP7D24ControlLevel,
    enumerate_p7_d24_control_levels,
    p7_factorization_audit,
)


P7_NORM_CLASSES = SILVER_P7_D24_SEED.canonical_mirror_norms
FIRST_PACKET_HEIGHT = Silver(5, 4)
SECOND_PACKET_HEIGHT = Silver(31, 22)

COMPLETE_P7_A15_D9_NORM_PROOF = (
    "Put pi=sqrt(2), epsilon=1+pi, delta=3+pi, and a=2*pi-1, so "
    "epsilon*a=delta and (delta) is a prime of norm seven. For primitive "
    "r=(x,y), reflection integrality says (q(r)) divides 2I(r), where "
    "I(r)=(a*x,(y,N(A15 D9))). Even unimodularity excludes every prime "
    "away from (pi) and (delta). At (delta), either y is primitive and the "
    "valuation is zero, or reflection integrality gives y in "
    "delta*N(A15 D9); then x is a unit and q=-a*x^2 mod delta^2 has "
    "valuation one. At (pi), if x is a unit then q=-a*x^2 mod pi^2 has "
    "valuation zero. Otherwise primitivity makes the spatial pairing ideal "
    "a unit; evenness gives valuation at least two and q|2I gives valuation "
    "exactly two, namely (2). Since Z[sqrt(2)] has class number one and "
    "every totally-positive unit is a square, the complete canonical norms "
    "are 1, 2, delta, and 2*delta."
)


def _has_delta(norm: Silver) -> bool:
    return norm in (P7_DELTA, 2 * P7_DELTA)


def _scaled_weight(
    weight: DominantWeight, scale: Silver, residue_function
) -> DominantWeight:
    weights = tuple(scale * value for value in weight.weights)
    return DominantWeight(
        scale * scale * weight.square,
        weights,
        tuple(scale * value for value in weight.coordinates),
        residue_function(weights),
    )


def p7_a15_d9_pair(left: SilverA15D9Candidate, right: SilverA15D9Candidate) -> Silver:
    return (
        -P7_TIMELIKE * left.control * right.control
        + sum(
            (
                left.a15.weights[index] * right.a15.coordinates[index]
                for index in range(15)
            ),
            start=SILVER_ZERO,
        )
        + sum(
            (
                left.d9.weights[index] * right.d9.coordinates[index]
                for index in range(9)
            ),
            start=SILVER_ZERO,
        )
    )


def _same_height_acute(roots: tuple[SilverA15D9Candidate, ...]) -> bool:
    """Exact batched acuteness check, using NumPy when it is available.

    All A15 and D9 inverse-Cartan denominators divide 16.  Multiplying each
    Lorentz pairing by 16 therefore gives an integer pair ``A+B*sqrt(2)``.
    The coefficients in the bounded packets here are tiny relative to int64.
    """

    if len(roots) < 2:
        return True
    try:
        import numpy as np
    except ImportError:
        return all(
            p7_a15_d9_pair(left, right).sign(0) <= 0
            for index, left in enumerate(roots)
            for right in roots[index + 1 :]
        )
    rows_a = []
    rows_b = []
    dual_a = []
    dual_b = []
    for root in roots:
        weights = (root.control, *root.a15.weights, *root.d9.weights)
        dual = (
            -P7_TIMELIKE * root.control,
            *root.a15.coordinates,
            *root.d9.coordinates,
        )
        scaled_dual = tuple(16 * value for value in dual)
        if not all(value.is_integral for value in (*weights, *scaled_dual)):
            raise AssertionError("prime-seven batched pairing left denominator 16")
        rows_a.append([int(value.a) for value in weights])
        rows_b.append([int(value.b) for value in weights])
        dual_a.append([int(value.a) for value in scaled_dual])
        dual_b.append([int(value.b) for value in scaled_dual])
    wa = np.asarray(rows_a, dtype=np.int64)
    wb = np.asarray(rows_b, dtype=np.int64)
    da = np.asarray(dual_a, dtype=np.int64)
    db = np.asarray(dual_b, dtype=np.int64)
    entry_count = int(wa.shape[1])
    max_wa = int(np.max(np.abs(wa)))
    max_wb = int(np.max(np.abs(wb)))
    max_da = int(np.max(np.abs(da)))
    max_db = int(np.max(np.abs(db)))
    rational_product_bound = entry_count * (
        max_wa * max_da + 2 * max_wb * max_db
    )
    radical_product_bound = entry_count * (
        max_wa * max_db + max_wb * max_da
    )
    if max(rational_product_bound, radical_product_bound) > np.iinfo(np.int64).max:
        raise OverflowError("prime-seven pairing matrix can exceed int64")
    rational = wa @ da.T + 2 * (wb @ db.T)
    radical = wa @ db.T + wb @ da.T
    if not np.array_equal(rational, rational.T) or not np.array_equal(
        radical, radical.T
    ):
        raise AssertionError("prime-seven batched pairing matrix is not symmetric")
    coefficient_bound = max(
        int(np.max(np.abs(rational))),
        int(np.max(np.abs(radical))),
    )
    if coefficient_bound > 1_000_000_000:
        raise OverflowError("prime-seven exact-sign squares exceed the int64 audit bound")
    for left in range(len(roots)):
        a_values = rational[left, left + 1 :]
        b_values = radical[left, left + 1 :]
        obvious_positive = (a_values > 0) & (b_values >= 0)
        mixed_positive = (a_values > 0) & (b_values < 0) & (
            a_values * a_values > 2 * b_values * b_values
        )
        radical_positive = (a_values <= 0) & (b_values > 0) & (
            a_values * a_values < 2 * b_values * b_values
        )
        if bool(np.any(obvious_positive | mixed_positive | radical_positive)):
            return False
    return True


def _primitive(candidate: SilverA15D9Candidate) -> bool:
    if candidate.norm not in (Silver(2), 2 * P7_DELTA):
        return True
    if not (candidate.control / SQRT2).is_integral:
        return True
    a = tuple(value / SQRT2 for value in candidate.a15.weights)
    d = tuple(value / SQRT2 for value in candidate.d9.weights)
    if not all(value.is_integral for value in (*a, *d)):
        return True
    return not a15_d9_glue_membership(a, d)


@lru_cache(maxsize=None)
def enumerate_p7_a15_d9_level(
    level: SilverP7D24ControlLevel,
    affine_caps_active: bool,
) -> tuple[SilverA15D9Candidate, ...]:
    """Enumerate one complete glue-aware prime-seven shell."""

    if level.norm not in P7_NORM_CLASSES or not p7_factorization_audit():
        raise ValueError("unsupported prime-seven A15D9 level")
    scale = P7_DELTA if _has_delta(level.norm) else Silver(1)
    reduced_target = level.target / (scale * scale)
    reduced_norm = level.norm / scale
    if not reduced_target.is_integral or any(
        reduced_target.sign(place) < 0 for place in (0, 1)
    ):
        return ()
    cap = P7_TIMELIKE * level.control / (2 * scale) if affine_caps_active else None
    d9 = enumerate_dominant_d9_weights(reduced_target, reduced_norm, cap)
    d9_vectors = tuple(
        vector
        for vector in d9.vectors
        if all((reduced_target - vector.square).sign(place) >= 0 for place in (0, 1))
    )
    targets = tuple(sorted({reduced_target - vector.square for vector in d9_vectors}))
    a15 = enumerate_dominant_a15_weights(targets, reduced_norm, cap)
    by_complement: dict[Silver, list[DominantWeight]] = {}
    for vector in d9_vectors:
        by_complement.setdefault(reduced_target - vector.square, []).append(vector)
    output = []
    for a_reduced in a15.vectors:
        for d_reduced in by_complement.get(a_reduced.square, ()):
            a = _scaled_weight(a_reduced, scale, a15_residue)
            d = _scaled_weight(d_reduced, scale, d9_residue)
            if not a15_d9_glue_membership(a.weights, d.weights):
                continue
            candidate = SilverA15D9Candidate(
                level.norm,
                level.height,
                level.control,
                level.target,
                a,
                d,
            )
            if _primitive(candidate):
                if p7_a15_d9_pair(candidate, candidate) != level.norm:
                    raise AssertionError("prime-seven A15D9 reconstruction changed norm")
                output.append(candidate)
    return tuple(
        sorted(output, key=lambda root: (root.control, root.a15.weights, root.d9.weights))
    )


@dataclass(frozen=True, slots=True)
class SilverP7A15D9Packet:
    height: Silver
    controls: tuple[SilverP7D24ControlLevel, ...]
    candidate_count: int
    rejected_count: int
    roots: tuple[SilverA15D9Candidate, ...]


@dataclass(frozen=True, slots=True)
class SilverP7A15D9Prefix:
    height_bound: int
    controls: tuple[SilverP7D24ControlLevel, ...]
    packets: tuple[SilverP7A15D9Packet, ...]
    accepted_roots: tuple[SilverA15D9Candidate, ...]


@lru_cache(maxsize=None)
def run_p7_a15_d9_prefix(height_bound: int = 63) -> SilverP7A15D9Prefix:
    controls = enumerate_p7_d24_control_levels(height_bound)
    accepted: list[SilverA15D9Candidate] = []
    packets = []
    for height, iterator in groupby(controls, key=lambda level: level.height):
        packet_controls = tuple(iterator)
        candidates = tuple(
            root
            for level in packet_controls
            for root in enumerate_p7_a15_d9_level(level, bool(accepted))
        )
        survivors = tuple(
            root
            for root in candidates
            if all(p7_a15_d9_pair(root, wall).sign(0) <= 0 for wall in accepted)
        )
        if not _same_height_acute(survivors):
            raise AssertionError("same-height prime-seven A15D9 walls are not acute")
        accepted.extend(survivors)
        packets.append(
            SilverP7A15D9Packet(
                height,
                packet_controls,
                len(candidates),
                len(candidates) - len(survivors),
                survivors,
            )
        )
    return SilverP7A15D9Prefix(height_bound, controls, tuple(packets), tuple(accepted))


def p7_d9_axis_witness(
    index: int = 6,
    weight: Fraction = Fraction(13, 25),
    height_bound: int = 30,
) -> SilverA15D9Witness:
    """Evaluate the frozen strict prime-seven D9 fundamental-axis ray."""

    if not 0 <= index < 9 or weight <= 0:
        raise ValueError("invalid D9 axis witness")
    a_weights = (SILVER_ZERO,) * 15
    d_weights = tuple(
        Silver(weight) if position == index else SILVER_ZERO
        for position in range(9)
    )
    inverse = _fraction_inverse(cartan_gram(D9_COMPONENT))
    norm = Silver(weight * weight * inverse[index][index]) - P7_TIMELIKE
    prefix = run_p7_a15_d9_prefix(height_bound)
    pairings = tuple(
        -P7_TIMELIKE * root.control
        + Silver(weight) * root.d9.coordinates[index]
        for root in prefix.accepted_roots
    )
    return SilverA15D9Witness(
        a_weights,
        d_weights,
        norm,
        pairings,
        True,
        all(pairing.sign(0) < 0 for pairing in pairings),
        norm.sign(0) > 0,
    )


__all__ = [
    "COMPLETE_P7_A15_D9_NORM_PROOF",
    "FIRST_PACKET_HEIGHT",
    "P7_NORM_CLASSES",
    "SECOND_PACKET_HEIGHT",
    "SilverP7A15D9Packet",
    "SilverP7A15D9Prefix",
    "enumerate_p7_a15_d9_level",
    "p7_a15_d9_pair",
    "p7_d9_axis_witness",
    "run_p7_a15_d9_prefix",
]
