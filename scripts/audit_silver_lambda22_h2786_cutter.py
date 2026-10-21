#!/usr/bin/env python3
"""Independent exact audit of the proposed fourth silver Lambda22 cutter.

This file is intentionally a work-in-progress until the chronological scan
and every surviving arithmetic extension have both been resolved.  The fixed
spatial vectors are ``d=e1`` and ``e=e2+e22`` in the catalogue basis, with
Gram matrix ``[[4,2],[2,10]]``.  The proposed norm-two root is

``(82+58*pi, (26d+e)+pi*(18d+e))``.

The scan below starts from the complete mirror-norm/control alphabet, applies
the binary minimum and local pairing gates, and then performs an exact 4-by-4
Schur extension.  It also imposes acuteness to the H14 wall based on ``e`` and
to the accepted H738 wall, both of which lie in the fixed two-plane.  The H164
and H1598 constraints and the final lattice-minimum resolutions are added in
the certificate stage below; no acceptance conclusion is drawn merely from
this necessary scan.
"""

from __future__ import annotations

import argparse
import itertools
import json
import shutil
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import isqrt, lcm
from pathlib import Path

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_decomp

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.cubic7_lambda22 import (
    lambda22_integer_gram_minimum_at_least,
)
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    enumerate_lambda22_controls,
    gated_relaxed_witness_bound,
    lambda22_gram,
    lambda22_pair_silver,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)


HEIGHT = Silver(2786, 1970)
CONTROL = Silver(82, 58)
REPLACEMENT_SCALE = Fraction(28, 87)
H14_CONTROL = Silver(6, 4)
H738_CONTROL = Silver(42, 30)
NORM_GATES = {
    Silver(1): (1, 1),
    Silver(2): (1, 1),
    Silver(4, 2): (2, 1),
    Silver(4): (2, 2),
    Silver(3): (3, 3),
    Silver(6): (3, 3),
    Silver(12, 6): (6, 3),
    Silver(12): (6, 6),
}
MINIMUM = {1: 4, 2: 12, 3: 24, 6: 96}
ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_lambda22_h2786_cutter.json"
REPORT_PATH = ROOT / "results/silver_lambda22_h2786_cutter.md"
D = (1,) + (0,) * 21
E = (0, 1) + (0,) * 19 + (1,)
U = tuple(26 * D[i] + E[i] for i in range(22))
V = tuple(18 * D[i] + E[i] for i in range(22))
SPATIAL = tuple(Silver(U[i], V[i]) for i in range(22))

# These disjoint slices were replayed independently after the exact hot loop
# had been changed from Fraction objects to integer pairs.  Their union is
# precisely the 374-control chronological scan through the candidate.
PINNED_SLICES = (
    (0, 276, 1067, 14, 17, 6),
    (276, 287, 67, 1, 1, 0),
    (287, 300, 112, 2, 2, 0),
    (300, 310, 114, 2, 2, 0),
    (310, 320, 100, 0, 0, 0),
    (320, 330, 86, 0, 0, 0),
    (330, 340, 96, 0, 0, 0),
    (340, 350, 127, 2, 2, 1),
    (350, 360, 68, 0, 0, 0),
    (360, 370, 86, 0, 0, 0),
    (370, 374, 74, 2, 2, 1),
)


def norm_gate(square: int, modulus: int) -> bool:
    return (
        (not square or square >= MINIMUM[modulus])
        and square % (2 * modulus) == 0
        and (modulus % 2 or square % 8 in (0, 4))
    )


@dataclass(frozen=True, slots=True)
class OneRow:
    candidate: tuple[int, int]
    h14: tuple[int, int]
    h738: tuple[int, int]
    d_dot: int
    e_dot: int
    schur: int


@lru_cache(maxsize=None)
def one_rows(square: int, modulus: int, kind: int) -> tuple[OneRow, ...]:
    """All one-vector extensions of Gram(d,e), sorted for the target wall."""

    rows: list[OneRow] = []
    for d_dot in range(-isqrt(4 * square), isqrt(4 * square) + 1):
        if d_dot % modulus:
            continue
        for e_dot in range(-isqrt(10 * square), isqrt(10 * square) + 1):
            if e_dot % modulus:
                continue
            # det Gram(d,e)=36 and adjugate [[10,-2],[-2,4]].
            schur = (
                36 * square
                - 10 * d_dot * d_dot
                + 4 * d_dot * e_dot
                - 4 * e_dot * e_dot
            )
            if schur < 0:
                continue
            if kind == 0:
                candidate = (26 * d_dot + e_dot, 18 * d_dot + e_dot)
                h14 = (e_dot, e_dot)
                h738 = (14 * d_dot - e_dot, 10 * d_dot - e_dot)
            else:
                candidate = (36 * d_dot + 2 * e_dot, 26 * d_dot + e_dot)
                h14 = (2 * e_dot, e_dot)
                h738 = (20 * d_dot - 2 * e_dot, 14 * d_dot - e_dot)
            rows.append(OneRow(candidate, h14, h738, d_dot, e_dot, schur))
    # Deliberately leave the exact rows unsorted.  At the larger controls the
    # repeated quadratic-field comparisons needed by sorting cost more than
    # the Schur enumeration itself.  The join below uses exact monotone bounds
    # in the remaining dot product instead.
    return tuple(rows)


@dataclass(frozen=True, slots=True)
class Extension:
    pairing: Silver
    d_dot_u: int
    d_dot_v: int
    e_dot_u: int
    e_dot_v: int
    schur_a: int
    schur_b: int
    schur_n: int
    h14_pairing: Silver
    h738_pairing: Silver


def cross_extensions(
    u_square: int,
    v_square: int,
    u_dot_v: int,
    u_modulus: int,
    v_modulus: int,
    other_control: Silver,
) -> tuple[Extension, ...]:
    """Exact nonnegative extensions acute to the fixed H14/H738 walls."""

    target_threshold = LAMBDA22_COEFFICIENT * CONTROL * other_control
    h14_threshold = LAMBDA22_COEFFICIENT * H14_CONTROL * other_control
    h738_threshold = LAMBDA22_COEFFICIENT * H738_CONTROL * other_control
    u_rows = one_rows(u_square, u_modulus, 0)
    v_rows = one_rows(v_square, v_modulus, 1)
    if not u_rows or not v_rows:
        return ()
    # All coefficients here are integral.  Inline pair arithmetic avoids the
    # Fraction normalisation overhead in millions of hot-loop comparisons.
    def add_pair(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        return left[0] + right[0], left[1] + right[1]

    def sub_pair(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        return left[0] - right[0], left[1] - right[1]

    def sign_pair(value: tuple[int, int]) -> int:
        a, b = value
        if not a:
            return (b > 0) - (b < 0)
        if not b or (a > 0) == (b > 0):
            return (a > 0) - (a < 0)
        return ((a > 0) - (a < 0)) if a * a > 2 * b * b else ((b > 0) - (b < 0))

    def ge_pair(left: tuple[int, int], right: tuple[int, int]) -> bool:
        return sign_pair(sub_pair(left, right)) >= 0

    def le_pair(left: tuple[int, int], right: tuple[int, int]) -> bool:
        return sign_pair(sub_pair(left, right)) <= 0

    def pair(value: Silver) -> tuple[int, int]:
        if value.a.denominator != 1 or value.b.denominator != 1:
            raise AssertionError("the chronological join ceased to be integral")
        return int(value.a), int(value.b)

    target_pair = pair(target_threshold)
    h14_pair = pair(h14_threshold)
    h738_pair = pair(h738_threshold)
    v_maximum = max(
        (row.candidate for row in v_rows),
        key=lambda value: float(value[0]) + 1.4142135623730951 * value[1],
    )

    def first_true(lower: int, upper: int, predicate: object) -> int:
        """First integer in a monotone-false/true interval; upper+1 if none."""

        lo, hi = lower, upper + 1
        while lo < hi:
            middle = (lo + hi) // 2
            if predicate(middle):  # type: ignore[operator]
                hi = middle
            else:
                lo = middle + 1
        return lo

    def last_true(lower: int, upper: int, predicate: object) -> int:
        """Last integer in a monotone-true/false interval; lower-1 if none."""

        lo, hi = lower - 1, upper
        while lo < hi:
            middle = (lo + hi + 1) // 2
            if predicate(middle):  # type: ignore[operator]
                lo = middle
            else:
                hi = middle - 1
        return lo

    output: list[Extension] = []
    for u_row in u_rows:
        u_candidate_pair = u_row.candidate
        u_h14_pair = u_row.h14
        u_h738_pair = u_row.h738
        if not ge_pair(add_pair(u_candidate_pair, v_maximum), target_pair):
            continue
        q_bound = isqrt(4 * v_square)
        e_bound = isqrt(10 * v_square)
        for d_dot_v in range(-q_bound, q_bound + 1):
            if d_dot_v % v_modulus:
                continue
            # The transparent floor forms are easier to audit than float
            # roots; all three following bounds are exact Silver comparisons.
            k_lower = -(e_bound // v_modulus)
            k_upper = e_bound // v_modulus

            def v_values(index: int) -> tuple[int, tuple[int, int], tuple[int, int], tuple[int, int]]:
                e_dot_v = v_modulus * index
                return (
                    e_dot_v,
                    (36 * d_dot_v + 2 * e_dot_v, 26 * d_dot_v + e_dot_v),
                    (2 * e_dot_v, e_dot_v),
                    (20 * d_dot_v - 2 * e_dot_v, 14 * d_dot_v - e_dot_v),
                )

            candidate_lower = first_true(
                k_lower,
                k_upper,
                lambda index: ge_pair(
                    add_pair(u_candidate_pair, v_values(index)[1]), target_pair
                ),
            )
            h738_lower = first_true(
                k_lower,
                k_upper,
                lambda index: le_pair(
                    add_pair(u_h738_pair, v_values(index)[3]), h738_pair
                ),
            )
            h14_upper = last_true(
                k_lower,
                k_upper,
                lambda index: le_pair(
                    add_pair(u_h14_pair, v_values(index)[2]), h14_pair
                ),
            )
            lower = max(k_lower, candidate_lower, h738_lower)
            upper = min(k_upper, h14_upper)
            for index in range(lower, upper + 1):
                e_dot_v, candidate_v, h14_v, h738_v = v_values(index)
                schur_b = (
                    36 * v_square
                    - 10 * d_dot_v * d_dot_v
                    + 4 * d_dot_v * e_dot_v
                    - 4 * e_dot_v * e_dot_v
                )
                if schur_b < 0:
                    continue
                spatial_pairing_pair = sub_pair(add_pair(u_candidate_pair, candidate_v), target_pair)
                h14_pairing_pair = sub_pair(add_pair(u_h14_pair, h14_v), h14_pair)
                h738_pairing_pair = sub_pair(add_pair(u_h738_pair, h738_v), h738_pair)
                # Redundant exact assertions guard the monotone interval code.
                if (
                    sign_pair(spatial_pairing_pair) < 0
                    or sign_pair(h14_pairing_pair) > 0
                    or sign_pair(h738_pairing_pair) > 0
                ):
                    raise AssertionError("an exact monotone join bound failed")
                schur_n = 36 * u_dot_v - (
                10 * u_row.d_dot * d_dot_v
                - 2
                * (
                    u_row.d_dot * e_dot_v
                    + u_row.e_dot * d_dot_v
                )
                + 4 * u_row.e_dot * e_dot_v
            )
                if u_row.schur * schur_b < schur_n * schur_n:
                    continue
                output.append(
                    Extension(
                        Silver(*spatial_pairing_pair),
                        u_row.d_dot,
                        d_dot_v,
                        u_row.e_dot,
                        e_dot_v,
                        u_row.schur,
                        schur_b,
                        schur_n,
                        Silver(*h14_pairing_pair),
                        Silver(*h738_pairing_pair),
                    )
                )
    return tuple(output)


def scan(start: int = 0, stop: int | None = None) -> None:
    levels = enumerate_lambda22_controls(HEIGHT)
    if stop is None:
        stop = len(levels)
    if not (0 <= start <= stop <= len(levels)):
        raise ValueError("invalid control-index slice")
    type_count = row_count = extension_count = positive_count = 0
    for level_index, level in enumerate(levels):
        if level_index < start or level_index >= stop:
            continue
        if level.height <= Silver(14, 10):
            continue
        u_modulus, v_modulus = NORM_GATES[level.norm]
        for item in candidate_leech_gram_triples(level.target):
            if has_short_binary_combination(item):
                continue
            gram = (item.u_square, item.v_square, item.u_dot_v)
            if (
                not norm_gate(gram[0], u_modulus)
                or not norm_gate(gram[1], v_modulus)
                or gram[2] % lcm(u_modulus, v_modulus)
            ):
                continue
            type_count += 1
            extensions = cross_extensions(
                *gram, u_modulus, v_modulus, level.control
            )
            if not extensions:
                continue
            row_count += 1
            extension_count += len(extensions)
            positives = tuple(row for row in extensions if row.pairing.sign(0) > 0)
            positive_count += len(positives)
            print(
                "ROW",
                str(level.height),
                str(level.norm),
                str(level.control),
                gram,
                "extensions",
                len(extensions),
                "positives",
                len(positives),
                flush=True,
            )
            for extension in positives:
                print("POSITIVE", str(level.height), gram, extension, flush=True)
        if level_index % 25 == 0:
            print(
                "PROGRESS",
                level_index,
                len(levels),
                "types",
                type_count,
                "rows",
                row_count,
                "positives",
                positive_count,
                flush=True,
            )
    print(
        "SUMMARY",
        f"slice={start}:{stop}",
        f"controls={len(levels)}",
        f"types={type_count}",
        f"rows={row_count}",
        f"extensions={extension_count}",
        f"positives={positive_count}",
    )


def value_json(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(0), value.sign(1)],
    }


def explicit_root_audit() -> dict[str, object]:
    gram = lambda22_gram()
    coefficient_gram = (
        lambda22_pair_silver(U, U),
        lambda22_pair_silver(V, V),
        lambda22_pair_silver(U, V),
    )
    norm = -LAMBDA22_COEFFICIENT * CONTROL * CONTROL + lambda22_pair_silver(
        SPATIAL, SPATIAL
    )
    spatial_pairings = tuple(
        sum(
            (Silver(gram[i][j]) * SPATIAL[j] for j in range(22)),
            start=Silver(),
        )
        for i in range(22)
    )
    reflective = (
        (2 * LAMBDA22_COEFFICIENT * CONTROL / norm).is_integral
        and all((2 * pairing / norm).is_integral for pairing in spatial_pairings)
    )
    old_ray_pairing = (
        -LAMBDA22_COEFFICIENT * CONTROL
        + Silver(10) * Silver(106, 74) / 31
    )
    if (
        (lambda22_pair_silver(D, D), lambda22_pair_silver(D, E), lambda22_pair_silver(E, E))
        != (4, 2, 10)
        or coefficient_gram != (2818, 1378, 1970)
        or norm != 2
        or LAMBDA22_COEFFICIENT * CONTROL * CONTROL / norm != HEIGHT
        or not reflective
        or Silver(1, 1) not in SPATIAL
        or old_ray_pairing != Silver(6, -4) / 31
        or old_ray_pairing.sign(0) <= 0
    ):
        raise AssertionError("the H2786 explicit-root arithmetic changed")
    return {
        "d": list(D),
        "e": list(E),
        "Gram_d_e": [[4, 2], [2, 10]],
        "u": list(U),
        "v": list(V),
        "coefficient_gram": [[2818, 1970], [1970, 1378]],
        "control": value_json(CONTROL),
        "norm": value_json(norm),
        "height": value_json(HEIGHT),
        "primitive": True,
        "reflection_integral": reflective,
        "old_10d_over_31_ray_pairing": value_json(old_ray_pairing),
    }


def accepted_wall_acuteness() -> dict[str, object]:
    """Complete H14 shell maxima and direct later-wall pairings."""

    packets = []
    for name, shell_square, control, modulus in (
        ("ordinary_q2", 10, Silver(6, 4), 1),
        ("dyadic_q4_S2", 20, Silver(8, 6), 2),
    ):
        best: tuple[Silver, int, int] | None = None
        psd = minimum_four = 0
        for d_dot in range(-isqrt(4 * shell_square), isqrt(4 * shell_square) + 1):
            if d_dot % modulus:
                continue
            for e_dot in range(-isqrt(10 * shell_square), isqrt(10 * shell_square) + 1):
                if e_dot % modulus:
                    continue
                determinant = (
                    36 * shell_square
                    - 10 * d_dot * d_dot
                    + 4 * d_dot * e_dot
                    - 4 * e_dot * e_dot
                )
                if determinant < 0:
                    continue
                psd += 1
                extension_gram = (
                    (4, 2, d_dot),
                    (2, 10, e_dot),
                    (d_dot, e_dot, shell_square),
                )
                if not lambda22_integer_gram_minimum_at_least(extension_gram, 4):
                    continue
                minimum_four += 1
                spatial = Silver(1, 1) * (
                    Silver(26, 18) * d_dot + Silver(1, 1) * e_dot
                )
                if best is None or spatial > best[0]:
                    best = spatial, d_dot, e_dot
        if best is None:
            raise AssertionError("an H14 shell lost every extension")
        margin = best[0] - LAMBDA22_COEFFICIENT * CONTROL * control
        packets.append(
            {
                "packet": name,
                "PSD_extension_count": psd,
                "minimum_four_extension_count": minimum_four,
                "maximizer": [best[1], best[2]],
                "maximum_pairing": value_json(margin),
            }
        )
    if [row["maximum_pairing"]["display"] for row in packets] != [
        "-65 - 46*sqrt(2)",
        "-40 - 28*sqrt(2)",
    ]:
        raise AssertionError("the complete H14 maxima changed")

    a164 = (-6, 31, 21, 15, 10, 20, 14, 7, 4, 2, 2, 1) + (0,) * 10
    f1598 = (0, 0, -1) + (0,) * 13 + (-1,) + (0,) * 5
    g738 = tuple(-entry for entry in E)
    prior_spatial = {
        "H164": (
            Silver(28, 20),
            tuple(Silver(8 * D[i] - a164[i], 6 * D[i] - a164[i]) for i in range(22)),
        ),
        "H738": (
            Silver(42, 30),
            tuple(Silver(14 * D[i] + g738[i], 10 * D[i] + g738[i]) for i in range(22)),
        ),
        "H1598": (
            Silver(62, 44),
            tuple(Silver(20 * D[i] + f1598[i], 14 * D[i] + f1598[i]) for i in range(22)),
        ),
    }
    pairings = {
        name: -LAMBDA22_COEFFICIENT * CONTROL * control
        + lambda22_pair_silver(SPATIAL, spatial)
        for name, (control, spatial) in prior_spatial.items()
    }
    if pairings != {
        "H164": Silver(-52, -36),
        "H738": Silver(-58, -40),
        "H1598": Silver(-19, -14),
    } or any(value.sign(0) >= 0 for value in pairings.values()):
        raise AssertionError("accepted later-wall acuteness changed")
    return {
        "complete_H14_packets": packets,
        "post_H14_pairings": {name: value_json(value) for name, value in pairings.items()},
        "all_strict": True,
    }


def h2404_embedding_exclusion() -> dict[str, object]:
    """Exhaust the sole new earlier positive row as a mod-nine coset."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required for the H2404 LLL replay")
    gram = Matrix(lambda22_gram())
    d = Matrix(D)
    e = Matrix(E)
    constraints = Matrix.vstack(d.T, e.T) * gram
    smith, left, right = smith_normal_decomp(constraints, domain=ZZ)
    if left * constraints * right != smith or tuple(smith.diagonal()) != (1, 1):
        raise AssertionError("the fixed-plane Smith decomposition changed")
    kernel = right[:, 2:]
    kernel_gram = kernel.T * gram * kernel
    if constraints * kernel != Matrix.zeros(2, 20) or int(kernel_gram.det()) != 432:
        raise AssertionError("the orthogonal kernel changed")

    gp_matrix = "[" + ";".join(
        ",".join(str(int(kernel_gram[i, j])) for j in range(20))
        for i in range(20)
    ) + "]"
    process = subprocess.run(
        [gp, "-q"],
        input=f"H={gp_matrix};T=qflllgram(H);print(T);quit()\n",
        text=True,
        capture_output=True,
        timeout=30,
        check=True,
    )
    raw = process.stdout.strip()
    transform = Matrix(
        [[int(entry.strip()) for entry in row.split(",")] for row in raw[1:-1].split(";")]
    )
    reduced = transform.T * kernel_gram * transform
    if int(transform.det()) not in (-1, 1) or tuple(reduced.diagonal()) != (4,) * 20:
        raise AssertionError("the deterministic LLL reduction changed")

    # The Schur row forces u-v=10d.  Clearing the projection denominator gives
    # w=9u-304d-13e, with w^2=18, w perpendicular to d,e, and
    # w congruent to 2d-4e modulo 9.  Express that congruence in the reduced
    # kernel coordinates.
    target = 2 * d - 4 * e
    target_coordinates = (right.inv() * target)[2:, 0]
    residues = tuple(int(value) % 9 for value in transform.inv() * target_coordinates)
    inverse = reduced.inv()
    bounds = []
    for index in range(20):
        cap = 18 * inverse[index, index]
        bounds.append(isqrt(int(cap.p) // int(cap.q)))
    domains = tuple(
        tuple(value for value in range(-bounds[i], bounds[i] + 1) if value % 9 == residues[i])
        for i in range(20)
    )
    candidate_count = 1
    for domain in domains:
        candidate_count *= len(domain)
    solutions = 0
    for vector in itertools.product(*domains):
        square = sum(
            int(reduced[i, j]) * vector[i] * vector[j]
            for i in range(20)
            for j in range(20)
        )
        solutions += square == 18
    if candidate_count != 13_824 or solutions:
        raise AssertionError(
            f"the H2404 coset exhaustion changed: {candidate_count}, {solutions}"
        )
    return {
        "height": value_json(Silver(2404, 1700)),
        "norm": value_json(Silver(4)),
        "control": value_json(Silver(108, 76)),
        "coefficient_gram": [4780, 2420, 3400],
        "candidate_pairing": value_json(Silver(2)),
        "Schur_data": [8, 8, 8],
        "forced_relation": "u-v=10d",
        "cleared_vector": "w=9u-304d-13e",
        "cleared_norm": 18,
        "cleared_congruence": "w=2d-4e mod 9",
        "orthogonal_kernel_rank": 20,
        "orthogonal_kernel_gram_determinant": 432,
        "reduced_gram_diagonal": [4] * 20,
        "dual_coordinate_bounds": bounds,
        "mod9_residues": list(residues),
        "finite_candidate_count": candidate_count,
        "solution_count": solutions,
        "conclusion": "The H2404 positive necessary row has no integral Lambda22 embedding.",
    }


def chronological_certificate() -> dict[str, object]:
    totals = [sum(row[index] for row in PINNED_SLICES) for index in range(2, 6)]
    if PINNED_SLICES[0][:2] != (0, 276) or PINNED_SLICES[-1][:2] != (370, 374):
        raise AssertionError("the chronological slices no longer cover the endpoints")
    if any(PINNED_SLICES[i][1] != PINNED_SLICES[i + 1][0] for i in range(len(PINNED_SLICES) - 1)):
        raise AssertionError("the chronological slices have a gap or overlap")
    if totals != [1997, 23, 26, 8]:
        raise AssertionError(f"the chronological totals changed: {totals}")
    return {
        "control_count": 374,
        "post_H14_control_count": 349,
        "gate_compatible_post_H14_type_count": totals[0],
        "nonnegative_row_count": totals[1],
        "nonnegative_extension_count": totals[2],
        "positive_extension_count": totals[3],
        "independent_slice_replays": [
            {
                "start": start,
                "stop": stop,
                "types": types,
                "rows": rows,
                "extensions": extensions,
                "positives": positives,
            }
            for start, stop, types, rows, extensions, positives in PINNED_SLICES
        ],
        "positive_resolutions": [
            "The two H21 q6 extensions have nonzero u-v in S3 of norm 6<24.",
            "The H42 q12 extension has nonzero u-v in S6 of norm 12<96.",
            "The H56 q4 extension has nonzero u-v in S2 of norm 4<12.",
            "The H274 q4 extension has nonzero u-12d in S2 of norm 4<12.",
            "The H656 q4 extension has nonzero 2u-3v in S2 of norm 4<12.",
            "The H2404 q4 extension is empty by the archived 13,824-vector mod-nine coset exhaustion.",
            "The sole remaining positive is the displayed H2786 norm-two root.",
        ],
        "H2404_exact_exclusion": h2404_embedding_exclusion(),
        "conclusion": "The displayed H2786 norm-two root is chronologically accepted.",
    }


REPLACEMENT_THREATS = (
    (
        "14 + 10*sqrt(2)", "4", "8 + 6*sqrt(2)",
        (28, 16, 20), (10, 8), Silver(Fraction(-68, 87), Fraction(50, 87)),
        "The relaxed d-pairing 10 is not attained in the complete H14 S2 norm-20 shell; its exact bound is |d.w|<=8.",
    ),
    (
        "21 + 15*sqrt(2)", "6", "12 + 9*sqrt(2)",
        (60, 36, 45), (15, 12), Silver(Fraction(-34, 29), Fraction(25, 29)),
        "u-v is a nonzero S3 vector of norm 6 below min(S3)=24.",
    ),
    (
        "42 + 30*sqrt(2)", "12", "24 + 18*sqrt(2)",
        (228, 144, 180), (30, 24), Silver(Fraction(-68, 29), Fraction(50, 29)),
        "u-v is a nonzero S6 vector of norm 12 below min(S6)=96.",
    ),
    (
        "56 + 40*sqrt(2)", "4", "16 + 12*sqrt(2)",
        (100, 64, 80), (20, 16), Silver(Fraction(-136, 87), Fraction(100, 87)),
        "u-v is a nonzero S2 vector of norm 4 below min(S2)=12.",
    ),
    (
        "1202 + 850*sqrt(2)", "4", "76 + 54*sqrt(2)",
        (2500, 1156, 1700), (100, 68), Silver(Fraction(16, 87), Fraction(-10, 87)),
        "Schur zero forces u=25d and v=17d; this violates the q4 spatial gate because Gd is not even.",
    ),
)


def full_replacement_scan() -> tuple[int, tuple[tuple[object, ...], ...]]:
    """Replay every gated necessary wall against the proposed 28d/87 ray."""

    retained = 0
    threats: list[tuple[object, ...]] = []
    scale = Silver(REPLACEMENT_SCALE)
    for level in enumerate_lambda22_controls(HEIGHT):
        u_modulus, v_modulus = NORM_GATES[level.norm]
        for item in candidate_leech_gram_triples(level.target):
            if has_short_binary_combination(item):
                continue
            gram = item.u_square, item.v_square, item.u_dot_v
            if (
                not norm_gate(gram[0], u_modulus)
                or not norm_gate(gram[1], v_modulus)
                or gram[2] % lcm(u_modulus, v_modulus)
            ):
                continue
            try:
                d_dot_u, d_dot_v, maximum = gated_relaxed_witness_bound(
                    *gram, u_modulus, v_modulus
                )
            except AssertionError:
                continue
            retained += 1
            margin = scale * maximum - LAMBDA22_COEFFICIENT * level.control
            if margin.sign(0) >= 0:
                threats.append(
                    (
                        str(level.height), str(level.norm), str(level.control),
                        gram, (d_dot_u, d_dot_v), margin,
                    )
                )
    return retained, tuple(threats)


def replacement_ray_audit(full: bool = False) -> dict[str, object]:
    expected = tuple(row[:6] for row in REPLACEMENT_THREATS)
    if full:
        retained, threats = full_replacement_scan()
        if retained != 2001 or threats != expected:
            raise AssertionError(
                f"the H2786 replacement frontier changed: {retained}, {threats!r}"
            )
    scale = Silver(REPLACEMENT_SCALE)
    norm = 4 * scale * scale - LAMBDA22_COEFFICIENT
    accepted_pairings = {
        "H164": scale * Silver(34, 26) - LAMBDA22_COEFFICIENT * Silver(28, 20),
        "H738": scale * Silver(54, 38) - LAMBDA22_COEFFICIENT * Silver(42, 30),
        "H1598": scale * Silver(80, 56) - LAMBDA22_COEFFICIENT * Silver(62, 44),
        "H2786": scale * Silver(106, 74) - LAMBDA22_COEFFICIENT * CONTROL,
    }
    h14_actual = scale * Silver(8, 8) - LAMBDA22_COEFFICIENT * Silver(8, 6)
    if (
        norm != Silver(Fraction(10705, 7569), -1)
        or norm.sign(0) <= 0
        or h14_actual != Silver(Fraction(-124, 87), Fraction(50, 87))
        or h14_actual.sign(0) >= 0
        or accepted_pairings
        != {
            "H164": Silver(Fraction(-92, 87), Fraction(32, 87)),
            "H738": Silver(Fraction(-18, 29), Fraction(20, 87)),
            "H1598": Silver(Fraction(-22, 87), Fraction(2, 87)),
            "H2786": Silver(Fraction(10, 87), Fraction(-16, 87)),
        }
        or any(value.sign(0) >= 0 for value in accepted_pairings.values())
    ):
        raise AssertionError("the 28d/87 replacement arithmetic changed")
    return {
        "vector": "z=(1,(28/87)d)",
        "scale": "28/87",
        "norm": value_json(norm),
        "control_count": 374,
        "retained_necessary_type_count": 2001,
        "relaxed_nonnegative_row_count": len(REPLACEMENT_THREATS),
        "relaxed_nonnegative_rows": [
            {
                "height": height,
                "norm": root_norm,
                "control": control,
                "gram": list(gram),
                "relaxed_pairings": list(pairings),
                "margin": value_json(margin),
                "resolution": resolution,
            }
            for height, root_norm, control, gram, pairings, margin, resolution
            in REPLACEMENT_THREATS
        ],
        "actual_H14_q4_margin": value_json(h14_actual),
        "accepted_wall_pairings": {
            name: value_json(value) for name, value in accepted_pairings.items()
        },
        "conclusion": (
            "All five relaxed threats are exactly excluded; the spacelike ray "
            "is strict against every possible wall through H2786."
        ),
    }


def report() -> str:
    return """# Silver Lambda22: fourth accepted cutter at H2786

For `<-(sqrt(2)-1)> + Lambda22` in H22, put `d=e1` and
`e=e2+e22` in the catalogue basis.  Their Gram matrix is `[[4,2],[2,10]]`.
The primitive reflection-integral norm-two root

`r=(82+58sqrt(2),(26d+e)+sqrt(2)(18d+e))`

has height `2786+1970sqrt(2)`.  It pairs positively with the former
`(1,10d/31)` hyperideal ray, so it genuinely cuts that ray.

The complete chronological necessary scan contains 374 controls, 1,997
post-H14 gate-compatible Gram types, 23 nonnegative rows, 26 extensions, and
eight positive extensions.  Six are the previously certified S2/S3/S6
minimum contradictions.  The only new earlier positive is a q4 row at
`H=2404+1700sqrt(2)`.  Its Schur equations force `u-v=10d`; then
`w=9u-304d-13e` would have norm 18, be orthogonal to `d,e`, and satisfy
`w=2d-4e mod 9`.  An exact Smith/LLL reduction leaves 13,824 vectors inside
the dual coordinate bounds, and none has norm 18.  Thus that row is empty.
The last positive is exactly the displayed root, which is therefore accepted.

Complete H14 shell maximization and direct calculation show strict acuteness
to both H14 packets and the accepted H164, H738, and H1598 cutters.

The old `10d/31` ray is replaced by `z=(1,(28/87)d)`.  Its exact norm is
`10705/7569-sqrt(2)>0`, and it lies strictly behind all four accepted later
walls.  The complete 374-control necessary-wall replay retains 2,001 gated
types and only five relaxed threats.  The H14 threat is removed by the exact
S2 shell bound; three are the frozen S3/S6/S2 minimum contradictions; and
the H1202 q4 threat violates its reflection gate.  Thus `z` is strict against
every possible wall through H2786, so this complete bounded prefix still has
an open hyperideal patch.

This is not a finite-volume example, a Vinberg-termination claim, or a global
nonreflectivity theorem.

```sh
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h2786_cutter.py --freeze
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h2786_cutter.py --start 340 --stop 350
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h2786_cutter.py
```
"""


def freeze() -> None:
    payload = {
        "schema_version": "0.2.0",
        "seed_id": "silver-H22-Lambda22-a-sqrt2-1",
        "scope": (
            "Exact fourth accepted cutter and strict replacement hyperideal ray "
            "through H2786; bounded Vinberg prefix only, with no finite-volume, "
            "termination, or global nonreflectivity claim."
        ),
        "explicit_cutter": explicit_root_audit(),
        "complete_acuteness": accepted_wall_acuteness(),
        "chronological_acceptance": chronological_certificate(),
        "replacement_hyperideal_ray": replacement_ray_audit(full=True),
        "conclusion": (
            "The H2786 root is a fourth genuine accepted cutter, but 28d/87 "
            "is a strict hyperideal replacement through the complete prefix."
        ),
        "discovery": None,
    }
    JSON_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    REPORT_PATH.write_text(report())
    print(JSON_PATH)
    print(REPORT_PATH)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--freeze", action="store_true")
    arguments = parser.parse_args()
    if arguments.freeze:
        freeze()
    else:
        scan(arguments.start, arguments.stop)


if __name__ == "__main__":
    main()
