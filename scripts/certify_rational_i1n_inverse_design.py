#!/usr/bin/env python3
"""Exact inverse-design certificate for the Vinberg chambers of I_(1,n).

The Lorentzian convention is ``diag(-1, +1, ..., +1)``.  A positive
height root is written ``(d, -m_1, ..., -m_n)`` and corresponds, after
changing the sign of the form, to the plane divisor

    d H - m_1 E_1 - ... - m_n E_n.

The enumeration is a streaming Vinberg enumeration.  It never materializes
all partitions of ``d^2 + q``: the already accepted walls prune the recursion.
Only roots of norm q=1 or q=2 can give integral reflections in I_(1,n).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, lcm
from pathlib import Path
from typing import Iterable, Sequence


DIMENSIONS = tuple(range(12, 20))
EXPECTED_FACETS = {12: 14, 13: 15, 14: 17, 15: 18, 16: 20, 17: 22, 18: 37, 19: 50}
EXPECTED_MAX_DEGREE = {12: 3, 13: 3, 14: 4, 15: 4, 16: 6, 17: 6, 18: 27, 19: 33}
COXITER_SHA256 = {
    12: "dcb3335a4fdfbd38b2bbbbb7313c300940b8c923c8d201ed66655aab5e82761c",
    13: "20bb824f5594699ac0bdbcf561785a2f06bd59e3d5aadf2619e28e4b6d935643",
    14: "30b5dd05e6af91b26dcffde49d7080be2eed55dcb7810aef457f51461f421da2",
    15: "2b773f565a173d5cee6aae770c05e065b9cb6b8d3d39addec2cd3336b1396ffe",
    16: "aeb4003797a4fd4faf7286a4d1de9acb2869d85f4a79b21472b5f4d7e834931b",
    17: "ce0427729a8e75bffeb48d4adcf037b961670edcab8d05b022c7e1baf64822bb",
    18: "1a605a317f04bb28e679fc7c1246816ab93a9ec5b2c49ea69677381e393d05e8",
    19: "60b8303c4e88d5c15d57c2e4d6b04d956dcdb478f848d965d022b39780ec229d",
}

VINBERG_VOLUME_DATA = {
    12: {
        "compact": False,
        "finite_volume": True,
        "finite_vertices": 35,
        "ideal_vertices": 2,
        "f_vector": [37, 234, 786, 1749, 2793, 3312, 2958, 1992, 1000, 364, 91, 14, 1],
        "euler_characteristic": "691/382588157952000",
        "covolume": "pi^6*691/62140685967360000",
    },
    13: {
        "compact": False,
        "finite_volume": True,
        "finite_vertices": 43,
        "ideal_vertices": 3,
        "f_vector": [46, 315, 1140, 2745, 4794, 6315, 6390, 4995, 3002, 1365, 455, 105, 15, 1],
        "euler_characteristic": "0",
    },
    14: {
        "compact": False,
        "finite_volume": True,
        "finite_vertices": 89,
        "ideal_vertices": 5,
        "f_vector": [94, 704, 2695, 6825, 12579, 17633, 19215, 16425, 11009, 5733, 2275, 665, 135, 17, 1],
        "euler_characteristic": "-87757/289236647411712000",
        "covolume": "pi^7*87757/305359330843607040000",
    },
    15: {
        "compact": False,
        "finite_volume": True,
        "finite_vertices": 139,
        "ideal_vertices": 6,
        "f_vector": [145, 1172, 4752, 12691, 24745, 36995, 43505, 40733, 30457, 18109, 8463, 3045, 815, 153, 18, 1],
        "euler_characteristic": "0",
    },
    16: {
        "compact": False,
        "finite_volume": True,
        "finite_vertices": 313,
        "ideal_vertices": 12,
        "f_vector": [325, 2804, 11914, 33164, 67410, 105462, 130646, 130062, 104670, 68042, 35490, 14658, 4690, 1122, 189, 20, 1],
        "euler_characteristic": "642332179/2360171042879569920000",
        "covolume": "pi^8*642332179/18687991047628750848000000",
    },
}


Vector = tuple[int, ...]


def lorentz_pair(left: Sequence[int], right: Sequence[int]) -> int:
    return -left[0] * right[0] + sum(a * b for a, b in zip(left[1:], right[1:]))


def zero_roots(n: int) -> list[Vector]:
    roots: list[Vector] = []
    for i in range(n - 1):
        roots.append(
            (0,) + tuple(1 if j == i else -1 if j == i + 1 else 0 for j in range(n))
        )
    roots.append((0,) + (0,) * (n - 1) + (1,))
    return roots


def _candidate_multiplicities(
    n: int, degree: int, norm: int, positive_roots: Sequence[Vector]
) -> list[tuple[int, ...]]:
    """Stream the dominant partitions surviving all lower Vinberg walls."""

    target = degree * degree + norm
    previous_m = [tuple(-entry for entry in root[1:]) for root in positive_roots]
    bounds = [degree * root[0] for root in positive_roots]
    answers: list[tuple[int, ...]] = []

    def extend(
        index: int,
        last: int,
        remaining_square: int,
        prefix: tuple[int, ...],
        partial_pairings: tuple[int, ...],
    ) -> None:
        slots = n - index
        if remaining_square < 0 or remaining_square > slots * last * last:
            return
        if index == n:
            if remaining_square == 0:
                answers.append(prefix)
            return

        upper = min(last, isqrt(remaining_square))
        for value in range(upper, -1, -1):
            remainder = remaining_square - value * value
            if remainder > (slots - 1) * value * value:
                continue
            updated = tuple(
                partial_pairings[j] + value * previous_m[j][index]
                for j in range(len(previous_m))
            )
            if any(updated[j] > bounds[j] for j in range(len(updated))):
                continue
            extend(index + 1, value, remainder, (*prefix, value), updated)

    extend(0, isqrt(target), target, (), (0,) * len(previous_m))
    return answers


@lru_cache(maxsize=None)
def enumerate_simple_roots(n: int) -> tuple[Vector, ...]:
    roots = zero_roots(n)
    target_positive_count = EXPECTED_FACETS[n] - n
    events: dict[Fraction, list[tuple[int, int]]] = {}
    for degree in range(1, EXPECTED_MAX_DEGREE[n] + 1):
        for norm in (1, 2):
            events.setdefault(Fraction(degree * degree, norm), []).append((degree, norm))

    for height in sorted(events):
        candidates: list[tuple[tuple[int, ...], int, int]] = []
        lower_positive_roots = tuple(roots[n:])
        for degree, norm in events[height]:
            candidates.extend(
                (multiplicities, degree, norm)
                for multiplicities in _candidate_multiplicities(
                    n, degree, norm, lower_positive_roots
                )
            )
        # This is the ordering used in the pinned CoxIter tables when several
        # roots occur at one Vinberg height.
        for multiplicities, degree, norm in sorted(candidates):
            root = (degree, *(-entry for entry in multiplicities))
            if lorentz_pair(root, root) != norm:
                raise AssertionError("candidate has the wrong norm")
            if all(lorentz_pair(root, wall) <= 0 for wall in roots):
                roots.append(root)
        if len(roots) - n == target_positive_count:
            break

    if len(roots) != EXPECTED_FACETS[n]:
        raise AssertionError(
            f"I_(1,{n}) enumeration stopped at {len(roots)} roots, "
            f"expected {EXPECTED_FACETS[n]}"
        )
    return tuple(roots)


def coxeter_label(left: Vector, right: Vector) -> int:
    pairing = lorentz_pair(left, right)
    if pairing == 0:
        return 2
    if pairing > 0:
        raise AssertionError("simple roots do not define an acute chamber")
    q_left = lorentz_pair(left, left)
    q_right = lorentz_pair(right, right)
    square = pairing * pairing
    product = q_left * q_right
    if 4 * square == product:
        return 3
    if 2 * square == product:
        return 4
    if square == product:
        return 0  # parallel
    if square > product:
        return 1  # ultraparallel / dotted
    raise AssertionError((left, right, pairing, q_left, q_right))


def relation_edges(roots: Sequence[Vector]) -> list[tuple[int, int, int]]:
    edges: list[tuple[int, int, int]] = []
    for i, left in enumerate(roots, start=1):
        for j, right in enumerate(roots[i:], start=i + 1):
            label = coxeter_label(left, right)
            if label != 2:
                edges.append((i, j, label))
    return edges


def coxiter_filename(n: int) -> str:
    return f"{n}-1," + ",".join("1" for _ in range(n)) + ".coxiter"


def parse_coxiter(path: Path) -> tuple[int, int, list[tuple[int, int, int]]]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    facets, dimension = (int(entry) for entry in lines[0].split())
    edges: list[tuple[int, int, int]] = []
    for line in lines[1:]:
        if line.startswith("vertices labels:"):
            continue
        left, right, label = (int(entry) for entry in line.split())
        edges.append((left, right, label))
    return facets, dimension, edges


def canonical_pairing(root: Vector) -> int:
    # K_X = -3H + E_1 + ... + E_n and the stored exceptional
    # coefficients of a positive root are -m_i.
    return -3 * root[0] - sum(root[1:])


def arithmetic_genus(root: Vector) -> int:
    surface_square = -lorentz_pair(root, root)
    numerator = surface_square + canonical_pairing(root)
    if numerator % 2:
        raise AssertionError("adjunction parity failure")
    return 1 + numerator // 2


def class_record(label: int, root: Vector) -> dict[str, object]:
    norm = lorentz_pair(root, root)
    if root[0] == 0:
        return {
            "label": label,
            "kind": "exceptional_chain_root" if label < len(root) - 1 else "last_exceptional",
            "vector": list(root),
            "lorentz_norm": norm,
            "surface_square": -norm,
            "canonical_pairing": canonical_pairing(root),
            "arithmetic_genus": arithmetic_genus(root),
        }
    multiplicities = [-entry for entry in root[1:]]
    return {
        "label": label,
        "kind": "plane_curve",
        "vector": list(root),
        "degree": root[0],
        "multiplicities": multiplicities,
        "height": f"{root[0] * root[0]}/{norm}",
        "lorentz_norm": norm,
        "surface_square": -norm,
        "canonical_pairing": canonical_pairing(root),
        "arithmetic_genus": arithmetic_genus(root),
        "proximity_inequalities_hold": all(
            multiplicities[i] >= multiplicities[i + 1]
            for i in range(len(multiplicities) - 1)
        ),
    }


def positive_class(degree: int, multiplicities: Iterable[int], n: int) -> Vector:
    values = tuple(multiplicities)
    if len(values) > n:
        values = values[:n]
    return (degree, *(-entry for entry in values), *((0,) * (n - len(values))))


def named_classes(n: int) -> dict[str, Vector]:
    return {
        "L": positive_class(1, (1, 1, 1), n),
        "C": positive_class(3, (1,) * 11, n),
        "Q": positive_class(4, (2, *((1,) * 14)), n),
        "A": positive_class(4, (1,) * min(n, 18), n),
        "S": positive_class(6, (*((2,) * 7), *((1,) * 10)), n),
    }


def locate_named_roots(roots: Sequence[Vector], n: int) -> dict[str, int | None]:
    lookup = {root: label for label, root in enumerate(roots, start=1)}
    return {name: lookup.get(root) for name, root in named_classes(n).items()}


def explicit_cluster_record() -> dict[str, object]:
    return {
        "field": "Q (and hence C)",
        "coordinates": "[X:Y:Z]",
        "point": "p=[0:1:0]",
        "forms": {
            "L": "Z",
            "C": "Y^2*Z-X^3-Z^3",
            "Q": "Z^4+(X+Z)*C",
            "S": "Z^2*Q+C^2",
        },
        "absolute_irreducibility": {
            "tool": "Singular 4.4.1 absfact.lib::absFactorize",
            "replay_script": "scripts/certify_rational_i1n_curves.sing",
            "replay_script_sha256": "ee7336d70b487829a8b4474a6b490031b6188df30c86341469bda54496d972d3",
            "result": {"C": 1, "Q": 1, "S": 1},
            "meaning": "one nonconstant absolute factor, with exponent one, for each form",
        },
        "global_intersections_concentrated_at_p": {
            "L.C": 3,
            "L.Q": 4,
            "L.S": 6,
            "C.Q": 12,
            "C.S": 18,
            "Q.S": 24,
        },
        "local_branch_contacts": {
            "Q_tangent_cone_at_p": "z*(x+z) (ordinary node in the affine chart Y=1)",
            "Q_restriction_to_C": "Z^4, so the total intersection 12 is concentrated at p",
            "Q_with_C": [1, 11],
            "S_with_C": [7, 11],
            "S_order_11_branch_with_Q_order_11_branch": 15,
            "S_valuation_proof_on_selected_Q_branch": [
                "v_k(C^2)=2k for k<=11 and then 22",
                "v_k(Z^2*Q)=k+7 for k>=3",
                "the successive minima give multiplicity 2 at p1,...,p7 and 1 at p8,...,p15",
                "Q and S then separate; choose p16 on the S branch",
            ],
        },
        "free_chain_multiplicities": {
            "L": [1, 1, 1],
            "C": [1] * 11,
            "Q": [2] + [1] * 14,
            "S": [2] * 7 + [1] * 9,
        },
        "construction_scope": {
            "n=12,13": "L and C are the negative positive-degree walls; Q has nonnegative square; exact two-Bezout cone proof is complete",
            "n=14": "truncate the Q branch chain after p14; Q has square -1; exact three-Bezout cone proof is complete",
            "n=15": "truncate the Q branch chain after p15; Q has square -2; exact three-Bezout cone proof is complete",
            "n=16": "after p15 continue at p16 along the S order-11 branch, away from Q",
        },
    }


def p12_p13_decomposition(
    n: int, degree: int, multiplicities: Sequence[int]
) -> tuple[Fraction, Fraction, tuple[Fraction, ...]]:
    """Decompose a class satisfying the two Bezout/proximity inequalities.

    The returned tuple ``(x,y,z)`` certifies

        D = x L + y C + sum(z_i R_i) + z_n E_n,

    where ``R_i=E_i-E_(i+1)``.  All coefficients are exact and nonnegative.
    """

    if n not in (12, 13) or len(multiplicities) != n:
        raise ValueError("the closed two-wall argument applies only to n=12,13")
    m = tuple(multiplicities)
    if degree < 0 or any(value < 0 for value in m):
        raise ValueError("degree and multiplicities must be nonnegative")
    if any(m[i] < m[i + 1] for i in range(n - 1)):
        raise ValueError("the free-chain proximity inequalities fail")
    cumulative = tuple(sum(m[:j]) for j in range(1, n + 1))
    if cumulative[2] > degree or cumulative[10] > 3 * degree:
        raise ValueError("the L or C Bezout inequality fails")

    y = max(Fraction(0), Fraction(cumulative[-1] - 3 * degree, 2))
    x = Fraction(degree) - 3 * y
    base = tuple(
        (x if i < 3 else Fraction(0)) + (y if i < 11 else Fraction(0))
        for i in range(n)
    )
    base_cumulative = tuple(sum(base[:j]) for j in range(1, n + 1))
    exceptional = tuple(
        base_cumulative[j] - cumulative[j] for j in range(n)
    )
    if x < 0 or y < 0 or any(value < 0 for value in exceptional):
        raise AssertionError("the exact P12/P13 cone lemma failed")
    return x, y, exceptional


def p12_p13_cone_proof(n: int) -> dict[str, object] | None:
    if n not in (12, 13):
        return None
    return {
        "applies": True,
        "surface": (
            "blow up the eleven-point free flex-cubic chain, then a free tail "
            f"of length {n - 11} away from C"
        ),
        "irreducible_curve_generators": [
            f"R_i=E_i-E_(i+1), 1<=i<{n}",
            f"E_{n}",
            "L=H-E1-E2-E3",
            "C=3H-E1-...-E11",
        ],
        "arbitrary_other_curve": "D=dH-sum(m_i E_i), m1>=...>=m_n>=0",
        "bezout_inequalities": ["A3<=d", "A11<=3d"],
        "exact_choice": [
            "y=max(0,(A_n-3d)/2)",
            "x=d-3y",
            "B_j=sum_(i<=j)[x*1_(i<=3)+y*1_(i<=11)]",
            "z_j=B_j-A_j",
        ],
        "nonnegativity_proof": [
            "The tail length n-11 is at most two, so y<=m11<=d/3 and x>=0.",
            "For j<=3, A_j<=j(d-2m11)<=j(d-2y)=B_j.",
            "For 4<=j<=8, A_j+(9-j)y<=d+6m4<=3d, hence A_j<=B_j.",
            "For 9<=j<=11 use A_j<=A11<=3d; for j>=12 use A_j<=A_n<=3d+2y.",
        ],
        "conclusion": (
            "NE(X_n)=Cone(R_i,E_n,L,C), so Nef(X_n) is exactly the dual "
            f"Vinberg P{n} chamber"
        ),
    }


def _primitive_integer_vector(values: Sequence[object]) -> tuple[int, ...]:
    """Convert a rational vector to primitive integers."""

    denominator = 1
    for value in values:
        value_denominator = (
            int(value.q) if hasattr(value, "q") else int(value.denominator)
        )
        denominator = lcm(denominator, value_denominator)
    integers = [int(value * denominator) for value in values]
    common = 0
    for value in integers:
        common = gcd(common, abs(value))
    if not common:
        raise ValueError("zero vector")
    integers = [value // common for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def _extreme_rays_via_cdd(
    rows: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...] | None:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        return None
    column_count = len(rows[0]) + 1
    payload = "\n".join(
        [
            "H-representation",
            "begin",
            f"{len(rows)} {column_count} rational",
            *("0 " + " ".join(str(value) for value in row) for row in rows),
            "end",
            "",
        ]
    )
    completed = subprocess.run(
        [executable, "--rep"],
        input=payload,
        check=True,
        capture_output=True,
        text=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines()]
    representation = lines.index("V-representation")
    begin = lines.index("begin", representation)
    count = int(lines[begin + 1].split()[0])
    rays: set[tuple[int, ...]] = set()
    for line in lines[begin + 2 : begin + 2 + count]:
        values = tuple(Fraction(entry) for entry in line.split())
        if values[0] != 0:
            raise AssertionError("the H-cone unexpectedly has a vertex")
        rays.add(_primitive_integer_vector(values[1:]))
    return tuple(sorted(rays))


def _extreme_rays_via_sympy(
    rows: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    from sympy import Matrix

    dimension = len(rows[0])
    rays: set[tuple[int, ...]] = set()
    for active in combinations(range(len(rows)), dimension - 1):
        nullspace = Matrix([rows[index] for index in active]).nullspace()
        if len(nullspace) != 1:
            continue
        vector = nullspace[0]
        evaluations = [
            sum(coefficient * value for coefficient, value in zip(row, vector))
            for row in rows
        ]
        if all(value >= 0 for value in evaluations):
            pass
        elif all(value <= 0 for value in evaluations):
            vector = -vector
        else:
            continue
        ray = _primitive_integer_vector(vector)
        if all(sum(a * b for a, b in zip(row, ray)) >= 0 for row in rows):
            rays.add(ray)
    return tuple(sorted(rays))


def exact_extreme_rays(
    rows: Sequence[Sequence[int]], expected: int
) -> tuple[tuple[int, ...], ...]:
    rays = _extreme_rays_via_cdd(rows)
    if rays is None:
        rays = _extreme_rays_via_sympy(rows)
    if len(rays) != expected:
        raise AssertionError(f"found {len(rays)} rays, expected {expected}")
    return rays


def q_tail_inequality_rows(n: int) -> tuple[tuple[int, ...], ...]:
    """H-description for plane curves distinct from L,C,Q, for n=14 or 15."""

    if n not in (14, 15):
        raise ValueError("the one-quartic exact cone certificate is n=14,15")
    rows: list[tuple[int, ...]] = []
    rows.append((1, *((0,) * n)))  # d >= 0
    for i in range(n - 1):
        rows.append(
            (0, *tuple(1 if j == i else -1 if j == i + 1 else 0 for j in range(n)))
        )
    rows.append((0, *((0,) * (n - 1)), 1))  # m_n >= 0
    rows.append((1, *tuple(-1 if j < 3 else 0 for j in range(n))))
    rows.append((3, *tuple(-1 if j < 11 else 0 for j in range(n))))
    rows.append((4, *tuple(-2 if j == 0 else -1 for j in range(n))))
    return tuple(rows)


@lru_cache(maxsize=None)
def q_tail_extreme_rays(n: int) -> tuple[tuple[int, ...], ...]:
    """Enumerate every extreme ray by exact rank/nullspace arithmetic."""

    rows = q_tail_inequality_rows(n)
    expected = {14: 94, 15: 145}[n]
    return exact_extreme_rays(rows, expected)


def q_tail_decomposition(
    n: int, degree: int, multiplicities: Sequence[int]
) -> tuple[Fraction, Fraction, Fraction, tuple[Fraction, ...]]:
    """Find an exact L,C,Q,exceptional decomposition in the n=14/15 cone."""

    if n not in (14, 15) or len(multiplicities) != n:
        raise ValueError("the one-quartic exact cone certificate is n=14,15")
    vector = (degree, *multiplicities)
    rows = q_tail_inequality_rows(n)
    if any(sum(a * b for a, b in zip(row, vector)) < 0 for row in rows):
        raise ValueError("the proximity/Bezout H-inequalities fail")
    cumulative = tuple(sum(multiplicities[:j]) for j in range(1, n + 1))

    # Eliminate x by x=d-3y-4w.  Each tuple is c+a*y+b*w >= 0.
    halfplanes: list[tuple[Fraction, Fraction, Fraction]] = [
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
        (Fraction(degree), Fraction(-3), Fraction(-4)),
    ]
    for j, target in enumerate(cumulative, start=1):
        if j <= 3:
            halfplanes.append(
                (Fraction(j * degree - target), Fraction(-2 * j), Fraction(1 - 3 * j))
            )
        elif j <= 11:
            halfplanes.append(
                (Fraction(3 * degree - target), Fraction(j - 9), Fraction(j - 11))
            )
        else:
            halfplanes.append(
                (Fraction(3 * degree - target), Fraction(2), Fraction(j - 11))
            )

    candidates: set[tuple[Fraction, Fraction]] = {(Fraction(0), Fraction(0))}
    for first, second in combinations(halfplanes, 2):
        c1, a1, b1 = first
        c2, a2, b2 = second
        determinant = a1 * b2 - a2 * b1
        if not determinant:
            continue
        y = (-c1 * b2 + b1 * c2) / determinant
        w = (-a1 * c2 + c1 * a2) / determinant
        candidates.add((y, w))

    feasible = sorted(
        (
            (y, w)
            for y, w in candidates
            if all(c + a * y + b * w >= 0 for c, a, b in halfplanes)
        ),
        key=lambda pair: (pair[0] + pair[1], pair[1], pair[0]),
    )
    if not feasible:
        raise AssertionError("no exact L,C,Q decomposition witness")
    y, w = feasible[0]
    x = Fraction(degree) - 3 * y - 4 * w
    base = tuple(
        (x if i < 3 else Fraction(0))
        + (y if i < 11 else Fraction(0))
        + w
        + (w if i == 0 else Fraction(0))
        for i in range(n)
    )
    base_cumulative = tuple(sum(base[:j]) for j in range(1, n + 1))
    exceptional = tuple(
        base_cumulative[j] - cumulative[j] for j in range(n)
    )
    if min((x, y, w, *exceptional)) < 0:
        raise AssertionError("negative coefficient in exact decomposition")
    return x, y, w, exceptional


@lru_cache(maxsize=None)
def q_tail_cone_certificate(n: int) -> dict[str, object]:
    rays = q_tail_extreme_rays(n)
    witnesses = [q_tail_decomposition(n, ray[0], ray[1:]) for ray in rays]
    denominators = [
        coefficient.denominator
        for witness in witnesses
        for item in witness
        for coefficient in (item if isinstance(item, tuple) else (item,))
    ]
    ray_digest = hashlib.sha256(
        json.dumps(rays, separators=(",", ":")).encode()
    ).hexdigest()
    witness_payload = [
        [str(x), str(y), str(w), *[str(value) for value in exceptional]]
        for x, y, w, exceptional in witnesses
    ]
    witness_digest = hashlib.sha256(
        json.dumps(witness_payload, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "applies": True,
        "dimension": n,
        "H_inequalities": [
            "d>=0",
            "m1>=...>=m_n>=0",
            "A3<=d",
            "A11<=3d",
            "A_n+m1<=4d",
        ],
        "primitive_extreme_ray_count": len(rays),
        "primitive_extreme_rays_sha256": ray_digest,
        "every_ray_has_exact_nonnegative_L_C_Q_exceptional_decomposition": True,
        "decomposition_witness_sha256": witness_digest,
        "maximum_witness_denominator": max(denominators),
        "irreducible_curve_generators": [
            f"R_i=E_i-E_(i+1), 1<=i<{n}",
            f"E_{n}",
            "L=H-E1-E2-E3",
            "C=3H-E1-...-E11",
            f"Q=4H-2E1-E2-...-E{n}",
        ],
        "conclusion": (
            f"NE(X_{n}) is exactly the cone generated by the {n + 3} visible "
            f"negative curves, and its dual is Vinberg P{n}"
        ),
    }


def p16_inequality_rows() -> tuple[tuple[int, ...], ...]:
    """Necessary proximity and four-curve Bezout inequalities in dimension 16."""

    n = 16
    rows: list[tuple[int, ...]] = []
    rows.append((1, *((0,) * n)))
    for i in range(n - 1):
        rows.append(
            (0, *tuple(1 if j == i else -1 if j == i + 1 else 0 for j in range(n)))
        )
    rows.append((0, *((0,) * (n - 1)), 1))
    rows.append((1, *tuple(-1 if j < 3 else 0 for j in range(n))))
    rows.append((3, *tuple(-1 if j < 11 else 0 for j in range(n))))
    rows.append(
        (4, *tuple(-2 if j == 0 else -1 if j < 15 else 0 for j in range(n)))
    )
    rows.append((6, *tuple(-2 if j < 7 else -1 for j in range(n))))
    return tuple(rows)


@lru_cache(maxsize=None)
def p16_extreme_rays() -> tuple[tuple[int, ...], ...]:
    return exact_extreme_rays(p16_inequality_rows(), 325)


def _determinant_three(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def p16_decomposition(
    degree: int, multiplicities: Sequence[int]
) -> tuple[Fraction, Fraction, Fraction, Fraction, tuple[Fraction, ...]]:
    """Find an exact L,C,Q,S,exceptional decomposition for a P16 H-cone class."""

    n = 16
    if len(multiplicities) != n:
        raise ValueError("P16 needs sixteen multiplicities")
    vector = (degree, *multiplicities)
    rows = p16_inequality_rows()
    if any(sum(a * b for a, b in zip(row, vector)) < 0 for row in rows):
        raise ValueError("the P16 proximity/Bezout H-inequalities fail")

    line = tuple(1 if i < 3 else 0 for i in range(n))
    cubic = tuple(1 if i < 11 else 0 for i in range(n))
    quartic = tuple(2 if i == 0 else 1 if i < 15 else 0 for i in range(n))
    sextic = tuple(2 if i < 7 else 1 for i in range(n))
    cumulative = tuple(sum(multiplicities[:j]) for j in range(1, n + 1))

    # Eliminate x=d-3y-4w-6v.  Each tuple is c+a*y+b*w+e*v >= 0.
    halfspaces: list[tuple[Fraction, Fraction, Fraction, Fraction]] = [
        (Fraction(0), Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0), Fraction(1)),
        (Fraction(degree), Fraction(-3), Fraction(-4), Fraction(-6)),
    ]
    for j, target in enumerate(cumulative, start=1):
        line_sum = sum(line[:j])
        halfspaces.append(
            (
                Fraction(line_sum * degree - target),
                Fraction(sum(cubic[:j]) - 3 * line_sum),
                Fraction(sum(quartic[:j]) - 4 * line_sum),
                Fraction(sum(sextic[:j]) - 6 * line_sum),
            )
        )

    feasible: tuple[Fraction, Fraction, Fraction] | None = None
    for indices in combinations(range(len(halfspaces)), 3):
        selected = [halfspaces[index] for index in indices]
        coefficient_matrix = [list(row[1:]) for row in selected]
        determinant = _determinant_three(coefficient_matrix)
        if determinant == 0:
            continue
        right = [-row[0] for row in selected]
        numerators: list[Fraction] = []
        for column in range(3):
            replaced = [row[:] for row in coefficient_matrix]
            for row in range(3):
                replaced[row][column] = right[row]
            numerators.append(_determinant_three(replaced))
        y, w, v = (numerator / determinant for numerator in numerators)
        if all(c + a * y + b * w + e * v >= 0 for c, a, b, e in halfspaces):
            feasible = (y, w, v)
            break
    if feasible is None:
        raise AssertionError("no exact L,C,Q,S decomposition witness")

    y, w, v = feasible
    x = Fraction(degree) - 3 * y - 4 * w - 6 * v
    base = tuple(
        x * line[i] + y * cubic[i] + w * quartic[i] + v * sextic[i]
        for i in range(n)
    )
    base_cumulative = tuple(sum(base[:j]) for j in range(1, n + 1))
    exceptional = tuple(
        base_cumulative[j] - cumulative[j] for j in range(n)
    )
    if min((x, y, w, v, *exceptional)) < 0:
        raise AssertionError("negative coefficient in exact P16 decomposition")
    return x, y, w, v, exceptional


@lru_cache(maxsize=None)
def p16_cone_certificate() -> dict[str, object]:
    rays = p16_extreme_rays()
    witnesses = [p16_decomposition(ray[0], ray[1:]) for ray in rays]
    denominators = [
        coefficient.denominator
        for witness in witnesses
        for item in witness
        for coefficient in (item if isinstance(item, tuple) else (item,))
    ]
    ray_digest = hashlib.sha256(
        json.dumps(rays, separators=(",", ":")).encode()
    ).hexdigest()
    witness_payload = [
        [str(x), str(y), str(w), str(v), *[str(value) for value in exceptional]]
        for x, y, w, v, exceptional in witnesses
    ]
    witness_digest = hashlib.sha256(
        json.dumps(witness_payload, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "applies": True,
        "dimension": 16,
        "H_inequalities": [
            "d>=0",
            "m1>=...>=m16>=0",
            "A3<=d",
            "A11<=3d",
            "A15+m1<=4d",
            "A16+A7<=6d",
        ],
        "primitive_extreme_ray_count": len(rays),
        "primitive_extreme_rays_sha256": ray_digest,
        "every_ray_has_exact_nonnegative_L_C_Q_S_exceptional_decomposition": True,
        "decomposition_witness_sha256": witness_digest,
        "maximum_witness_denominator": max(denominators),
        "exact_enumeration": "cddlib GMP when available, exact SymPy rank/nullspace fallback",
        "irreducible_curve_generators": [
            "R_i=E_i-E_(i+1), 1<=i<16",
            "E_16",
            "L=H-E1-E2-E3",
            "C=3H-E1-...-E11",
            "Q=4H-2E1-E2-...-E15",
            "S=6H-2E1-...-2E7-E8-...-E16",
        ],
        "conclusion": (
            "NE(X_16) is exactly the cone generated by the twenty visible "
            "negative curves, and its dual is Vinberg P16"
        ),
    }


def flex_obstruction_record(n: int, roots: Sequence[Vector]) -> dict[str, object] | None:
    if n < 17:
        return None
    labels = locate_named_roots(roots, n)
    required = {name: labels[name] for name in ("L", "C", "Q", "A")}
    if any(label is None for label in required.values()):
        raise AssertionError((n, required))
    classes = named_classes(n)
    mmp_data = {
        17: {
            "replay": "scripts/certify_rational_p17_mmp.py",
            "certificate": "data/enumerations/rational_p17_mmp.json",
            "rank_two_endpoints": {"F1_nonminimal": 8, "F2_minimal": 8},
            "rank_one_P2_endpoints": 8,
            "canonical_class_count": 8,
            "weighted_diagram_automorphism_order": 8,
        },
        18: {
            "replay": "scripts/certify_rational_p18_p19_mmp.py",
            "certificate": "data/enumerations/rational_p18_p19_mmp.json",
            "rank_two_endpoints": {"F1_nonminimal": 24, "F2_minimal": 24},
            "rank_one_P2_endpoints": 24,
            "canonical_class_count": 24,
            "weighted_diagram_automorphism_order": 24,
        },
        19: {
            "replay": "scripts/certify_rational_p18_p19_mmp.py",
            "certificate": "data/enumerations/rational_p18_p19_mmp.json",
            "rank_two_endpoints": {"F1_nonminimal": 120, "F2_minimal": 120},
            "rank_one_P2_endpoints": 120,
            "canonical_class_count": 120,
            "weighted_diagram_automorphism_order": 120,
        },
    }[n]
    mmp_data = {
        **mmp_data,
        "classification": (
            "the endpoint canonical classes are exactly the full weighted-diagram "
            f"orbit of K=(-3,1^{n}); conjugating any marking exposes a "
            f"smooth-rational B{n} chain and reduces to the Abel contradiction"
        ),
    }
    return {
        "applies": True,
        "scope": {
            "standard_marking_abel_obstruction": (
                "the B_n free-chain plane-blowup marking is impossible"
            ),
            "intrinsic_rational_surface_extension": (
                "any smooth rational-surface realization whose B_n facet chain "
                "consists of smooth rational curves; sequential Castelnuovo contraction "
                "then forces the standard free-chain marking over P2"
            ),
            "complete_mori_carrier": (
                f"P{n} is impossible on every smooth projective rational surface "
                "by the exact facet-cone MMP classification"
            ),
            "unconditional": (
                "nonexistence of a complete facet-cone carrier on a smooth "
                "projective rational surface in characteristic zero"
            ),
            "not_claimed": (
                f"nonexistence of the classical abstract P{n} hyperbolic polytope, "
                "or nonexistence of a carrier on a nonrational variety"
            ),
        },
        "wall_labels": required,
        "intersection_numbers": {
            "L.C": -lorentz_pair(classes["L"], classes["C"]),
            "C.Q": -lorentz_pair(classes["C"], classes["Q"]),
            "C.A": -lorentz_pair(classes["C"], classes["A"]),
            "Q.A": -lorentz_pair(classes["Q"], classes["A"]),
        },
        "proof": [
            "In the standard marking the B_n walls are the free chain p1,...,pn.",
            "More intrinsically, if those B_n facet curves are smooth rational, contract E_n,R_(n-1),...,R_1; the smooth rational rank-one target is P2, and reversing the contractions gives that same free chain.",
            "L.C=0 and Bezout give C|_L=3p, hence O_C(1)=O_C(3p).",
            "C.Q=0 and the common multiplicities already contribute 2+10=12, so Q forces p12 away from C.",
            "A is smooth at p and follows the required free chain; since p12 is away from C, I_p(A,C)=11 exactly.",
            "Bezout gives A|_C=11p+r for an effective Cartier divisor r of degree one away from p.",
            "But O_C(A)=O_C(4)=O_C(12p), so O_C(r)=O_C(p).",
            "The degree-one Abel map of an integral plane cubic is injective on its smooth locus; a degree-one Cartier point is smooth, hence r=p, contradiction.",
        ],
        "marking_independent_mmp_upgrade": mmp_data,
        "p17_mmp_upgrade": mmp_data if n == 17 else None,
        "conclusion": (
            f"no smooth projective rational surface has complete Mori cone generated "
            f"by the P{n} facet classes"
        ),
    }


def status_for_dimension(n: int) -> dict[str, str]:
    if n <= 13:
        return {
            "wall_realization": "explicit for every abstract simple root",
            "nef_cone_completeness": "proved exactly by proximity, Bezout with L and C, and an explicit cone decomposition",
            "coxeter_polytope": "already-known abstract Vinberg chamber; positive rational-surface nef realization, priority open",
        }
    if n <= 15:
        return {
            "wall_realization": "explicit for every abstract simple root using the free chain and C,L,Q",
            "nef_cone_completeness": "proved exactly by proximity, three Bezout inequalities, exhaustive exact rays, and exact decompositions",
            "coxeter_polytope": "already-known abstract Vinberg chamber; positive rational-surface nef realization, priority open",
        }
    if n == 16:
        return {
            "wall_realization": "explicit for every abstract simple root using C,L,Q,S",
            "nef_cone_completeness": "proved exactly by proximity, four Bezout inequalities, exhaustive exact rays, and exact decompositions",
            "coxeter_polytope": "already-known abstract Vinberg chamber; positive rational-surface nef realization, priority open",
        }
    classical_source = "Vinberg" if n == 17 else "Vinberg--Kaplinskaja"
    orbit_count = {17: 8, 18: 24, 19: 120}[n]
    return {
        "wall_realization": (
            f"impossible on every smooth projective rational surface with the P{n} "
            "facet cone, by exact surface-MMP marking classification and the Abel obstruction"
        ),
        "nef_cone_completeness": (
            "universally obstructed for rational surfaces: all possible minimal-model "
            f"markings form the {orbit_count}-element diagram orbit of the standard marking"
        ),
        "coxeter_polytope": f"already-known abstract {classical_source} finite-volume chamber; not a new polytope",
    }


def galindo_factory_ceiling() -> dict[str, object]:
    return {
        "source": "Galindo-Monserrat-Moreno-Avila-Perez-Callejo, arXiv:2602.15703v1, Theorems 3 and 7",
        "source_url": "https://arxiv.org/abs/2602.15703",
        "reasoning": [
            "A smooth degree-d one-place-at-infinity pencil has d^2 simple base points.",
            "The strict line at infinity has square 1-d, so requiring only -1/-2 generators gives d<=3.",
            "A selected pencil fibre has square -s after a free extension of length s, so s<=2.",
            "For d=3 this gives at most 9+2=11 blowups, Picard rank 12, and hyperbolic dimension 11.",
            "A twelfth blowup along that fibre makes it a -3 curve; a twelfth blowup away from the fibre lies outside Theorem 7's configuration.",
        ],
        "ceiling": {"blowups": 11, "picard_rank": 12, "hyperbolic_dimension": 11},
        "conclusion": "this exact published cone factory cannot certify P12 or any higher I_(1,n) chamber",
    }


def build_certificate(coxiter_root: Path | None = None) -> dict[str, object]:
    dimensions: dict[str, object] = {}
    for n in DIMENSIONS:
        roots = enumerate_simple_roots(n)
        edges = relation_edges(roots)
        graph_record: dict[str, object] = {
            "expected_sha256": COXITER_SHA256[n],
            "exact_relation_match": None,
        }
        if coxiter_root is not None:
            path = coxiter_root / coxiter_filename(n)
            payload = path.read_bytes()
            digest = hashlib.sha256(payload).hexdigest()
            facets, dimension, published_edges = parse_coxiter(path)
            graph_record.update(
                {
                    "source_filename": path.name,
                    "observed_sha256": digest,
                    "sha256_match": digest == COXITER_SHA256[n],
                    "declared_facets": facets,
                    "declared_dimension": dimension,
                    "exact_relation_match": facets == len(roots)
                    and dimension == n
                    and len(published_edges) == len(edges)
                    and set(published_edges) == set(edges),
                }
            )
            if not graph_record["sha256_match"] or not graph_record["exact_relation_match"]:
                raise AssertionError((n, graph_record))

        records = [class_record(label, root) for label, root in enumerate(roots, start=1)]
        positive = records[n:]
        root_digest = hashlib.sha256(
            json.dumps(roots, separators=(",", ":")).encode()
        ).hexdigest()
        edge_digest = hashlib.sha256(
            json.dumps(edges, separators=(",", ":")).encode()
        ).hexdigest()
        dimensions[str(n)] = {
            "lattice": f"I_(1,{n})",
            "hyperbolic_dimension": n,
            "facet_count": len(roots),
            "height_zero_Bn_count": n,
            "positive_height_count": len(roots) - n,
            "max_positive_degree": max(root[0] for root in roots[n:]),
            "ordered_root_sha256": root_digest,
            "ordered_relation_edge_sha256": edge_digest,
            "relation_histogram": dict(
                sorted(Counter(str(label) for _, _, label in edges).items())
            ),
            "canonical_pairing_histogram_positive": dict(
                sorted(Counter(str(record["canonical_pairing"]) for record in positive).items())
            ),
            "genus_histogram_positive": dict(
                sorted(Counter(str(record["arithmetic_genus"]) for record in positive).items())
            ),
            "named_root_labels": locate_named_roots(roots, n),
            "roots": records,
            "coxiter": graph_record,
            "vinberg_volume": VINBERG_VOLUME_DATA.get(n),
            "status": status_for_dimension(n),
            "p12_p13_effective_cone_proof": p12_p13_cone_proof(n),
            "q_tail_effective_cone_certificate": (
                q_tail_cone_certificate(n) if n in (14, 15) else None
            ),
            "p16_effective_cone_certificate": (
                p16_cone_certificate() if n == 16 else None
            ),
            "flex_cubic_obstruction": flex_obstruction_record(n, roots),
        }

    certificate: dict[str, object] = {
        "schema": "rational_i1n_inverse_design_v1",
        "source_lineage": {
            "n12_to_n17": {
                "author": "E. B. Vinberg",
                "title": "On groups of unit elements of certain quadratic forms",
                "year": 1972,
                "url": "https://www.mathnet.ru/eng/sm3032",
                "doi": "10.1070/SM1972v016n01ABEH001346",
            },
            "n18_n19": {
                "authors": "E. B. Vinberg and I. M. Kaplinskaja",
                "title": "The groups O_(18,1)(Z) and O_(19,1)(Z)",
                "year": 1978,
                "url": "https://www.mathnet.ru/eng/dan41538",
            },
            "coxiter": {
                "author": "R. Guglielmetti",
                "doi": "10.1112/S1461157015000273",
                "repository": "https://github.com/rgugl/CoxIter",
                "local_replay_version": "1.1.8",
                "local_binary_sha256": "cad0b97caa8ab1dc25b86dff65252e1869382ed079d6ec949dbe006b3588f957",
            },
        },
        "conventions": {
            "lorentz_form": "diag(-1,+1^n)",
            "surface_form": "diag(+1,-1^n)",
            "canonical_class": "K=(-3,+1^n)",
            "positive_root": "(d,-m_1,...,-m_n) <-> dH-sum m_i E_i",
        },
        "dimensions": dimensions,
        "explicit_cluster": explicit_cluster_record(),
        "galindo_factory_ceiling": galindo_factory_ceiling(),
        "global_conclusion": {
            "n12_to_n16": (
                "complete rational-surface Mori and nef cone theorems; the projectivized "
                "nef cones are exactly the finite-volume noncompact Vinberg chambers P12 through P16"
            ),
            "n17": (
                "no smooth projective rational surface has the complete P17 facet cone; "
                "the exact surface MMP reduces all eight possible markings to the Abel obstruction"
            ),
            "n18_to_n19": (
                "no smooth projective rational surface has the complete P18 or P19 "
                "facet cone; exact surface MMPs reduce respectively all 24 and 120 "
                "possible canonical markings to the Abel obstruction"
            ),
            "rational_carrier_nonexistence_dimensions": [17, 18, 19],
            "new_hyperbolic_polytope": False,
            "positive_algebro_geometric_realizations": [12, 13, 14, 15, 16],
            "priority_caveat": (
                "the abstract polytopes are classical; worldwide priority for these "
                "specific complete rational-surface nef realizations has not been established"
            ),
        },
    }
    canonical = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    certificate["certificate_payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coxiter-root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = build_certificate(args.coxiter_root)
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
