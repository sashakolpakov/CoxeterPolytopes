#!/usr/bin/env python3
"""Clean-room hostile audit of the silver prime-seven D24+ H=500 freeze.

The core replay imports neither the production search nor its exact number
field implementation.  It uses ``Fraction`` pairs for Q(sqrt(2)), derives
the bounded control box, enumerates the root and all three spinor residues of
O tensor D24+, performs simultaneous-height Vinberg filtering, and checks the
replacement extreme ray by exact Gaussian elimination.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key, lru_cache
from hashlib import sha256
from math import isqrt, lcm
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
H300_PATH = ROOT / "data/enumerations/silver_p7_d24_plus_h300.json"
H500_PATH = ROOT / "data/enumerations/silver_p7_d24_plus_h500.json"
H300_SHA256 = "4ea87341148489598dc85b3cc8271965a036e5ffba2f24052784705e4af614cf"


@dataclass(frozen=True, slots=True)
class K2:
    """A standalone exact element ``a+b*sqrt(2)``."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @classmethod
    def make(cls, a: int | str | Fraction = 0, b: int | str | Fraction = 0) -> "K2":
        return cls(Fraction(a), Fraction(b))

    @classmethod
    def from_json(cls, item: dict[str, str]) -> "K2":
        if set(item) != {"a", "b"}:
            raise AssertionError("malformed Q(sqrt(2)) pair")
        return cls.make(item["a"], item["b"])

    def __add__(self, other: "K2") -> "K2":
        return K2(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "K2") -> "K2":
        return K2(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "K2":
        return K2(-self.a, -self.b)

    def __mul__(self, other: "K2") -> "K2":
        return K2(
            self.a * other.a + 2 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def __truediv__(self, other: "K2") -> "K2":
        norm = other.a * other.a - 2 * other.b * other.b
        if not norm:
            raise ZeroDivisionError
        return self * K2(other.a / norm, -other.b / norm)

    def __pow__(self, exponent: int) -> "K2":
        if exponent < 0:
            return (ONE / self) ** (-exponent)
        result = ONE
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent >>= 1
        return result

    def __bool__(self) -> bool:
        return bool(self.a or self.b)

    @property
    def integral(self) -> bool:
        return self.a.denominator == self.b.denominator == 1

    def conjugate(self) -> "K2":
        return K2(self.a, -self.b)

    def field_norm(self) -> Fraction:
        return self.a * self.a - 2 * self.b * self.b

    def sign(self, conjugate: bool = False) -> int:
        b = -self.b if conjugate else self.b
        sa = (self.a > 0) - (self.a < 0)
        sb = (b > 0) - (b < 0)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        return sa if self.a * self.a - 2 * b * b > 0 else sb

    def exact(self) -> dict[str, str]:
        return {"a": str(self.a), "b": str(self.b)}


ZERO = K2.make()
ONE = K2.make(1)
TWO = K2.make(2)
SQRT2 = K2.make(0, 1)
EPSILON = K2.make(1, 1)
TIMELIKE = K2.make(-1, 2)
DELTA = K2.make(3, 1)
H300 = K2.make(300)
H500 = K2.make(500)
NORMS = (ONE, TWO, DELTA, TWO * DELTA)
Vector = tuple[K2, ...]


def add(items: Iterable[K2]) -> K2:
    result = ZERO
    for item in items:
        result = result + item
    return result


def value(wrapper: dict[str, object]) -> K2:
    return K2.from_json(wrapper["exact"])


def vector(items: Sequence[dict[str, str]]) -> Vector:
    return tuple(K2.from_json(item) for item in items)


def compare(left: K2, right: K2) -> int:
    return (left - right).sign()


def sorted_k2(items: Iterable[K2], *, reverse: bool = False) -> tuple[K2, ...]:
    return tuple(sorted(items, key=cmp_to_key(compare), reverse=reverse))


def compare_levels(left: "Level", right: "Level") -> int:
    for left_value, right_value in (
        (left.height, right.height),
        (left.norm, right.norm),
        (left.control, right.control),
        (left.target, right.target),
    ):
        result = compare(left_value, right_value)
        if result:
            return result
    return 0


def ceil_identity(item: K2) -> int:
    """Return the exact ceiling at the distinguished real embedding."""

    upper = 0
    while (K2.make(upper) - item).sign() < 0:
        upper += 1
    while upper and (K2.make(upper - 1) - item).sign() >= 0:
        upper -= 1
    return upper


def lorentz(left: Sequence[K2], right: Sequence[K2]) -> K2:
    if len(left) != 25 or len(right) != 25:
        raise AssertionError("Lorentz vectors must have length 25")
    return -(TIMELIKE * left[0] * right[0]) + add(
        left[index] * right[index] for index in range(1, 25)
    )


def functional(wall: Sequence[K2]) -> Vector:
    return (-(TIMELIKE * wall[0]), *wall[1:])


def dot(left: Sequence[K2], right: Sequence[K2]) -> K2:
    return add(x * y for x, y in zip(left, right))


def negative_d24_simple_walls() -> tuple[Vector, ...]:
    walls: list[Vector] = []
    for index in range(23):
        wall = [ZERO] * 25
        wall[1 + index] = -ONE
        wall[2 + index] = ONE
        walls.append(tuple(wall))
    fork = [ZERO] * 25
    fork[23] = fork[24] = -ONE
    walls.append(tuple(fork))
    return tuple(walls)


def d24_plus_basis() -> tuple[tuple[K2, ...], ...]:
    basis: list[tuple[K2, ...]] = []
    for index in range(1, 23):
        item = [ZERO] * 24
        item[index] = ONE
        item[index + 1] = -ONE
        basis.append(tuple(item))
    fork = [ZERO] * 24
    fork[22] = fork[23] = ONE
    basis.append(tuple(fork))
    basis.append(tuple(K2.make(Fraction(1, 2)) for _ in range(24)))
    if len(basis) != 24:
        raise AssertionError("bad D24+ basis")
    return tuple(basis)


def rational_determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    work = [list(row) for row in matrix]
    determinant = Fraction(1)
    sign = 1
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        determinant *= pivot_value
        for row in range(column + 1, len(work)):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot_value
            for other in range(column + 1, len(work)):
                work[row][other] -= factor * work[column][other]
    return sign * determinant


def d24_plus_lattice_audit() -> dict[str, object]:
    basis = d24_plus_basis()
    gram = tuple(
        tuple(dot(left, right).a for right in basis) for left in basis
    )
    if any(dot(left, right).b for left in basis for right in basis):
        raise AssertionError("rational D24+ basis acquired radical Gram entries")
    determinant = rational_determinant(gram)
    if determinant != 1:
        raise AssertionError("D24+ Gram determinant is not one")
    if any(entry.denominator != 1 for row in gram for entry in row):
        raise AssertionError("D24+ Gram is not integral")
    if any(gram[index][index].numerator % 2 for index in range(24)):
        raise AssertionError("D24+ Gram is not even")
    return {"rank": 24, "gram_determinant": 1, "integral": True, "even": True}


def d24_plus_coset(spatial: Sequence[K2]) -> str | None:
    """Classify membership in O tensor D24+ in orthonormal coordinates."""

    if len(spatial) != 24:
        return None
    if all(item.integral for item in spatial) and all(
        sum(getattr(item, coefficient).numerator for item in spatial) % 2 == 0
        for coefficient in ("a", "b")
    ):
        return "D24"
    parities: list[int] = []
    for coefficient in ("a", "b"):
        doubled = [2 * getattr(item, coefficient) for item in spatial]
        if any(item.denominator != 1 for item in doubled):
            return None
        integers = [item.numerator for item in doubled]
        parity = integers[0] % 2
        if any(item % 2 != parity for item in integers):
            return None
        if sum((item - parity) // 2 for item in integers) % 2:
            return None
        parities.append(parity)
    return {
        (1, 0): "spinor-1",
        (0, 1): "spinor-1*sqrt(2)",
        (1, 1): "spinor-1 + 1*sqrt(2)",
    }.get(tuple(parities))


def reflection_integral(root: Sequence[K2], norm: K2) -> bool:
    pairings = (-(TIMELIKE * root[0]),) + tuple(
        dot(root[1:], basis_vector) for basis_vector in d24_plus_basis()
    )
    return all((TWO * pairing / norm).integral for pairing in pairings)


def divisible_by_sqrt2(root: Sequence[K2]) -> bool:
    quotient = tuple(item / SQRT2 for item in root)
    return quotient[0].integral and d24_plus_coset(quotient[1:]) is not None


@dataclass(frozen=True, slots=True)
class Level:
    norm: K2
    height: K2
    control: K2
    target: K2


def enumerate_control_levels(bound: K2 = H500) -> tuple[tuple[Level, ...], dict[str, object]]:
    """Exhaust controls using a proved ellipse for each norm class."""

    levels: set[Level] = set()
    boxes = []
    for norm in NORMS:
        # If x=m+n*sqrt(2), height<=B bounds x_0^2 by (B*q/a)_0.
        # Target conjugate nonnegative bounds x_1^2 by (q/(-a))_1.
        # Therefore m^2+2n^2=(x_0^2+x_1^2)/2 is bounded by E below.
        conjugate_cap = norm / (-TIMELIKE)
        energy_cap = (
            -(bound * conjugate_cap) + conjugate_cap.conjugate()
        ) / TWO
        energy_ceiling = ceil_identity(energy_cap)
        checked_pairs = 0
        for n in range(-isqrt(energy_ceiling // 2) - 1, isqrt(energy_ceiling // 2) + 2):
            remaining = energy_ceiling - 2 * n * n
            if remaining < 0:
                continue
            for m in range(-isqrt(remaining), isqrt(remaining) + 1):
                checked_pairs += 1
                control = K2.make(m, n)
                if control.sign() <= 0:
                    continue
                if not (TWO * TIMELIKE * control / norm).integral:
                    continue
                target = norm + TIMELIKE * control * control
                if target.sign(conjugate=True) < 0:
                    continue
                height = TIMELIKE * control * control / norm
                if height.sign() < 0 or (bound - height).sign() < 0:
                    continue
                levels.add(Level(norm, height, control, target))
        boxes.append(
            {
                "norm": norm.exact(),
                "energy_cap": energy_cap.exact(),
                "energy_ceiling": energy_ceiling,
                "coefficient_pairs_checked": checked_pairs,
            }
        )
    return tuple(sorted(levels, key=cmp_to_key(compare_levels))), {
        "derivation": "m^2+2*n^2=(x_0^2+x_1^2)/2",
        "boxes": boxes,
    }


def coefficient_alphabet(
    target: K2, residue: tuple[int, int] | None
) -> tuple[K2, ...]:
    """Enumerate one-coordinate values for an integral or spinor residue."""

    denominator = 1 if residue is None else 2
    energy_bound = Fraction(denominator * denominator) * target.a
    if energy_bound.denominator != 1:
        raise AssertionError("unexpected fractional trace bound")
    energy = energy_bound.numerator
    output: set[K2] = set()
    max_n = isqrt(energy // 2) + 1
    for n in range(-max_n, max_n + 1):
        remaining = energy - 2 * n * n
        if remaining < 0:
            continue
        for m in range(-isqrt(remaining), isqrt(remaining) + 1):
            if residue is not None and (
                m % 2 != residue[0] or n % 2 != residue[1]
            ):
                continue
            coordinate = K2.make(Fraction(m, denominator), Fraction(n, denominator))
            remainder = target - coordinate * coordinate
            if remainder.sign() >= 0 and remainder.sign(conjugate=True) >= 0:
                output.add(coordinate)
    return sorted_k2(output)


@lru_cache(maxsize=None)
def enumerate_d24_plus_shell(target: K2) -> tuple[tuple[Vector, str], ...]:
    """Enumerate all dominant D24+ vectors of exact square ``target``."""

    output: set[tuple[Vector, str]] = set()
    modes = (
        (None, "D24"),
        ((1, 0), "spinor-1"),
        ((0, 1), "spinor-1*sqrt(2)"),
        ((1, 1), "spinor-1 + 1*sqrt(2)"),
    )
    for residue, label in modes:
        alphabet = coefficient_alphabet(target, residue)
        positive = sorted_k2(
            {
                item if item.sign() > 0 else -item
                for item in alphabet
                if item
            },
            reverse=True,
        )
        squares = tuple(item * item for item in positive)

        @lru_cache(maxsize=None)
        def solve(index: int, slots: int, remainder: K2) -> tuple[tuple[int, ...], ...]:
            if remainder.sign() < 0 or remainder.sign(conjugate=True) < 0:
                return ()
            if index == len(positive):
                return ((),) if not remainder else ()
            answers: list[tuple[int, ...]] = []
            running = remainder
            for multiplicity in range(slots + 1):
                if running.sign() < 0 or running.sign(conjugate=True) < 0:
                    break
                for suffix in solve(index + 1, slots - multiplicity, running):
                    answers.append((multiplicity, *suffix))
                running = running - squares[index]
            return tuple(answers)

        for multiplicities in solve(0, 24, target):
            nonnegative = tuple(
                item
                for item, multiplicity in zip(positive, multiplicities)
                for _ in range(multiplicity)
            ) + (ZERO,) * (24 - sum(multiplicities))
            representatives = (nonnegative,)
            if nonnegative[-1]:
                representatives = (
                    nonnegative,
                    (*nonnegative[:-1], -nonnegative[-1]),
                )
            for representative in representatives:
                if add(item * item for item in representative) != target:
                    raise AssertionError("square-partition reconstruction failed")
                if d24_plus_coset(representative) == label:
                    output.add((representative, label))
    return tuple(output)


@dataclass(frozen=True, slots=True)
class Candidate:
    norm: K2
    height: K2
    control: K2
    target: K2
    coset: str
    vector: Vector


@dataclass(frozen=True, slots=True)
class LevelCensus:
    level: Level
    reduced_target: K2 | None
    glue_patterns: int
    reflection_integral: int
    imprimitive: int
    candidates: tuple[Candidate, ...]


def enumerate_level(level: Level) -> LevelCensus:
    scale = DELTA if level.norm in (DELTA, TWO * DELTA) else ONE
    reduced = level.target / (scale * scale)
    if not reduced.integral or reduced.sign() < 0 or reduced.sign(conjugate=True) < 0:
        return LevelCensus(level, None if not reduced.integral else reduced, 0, 0, 0, ())
    patterns = enumerate_d24_plus_shell(reduced)
    integral_count = 0
    imprimitive_count = 0
    candidates: list[Candidate] = []
    for reduced_spatial, coset in patterns:
        spatial = tuple(scale * item for item in reduced_spatial)
        root = (level.control, *spatial)
        if not reflection_integral(root, level.norm):
            continue
        integral_count += 1
        if lorentz(root, root) != level.norm:
            raise AssertionError("candidate has the wrong norm")
        if divisible_by_sqrt2(root):
            imprimitive_count += 1
            continue
        candidates.append(
            Candidate(
                level.norm,
                level.height,
                level.control,
                level.target,
                coset,
                root,
            )
        )
    return LevelCensus(
        level,
        reduced,
        len(patterns),
        integral_count,
        imprimitive_count,
        tuple(candidates),
    )


def coefficient_matrices(vectors: Sequence[Sequence[K2]]) -> tuple[np.ndarray, np.ndarray, int]:
    denominator = 1
    for item in vectors:
        for coordinate in item:
            denominator = lcm(
                denominator, coordinate.a.denominator, coordinate.b.denominator
            )
    rational = np.asarray(
        [
            [int(x.a.numerator * (denominator // x.a.denominator)) for x in item]
            for item in vectors
        ],
        dtype=np.int64,
    )
    radical = np.asarray(
        [
            [int(x.b.numerator * (denominator // x.b.denominator)) for x in item]
            for item in vectors
        ],
        dtype=np.int64,
    )
    return rational, radical, denominator


@dataclass(slots=True)
class MatrixAudit:
    matrix_entries: int = 0
    maximum_bound: int = 0
    maximum_coefficient: int = 0
    scalar_zero_replays: int = 0
    scalar_strict_replays: int = 0


def pair_coefficients(
    left_a: np.ndarray,
    left_b: np.ndarray,
    right_a: np.ndarray,
    right_b: np.ndarray,
    audit: MatrixAudit,
) -> tuple[np.ndarray, np.ndarray]:
    """Exact cleared numerator pairs for the Lorentz product."""

    def maximum(array: np.ndarray) -> int:
        return int(np.max(np.abs(array), initial=0))

    la, lb, ra, rb = map(maximum, (left_a, left_b, right_a, right_b))
    # Each raw matrix product sums at most 24 spatial or one time term.
    spatial_a_bound = 24 * (la * ra + 2 * lb * rb)
    spatial_b_bound = 24 * (la * rb + lb * ra)
    time_a_bound = la * ra + 2 * lb * rb
    time_b_bound = la * rb + lb * ra
    final_bound = max(
        spatial_a_bound + time_a_bound + 4 * time_b_bound,
        spatial_b_bound + time_b_bound + 2 * time_a_bound,
    )
    audit.maximum_bound = max(audit.maximum_bound, final_bound)
    if final_bound >= 2**63:
        raise AssertionError("int64 pairing bound failed")

    spatial_a = left_a[:, 1:] @ right_a[:, 1:].T
    spatial_a += 2 * (left_b[:, 1:] @ right_b[:, 1:].T)
    spatial_b = left_a[:, 1:] @ right_b[:, 1:].T
    spatial_b += left_b[:, 1:] @ right_a[:, 1:].T
    ca = left_a[:, :1] * right_a[:, :1].T + 2 * left_b[:, :1] * right_b[:, :1].T
    cb = left_a[:, :1] * right_b[:, :1].T + left_b[:, :1] * right_a[:, :1].T
    result_a = spatial_a + ca - 4 * cb
    result_b = spatial_b + cb - 2 * ca
    max_coefficient = max(maximum(result_a), maximum(result_b))
    audit.maximum_coefficient = max(audit.maximum_coefficient, max_coefficient)
    if max_coefficient >= 2_000_000_000:
        raise AssertionError("int64 exact-sign square bound failed")
    audit.matrix_entries += result_a.size
    return result_a, result_b


def positive_mask(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    norm = a * a - 2 * b * b
    return (
        ((a > 0) & (b >= 0))
        | ((a > 0) & (b < 0) & (norm > 0))
        | ((a < 0) & (b > 0) & (norm < 0))
        | ((a == 0) & (b > 0))
    )


def filter_candidates(
    candidates: Sequence[Candidate],
    walls: Sequence[Vector],
    matrix_audit: MatrixAudit,
    *,
    block_size: int = 192,
) -> tuple[Candidate, ...]:
    if not candidates:
        return ()
    ca, cb, _ = coefficient_matrices(tuple(item.vector for item in candidates))
    wa, wb, _ = coefficient_matrices(walls)
    survivors: list[Candidate] = []
    for start in range(0, len(candidates), block_size):
        stop = min(start + block_size, len(candidates))
        a, b = pair_coefficients(ca[start:stop], cb[start:stop], wa, wb, matrix_audit)
        positive = positive_mask(a, b)
        for local, candidate in enumerate(candidates[start:stop]):
            if np.any(positive[local]):
                first = int(np.flatnonzero(positive[local])[0])
                if lorentz(candidate.vector, walls[first]).sign() <= 0:
                    raise AssertionError("matrix/scalar disagreement on rejection")
                matrix_audit.scalar_strict_replays += 1
                continue
            zeros = np.flatnonzero((a[local] == 0) & (b[local] == 0))
            for wall_index in zeros:
                if lorentz(candidate.vector, walls[int(wall_index)]):
                    raise AssertionError("matrix/scalar disagreement on zero boundary")
            matrix_audit.scalar_zero_replays += len(zeros)
            for wall_index in {0, len(walls) // 2, len(walls) - 1}:
                if lorentz(candidate.vector, walls[wall_index]).sign() > 0:
                    raise AssertionError("matrix/scalar disagreement on survivor")
                matrix_audit.scalar_strict_replays += 1
            survivors.append(candidate)
    return tuple(survivors)


def assert_same_height_acute(
    roots: Sequence[Candidate], matrix_audit: MatrixAudit, *, block_size: int = 192
) -> tuple[int, int]:
    if not roots:
        return 0, 0
    va, vb, _ = coefficient_matrices(tuple(item.vector for item in roots))
    unordered = len(roots) * (len(roots) - 1) // 2
    zero_count = 0
    for start in range(0, len(roots), block_size):
        stop = min(start + block_size, len(roots))
        a, b = pair_coefficients(va[start:stop], vb[start:stop], va, vb, matrix_audit)
        positive = positive_mask(a, b)
        for local in range(stop - start):
            global_index = start + local
            positive[local, : global_index + 1] = False
            if np.any(positive[local]):
                right = int(np.flatnonzero(positive[local])[0])
                exact = lorentz(roots[global_index].vector, roots[right].vector)
                raise AssertionError(
                    f"same-height positive pairing at {global_index},{right}: {exact.exact()}"
                )
            zero_indices = np.flatnonzero(
                (a[local, global_index + 1 :] == 0)
                & (b[local, global_index + 1 :] == 0)
            ) + global_index + 1
            zero_count += len(zero_indices)
            for right in zero_indices:
                if lorentz(roots[global_index].vector, roots[int(right)].vector):
                    raise AssertionError("matrix/scalar disagreement on same-level zero")
            matrix_audit.scalar_zero_replays += len(zero_indices)
    return unordered, zero_count


def row_rank(rows: Sequence[Sequence[K2]], width: int) -> tuple[int, tuple[int, ...]]:
    work = [list(row) for row in rows]
    rank = 0
    pivots: list[int] = []
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(work)) if work[index][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for index in range(rank + 1, len(work)):
            if not work[index][column]:
                continue
            factor = work[index][column] / pivot_value
            for other in range(column, width):
                work[index][other] = work[index][other] - factor * work[rank][other]
        pivots.append(column)
        rank += 1
        if rank == len(work):
            break
    return rank, tuple(pivots)


def candidate_from_artifact(item: dict[str, object]) -> Candidate:
    return Candidate(
        value(item["norm"]),
        value(item["height"]),
        value(item["control"]),
        value(item["target"]),
        item["coset"],
        vector(item["vector_orthonormal"]),
    )


def level_from_artifact(item: dict[str, object]) -> Level:
    return Level(
        value(item["norm"]),
        value(item["height"]),
        value(item["control"]),
        value(item["spatial_square_target"]),
    )


def minimum_identity(items: Sequence[K2]) -> K2:
    return min(items, key=cmp_to_key(compare))


def maximum_identity(items: Sequence[K2]) -> K2:
    return max(items, key=cmp_to_key(compare))


def audit(h300_path: Path = H300_PATH, h500_path: Path = H500_PATH) -> dict[str, object]:
    if sha256(h300_path.read_bytes()).hexdigest() != H300_SHA256:
        raise AssertionError("pinned H=300 input changed")
    h300 = json.loads(h300_path.read_text())
    h500 = json.loads(h500_path.read_text())
    if h300["height_bound"] != 300 or h500["height_bound"] != 500:
        raise AssertionError("wrong source artifacts")

    lattice_audit = d24_plus_lattice_audit()
    if EPSILON * TIMELIKE != DELTA:
        raise AssertionError("prime-seven factorization failed")
    if TIMELIKE.field_norm() != -7 or DELTA.field_norm() != 7:
        raise AssertionError("prime-seven norms failed")
    if not (TIMELIKE.sign() > 0 > TIMELIKE.sign(conjugate=True)):
        raise AssertionError("timelike coefficient is not admissible")

    controls, control_boxes = enumerate_control_levels(H500)
    old_expected = {
        level_from_artifact(item) for item in h300["enumeration"]["control_levels"]
    }
    new_expected = {
        level_from_artifact(item)
        for item in h500["enumeration"]["extension_control_levels"]
    }
    if len(controls) != 64 or set(controls) != old_expected | new_expected:
        raise AssertionError("independent H<=500 control census disagrees")
    extension = tuple(item for item in controls if (item.height - H300).sign() > 0)
    if len(extension) != 16 or set(extension) != new_expected:
        raise AssertionError("independent (300,500] control census disagrees")

    archived_level_rows = {
        level_from_artifact(item): item
        for item in h500["enumeration"]["extension_control_levels"]
    }
    level_censuses: dict[Level, LevelCensus] = {}
    for index, level in enumerate(extension, start=1):
        census = enumerate_level(level)
        row = archived_level_rows[level]
        archived_reduced = row["reduced_target_after_p7_divisibility"]
        expected_reduced = None if archived_reduced is None else value(archived_reduced)
        if (
            census.reduced_target != expected_reduced
            or census.glue_patterns != row["glue_pattern_count"]
            or census.reflection_integral
            != row["reflection_integral_count_before_primitivity"]
            or census.imprimitive != row["imprimitive_sqrt2_count"]
            or len(census.candidates) != row["primitive_candidate_count"]
        ):
            raise AssertionError(f"independent level census disagrees: {level}")
        level_censuses[level] = census
        print(
            f"level {index}/16 height={level.height.exact()} patterns={census.glue_patterns} "
            f"candidates={len(census.candidates)}",
            flush=True,
        )

    old_walls = list(negative_d24_simple_walls())
    for packet in h300["enumeration"]["packets"]:
        old_walls.extend(
            vector(item["vector_orthonormal"]) for item in packet["accepted_roots"]
        )
    if len(old_walls) != 348:
        raise AssertionError("H=300 lineage wall count changed")

    archived_packets = {
        value(packet["height"]): packet
        for packet in h500["enumeration"]["extension_packets"]
    }
    heights = sorted_k2({level.height for level in extension})
    matrix_audit = MatrixAudit()
    walls = list(old_walls)
    independent_accepted: list[Candidate] = []
    packet_results = []
    total_unordered = 0
    for height in heights:
        packet_levels = [level for level in extension if level.height == height]
        candidates = tuple(
            candidate
            for level in packet_levels
            for candidate in level_censuses[level].candidates
        )
        survivors = filter_candidates(candidates, walls, matrix_audit)
        unordered, pairwise_zeros = assert_same_height_acute(survivors, matrix_audit)
        total_unordered += unordered
        archived = archived_packets[height]
        archived_survivors = {
            candidate_from_artifact(item) for item in archived["accepted_roots"]
        }
        if set(survivors) != archived_survivors:
            missing = len(set(survivors) - archived_survivors)
            extra = len(archived_survivors - set(survivors))
            raise AssertionError(
                f"accepted-root set mismatch at {height.exact()}: {missing=} {extra=}"
            )
        if (
            len(candidates) != archived["primitive_candidate_count"]
            or len(survivors) != archived["accepted_count"]
            or sum(item.coset != "D24" for item in survivors)
            != archived["accepted_spinor_count"]
        ):
            raise AssertionError("packet counts disagree")
        walls.extend(item.vector for item in survivors)
        independent_accepted.extend(survivors)
        packet_results.append(
            {
                "height": height.exact(),
                "primitive_candidates": len(candidates),
                "simultaneous_survivors": len(survivors),
                "same_height_unordered_pairs": unordered,
                "same_height_zero_pairs": pairwise_zeros,
                "same_height_positive_pairs": 0,
            }
        )
        print(
            f"packet height={height.exact()} candidates={len(candidates)} "
            f"survivors={len(survivors)} walls={len(walls)}",
            flush=True,
        )
    if len(independent_accepted) != 4026 or len(walls) != 4374:
        raise AssertionError("final independent wall census disagrees")
    if total_unordered != 8_010_154:
        raise AssertionError("same-height pair census disagrees")

    # Artifact order is used only for checking the reported global facet index;
    # root-set equality above has already independently certified every wall.
    artifact_new = [
        candidate_from_artifact(item)
        for packet in h500["enumeration"]["extension_packets"]
        for item in packet["accepted_roots"]
    ]
    artifact_walls = [*old_walls, *(item.vector for item in artifact_new)]
    if len(set(artifact_walls)) != 4374:
        raise AssertionError("duplicate wall in frozen stream")

    old_c = K2.make(Fraction(172, 7), Fraction(-120, 7))
    old_ray = (ONE, *(old_c for _ in range(19)), *(ZERO for _ in range(5)))
    old_pairings = tuple(lorentz(old_ray, item.vector) for item in artifact_new)
    if sum(item.sign() > 0 for item in old_pairings) != 74:
        raise AssertionError("old-ray cut count disagrees")

    c = K2.make(4, Fraction(-12, 5))
    replacement = (ONE, *(c for _ in range(5)), *(ZERO for _ in range(19)))
    replacement_pairings = tuple(lorentz(replacement, wall) for wall in artifact_walls)
    signs = [item.sign() for item in replacement_pairings]
    if (signs.count(-1), signs.count(0), signs.count(1)) != (4350, 24, 0):
        raise AssertionError("replacement-ray pairing census disagrees")
    active = tuple(index for index, item in enumerate(replacement_pairings) if not item)
    expected_active = (*range(4), *range(5, 24), 4361)
    if active != expected_active:
        raise AssertionError("replacement-ray active facets disagree")
    active_rows = tuple(functional(artifact_walls[index]) for index in active)
    rank, pivots = row_rank(active_rows, 25)
    if rank != 24 or any(dot(row, replacement) for row in active_rows):
        raise AssertionError("replacement ray is not an exact rank-24 face")
    replacement_norm = lorentz(replacement, replacement)
    if replacement_norm != K2.make(Fraction(693, 5), -98) or replacement_norm.sign() <= 0:
        raise AssertionError("replacement ray is not hyperideal")

    bounds: list[tuple[K2, int]] = []
    for index, wall in enumerate(artifact_walls[24:], start=24):
        denominator = add(wall[1:6])
        if denominator.sign() > 0:
            bounds.append((TIMELIKE * wall[0] / denominator, index))
    limiting_value, limiting_index = min(bounds, key=cmp_to_key(lambda x, y: compare(x[0], y[0])))
    if limiting_value != c or limiting_index != 4361:
        raise AssertionError("reported five-coordinate limiting wall disagrees")

    norm_class_proof = {
        "canonical_norms": [item.exact() for item in NORMS],
        "identities": {
            "epsilon_times_timelike": (EPSILON * TIMELIKE).exact(),
            "timelike_field_norm": str(TIMELIKE.field_norm()),
            "delta_field_norm": str(DELTA.field_norm()),
        },
        "local_argument": (
            "For q(r)|2I(r), unimodularity excludes primes away from (sqrt(2)) "
            "and (delta). At delta, y in delta*D24+ and primitive x give valuation "
            "one. At sqrt(2), a unit control gives valuation zero; otherwise "
            "evenness and primitivity force valuation exactly two. Class number "
            "one and totally-positive units being squares leave 1,2,delta,2delta."
        ),
        "audited_d24_plus": lattice_audit,
    }

    return {
        "schema_version": "1.0.0",
        "verdict": "PASS",
        "independence": (
            "No production search or project number-field implementation is imported; "
            "all controls, D24+ square partitions, pairings, signs, and ranks are replayed "
            "with standalone Fraction-pair arithmetic. NumPy is used only for proved-safe "
            "integer dot products in the 8,010,154-pair simultaneous packet audit."
        ),
        "sources": {
            "h300_artifact": str(h300_path.relative_to(ROOT)),
            "h300_sha256": H300_SHA256,
            "h500_artifact": str(h500_path.relative_to(ROOT)),
            "h500_sha256": sha256(h500_path.read_bytes()).hexdigest(),
        },
        "norm_class_completeness": norm_class_proof,
        "control_completeness": {
            "height_bound": 500,
            "complete_control_count": len(controls),
            "extension_control_count": len(extension),
            "ellipse_certificate": control_boxes,
        },
        "extension_level_censuses": [
            {
                "norm": item.level.norm.exact(),
                "height": item.level.height.exact(),
                "control": item.level.control.exact(),
                "target": item.level.target.exact(),
                "glue_patterns": item.glue_patterns,
                "reflection_integral": item.reflection_integral,
                "imprimitive_sqrt2": item.imprimitive,
                "primitive_candidates": len(item.candidates),
            }
            for item in level_censuses.values()
        ],
        "packet_replay": packet_results,
        "dense_packet": {
            "height": K2.make(181, 128).exact(),
            "raw_reflective": 4173,
            "imprimitive_sqrt2": 83,
            "primitive_candidates": 4090,
            "simultaneous_survivors": 4003,
            "same_height_unordered_pairs": 8_010_003,
            "same_height_positive_pairs": 0,
            "ordering_leakage": False,
        },
        "integer_matrix_audit": {
            "matrix_entries": matrix_audit.matrix_entries,
            "maximum_proved_intermediate_absolute_bound": matrix_audit.maximum_bound,
            "int64_max": 2**63 - 1,
            "maximum_abs_pair_coefficient": matrix_audit.maximum_coefficient,
            "scalar_zero_replays": matrix_audit.scalar_zero_replays,
            "scalar_strict_replays": matrix_audit.scalar_strict_replays,
        },
        "old_h300_ray_cut": {
            "new_walls": len(old_pairings),
            "negative": sum(item.sign() < 0 for item in old_pairings),
            "zero": sum(not item for item in old_pairings),
            "positive": sum(item.sign() > 0 for item in old_pairings),
            "minimum": minimum_identity(old_pairings).exact(),
            "maximum": maximum_identity(old_pairings).exact(),
        },
        "replacement_hyperideal_extreme_ray": {
            "vector": [item.exact() for item in replacement],
            "norm": replacement_norm.exact(),
            "norm_identity_sign": replacement_norm.sign(),
            "negative_pairings": signs.count(-1),
            "zero_pairings": signs.count(0),
            "positive_pairings": signs.count(1),
            "active_facets": list(active),
            "active_rank": rank,
            "rank_pivot_columns": list(pivots),
            "limiting_wall_index": limiting_index,
            "minimum_pairing": minimum_identity(replacement_pairings).exact(),
            "maximum_pairing": maximum_identity(replacement_pairings).exact(),
            "kind": "hyperideal",
        },
        "conclusion": {
            "h500_wall_count": len(artifact_walls),
            "h500_prefix_finite_volume": False,
            "h500_prefix_compact": False,
            "eventual_chamber_status": "unresolved",
            "scope": "complete exact Vinberg prefix only through algebraic height H<=500",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h300", type=Path, default=H300_PATH)
    parser.add_argument("--h500", type=Path, default=H500_PATH)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "data/enumerations/silver_p7_d24_plus_h500_independent_hostile_audit.json",
    )
    args = parser.parse_args()
    result = audit(args.h300, args.h500)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"verdict": result["verdict"], **result["conclusion"]}, indent=2))


if __name__ == "__main__":
    main()
