#!/usr/bin/env python3
"""Clean-room exact audit of the silver prime-seven D24+ H=300 witness.

This script deliberately imports no project arithmetic or search code.  It
reconstructs the 24 negative D24 simple walls, reads the 324 positive-height
walls from the frozen schema-0.3 archive, and checks the replacement ray with
``fractions.Fraction`` arithmetic in Q(sqrt(2)).
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence


@dataclass(frozen=True, slots=True)
class K2:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @classmethod
    def from_json(cls, value: dict[str, str]) -> "K2":
        if set(value) != {"a", "b"}:
            raise AssertionError("malformed Q(sqrt(2)) value")
        return cls(Fraction(value["a"]), Fraction(value["b"]))

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

    def __bool__(self) -> bool:
        return bool(self.a or self.b)

    def sign(self, conjugate: bool = False) -> int:
        b = -self.b if conjugate else self.b
        sa = (self.a > 0) - (self.a < 0)
        sb = (b > 0) - (b < 0)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        rational_norm = self.a * self.a - 2 * b * b
        return sa if rational_norm > 0 else sb

    @property
    def integral(self) -> bool:
        return self.a.denominator == self.b.denominator == 1

    def exact(self) -> tuple[str, str]:
        return str(self.a), str(self.b)


ZERO = K2()
ONE = K2(Fraction(1))
TWO = K2(Fraction(2))
TIMELIKE = K2(Fraction(-1), Fraction(2))
DELTA = K2(Fraction(3), Fraction(1))
Vector = tuple[K2, ...]


def value(wrapper: dict[str, object]) -> K2:
    return K2.from_json(wrapper["exact"])


def vector(items: Sequence[dict[str, str]]) -> Vector:
    return tuple(K2.from_json(item) for item in items)


def add(values: Iterable[K2]) -> K2:
    total = ZERO
    for item in values:
        total = total + item
    return total


def lorentz(left: Sequence[K2], right: Sequence[K2]) -> K2:
    assert len(left) == len(right) == 25
    return -(TIMELIKE * left[0] * right[0]) + add(
        left[index] * right[index] for index in range(1, 25)
    )


def functional(wall: Sequence[K2]) -> Vector:
    return (-(TIMELIKE * wall[0]), *wall[1:])


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


def d24_plus_coset(spatial: Sequence[K2]) -> str | None:
    """Classify an orthonormal vector in O tensor D24+ from first principles."""

    if len(spatial) != 24:
        return None
    if all(item.integral for item in spatial):
        if all(
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


def d24_plus_basis() -> tuple[tuple[K2, ...], ...]:
    basis: list[tuple[K2, ...]] = []
    # alpha_2,...,alpha_23 in the standard D24 convention.
    for index in range(1, 23):
        item = [ZERO] * 24
        item[index] = ONE
        item[index + 1] = -ONE
        basis.append(tuple(item))
    fork = [ZERO] * 24
    fork[22] = fork[23] = ONE
    basis.append(tuple(fork))
    basis.append(tuple(K2(Fraction(1, 2)) for _ in range(24)))
    assert len(basis) == 24
    return tuple(basis)


def dot(left: Sequence[K2], right: Sequence[K2]) -> K2:
    return add(x * y for x, y in zip(left, right))


def reflection_integral(root: Vector, norm: K2) -> bool:
    pairings = (-(TIMELIKE * root[0]),) + tuple(
        dot(root[1:], basis_vector) for basis_vector in d24_plus_basis()
    )
    return all((TWO * item / norm).integral for item in pairings)


def row_rank(rows: Sequence[Sequence[K2]], width: int) -> tuple[int, tuple[int, ...]]:
    work = [list(row) for row in rows]
    rank = 0
    pivots: list[int] = []
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(work)) if work[index][column]),
            None,
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
                work[index][other] = (
                    work[index][other] - factor * work[rank][other]
                )
        pivots.append(column)
        rank += 1
        if rank == len(work):
            break
    return rank, tuple(pivots)


def audit(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    assert payload["schema_version"] == "0.3.0"
    assert payload["height_bound"] == 300
    enumeration = payload["enumeration"]
    assert enumeration["control_level_count"] == len(enumeration["control_levels"]) == 48
    assert enumeration["height_packet_count"] == len(enumeration["packets"]) == 33

    walls = list(negative_d24_simple_walls())
    assert all(lorentz(wall, wall) == TWO for wall in walls)
    assert all(reflection_integral(wall, TWO) for wall in walls)
    assert all(
        lorentz(left, right).sign() <= 0
        for index, left in enumerate(walls)
        for right in walls[index + 1 :]
    )
    accepted_count = 0
    packet_by_global_wall: dict[int, tuple[K2, int]] = {}
    for packet in enumeration["packets"]:
        height = value(packet["height"])
        roots = packet["accepted_roots"]
        assert packet["accepted_count"] == len(roots)
        for packet_index, item in enumerate(roots):
            root = vector(item["vector_orthonormal"])
            norm = value(item["norm"])
            target = value(item["target"])
            assert len(root) == 25 and root[0] == value(item["control"])
            assert root[0].sign() > 0
            assert norm.sign() > 0 and norm.sign(conjugate=True) > 0
            assert lorentz(root, root) == norm
            assert add(coordinate * coordinate for coordinate in root[1:]) == target
            assert TIMELIKE * root[0] * root[0] / norm == height
            assert (K2(Fraction(300)) - height).sign() >= 0
            # For prime-seven norms the archive's coset label belongs to the
            # reduced vector before multiplication by delta.  Multiplication
            # by delta permutes the three spinor residues modulo 2.
            reduced_spatial = (
                tuple(coordinate / DELTA for coordinate in root[1:])
                if norm in (DELTA, TWO * DELTA)
                else root[1:]
            )
            assert d24_plus_coset(root[1:]) is not None
            assert d24_plus_coset(reduced_spatial) == item["coset"]
            assert reflection_integral(root, norm)
            # This verifies the dominant-chamber orientation and every later
            # simultaneous-height orientation, not just the witness pairing.
            assert all(lorentz(root, previous).sign() <= 0 for previous in walls)
            packet_by_global_wall[len(walls)] = (height, packet_index)
            walls.append(root)
            accepted_count += 1

    assert accepted_count == enumeration["accepted_root_count"] == 324
    assert len(walls) == enumeration["final_wall_count"] == 348
    assert len(set(walls)) == 348

    c = K2(Fraction(172, 7), Fraction(-120, 7))
    witness = (ONE, *(c for _ in range(19)), *(ZERO for _ in range(5)))
    frozen = payload["h300_replacement_hyperideal_ray"]
    assert vector(frozen["vector"]) == witness
    assert witness[0].sign() > 0
    pairings = tuple(lorentz(witness, wall) for wall in walls)
    assert all(item.sign() <= 0 for item in pairings)
    active = tuple(index for index, item in enumerate(pairings) if not item)
    expected_active = (*range(18), *range(19, 24), 59)
    assert active == expected_active
    assert len(active) == frozen["zero_pairing_count"] == 24
    assert list(active) == frozen["active_facets"]

    active_rows = tuple(functional(walls[index]) for index in active)
    assert all(dot(row, witness) == ZERO for row in active_rows)
    rank, pivots = row_rank(active_rows, 25)
    assert rank == frozen["active_rank"] == 24

    witness_norm = lorentz(witness, witness)
    expected_norm = K2(Fraction(1109345, 49), Fraction(-784418, 49))
    assert witness_norm == expected_norm == value(frozen["norm"])
    assert witness_norm.sign() > 0
    assert witness_norm.a * witness_norm.a - 2 * witness_norm.b * witness_norm.b == Fraction(67439, 7)
    identity_min = pairings[0]
    identity_max = pairings[0]
    for item in pairings[1:]:
        if (item - identity_min).sign() < 0:
            identity_min = item
        if (item - identity_max).sign() > 0:
            identity_max = item
    assert identity_min == value(frozen["minimum_pairing"])
    assert identity_max == value(frozen["maximum_pairing"]) == ZERO
    assert frozen["all_wall_pairings_nonpositive"] is True
    assert frozen["finite_volume"] is False and frozen["compact"] is False

    active_height, active_packet_index = packet_by_global_wall[59]
    assert active_height == K2(Fraction(62), Fraction(44))
    assert active_packet_index == 3
    assert payload["geometry_status"]["eventual_chamber_status"] == "unresolved"

    return {
        "schema_version": payload["schema_version"],
        "wall_count": len(walls),
        "accepted_positive_height_wall_count": accepted_count,
        "all_walls_valid_and_chronologically_oriented": True,
        "all_witness_pairings_nonpositive": True,
        "zero_pairing_count": len(active),
        "active_facets": list(active),
        "active_rank_over_Qsqrt2": rank,
        "rank_pivot_columns": list(pivots),
        "witness_norm_exact": {
            "a": witness_norm.exact()[0],
            "b": witness_norm.exact()[1],
        },
        "witness_norm_identity_sign": witness_norm.sign(),
        "active_positive_height_wall": {
            "global_index": 59,
            "packet_height": {"a": str(active_height.a), "b": str(active_height.b)},
            "zero_based_index_within_packet": active_packet_index,
        },
        "h300_prefix_finite_volume": False,
        "eventual_chamber_status": "unresolved",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "artifact",
        nargs="?",
        type=Path,
        default=Path("data/enumerations/silver_p7_d24_plus_h300.json"),
    )
    args = parser.parse_args()
    print(json.dumps(audit(args.artifact), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
