#!/usr/bin/env python3
"""Clean-room exact replay of the golden D24/D24+ prefix through x=phi^2.

This certificate intentionally does not import the production ``golden_dn``
search or its exact-arithmetic implementation.  It enumerates every dominant
norm-two spatial shell in the four residue modes of D24+ over Z[phi], then
applies the Vinberg inequalities against the D24 simple roots and previously
accepted positive-height walls.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key, lru_cache
from pathlib import Path
from typing import Iterable, Sequence


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True, slots=True)
class QPhi:
    """An independent exact ``a+b*phi`` value, with ``phi^2=phi+1``."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", Fraction(self.a))
        object.__setattr__(self, "b", Fraction(self.b))

    @classmethod
    def coerce(cls, value: QPhi | int | Fraction) -> QPhi:
        return value if isinstance(value, cls) else cls(Fraction(value))

    def __add__(self, other: QPhi | int | Fraction) -> QPhi:
        other = self.coerce(other)
        return QPhi(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> QPhi:
        return QPhi(-self.a, -self.b)

    def __sub__(self, other: QPhi | int | Fraction) -> QPhi:
        return self + -self.coerce(other)

    def __rsub__(self, other: QPhi | int | Fraction) -> QPhi:
        return self.coerce(other) - self

    def __mul__(self, other: QPhi | int | Fraction) -> QPhi:
        other = self.coerce(other)
        return QPhi(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    __rmul__ = __mul__

    def __truediv__(self, denominator: int) -> QPhi:
        if not isinstance(denominator, int):
            return NotImplemented
        return QPhi(self.a / denominator, self.b / denominator)

    def __pow__(self, exponent: int) -> QPhi:
        if exponent < 0:
            raise ValueError("negative powers are not needed in this audit")
        result = ONE
        base = self
        while exponent:
            if exponent & 1:
                result *= base
            base *= base
            exponent >>= 1
        return result

    def conjugate(self) -> QPhi:
        return QPhi(self.a + self.b, -self.b)

    def sign(self, place: int = 0) -> int:
        if place == 1:
            return self.conjugate().sign(0)
        if place != 0:
            raise ValueError("place must be 0 or 1")
        # 2(a+b phi)=(2a+b)+b sqrt(5).  Compare squares if the two
        # rational summands have opposite signs.
        rational = 2 * self.a + self.b
        radical = self.b
        rs = (rational > 0) - (rational < 0)
        bs = (radical > 0) - (radical < 0)
        if rs == 0:
            return bs
        if bs == 0 or rs == bs:
            return rs
        return rs if rational * rational > 5 * radical * radical else bs

    def __bool__(self) -> bool:
        return bool(self.a or self.b)

    def to_json(self) -> dict[str, str]:
        return {"a": _fraction_text(self.a), "b": _fraction_text(self.b)}

    def __str__(self) -> str:
        if not self.b:
            return _fraction_text(self.a)
        if not self.a:
            return f"{_fraction_text(self.b)}*phi"
        operator = "+" if self.b > 0 else "-"
        return f"{_fraction_text(self.a)} {operator} {_fraction_text(abs(self.b))}*phi"


ZERO = QPhi(0)
ONE = QPhi(1)
PHI = QPhi(0, 1)

Vector = tuple[QPhi, ...]
Residue = tuple[int, int]
RESIDUES: tuple[Residue, ...] = ((0, 0), (1, 0), (0, 1), (1, 1))
RESIDUE_NAMES = {
    (0, 0): "D24-root",
    (1, 0): "D24+-spinor-1",
    (0, 1): "D24+-spinor-phi",
    (1, 1): "D24+-spinor-1+phi",
}


def _cmp_at(left: QPhi, right: QPhi, place: int = 0) -> int:
    return (left - right).sign(place)


def _descending(values: Iterable[QPhi]) -> tuple[QPhi, ...]:
    return tuple(sorted(values, key=cmp_to_key(lambda left, right: -_cmp_at(left, right))))


def _square_alphabet(target4: QPhi, residue: Residue) -> tuple[QPhi, ...]:
    """Exhaustive doubled-coordinate alphabet for one O/2O residue.

    For every level here, ``z^2 <= target4`` at both real places gives
    ``|z| < 12`` and ``|z'| < 4``.  Writing ``z=A+B*phi`` then gives
    ``|B| <= 7`` and ``|A| <= 11``.  The slightly wider A interval below
    makes the boundary coverage visible without relying on approximation.
    """

    if (QPhi(144) - target4).sign(0) <= 0 or (QPhi(16) - target4).sign(1) <= 0:
        raise AssertionError("the documented coefficient box does not cover this target")
    pa, pb = residue
    values = []
    for a in range(-12, 13):
        for b in range(-7, 8):
            if a % 2 != pa or b % 2 != pb:
                continue
            value = QPhi(a, b)
            if all((target4 - value * value).sign(place) >= 0 for place in (0, 1)):
                values.append(value)
    return _descending(values)


def _extreme_square(
    values: Sequence[QPhi], start: int, place: int, *, maximum: bool
) -> QPhi:
    squares = (value * value for value in values[start:])
    comparison = (lambda left, right: _cmp_at(left, right, place))
    return (max if maximum else min)(squares, key=cmp_to_key(comparison))


def _enumerate_dominant_shell(
    dimension: int, target: QPhi, residue: Residue
) -> tuple[Vector, ...]:
    """Solve ``sum y_i^2=target`` in one D24+ residue exactly.

    Coordinates are represented as ``y_i=z_i/2``.  Nonincreasing identity
    values impose the D_n permutation chamber; the final inequality imposes
    ``y_(n-1)>=-y_n``.  The two parity accumulators impose the even D24
    translate sum coefficientwise.
    """

    target4 = 4 * target
    alphabet = _square_alphabet(target4, residue)
    pa, pb = residue
    suffix_extrema = {
        (start, place, maximum): _extreme_square(
            alphabet, start, place, maximum=maximum
        )
        for start in range(len(alphabet))
        for place in (0, 1)
        for maximum in (False, True)
    }

    @lru_cache(maxsize=None)
    def extend(
        slots: int,
        start: int,
        remainder: QPhi,
        translated_a_parity: int,
        translated_b_parity: int,
    ) -> tuple[Vector, ...]:
        if slots == 0:
            if remainder or translated_a_parity or translated_b_parity:
                return ()
            return ((),)
        output = []
        for index in range(start, len(alphabet)):
            value = alphabet[index]
            next_remainder = remainder - value * value
            if any(next_remainder.sign(place) < 0 for place in (0, 1)):
                continue
            remaining_slots = slots - 1
            if remaining_slots:
                feasible = True
                for place in (0, 1):
                    minimum = suffix_extrema[(index, place, False)]
                    maximum = suffix_extrema[(index, place, True)]
                    if (next_remainder - remaining_slots * minimum).sign(place) < 0:
                        feasible = False
                    if (remaining_slots * maximum - next_remainder).sign(place) < 0:
                        feasible = False
                if not feasible:
                    continue
            elif next_remainder:
                continue
            a_translate = (int(value.a) - pa) // 2
            b_translate = (int(value.b) - pb) // 2
            for tail in extend(
                remaining_slots,
                index,
                next_remainder,
                (translated_a_parity + a_translate) % 2,
                (translated_b_parity + b_translate) % 2,
            ):
                output.append((value, *tail))
        return tuple(output)

    doubled = extend(dimension, 0, target4, 0, 0)
    output = []
    for vector in doubled:
        if (vector[-2] + vector[-1]).sign(0) < 0:
            continue
        output.append(tuple(value / 2 for value in vector))
    return tuple(output)


def _d24_membership(vector: Sequence[QPhi]) -> bool:
    if any(value.a.denominator != 1 or value.b.denominator != 1 for value in vector):
        return False
    return sum(int(value.a) for value in vector) % 2 == 0 and sum(
        int(value.b) for value in vector
    ) % 2 == 0


def _spinor_residue(vector: Sequence[QPhi]) -> Residue | None:
    doubled = [2 * value for value in vector]
    if any(value.a.denominator != 1 or value.b.denominator != 1 for value in doubled):
        return None
    residues = {(int(value.a) % 2, int(value.b) % 2) for value in doubled}
    if len(residues) != 1:
        return None
    residue = residues.pop()
    if residue == (0, 0):
        return None
    pa, pb = residue
    translated_a = sum((int(value.a) - pa) // 2 for value in doubled)
    translated_b = sum((int(value.b) - pb) // 2 for value in doubled)
    return residue if translated_a % 2 == translated_b % 2 == 0 else None


def _pair(left: Sequence[QPhi], right: Sequence[QPhi]) -> QPhi:
    return -2 * PHI * left[0] * right[0] + sum(
        (x * y for x, y in zip(left[1:], right[1:])), start=ZERO
    )


def _height_zero_walls(dimension: int) -> tuple[Vector, ...]:
    walls = []
    for index in range(dimension - 1):
        wall = [ZERO] * (dimension + 1)
        wall[index + 1] = -ONE
        wall[index + 2] = ONE
        walls.append(tuple(wall))
    final = [ZERO] * (dimension + 1)
    final[-2] = final[-1] = -ONE
    walls.append(tuple(final))
    return tuple(walls)


def _vector_json(vector: Sequence[QPhi]) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def _support_json(vector: Sequence[QPhi]) -> list[dict[str, object]]:
    return [
        {"coordinate": index + 1, "value": value.to_json(), "display": str(value)}
        for index, value in enumerate(vector)
        if value
    ]


def _candidate_record(
    control: QPhi,
    spatial: Vector,
    residue: Residue,
    prior_walls: Sequence[Vector],
) -> dict[str, object]:
    vector = (control, *spatial)
    pairings = tuple(_pair(vector, wall) for wall in prior_walls)
    positives = [
        {"wall_index": index, "pairing": value.to_json(), "display": str(value)}
        for index, value in enumerate(pairings)
        if value.sign(0) > 0
    ]
    norm = _pair(vector, vector)
    spinor = _spinor_residue(spatial)
    membership_ok = _d24_membership(spatial) if residue == (0, 0) else spinor == residue
    return {
        "residue": list(residue),
        "coset": RESIDUE_NAMES[residue],
        "spatial_vector": _vector_json(spatial),
        "nonzero_support": _support_json(spatial),
        "norm": norm.to_json(),
        "norm_display": str(norm),
        "d24_root_member": _d24_membership(spatial),
        "d24_plus_spinor_residue": list(spinor) if spinor is not None else None,
        "membership_pass": membership_ok,
        "prior_wall_pairings": [_value_record(value) for value in pairings],
        "positive_pairing_witnesses": positives,
        "accepted": not positives,
        "_vector": vector,
    }


def _value_record(value: QPhi) -> dict[str, object]:
    return {"exact": value.to_json(), "display": str(value), "identity_sign": value.sign(0)}


def _strip_private(value: object) -> object:
    if isinstance(value, dict):
        return {key: _strip_private(item) for key, item in value.items() if not key.startswith("_")}
    if isinstance(value, list):
        return [_strip_private(item) for item in value]
    return value


def build_audit() -> dict[str, object]:
    dimension = 24
    controls = (ONE, PHI, PHI * PHI)
    height_zero = _height_zero_walls(dimension)
    shells: list[dict[str, object]] = []
    shell_vectors: dict[QPhi, dict[Residue, tuple[Vector, ...]]] = {}
    alphabet_sizes: dict[QPhi, dict[Residue, int]] = {}
    for control in controls:
        target = 2 + 2 * PHI * control * control
        target4 = 4 * target
        modes = {}
        sizes = {}
        for residue in RESIDUES:
            modes[residue] = _enumerate_dominant_shell(dimension, target, residue)
            sizes[residue] = len(_square_alphabet(target4, residue))
        shell_vectors[control] = modes
        alphabet_sizes[control] = sizes

    runs = {}
    accepted_d24_plus: list[Vector] = []
    accepted_records_d24_plus: list[dict[str, object]] = []
    for run_name, residues in (
        ("D24-root-lattice", ((0, 0),)),
        ("D24-plus", RESIDUES),
    ):
        walls = list(height_zero)
        levels = []
        run_accepted: list[Vector] = []
        for control in controls:
            target = 2 + 2 * PHI * control * control
            records = [
                _candidate_record(control, spatial, residue, walls)
                for residue in residues
                for spatial in shell_vectors[control][residue]
            ]
            survivors = [record for record in records if record["accepted"]]
            survivor_vectors = [record["_vector"] for record in survivors]
            same_level = [
                {
                    "left": left,
                    "right": right,
                    "pairing": _value_record(_pair(survivor_vectors[left], survivor_vectors[right])),
                }
                for left in range(len(survivor_vectors))
                for right in range(left + 1, len(survivor_vectors))
            ]
            if any(item["pairing"]["identity_sign"] > 0 for item in same_level):
                raise AssertionError("same-level survivors are not mutually acute")
            if any(record["norm"] != QPhi(2).to_json() or not record["membership_pass"] for record in records):
                raise AssertionError("candidate norm or lattice membership failed")
            walls.extend(survivor_vectors)
            run_accepted.extend(survivor_vectors)
            levels.append(
                {
                    "control": control.to_json(),
                    "control_display": str(control),
                    "height": (PHI * control * control).to_json(),
                    "target": target.to_json(),
                    "alphabet_sizes_by_residue": {
                        RESIDUE_NAMES[residue]: alphabet_sizes[control][residue]
                        for residue in residues
                    },
                    "candidate_counts_by_residue": {
                        RESIDUE_NAMES[residue]: len(shell_vectors[control][residue])
                        for residue in residues
                    },
                    "candidate_count": len(records),
                    "accepted_count": len(survivors),
                    "candidates": records,
                    "same_level_survivor_pairings": same_level,
                }
            )
        runs[run_name] = {
            "levels": levels,
            "candidate_counts": [level["candidate_count"] for level in levels],
            "accepted_counts": [level["accepted_count"] for level in levels],
            "final_wall_count": len(walls),
        }
        if run_name == "D24-plus":
            accepted_d24_plus = run_accepted
            accepted_records_d24_plus = [
                record
                for level in levels
                for record in level["candidates"]
                if record["accepted"]
            ]

    # Canonical presentation order: control, then wider support first.  This
    # fixes the two-root x=phi packet as phi^4 before 2phi.
    accepted_order = sorted(
        zip(accepted_d24_plus, accepted_records_d24_plus),
        key=lambda item: (
            controls.index(item[0][0]),
            -sum(bool(value) for value in item[0][1:]),
        ),
    )
    ordered_vectors = [item[0] for item in accepted_order]
    ordered_records = [item[1] for item in accepted_order]
    accepted_pairing_matrix = [
        [_value_record(_pair(left, right)) for right in ordered_vectors]
        for left in ordered_vectors
    ]
    height_zero_pairings = [
        [_value_record(_pair(vector, wall)) for wall in height_zero]
        for vector in ordered_vectors
    ]

    expected_candidate_counts = [1, 3, 4]
    expected_accepted_counts = [1, 2, 1]
    for run in runs.values():
        if run["candidate_counts"] != expected_candidate_counts:
            raise AssertionError("unexpected complete candidate counts")
        if run["accepted_counts"] != expected_accepted_counts:
            raise AssertionError("unexpected accepted counts")
    if any(
        shell_vectors[control][residue]
        for control in controls
        for residue in RESIDUES[1:]
    ):
        raise AssertionError("a half-spinor shell unexpectedly appeared in this prefix")
    if any(
        entry["identity_sign"] > 0
        for row in accepted_pairing_matrix
        for entry in row
        if entry not in ()
    ):
        # Diagonal entries are +2, so only inspect off-diagonal below instead.
        pass
    if any(
        accepted_pairing_matrix[left][right]["identity_sign"] > 0
        for left in range(len(ordered_vectors))
        for right in range(left + 1, len(ordered_vectors))
    ):
        raise AssertionError("accepted roots have a positive mutual pairing")
    if any(entry["identity_sign"] > 0 for row in height_zero_pairings for entry in row):
        raise AssertionError("an accepted root violates a height-zero wall")

    payload = {
        "schema": "golden-d24-prefix-independent-audit-0.1.0",
        "scope": "Norm-two Vinberg prefix through controls x=1,phi,phi^2; no later-control or other-norm completeness claim.",
        "independence": {
            "production_search_imported": False,
            "production_exact_arithmetic_imported": False,
            "implementation": "standalone Fraction arithmetic and exhaustive four-residue enumeration",
        },
        "field": "Q(phi), phi^2=phi+1",
        "host": "<-2phi> direct-sum D24 or D24+",
        "dimension": 24,
        "coefficient_box_proof": {
            "doubled_coordinate": "z=A+B*phi=2y",
            "equation": "sum(z_i^2)=4*(2+2*phi*x^2)",
            "place_bounds": ["|z|<12", "|z'|<4"],
            "deduction_b": "|B|*sqrt(5)=|z-z'|<16, hence |B|<=7",
            "deduction_a": "|2A+B|=|z+z'|<16 and |B|<=7, hence |A|<=11",
            "enumerated_ranges": {"A": [-12, 12], "B": [-7, 7]},
        },
        "residue_modes": [
            {"residue": list(residue), "name": RESIDUE_NAMES[residue]}
            for residue in RESIDUES
        ],
        "runs": runs,
        "spinor_shell_counts": {
            str(control): {
                RESIDUE_NAMES[residue]: len(shell_vectors[control][residue])
                for residue in RESIDUES[1:]
            }
            for control in controls
        },
        "accepted_roots_canonical_order": ordered_records,
        "accepted_root_pairing_matrix": accepted_pairing_matrix,
        "accepted_to_height_zero_pairings": height_zero_pairings,
        "conclusion": {
            "root_lattice_candidate_counts": expected_candidate_counts,
            "root_lattice_accepted_counts": expected_accepted_counts,
            "d24_plus_candidate_counts": expected_candidate_counts,
            "d24_plus_accepted_counts": expected_accepted_counts,
            "spinor_candidates_through_phi_squared": 0,
            "accepted_root_count": 4,
            "accepted_roots_all_in_d24_root_coset": True,
            "all_accepted_norm_two": True,
            "all_accepted_pairings_acute": True,
        },
    }
    return _strip_private(payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/golden_d24_prefix_independent_audit.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload["conclusion"], sort_keys=True))


if __name__ == "__main__":
    main()
