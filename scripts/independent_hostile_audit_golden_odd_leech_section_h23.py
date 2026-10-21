#!/usr/bin/env python3
"""Clean-room audit of the two golden odd-Leech-section H23 prefixes.

No project module is imported.  The section is reconstructed directly from
the twelve catalogue Golay generator masks and Construction B.  Standalone
``Q(phi)`` arithmetic then replays the norm alphabet consequences, height
zero, controls, coefficient Gram types, shell counts, packet acuteness,
reflection integrality, primitivity, and strict hyperideal witnesses.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key, lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_PATH = ROOT / "data/enumerations/golden_odd_leech_section_h23_first_packets.json"
OUTPUT_PATH = (
    ROOT
    / "data/enumerations/golden_odd_leech_section_h23_first_packets_independent_hostile_audit.json"
)

# These are the twelve binary rows encoded by the Nebe--Sloane catalogue
# basis: eleven displayed octad rows and the all-ones word.  Keeping the
# masks here avoids importing either the golden owner or its lattice helpers.
GOLAY_GENERATORS = (
    0x0000FF,
    0x000F0F,
    0x003333,
    0x005555,
    0x009999,
    0x030395,
    0x050539,
    0x090953,
    0x11111E,
    0x333300,
    0x555500,
    0xFFFFFF,
)


@dataclass(frozen=True, slots=True)
class K5:
    """Exact ``a+b*phi``, where ``phi^2=phi+1``."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @classmethod
    def make(
        cls, a: int | str | Fraction = 0, b: int | str | Fraction = 0
    ) -> "K5":
        return cls(Fraction(a), Fraction(b))

    @classmethod
    def from_exact(cls, value: dict[str, str]) -> "K5":
        return cls.make(value["a"], value["b"])

    def __add__(self, other: "K5") -> "K5":
        return K5(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "K5") -> "K5":
        return K5(self.a - other.a, self.b - other.b)

    def __neg__(self) -> "K5":
        return K5(-self.a, -self.b)

    def __mul__(self, other: "K5") -> "K5":
        return K5(
            self.a * other.a + self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    def __truediv__(self, other: "K5") -> "K5":
        norm = other.a * other.a + other.a * other.b - other.b * other.b
        if not norm:
            raise ZeroDivisionError
        return self * K5((other.a + other.b) / norm, -other.b / norm)

    def __pow__(self, exponent: int) -> "K5":
        if exponent < 0:
            return (ONE / self) ** (-exponent)
        output = ONE
        base = self
        while exponent:
            if exponent & 1:
                output = output * base
            base = base * base
            exponent >>= 1
        return output

    def norm(self) -> Fraction:
        return self.a * self.a + self.a * self.b - self.b * self.b

    def sign(self, embedding: int = 0) -> int:
        if embedding == 1:
            return K5(self.a + self.b, -self.b).sign(0)
        if embedding != 0:
            raise ValueError("embedding must be zero or one")
        rational = 2 * self.a + self.b
        radical = self.b
        sr = (rational > 0) - (rational < 0)
        sb = (radical > 0) - (radical < 0)
        if not sr:
            return sb
        if not sb or sr == sb:
            return sr
        return sr if rational * rational - 5 * radical * radical > 0 else sb

    @property
    def integral(self) -> bool:
        return self.a.denominator == self.b.denominator == 1

    def exact(self) -> dict[str, str]:
        return {"a": str(self.a), "b": str(self.b)}

    def display(self) -> str:
        if not self.b:
            return str(self.a)
        if not self.a:
            return f"{self.b}*phi"
        operator = "+" if self.b > 0 else "-"
        return f"{self.a} {operator} {abs(self.b)}*phi"


ZERO = K5.make()
ONE = K5.make(1)
TWO = K5.make(2)
PHI = K5.make(0, 1)
NORMS = tuple(K5.make(value) for value in (1, 2, 3, 6))
THREE_GATED = frozenset((K5.make(3), K5.make(6)))
COEFFICIENT_BOX = 12


@lru_cache(maxsize=1)
def golay_words() -> tuple[int, ...]:
    words = []
    for message in range(1 << 12):
        word = 0
        for bit, generator in enumerate(GOLAY_GENERATORS):
            if message & (1 << bit):
                word ^= generator
        words.append(word)
    words = tuple(sorted(set(words)))
    distribution = Counter(word.bit_count() for word in words)
    expected = Counter({0: 1, 8: 759, 12: 2576, 16: 759, 24: 1})
    if len(words) != 4096 or distribution != expected:
        raise AssertionError("standalone Golay reconstruction failed")
    if any(
        (left & generator).bit_count() % 2
        for left in words
        for generator in GOLAY_GENERATORS
    ):
        raise AssertionError("the reconstructed code is not self-orthogonal")
    return words


def rank_mod_prime(rows: list[list[int]], prime: int) -> int:
    work = [[value % prime for value in row] for row in rows]
    rank = 0
    for column in range(len(work[0])):
        pivot = next(
            (index for index in range(rank, len(work)) if work[index][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(inverse * value) % prime for value in work[rank]]
        for index in range(len(work)):
            if index == rank or not work[index][column]:
                continue
            factor = work[index][column]
            work[index] = [
                (value - factor * other) % prime
                for value, other in zip(work[index], work[rank])
            ]
        rank += 1
    return rank


def null_vector_mod_prime(matrix: Matrix, prime: int) -> tuple[int, ...]:
    size = matrix.rows
    work = [[int(matrix[i, j]) % prime for j in range(size)] for i in range(size)]
    pivots = []
    rank = 0
    for column in range(size):
        pivot = next(
            (index for index in range(rank, size) if work[index][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(inverse * value) % prime for value in work[rank]]
        for index in range(size):
            if index == rank or not work[index][column]:
                continue
            factor = work[index][column]
            work[index] = [
                (value - factor * other) % prime
                for value, other in zip(work[index], work[rank])
            ]
        pivots.append(column)
        rank += 1
    free = tuple(column for column in range(size) if column not in pivots)
    if len(free) != 1:
        raise AssertionError("the mod-prime kernel is not one-dimensional")
    vector = [0] * size
    vector[free[0]] = 1
    for row, pivot in enumerate(pivots):
        vector[pivot] = (-work[row][free[0]]) % prime
    return tuple(vector)


def section_model() -> dict[str, object]:
    """Reconstruct the integral section from Construction B generators."""

    columns = [
        Matrix([2 * ((generator >> index) & 1) for index in range(24)])
        for generator in GOLAY_GENERATORS
    ]
    # 4*Z_+^24 is generated by 4(e_i+e_23) and 8e_23.
    for index in range(23):
        vector = [0] * 24
        vector[index] = vector[23] = 4
        columns.append(Matrix(vector))
    vector = [0] * 24
    vector[23] = 8
    columns.append(Matrix(vector))
    columns.append(Matrix([1] * 24))
    ambient = hermite_normal_form(Matrix.hstack(*columns))
    if ambient.shape != (24, 24):
        raise AssertionError("Construction B did not produce a full lattice")
    ambient_gram = ambient.T * ambient / 8
    if any(value.q != 1 for value in ambient_gram) or int(ambient_gram.det()) != 1:
        raise AssertionError("the reconstructed odd-Leech lattice is not integral unimodular")

    pairings = [int(sum(ambient[:, index]) / 8) for index in range(24)]
    pivot = next((index for index, value in enumerate(pairings) if abs(value) == 1), None)
    if pivot is None:
        raise AssertionError("the section normal ceased to be primitive")
    pivot_value = pairings[pivot]
    section = Matrix.hstack(
        *(
            ambient[:, index]
            - pairings[index] * pivot_value * ambient[:, pivot]
            for index in range(24)
            if index != pivot
        )
    )
    if any(sum(section[:, index]) for index in range(23)):
        raise AssertionError("a section basis vector is not zero-sum")
    gram = section.T * section / 8
    if any(value.q != 1 for value in gram):
        raise AssertionError("the section Gram matrix is not integral")
    gram = Matrix([[int(value) for value in gram.row(index)] for index in range(23)])
    smith = smith_normal_form(gram, domain=ZZ)
    invariants = tuple(abs(int(smith[i, i])) for i in range(23))
    if int(gram.det()) != 3 or invariants != (1,) * 22 + (3,):
        raise AssertionError("the section determinant/Smith form drifted")
    radical = null_vector_mod_prime(gram, 3)
    radical_raw = section * Matrix(radical)
    residue = tuple(int(value) % 3 for value in radical_raw)
    if len(set(residue)) != 1 or residue[0] == 0:
        raise AssertionError("the section radical lost its common raw residue")

    minimal = [
        [(-1 if word & (1 << index) else 1) for index in range(24)]
        for word in golay_words()
        if word.bit_count() == 12
    ]
    if len(minimal) != 2576 or rank_mod_prime(minimal, 1_000_003) != 23:
        raise AssertionError("the section minimal shell failed to span")
    digest = sha256()
    for row in minimal:
        digest.update(bytes(value + 1 for value in row))
    return {
        "ambient_raw_index": abs(int(ambient.det())),
        "ambient_gram_determinant": int(ambient_gram.det()),
        "rank": 23,
        "determinant": int(gram.det()),
        "smith_invariants": list(invariants),
        "gram_rank_mod_three": rank_mod_prime(
            [[int(gram[i, j]) for j in range(23)] for i in range(23)], 3
        ),
        "pairing_radical_dimension_mod_three": 1,
        "radical_raw_residue_mod_three": list(residue),
        "minimum": 3,
        "minimum_proof": (
            "The shifted Construction-B coset has raw square at least 24; "
            "the even coset has raw square at least 32 because a nonzero "
            "Golay parity support has weight at least eight and the zero "
            "support requires at least two nonzero 4-multiples by its even "
            "residual-sum condition."
        ),
        "minimal_shell_count": len(minimal),
        "minimal_shell_span_rank": 23,
        "minimal_shell_sha256": digest.hexdigest(),
    }


def construction_b_count(
    target_raw_square: int,
    *,
    require_common_residue_mod_three: bool = False,
) -> tuple[int, int]:
    """Count all zero-sum Construction-B vectors at an exact raw square."""

    maximum = int(target_raw_square**0.5) + 1
    counts = []
    for shift in (0, 1):
        total = 0
        residues = range(3) if require_common_residue_mod_three else (None,)
        for residue in residues:
            for word in golay_words():
                states: dict[tuple[int, int, int], int] = {(0, 0, 0): 1}
                for index in range(24):
                    bit = (word >> index) & 1
                    values = []
                    for value in range(-maximum, maximum + 1):
                        if value % 2 != shift:
                            continue
                        if ((value - shift) // 2) % 2 != bit:
                            continue
                        if residue is not None and value % 3 != residue:
                            continue
                        residual = (value - shift - 2 * bit) // 4
                        values.append((value, residual % 2))
                    new_states: dict[tuple[int, int, int], int] = {}
                    for (coordinate_sum, square, parity), multiplicity in states.items():
                        for value, residual_parity in values:
                            new_square = square + value * value
                            if new_square > target_raw_square:
                                continue
                            key = (
                                coordinate_sum + value,
                                new_square,
                                parity ^ residual_parity,
                            )
                            new_states[key] = new_states.get(key, 0) + multiplicity
                    states = new_states
                    if not states:
                        break
                total += states.get((0, target_raw_square, 0), 0)
        counts.append(total)
    return counts[0], counts[1]


def radical_small_vector_count() -> int:
    """Count radical vectors of every positive raw square at most 48."""

    total = 0
    for raw_square in range(1, 49):
        total += sum(
            construction_b_count(
                raw_square, require_common_residue_mod_three=True
            )
        )
    return total


def binary_minimum(gram: tuple[int, int, int]) -> int:
    aa, cc, bb = gram
    if aa < 0 or cc < 0 or bb * bb > aa * cc:
        raise ValueError("expected a positive-semidefinite Gram")
    if aa == 0:
        return cc
    if cc == 0:
        return aa
    while True:
        if aa > cc:
            aa, cc = cc, aa
        if aa == 0:
            return cc
        quotient = (bb + aa // 2) // aa
        if quotient:
            cc, bb = (
                cc - 2 * quotient * bb + quotient * quotient * aa,
                bb - quotient * aa,
            )
            continue
        if abs(2 * bb) > aa:
            quotient = 1 if bb > 0 else -1
            cc, bb = (
                cc - 2 * quotient * bb + quotient * quotient * aa,
                bb - quotient * aa,
            )
            continue
        return aa


def gram_types(
    target: K5, *, pairing_modulus: int = 1, coefficient_minimum: int = 3
) -> tuple[tuple[int, int, int], ...]:
    if not target.integral or target.a < 0:
        return ()
    rational_target, phi_target = int(target.a), int(target.b)
    output = []
    for v_square in range(rational_target + 1):
        u_square = rational_target - v_square
        numerator = phi_target - v_square
        if numerator % 2:
            continue
        dot = numerator // 2
        if u_square and u_square < 3:
            continue
        if v_square and v_square < 3:
            continue
        if dot * dot > u_square * v_square:
            continue
        gram = (u_square, v_square, dot)
        if any(value % pairing_modulus for value in gram):
            continue
        if binary_minimum(gram) < coefficient_minimum:
            continue
        output.append(gram)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class Control:
    height: K5
    norm: K5
    x: K5
    target: K5
    modulus: int
    minimum: int
    grams: tuple[tuple[int, int, int], ...]


def compare_k5(left: K5, right: K5) -> int:
    return (left - right).sign(0)


def compare_controls(left: Control, right: Control) -> int:
    for first, second in (
        (left.height, right.height),
        (left.norm, right.norm),
        (left.x, right.x),
    ):
        comparison = compare_k5(first, second)
        if comparison:
            return comparison
    return 0


def enumerate_controls(coefficient: K5, endpoint: K5) -> tuple[Control, ...]:
    """Exhaust a proved coefficient box through the selected endpoint."""

    # For both coefficients, height and conjugate-target inequalities imply
    # |x_0|,|x_1|<5.  Thus |n|<=4 and |m|<=12 for x=m+n*phi.
    if (K5.make(25) * coefficient - K5.make(6) * endpoint).sign(0) <= 0:
        raise AssertionError("the identity control bound five was not proved")
    if (K5.make(6) + K5.make(25) * K5(coefficient.a + coefficient.b, -coefficient.b)).sign(0) >= 0:
        raise AssertionError("the conjugate control bound five was not proved")
    output = set()
    for norm in NORMS:
        for rational in range(-COEFFICIENT_BOX, COEFFICIENT_BOX + 1):
            for radical in range(-4, 5):
                x = K5.make(rational, radical)
                if x.sign(0) <= 0:
                    continue
                if not (TWO * coefficient * x / norm).integral:
                    continue
                target = norm + coefficient * x * x
                if target.sign(1) < 0:
                    continue
                height = coefficient * x * x / norm
                if (endpoint - height).sign(0) < 0:
                    continue
                gated = norm in THREE_GATED
                modulus = 3 if gated else 1
                minimum = 7 if gated else 3
                output.add(
                    Control(
                        height,
                        norm,
                        x,
                        target,
                        modulus,
                        minimum,
                        gram_types(
                            target,
                            pairing_modulus=modulus,
                            coefficient_minimum=minimum,
                        ),
                    )
                )
    return tuple(sorted(output, key=cmp_to_key(compare_controls)))


def control_payload(control: Control) -> dict[str, object]:
    return {
        "height": control.height.exact(),
        "norm": control.norm.exact(),
        "control": control.x.exact(),
        "target": control.target.exact(),
        "pairing_modulus": control.modulus,
        "coefficient_minimum": control.minimum,
        "necessary_gram_triples": [
            {"u_square": row[0], "v_square": row[1], "u_dot_v": row[2]}
            for row in control.grams
        ],
    }


def compare_production(prefixes: dict[str, object], shells: dict[str, object]) -> None:
    production = json.loads(PRODUCTION_PATH.read_text())
    if production["section_shells"]["norm_4"]["total_count"] != shells["norm_4_total"]:
        raise AssertionError("production norm-four shell count disagrees")
    if production["section_shells"]["norm_5"]["total_count"] != shells["norm_5_total"]:
        raise AssertionError("production norm-five shell count disagrees")
    for label in ("a=phi", "a=2phi"):
        frozen = production["prefixes"][label]
        replay = prefixes[label]
        if frozen["complete_controls"] != replay["complete_controls"]:
            raise AssertionError(f"production {label} control census disagrees")
        if frozen["first_packet"] != replay["first_packet"]:
            raise AssertionError(f"production {label} packet disagrees")
        if frozen["witness"]["square"] != replay["witness_square"]:
            raise AssertionError(f"production {label} witness square disagrees")
        if frozen["witness"]["closest_wall_pairing"] != replay["closest_wall_pairing"]:
            raise AssertionError(f"production {label} witness pairing disagrees")


def packet_audits() -> dict[str, object]:
    specifications = (
        (
            "a=phi",
            PHI,
            K5.make(Fraction(3, 2), Fraction(5, 2)),
            K5.make(1, 1),
            5,
            553_104,
            3,
            4,
            Fraction(2, 3),
        ),
        (
            "a=2phi",
            K5.make(0, 2),
            K5.make(1, 2),
            PHI,
            4,
            53_682,
            2,
            3,
            Fraction(4, 3),
        ),
    )
    output = {}
    for (
        label,
        coefficient,
        endpoint,
        x,
        spatial_norm,
        shell_count,
        distinct_dot_bound,
        witness_norm,
        witness_scale,
    ) in specifications:
        controls = enumerate_controls(coefficient, endpoint)
        nonempty = tuple(control for control in controls if control.grams)
        expected_count = 7 if label == "a=phi" else 2
        expected_gram = (0, spatial_norm, 0)
        if (
            len(controls) != expected_count
            or len(nonempty) != 1
            or nonempty[0].height != endpoint
            or nonempty[0].norm != TWO
            or nonempty[0].x != x
            or nonempty[0].grams != (expected_gram,)
        ):
            raise AssertionError(f"the {label} first prefix drifted")
        scale = PHI
        root_square = -coefficient * x * x + K5.make(spatial_norm) * scale * scale
        time_reflection_coefficient = TWO * coefficient * x / TWO
        spatial_reflection_coefficient = TWO * scale / TWO
        if (
            root_square != TWO
            or not time_reflection_coefficient.integral
            or not spatial_reflection_coefficient.integral
            or abs(x.norm()) != 1
        ):
            raise AssertionError(f"the {label} packet roots lost validity")
        pairing_bound = (
            -coefficient * x * x
            + K5.make(distinct_dot_bound) * scale * scale
        )
        witness_square = (
            K5.make(Fraction(witness_scale * witness_scale) * witness_norm)
            - coefficient
        )
        closest = (
            -coefficient * x
            + K5.make(Fraction(witness_scale) * distinct_dot_bound) * scale
        )
        if pairing_bound.sign(0) >= 0 or witness_square.sign(0) <= 0 or closest.sign(0) >= 0:
            raise AssertionError(f"the {label} acuteness/witness audit failed")
        first_packet = {
            "height": endpoint.exact(),
            "norm": TWO.exact(),
            "control": x.exact(),
            "target": nonempty[0].target.exact(),
            "spatial_scale": scale.exact(),
            "spatial_norm": spatial_norm,
            "gram_type": {
                "u_square": 0,
                "v_square": spatial_norm,
                "u_dot_v": 0,
            },
            "root_count": shell_count,
            "pairing_upper_bound": pairing_bound.exact(),
        }
        output[label] = {
            "complete_controls": [control_payload(control) for control in controls],
            "first_packet": first_packet,
            "packet_logic": {
                "control_is_unit": True,
                "primitive_because_control_is_unit": True,
                "root_square": root_square.exact(),
                "time_reflection_coefficient": time_reflection_coefficient.exact(),
                "spatial_reflection_coefficient": spatial_reflection_coefficient.exact(),
                "distinct_spatial_dot_upper_bound": distinct_dot_bound,
                "mutually_acute": True,
            },
            "witness_square": witness_square.exact(),
            "closest_wall_pairing": closest.exact(),
            "strict_hyperideal": True,
        }
    return output


def audit() -> dict[str, object]:
    model = section_model()
    if model["minimal_shell_sha256"] != "2860f05b59483cf95a28433e76915d6862cf05fd1348303e938233033e2923b9":
        raise AssertionError("the independently reconstructed minimal-shell digest changed")
    radical_count = radical_small_vector_count()
    if radical_count:
        raise AssertionError("the mod-three radical has a vector of norm at most six")
    shell4 = construction_b_count(32)
    shell5 = construction_b_count(40)
    shells = {
        "norm_4_by_coset": list(shell4),
        "norm_4_total": sum(shell4),
        "norm_5_by_coset": list(shell5),
        "norm_5_total": sum(shell5),
        "method": (
            "Full Construction-B dynamic program over all 4096 reconstructed "
            "Golay words, enforcing coordinate sum zero, exact raw square, "
            "code congruence, and even residual sum."
        ),
    }
    if shell4 != (53_682, 0) or shell5 != (0, 553_104):
        raise AssertionError("the standalone section shell census drifted")

    height_zero_ambient = [len(gram_types(norm)) for norm in NORMS]
    height_zero_gated = [
        len(
            gram_types(
                norm,
                pairing_modulus=3 if norm in THREE_GATED else 1,
                coefficient_minimum=7 if norm in THREE_GATED else 3,
            )
        )
        for norm in NORMS
    ]
    if height_zero_ambient != [0, 0, 1, 1] or height_zero_gated != [0, 0, 0, 0]:
        raise AssertionError("the standalone height-zero census drifted")

    prefixes = packet_audits()
    compare_production(prefixes, shells)
    return {
        "schema_version": "0.1.0",
        "verdict": "PASS_WITH_SCHEMA_FIX",
        "scope": (
            "Independent exact audit of the a=phi and a=2phi golden S23 "
            "prefixes through their first nonempty packet; bounded prefixes only."
        ),
        "independence": (
            "No project imports.  The odd-Leech section is reconstructed from "
            "standalone Golay masks and Construction B; Q(phi), controls, "
            "Gram types, and geometric inequalities are independently replayed."
        ),
        "source": {
            "production_artifact": str(PRODUCTION_PATH.relative_to(ROOT)),
            "production_sha256": sha256(PRODUCTION_PATH.read_bytes()).hexdigest(),
        },
        "section_model": model,
        "radical_minimum": {
            "nonzero_vectors_of_raw_square_at_most_48": radical_count,
            "conclusion": "minimum at least 7",
        },
        "norm_alphabet": {
            "canonical": [norm.exact() for norm in NORMS],
            "prime_two_inert": all((x * x - x - 1) % 2 for x in range(2)),
            "prime_three_inert": all((x * x - x - 1) % 3 for x in range(3)),
            "local_exponents": {"at_(2)": [0, 1], "at_(3)": [0, 1]},
            "normalization": (
                "Z[phi] is norm-Euclidean; its units are +/-phi^n, and its "
                "totally-positive units phi^(2n) are unit squares.  Combining "
                "the independently recovered Smith support {3} with q|2I(r) "
                "gives exactly the four products 1,2,3,6."
            ),
        },
        "height_zero": {
            "ambient_type_counts": height_zero_ambient,
            "gated_type_counts": height_zero_gated,
            "wall_count": 0,
        },
        "section_shells": shells,
        "prefixes": prefixes,
        "schema_defect": {
            "status": "FIXED_IN_PRODUCTION_GENERATOR",
            "old_key": "section_model.radical_mod_three_rank=22",
            "reason": (
                "22 is the Gram rank modulo three, not the radical rank; "
                "the pairing radical has dimension one."
            ),
            "replacement_keys": {
                "gram_rank_mod_three": 22,
                "pairing_radical_dimension_mod_three": 1,
            },
        },
        "conclusion": (
            "Both production bounded-prefix results pass.  Each first packet "
            "is complete, primitive, reflective, and mutually acute, but its "
            "prefix cone contains the certified strict spacelike witness.  "
            "Neither is finite-volume or compact at this prefix; the eventual "
            "chambers remain unresolved."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    arguments = parser.parse_args()
    payload = audit()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "verdict": payload["verdict"],
                "norm_4": payload["section_shells"]["norm_4_total"],
                "norm_5": payload["section_shells"]["norm_5_total"],
                "phi_controls": len(payload["prefixes"]["a=phi"]["complete_controls"]),
                "twice_phi_controls": len(payload["prefixes"]["a=2phi"]["complete_controls"]),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
