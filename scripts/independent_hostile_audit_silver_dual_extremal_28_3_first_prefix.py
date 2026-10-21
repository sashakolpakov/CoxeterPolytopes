#!/usr/bin/env python3
"""Owner-import-free hostile audit of the silver rank-28 first prefix.

The production module is parsed only as inert source text to recover its
literal 28-by-28 catalogue Gram matrix.  This replay supplies its own
Q(sqrt(2)) arithmetic, mod-three kernel, local norm alphabet, control
enumerator, binary Gram gate, augmented Gram maximizer, and witness checks.
It never imports ``silver_dual_extremal_28_3``.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import isqrt
from pathlib import Path
from typing import Any

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
OWNER_SOURCE = (
    ROOT / "src/coxeter_polytopes/search/silver_dual_extremal_28_3.py"
)
PRODUCTION_ARTIFACT = (
    ROOT / "data/enumerations/silver_dual_extremal_28_3_first_prefix.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "silver_dual_extremal_28_3_first_prefix_independent_hostile_audit.json"
)
REPORT = (
    ROOT
    / "results/"
    "silver_dual_extremal_28_3_first_prefix_independent_hostile_audit.md"
)

CATALOGUE_URL = (
    "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/"
    "DualExtremal_28_3.html"
)
EXPECTED_GRAM_STREAM_SHA256 = (
    "d61f9a887f3b31db364b6c38c9b0e59b615c4cd178201ae3e8d0036a87629c54"
)
DIMENSION = 28
EXPECTED_DETERMINANT = 3**26
EXPECTED_SMITH = (1, 1) + (3,) * 26


@dataclass(frozen=True, slots=True)
class K2:
    """Standalone exact element ``a+b*sqrt(2)``."""

    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @classmethod
    def make(
        cls, a: int | str | Fraction = 0, b: int | str | Fraction = 0
    ) -> "K2":
        return cls(Fraction(a), Fraction(b))

    @classmethod
    def coerce(cls, value: "K2 | int | Fraction") -> "K2":
        return value if isinstance(value, cls) else cls.make(value)

    def __add__(self, other: "K2 | int | Fraction") -> "K2":
        other = K2.coerce(other)
        return K2(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __sub__(self, other: "K2 | int | Fraction") -> "K2":
        other = K2.coerce(other)
        return K2(self.a - other.a, self.b - other.b)

    def __rsub__(self, other: "K2 | int | Fraction") -> "K2":
        return K2.coerce(other) - self

    def __neg__(self) -> "K2":
        return K2(-self.a, -self.b)

    def __mul__(self, other: "K2 | int | Fraction") -> "K2":
        other = K2.coerce(other)
        return K2(
            self.a * other.a + 2 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    __rmul__ = __mul__

    def __truediv__(self, other: "K2 | int | Fraction") -> "K2":
        other = K2.coerce(other)
        norm = other.norm()
        if not norm:
            raise ZeroDivisionError
        return self * K2(other.a / norm, -other.b / norm)

    def __rtruediv__(self, other: "K2 | int | Fraction") -> "K2":
        return K2.coerce(other) / self

    def __pow__(self, exponent: int) -> "K2":
        if exponent < 0:
            return (ONE / self) ** (-exponent)
        output = ONE
        base = self
        while exponent:
            if exponent & 1:
                output *= base
            base *= base
            exponent >>= 1
        return output

    def norm(self) -> Fraction:
        return self.a * self.a - 2 * self.b * self.b

    def conjugate(self) -> "K2":
        return K2(self.a, -self.b)

    def sign(self, place: int = 0) -> int:
        if place == 1:
            return self.conjugate().sign()
        if place != 0:
            raise ValueError("place must be zero or one")
        sa = (self.a > 0) - (self.a < 0)
        sb = (self.b > 0) - (self.b < 0)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        return sa if self.norm() > 0 else sb

    @property
    def is_integral(self) -> bool:
        return self.a.denominator == self.b.denominator == 1

    def __lt__(self, other: "K2 | int | Fraction") -> bool:
        return (self - other).sign() < 0

    def __le__(self, other: "K2 | int | Fraction") -> bool:
        return (self - other).sign() <= 0

    def __str__(self) -> str:
        def part(value: Fraction) -> str:
            return str(value.numerator) if value.denominator == 1 else str(value)

        if not self.b:
            return part(self.a)
        if not self.a:
            return f"{part(self.b)}*sqrt(2)"
        operator = "+" if self.b > 0 else "-"
        return f"{part(self.a)} {operator} {part(abs(self.b))}*sqrt(2)"

    def exact(self) -> dict[str, str]:
        return {"a": str(self.a), "b": str(self.b)}


ZERO = K2.make()
ONE = K2.make(1)
SQRT2 = K2.make(0, 1)
RHO = 2 * SQRT2 - ONE
EPSILON = ONE + SQRT2
P7 = EPSILON * RHO
HEIGHT_BOUND = K2.make(93, 66)


def parse_literal_gram() -> tuple[Matrix, str, str]:
    """Parse the owner's literal data without executing or importing it."""

    source_text = OWNER_SOURCE.read_text()
    syntax = ast.parse(source_text, filename=str(OWNER_SOURCE))
    gram_text = None
    source_url = None
    for node in syntax.body:
        if not isinstance(node, ast.Assign):
            continue
        names = {target.id for target in node.targets if isinstance(target, ast.Name)}
        if "_GRAM_TEXT" in names:
            gram_text = ast.literal_eval(node.value)
        if "SOURCE_URL" in names:
            source_url = ast.literal_eval(node.value)
    if not isinstance(gram_text, str) or source_url != CATALOGUE_URL:
        raise AssertionError("owner source lost the literal primary catalogue data")
    rows = [list(map(int, line.split())) for line in gram_text.splitlines() if line]
    if len(rows) != DIMENSION or any(len(row) != DIMENSION for row in rows):
        raise AssertionError("literal Gram has the wrong shape")
    normalized = "".join(" ".join(map(str, row)) + "\n" for row in rows)
    digest = hashlib.sha256(normalized.encode("ascii")).hexdigest()
    if digest != EXPECTED_GRAM_STREAM_SHA256:
        raise AssertionError(f"literal primary Gram digest changed: {digest}")
    return Matrix(rows), digest, hashlib.sha256(source_text.encode()).hexdigest()


def rank_mod_prime(matrix: Matrix, prime: int) -> int:
    work = [
        [int(matrix[i, j]) % prime for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]
    rank = 0
    for column in range(matrix.cols):
        pivot = next(
            (row for row in range(rank, matrix.rows) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(inverse * entry) % prime for entry in work[rank]]
        for row in range(matrix.rows):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (left - factor * right) % prime
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
    return rank


def pairing_radical_basis_mod_three(gram: Matrix) -> Matrix:
    """Independent HNF basis of ``{x:Gx=0 mod 3}``."""

    prime = 3
    reduced = [
        [int(gram[i, j]) % prime for j in range(DIMENSION)]
        for i in range(DIMENSION)
    ]
    pivots: list[int] = []
    row = 0
    for column in range(DIMENSION):
        pivot = next(
            (index for index in range(row, DIMENSION) if reduced[index][column]),
            None,
        )
        if pivot is None:
            continue
        reduced[row], reduced[pivot] = reduced[pivot], reduced[row]
        inverse = pow(reduced[row][column], -1, prime)
        reduced[row] = [(inverse * entry) % prime for entry in reduced[row]]
        for index in range(DIMENSION):
            if index == row or not reduced[index][column]:
                continue
            factor = reduced[index][column]
            reduced[index] = [
                (left - factor * right) % prime
                for left, right in zip(reduced[index], reduced[row])
            ]
        pivots.append(column)
        row += 1
    radical = []
    for free in (column for column in range(DIMENSION) if column not in pivots):
        vector = [0] * DIMENSION
        vector[free] = 1
        for index, pivot in enumerate(pivots):
            vector[pivot] = (-reduced[index][free]) % prime
        radical.append(Matrix(vector))
    basis = hermite_normal_form(
        Matrix.hstack(*(prime * Matrix.eye(DIMENSION)).columnspace(), *radical)
    )
    if basis.shape != (DIMENSION, DIMENSION):
        raise AssertionError("HNF did not return a full pairing-radical basis")
    if any(int(value) % prime for value in gram * basis):
        raise AssertionError("pairing-radical basis violates Gx=0 mod 3")
    return basis


def gp_matrix(matrix: Matrix) -> str:
    return "[" + ";".join(
        ",".join(str(int(matrix[i, j])) for j in range(matrix.cols))
        for i in range(matrix.rows)
    ) + "]"


def pari_shortest_audit(gram: Matrix, s3_basis: Matrix) -> dict[str, Any]:
    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required for the hostile replay")
    s3_gram = s3_basis.T * gram * s3_basis
    program = ";".join(
        (
            f"G={gp_matrix(gram)}",
            "U=qflllgram(G)",
            "Q=U~*G*U",
            "R=qfminim(Q,,1)",
            f"H={gp_matrix(s3_gram)}",
            "V=qflllgram(H)",
            "S=V~*H*V",
            "T=qfminim(S,,1)",
            'print("HOST ",R[1]," ",floor(R[2]+1/2))',
            'print("S3 ",T[1]," ",floor(T[2]+1/2))',
            "quit()",
        )
    )
    completed = subprocess.run(
        [gp, "-q", "-s", "4000000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    parsed: dict[str, tuple[int, int]] = {}
    for line in completed.stdout.splitlines():
        fields = line.split()
        if len(fields) == 3 and fields[0] in ("HOST", "S3"):
            parsed[fields[0]] = (int(fields[1]), int(fields[2]))
    if parsed != {"HOST": (1_260, 8), "S3": (65_520, 12)}:
        raise AssertionError(f"independent shortest-vector census changed: {parsed}")
    return {
        "engine": "independent HNF bases; PARI qflllgram and qfminim(Q,,m=1)",
        "host_signed_shortest_count": parsed["HOST"][0],
        "host_minimum": parsed["HOST"][1],
        "S3_signed_shortest_count": parsed["S3"][0],
        "S3_minimum": parsed["S3"][1],
    }


@dataclass(frozen=True, slots=True)
class NormGate:
    norm: K2
    dyadic: bool
    p7: bool
    discriminant_three: bool


NORM_GATES = tuple(
    NormGate(
        K2.make(2 if dyadic else 1)
        * (P7 if p7 else ONE)
        * (3 if discriminant_three else 1),
        dyadic,
        p7,
        discriminant_three,
    )
    for dyadic, p7, discriminant_three in product((False, True), repeat=3)
)


@dataclass(frozen=True, slots=True)
class Control:
    gate: NormGate
    height: K2
    time: K2
    target: K2


def integer_sqrt_ceiling(value: int) -> int:
    root = isqrt(value)
    return root if root * root == value else root + 1


def embedding_integer_upper(value: K2, place: int) -> int:
    if value.sign(place) < 0:
        raise ValueError("negative embedding cap")
    upper = 0 if value == ZERO else 1
    while (K2.make(upper) - value).sign(place) < 0:
        upper *= 2
    return upper


def bounded_time_controls(gate: NormGate) -> tuple[K2, ...]:
    """Independent Minkowski-box exhaustion for one mirror norm."""

    q = gate.norm
    cap0 = q * HEIGHT_BOUND / RHO
    cap1 = q / (-RHO)
    upper0 = embedding_integer_upper(cap0, 0)
    upper1 = embedding_integer_upper(cap1, 1)
    coefficient_bound = (
        integer_sqrt_ceiling(upper0) + integer_sqrt_ceiling(upper1) + 1
    ) // 2
    output = []
    for a in range(-coefficient_bound, coefficient_bound + 1):
        for b in range(-coefficient_bound, coefficient_bound + 1):
            value = K2.make(a, b)
            square = value * value
            if (square - cap0).sign(0) > 0 or (square - cap1).sign(1) > 0:
                continue
            output.append(value)
    return tuple(output)


def enumerate_controls() -> tuple[Control, ...]:
    output: set[Control] = set()
    for gate in NORM_GATES:
        q = gate.norm
        for time in bounded_time_controls(gate):
            if time.sign(0) <= 0 or not (2 * RHO * time / q).is_integral:
                continue
            target = q + RHO * time * time
            height = RHO * time * time / q
            if target.sign(1) < 0 or (height - HEIGHT_BOUND).sign(0) > 0:
                continue
            output.add(Control(gate, height, time, target))
    return tuple(
        sorted(output, key=lambda row: (row.height, row.gate.norm, row.time))
    )


@dataclass(frozen=True, slots=True)
class GramTriple:
    u2: int
    v2: int
    uv: int

    def square(self, a: int, b: int) -> int:
        return a * a * self.u2 + b * b * self.v2 + 2 * a * b * self.uv


def reduced_target(control: Control) -> K2:
    return control.target / (RHO * RHO) if control.gate.p7 else control.target


def necessary_triples(control: Control) -> tuple[GramTriple, ...]:
    target = reduced_target(control)
    if not target.is_integral or (target.b / 2).denominator != 1:
        return ()
    rational = int(target.a)
    dot = int(target.b / 2)
    minimum = 12 if control.gate.discriminant_three else 8
    output = []
    for v2 in range(rational // 2 + 1):
        u2 = rational - 2 * v2
        if u2 < 0 or u2 % 2 or v2 % 2:
            continue
        if (u2 and u2 < minimum) or (v2 and v2 < minimum):
            continue
        if dot * dot > u2 * v2:
            continue
        if control.gate.discriminant_three and any(
            value % 3 for value in (u2, v2, dot)
        ):
            continue
        triple = GramTriple(u2, v2, dot)
        if any(
            0 < triple.square(a, b) < minimum
            for a, b in product(range(-3, 4), repeat=2)
            if a or b
        ):
            continue
        output.append(triple)
    return tuple(output)


def height_zero_triples(gate: NormGate) -> tuple[GramTriple, ...]:
    row = Control(gate, ZERO, ZERO, gate.norm)
    return necessary_triples(row)


@dataclass(frozen=True, slots=True)
class AugmentedMaximum:
    maximum: K2
    pair: tuple[int, int]
    determinant: int
    retained: int


def augmented_maximum(triple: GramTriple) -> AugmentedMaximum:
    """Independent exhaustive integral extension by a norm-eight vector."""

    best: K2 | None = None
    best_pair = (0, 0)
    best_determinant = 0
    retained = 0
    p_bound = isqrt(8 * triple.u2)
    q_bound = isqrt(8 * triple.v2)
    for p in range(-p_bound, p_bound + 1):
        for q in range(-q_bound, q_bound + 1):
            determinant = (
                8 * (triple.u2 * triple.v2 - triple.uv * triple.uv)
                - triple.v2 * p * p
                - triple.u2 * q * q
                + 2 * triple.uv * p * q
            )
            if determinant < 0:
                continue
            augmented = (
                (8, p, q),
                (p, triple.u2, triple.uv),
                (q, triple.uv, triple.v2),
            )
            if any(
                0
                < sum(
                    coeff[i] * augmented[i][j] * coeff[j]
                    for i in range(3)
                    for j in range(3)
                )
                < 8
                for coeff in product(range(-3, 4), repeat=3)
                if coeff != (0, 0, 0)
            ):
                continue
            retained += 1
            value = K2.make(p, q)
            if best is None or best < value:
                best = value
                best_pair = (p, q)
                best_determinant = determinant
    if best is None:
        raise AssertionError("augmented extension set is empty")
    return AugmentedMaximum(best, best_pair, best_determinant, retained)


def chronology_digest(rows: tuple[Control, ...]) -> str:
    stream = "".join(
        "|".join(
            (
                str(row.height.a),
                str(row.height.b),
                str(row.gate.norm.a),
                str(row.gate.norm.b),
                str(row.time.a),
                str(row.time.b),
                str(row.target.a),
                str(row.target.b),
                str(int(row.gate.dyadic)),
                str(int(row.gate.p7)),
                str(int(row.gate.discriminant_three)),
            )
        )
        + "\n"
        for row in rows
    )
    return hashlib.sha256(stream.encode("ascii")).hexdigest()


def audit(production_path: Path = PRODUCTION_ARTIFACT) -> dict[str, Any]:
    gram, gram_digest, owner_source_digest = parse_literal_gram()
    if gram != gram.T or any(int(gram[i, i]) % 2 for i in range(DIMENSION)):
        raise AssertionError("literal Gram is not even integral symmetric")
    determinant = int(gram.det())
    smith = tuple(
        abs(int(value))
        for value in smith_normal_form(gram, domain=ZZ).diagonal()
    )
    rank3 = rank_mod_prime(gram, 3)
    if (
        determinant != EXPECTED_DETERMINANT
        or smith != EXPECTED_SMITH
        or rank3 != 2
    ):
        raise AssertionError("literal Gram determinant/Smith/rank audit failed")

    s3_basis = pairing_radical_basis_mod_three(gram)
    s3_index = abs(int(s3_basis.det()))
    if s3_index != 9:
        raise AssertionError(f"S3 preimage index is {s3_index}, not 9")
    pari = pari_shortest_audit(gram, s3_basis)

    if RHO.norm() != -7 or P7 != K2.make(3, 1) or P7.norm() != 7:
        raise AssertionError("prime-seven factorization failed")
    if len(NORM_GATES) != 8 or len({gate.norm for gate in NORM_GATES}) != 8:
        raise AssertionError("the local norm alphabet is not eight distinct classes")
    expected_norms = (
        K2.make(1),
        K2.make(3),
        K2.make(3, 1),
        K2.make(9, 3),
        K2.make(2),
        K2.make(6),
        K2.make(6, 2),
        K2.make(18, 6),
    )
    if tuple(gate.norm for gate in NORM_GATES) != expected_norms:
        raise AssertionError("the complete norm representatives changed")
    if any(height_zero_triples(gate) for gate in NORM_GATES):
        raise AssertionError("a height-zero necessary Gram type appeared")

    rows = enumerate_controls()
    digest = chronology_digest(rows)
    # This digest is over all 51 rows, not a production count or endpoint list.
    expected_digest = (
        "ff6d4b175d8c9268ec565e4ccf06e979342b6cf55dc325c91bd6ba8574bdcedd"
    )
    if digest != expected_digest:
        raise AssertionError(f"full control chronology digest changed: {digest}")
    if len(rows) != 51 or len({row.height for row in rows}) != 35:
        raise AssertionError("the independent control chronology changed")
    # The time-pairing gate must independently force x in 3O for every q3 row.
    if any(
        (row.time.a.denominator != 1 or row.time.b.denominator != 1)
        or int(row.time.a) % 3
        or int(row.time.b) % 3
        for row in rows
        if row.gate.discriminant_three
    ):
        raise AssertionError("a discriminant-three control escaped x in 3O")

    nonempty = tuple((row, necessary_triples(row)) for row in rows if necessary_triples(row))
    earlier = tuple(row for row, _ in nonempty if row.height < HEIGHT_BOUND)
    endpoint = tuple((row, triples) for row, triples in nonempty if row.height == HEIGHT_BOUND)
    if earlier:
        raise AssertionError(f"an earlier necessary type survived: {earlier!r}")
    expected_endpoint = (
        (K2.make(3), K2.make(9, 6), (GramTriple(138, 72, 99),)),
        (
            K2.make(6),
            K2.make(12, 9),
            (
                GramTriple(300, 132, 198),
                GramTriple(288, 138, 198),
                GramTriple(276, 144, 198),
                GramTriple(264, 150, 198),
                GramTriple(252, 156, 198),
            ),
        ),
    )
    observed_endpoint = tuple(
        (row.gate.norm, row.time, triples) for row, triples in endpoint
    )
    if observed_endpoint != expected_endpoint:
        raise AssertionError(f"six endpoint types changed: {observed_endpoint!r}")

    maxima = tuple(
        augmented_maximum(triple)
        for _row, triples in endpoint
        for triple in triples
    )
    expected_maxima = (
        AugmentedMaximum(K2.make(24, 18), (24, 18), 432, 169),
        AugmentedMaximum(K2.make(45, 30), (45, 30), 468, 361),
        AugmentedMaximum(K2.make(44, 31), (44, 31), 528, 481),
        AugmentedMaximum(K2.make(44, 32), (44, 32), 480, 529),
        AugmentedMaximum(K2.make(40, 31), (40, 31), 504, 417),
        AugmentedMaximum(K2.make(33, 27), (33, 27), 108, 199),
    )
    if maxima != expected_maxima:
        raise AssertionError(f"augmented maxima changed: {maxima!r}")

    scale = K2.make(Fraction(1, 2))
    witness_norm = K2.make(8) * scale * scale - RHO
    margins = tuple(
        RHO * row.time - scale * augmented_maximum(triple).maximum
        for row, triples in endpoint
        for triple in triples
    )
    expected_margins = (
        K2.make(3, 3),
        K2.make(Fraction(3, 2)),
        K2.make(2, Fraction(-1, 2)),
        K2.make(2, -1),
        K2.make(4, Fraction(-1, 2)),
        K2.make(Fraction(15, 2), Fraction(3, 2)),
    )
    if witness_norm != K2.make(3, -2) or witness_norm.sign() <= 0:
        raise AssertionError("d/2 ceased to be hyperideal")
    if margins != expected_margins or any(margin.sign() <= 0 for margin in margins):
        raise AssertionError(f"a d/2 margin ceased to be strict: {margins!r}")
    # Use the literal minimum vector d=e0, not an existential shell choice.
    d = Matrix([1] + [0] * 27)
    if int((d.T * gram * d)[0]) != 8:
        raise AssertionError("the explicit witness direction d=e0 lost norm eight")

    production = json.loads(production_path.read_text())
    production_components = production["controls"]["endpoint_components"]
    flattened = [row for component in production_components for row in component["necessary_gram_types"]]
    if production["source"] != CATALOGUE_URL:
        raise AssertionError("production primary source URL drifted")
    if production["spatial_lattice"]["smith_invariants"] != list(smith):
        raise AssertionError("production Smith form disagrees with replay")
    if production["mirror_norms"]["complete_classes"] != [str(value) for value in expected_norms]:
        raise AssertionError("production norm alphabet disagrees with replay")
    if production["controls"]["count"] != len(rows) or production["controls"]["distinct_height_count"] != 35:
        raise AssertionError("production control counts disagree with replay")
    if [row["gram"] for row in flattened] != [
        [triple.u2, triple.v2, triple.uv]
        for _row, triples in endpoint
        for triple in triples
    ]:
        raise AssertionError("production endpoint triples disagree with replay")
    if [row["maximum_d_dot_y"] for row in flattened] != [str(row.maximum) for row in maxima]:
        raise AssertionError("production augmented maxima disagree with replay")
    if [row["witness_margin"] for row in flattened] != [str(value) for value in margins]:
        raise AssertionError("production witness margins disagree with replay")
    if production["witness"]["lorentz_norm"] != str(witness_norm):
        raise AssertionError("production witness norm disagrees with replay")

    return {
        "schema_version": "0.1.0",
        "verdict": "PASS",
        "scope": (
            "Owner-import-free exact hostile replay of the complete necessary "
            "prefix through H=93+66sqrt(2); no later-wall or termination claim."
        ),
        "independence": (
            "The owner module is never imported. Its Gram literal is parsed as "
            "inert AST data; all arithmetic, modular kernels, controls, gates, "
            "Gram screens, augmented maxima, and witness inequalities are "
            "implemented independently in this script."
        ),
        "sources": {
            "catalogue_url": CATALOGUE_URL,
            "owner_source": str(OWNER_SOURCE.relative_to(ROOT)),
            "owner_source_sha256": owner_source_digest,
            "production_artifact": str(production_path.relative_to(ROOT)),
            "production_artifact_sha256": hashlib.sha256(
                production_path.read_bytes()
            ).hexdigest(),
        },
        "literal_gram": {
            "normalized_stream_sha256": gram_digest,
            "dimension": DIMENSION,
            "even_integral_symmetric": True,
            "determinant": determinant,
            "smith_invariants": list(smith),
            "rank_mod_three": rank3,
            "minimum": pari["host_minimum"],
            "minimum_shell_signed_count": pari["host_signed_shortest_count"],
        },
        "pairing_radical_S3": {
            "construction": "independent GF(3) row reduction, lifts plus 3I, column HNF",
            "preimage_index": s3_index,
            "minimum": pari["S3_minimum"],
            "minimum_shell_signed_count": pari["S3_signed_shortest_count"],
            "gram_pairings_divisible_by_three": True,
            "pari": pari,
        },
        "local_norm_gates": {
            "class_count": len(NORM_GATES),
            "classes": [
                {
                    "norm": str(gate.norm),
                    "dyadic_valuation": 2 if gate.dyadic else 0,
                    "p7_spatial_gate": "y in rho*L" if gate.p7 else None,
                    "discriminant_three_gate": (
                        "x in 3O and coefficient vectors u,v in S3"
                        if gate.discriminant_three
                        else None
                    ),
                }
                for gate in NORM_GATES
            ],
            "height_zero_necessary_type_count": 0,
            "completeness": (
                "Odd determinant makes L even-unimodular at pi, giving dyadic "
                "valuations 0 or 2. Smith exponents at 3 are all at most one, "
                "and rho is squarefree; their valuations are independently 0 "
                "or 1. Local unimodularity and primitivity exclude every other prime."
            ),
        },
        "chronology": {
            "inclusive_height_bound": str(HEIGHT_BOUND),
            "control_count": len(rows),
            "distinct_height_count": len({row.height for row in rows}),
            "full_control_stream_sha256": digest,
            "earlier_nonempty_necessary_controls": len(earlier),
            "endpoint_component_count": len(endpoint),
            "endpoint_type_count": sum(len(triples) for _row, triples in endpoint),
            "endpoint": [
                {
                    "norm": str(row.gate.norm),
                    "time": str(row.time),
                    "types": [
                        {
                            "gram": [triple.u2, triple.v2, triple.uv],
                            "maximum": str(maximum.maximum),
                            "maximizing_pair": list(maximum.pair),
                            "augmented_determinant": maximum.determinant,
                            "retained_extensions": maximum.retained,
                            "strict_margin": str(margin),
                            "strict_margin_exact": margin.exact(),
                        }
                        for triple, maximum, margin in zip(
                            triples,
                            maxima[
                                sum(len(packet) for _prior, packet in endpoint[:index]) :
                                sum(len(packet) for _prior, packet in endpoint[: index + 1])
                            ],
                            margins[
                                sum(len(packet) for _prior, packet in endpoint[:index]) :
                                sum(len(packet) for _prior, packet in endpoint[: index + 1])
                            ],
                            strict=True,
                        )
                    ],
                }
                for index, (row, triples) in enumerate(endpoint)
            ],
        },
        "strict_hyperideal_witness": {
            "explicit_direction": "d=e0 in the literal Gram basis, d^2=8",
            "section_vector": "z=(1,d/2)",
            "norm": str(witness_norm),
            "norm_exact": witness_norm.exact(),
            "all_six_margins_strict": True,
            "bounded_prefix_finite_volume": False,
        },
        "production_comparison": {
            "all_requested_claims_match": True,
            "defects_found": [],
        },
        "conclusion": {
            "bounded_prefix_has_open_hyperideal_patch": True,
            "compact_or_finite_volume_example": False,
            "eventual_chamber": "UNRESOLVED",
        },
    }


def render_report(payload: dict[str, Any]) -> str:
    chronology = payload["chronology"]
    rows = []
    for component in chronology["endpoint"]:
        for row in component["types"]:
            rows.append(
                f"| {component['norm']} | {component['time']} | "
                f"{tuple(row['gram'])} | {row['maximum']} | "
                f"{row['strict_margin']} |"
            )
    return f"""# Independent hostile audit: silver `DualExtremal(28,3)`

**PASS.** This replay never imports the production owner module. It parses the
Gram literal as inert source data and independently reconstructs every claim
needed for the bounded prefix through `H=93+66sqrt(2)`.

The 28-dimensional Gram is even integral symmetric, has determinant `3^26`,
Smith form `1^2,3^26`, rank two modulo three, minimum 8, and 1,260 signed
minimal vectors. An independent GF(3)-kernel lift and column-HNF computation
gives the pairing-radical preimage index 9; PARI finds its exact minimum 12
and 65,520 signed minimal vectors.

The local replay gives exactly eight dyadic/prime-seven/discriminant-three
norm classes. It enumerates 51 controls on 35 heights. Every earlier control
has an empty rigorous necessary Gram superset. The endpoint contains exactly
six types:

| norm | time | Gram `(u²,v²,u.v)` | augmented maximum | `d/2` margin |
|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

Here the witness is explicit: `d=e0`, so `d²=8`, and `z=(1,d/2)` has
`z²=3-2sqrt(2)>0`. All six margins are strictly positive. Therefore the
complete necessary prefix has an open hyperideal patch and infinite volume.

No production defect was found. This remains a bounded negative result: it
does not determine later walls, Vinberg termination, or the eventual chamber.
"""


def encoded(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--production", type=Path, default=PRODUCTION_ARTIFACT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--print-digest", action="store_true")
    args = parser.parse_args()
    payload = audit(args.production)
    if args.print_digest:
        print(payload["chronology"]["full_control_stream_sha256"])
        return
    artifact = encoded(payload)
    report = render_report(payload).encode()
    if args.verify:
        if args.output.read_bytes() != artifact or args.report.read_bytes() != report:
            raise SystemExit("frozen hostile audit differs from exact replay")
        print(f"PASS {hashlib.sha256(artifact).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(artifact)
    args.report.write_bytes(report)
    print(f"PASS {hashlib.sha256(artifact).hexdigest()}")


if __name__ == "__main__":
    main()
