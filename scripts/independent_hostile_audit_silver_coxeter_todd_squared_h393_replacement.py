#!/usr/bin/env python3
"""Clean-room hostile replay of the silver K12^2 H393 certificate.

This file deliberately imports no project code.  It reconstructs arithmetic
in ``Z[sqrt(2)]``, the catalogue K12 Gram matrix, all sixteen local norm
classes, the complete control/necessary-Gram chronology, the dual-coordinate
gate, the explicit cutter, every required cross-Gram acuteness bound, the
same-height norm-two lemma, and the 7/4 replacement witness.

The binary and cross-Gram gates are necessary supersets of actual K12 data.
Therefore their maxima are safe upper bounds: nonpositive cutter pairings and
strictly negative replacement pairings proved here apply to every actual wall.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd, isqrt
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
OWNER_MODULE = ROOT / "src/coxeter_polytopes/search/silver_coxeter_todd_squared.py"
OWNER_SCRIPT = ROOT / "scripts/audit_silver_coxeter_todd_squared_h393_replacement.py"
OWNER_JSON = ROOT / "data/enumerations/silver_coxeter_todd_squared_h393_replacement.json"
AUDIT_JSON = (
    ROOT
    / "data/enumerations/silver_coxeter_todd_squared_h393_independent_hostile_audit.json"
)
AUDIT_REPORT = (
    ROOT / "results/silver_coxeter_todd_squared_h393_independent_hostile_audit.md"
)
AUDIT_SHA = (
    ROOT
    / "data/enumerations/silver_coxeter_todd_squared_h393_independent_hostile_audit.sha256"
)


def _f(value: int | Fraction | str) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


def _sgn(value: Fraction) -> int:
    return (value > 0) - (value < 0)


@dataclass(frozen=True, slots=True)
class Q2:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", _f(self.a))
        object.__setattr__(self, "b", _f(self.b))

    @staticmethod
    def coerce(value: "Q2 | int | Fraction") -> "Q2":
        return value if isinstance(value, Q2) else Q2(value)

    def __add__(self, other: "Q2 | int | Fraction") -> "Q2":
        other = Q2.coerce(other)
        return Q2(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> "Q2":
        return Q2(-self.a, -self.b)

    def __sub__(self, other: "Q2 | int | Fraction") -> "Q2":
        return self + (-Q2.coerce(other))

    def __rsub__(self, other: "Q2 | int | Fraction") -> "Q2":
        return Q2.coerce(other) - self

    def __mul__(self, other: "Q2 | int | Fraction") -> "Q2":
        other = Q2.coerce(other)
        return Q2(
            self.a * other.a + 2 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    __rmul__ = __mul__

    def inverse(self) -> "Q2":
        norm = self.a * self.a - 2 * self.b * self.b
        if not norm:
            raise ZeroDivisionError
        return Q2(self.a / norm, -self.b / norm)

    def __truediv__(self, other: "Q2 | int | Fraction") -> "Q2":
        return self * Q2.coerce(other).inverse()

    def __rtruediv__(self, other: "Q2 | int | Fraction") -> "Q2":
        return Q2.coerce(other) / self

    def __pow__(self, exponent: int) -> "Q2":
        if exponent < 0:
            return self.inverse() ** (-exponent)
        result = ONE
        base = self
        power = exponent
        while power:
            if power & 1:
                result *= base
            base *= base
            power >>= 1
        return result

    def conjugate(self) -> "Q2":
        return Q2(self.a, -self.b)

    def sign(self, place: int = 0) -> int:
        if place == 1:
            return self.conjugate().sign(0)
        if place != 0:
            raise ValueError("place must be zero or one")
        sa, sb = _sgn(self.a), _sgn(self.b)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        return sa if self.a * self.a - 2 * self.b * self.b > 0 else sb

    def __lt__(self, other: "Q2 | int | Fraction") -> bool:
        return (self - other).sign() < 0

    def __le__(self, other: "Q2 | int | Fraction") -> bool:
        return (self - other).sign() <= 0

    def __gt__(self, other: "Q2 | int | Fraction") -> bool:
        return (self - other).sign() > 0

    def __ge__(self, other: "Q2 | int | Fraction") -> bool:
        return (self - other).sign() >= 0

    @property
    def integral(self) -> bool:
        return self.a.denominator == self.b.denominator == 1

    def display(self) -> str:
        def q(value: Fraction) -> str:
            return str(value.numerator) if value.denominator == 1 else str(value)

        if not self.b:
            return q(self.a)
        if not self.a:
            return f"{q(self.b)}*sqrt(2)"
        sign = "+" if self.b > 0 else "-"
        return f"{q(self.a)} {sign} {q(abs(self.b))}*sqrt(2)"

    def json(self) -> dict[str, object]:
        return {
            "exact": {"a": str(self.a), "b": str(self.b)},
            "display": self.display(),
            "signs": [self.sign(0), self.sign(1)],
        }


ZERO = Q2(0)
ONE = Q2(1)
R = Q2(0, 1)
A = Q2(-13, 12)
P7 = Q2(3, 1)
P17 = Q2(5, -2)
FIRST_HEIGHT = Q2(67, 48)
CUT_HEIGHT = Q2(393, 278)
CUT_CONTROL = Q2(10, 7)
OLD_SCALE = Fraction(9, 5)
NEW_SCALE = Fraction(7, 4)


K12_GRAM: tuple[tuple[int, ...], ...] = (
    (4, 0, 0, -2, 0, 0, 2, -1, -1, -1, 2, -1),
    (0, 4, 0, 0, -2, 0, 2, -1, -1, -1, -1, 2),
    (0, 0, 4, 0, 0, -2, 2, 2, 2, -1, -1, -1),
    (-2, 0, 0, 4, 0, 0, -1, -1, 2, 2, -1, -1),
    (0, -2, 0, 0, 4, 0, -1, 2, -1, 2, -1, -1),
    (0, 0, -2, 0, 0, 4, -1, -1, -1, 2, 2, 2),
    (2, 2, 2, -1, -1, -1, 4, 0, 0, -2, 0, 0),
    (-1, -1, 2, -1, 2, -1, 0, 4, 0, 0, -2, 0),
    (-1, -1, 2, 2, -1, -1, 0, 0, 4, 0, 0, -2),
    (-1, -1, -1, 2, 2, 2, -2, 0, 0, 4, 0, 0),
    (2, -1, -1, -1, -1, 2, 0, -2, 0, 0, 4, 0),
    (-1, 2, -1, -1, -1, 2, 0, 0, -2, 0, 0, 4),
)

H: tuple[Fraction, ...] = (
    Fraction(4, 3), 0, 0, Fraction(2, 3), 0, 0,
    Fraction(-2, 3), Fraction(1, 3), Fraction(1, 3),
    Fraction(-1, 3), Fraction(-1, 3), Fraction(2, 3),
)

D_VECTOR = (6, 0, 0, 2, 0, 0, -3, 1, 2, -1, -2, 3)
E_VECTOR = (10, 1, 0, 4, -1, 0, -6, 3, 2, -2, -3, 3)
U_VECTOR = tuple(2 * d + e for d, e in zip(D_VECTOR, E_VECTOR, strict=True))
V_VECTOR = tuple(d + e for d, e in zip(D_VECTOR, E_VECTOR, strict=True))


def pair(left: Iterable[int | Fraction], right: Iterable[int | Fraction]) -> Fraction:
    left, right = tuple(map(Fraction, left)), tuple(map(Fraction, right))
    if len(left) != 12 or len(right) != 12:
        raise ValueError("K12 vectors must have length twelve")
    return sum(
        left[i] * K12_GRAM[i][j] * right[j]
        for i in range(12)
        for j in range(12)
    )


def bareiss(matrix: tuple[tuple[int, ...], ...]) -> int:
    work = [list(row) for row in matrix]
    previous = 1
    sign = 1
    for col in range(len(work) - 1):
        pivot = next((i for i in range(col, len(work)) if work[i][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            sign = -sign
        value = work[col][col]
        for i in range(col + 1, len(work)):
            for j in range(col + 1, len(work)):
                numerator = work[i][j] * value - work[i][col] * work[col][j]
                if numerator % previous:
                    raise AssertionError("nonexact Bareiss step")
                work[i][j] = numerator // previous
        previous = value
        for i in range(col + 1, len(work)):
            work[i][col] = 0
    return sign * work[-1][-1]


@dataclass(frozen=True, slots=True)
class NormClass:
    norm: Q2
    ternary: bool
    p7: bool
    p17: bool

    @property
    def odd_divisor(self) -> Q2:
        return (P7 if self.p7 else ONE) * (P17 if self.p17 else ONE)


NORM_CLASSES = tuple(
    NormClass(
        Q2(2 if dyadic else 1)
        * (3 if ternary else 1)
        * (P7 if p7 else ONE)
        * (P17 if p17 else ONE),
        ternary,
        p7,
        p17,
    )
    for dyadic, ternary, p7, p17 in product((False, True), repeat=4)
)


def ceil_sqrt_at(value: Q2, place: int) -> int:
    if value.sign(place) < 0:
        raise ValueError("negative square cap")
    high = 1
    while (Q2(high * high) - value).sign(place) < 0:
        high *= 2
    low = 0
    while low + 1 < high:
        mid = (low + high) // 2
        if (Q2(mid * mid) - value).sign(place) >= 0:
            high = mid
        else:
            low = mid
    return high


def bounded_ring_integers(cap0: Q2, cap1: Q2) -> tuple[Q2, ...]:
    b0 = ceil_sqrt_at(cap0, 0)
    b1 = ceil_sqrt_at(cap1, 1)
    bound = (b0 + b1 + 1) // 2 + 1
    output = []
    for m in range(-bound, bound + 1):
        for n in range(-bound, bound + 1):
            x = Q2(m, n)
            square = x * x
            if (cap0 - square).sign(0) >= 0 and (cap1 - square).sign(1) >= 0:
                output.append(x)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class Control:
    norm_class: NormClass
    height: Q2
    x: Q2
    target: Q2


@lru_cache(maxsize=1)
def controls() -> tuple[Control, ...]:
    output: set[Control] = set()
    for norm_class in NORM_CLASSES:
        q = norm_class.norm
        for x in bounded_ring_integers(q * CUT_HEIGHT / A, q / (-A)):
            if x.sign(0) <= 0 or not (2 * A * x / q).integral:
                continue
            target = q + A * x * x
            height = A * x * x / q
            if target.sign(1) < 0 or height > CUT_HEIGHT:
                continue
            output.add(Control(norm_class, height, x, target))
    return tuple(sorted(output, key=lambda row: (row.height, row.norm_class.norm, row.x)))


def has_combination_below(triple: tuple[int, int, int], minimum: int) -> bool:
    s, t, dot = triple
    determinant = s * t - dot * dot
    if determinant < 0:
        return True
    if determinant == 0:
        return gcd(gcd(s, t), abs(dot)) < minimum
    ib = isqrt(((minimum - 1) * t) // determinant) + 1
    jb = isqrt(((minimum - 1) * s) // determinant) + 1
    return any(
        0 < i * i * s + j * j * t + 2 * i * j * dot < minimum
        for i in range(-ib, ib + 1)
        for j in range(-jb, jb + 1)
        if i or j
    )


@lru_cache(maxsize=None)
def gram_types(control: Control) -> tuple[tuple[int, int, int], ...]:
    divisor = control.norm_class.odd_divisor
    target = control.target / (divisor * divisor)
    if not target.integral or (target.b / 2).denominator != 1:
        return ()
    rational = int(target.a)
    dot = int(target.b / 2)
    minimum = 12 if control.norm_class.ternary else 4
    output = []
    for v_square in range(0, rational // 2 + 1, 2):
        u_square = rational - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if u_square not in (0,) and u_square < 4:
            continue
        if v_square not in (0,) and v_square < 4:
            continue
        if dot * dot > u_square * v_square:
            continue
        triple = (u_square, v_square, dot)
        if control.norm_class.ternary and any(value % 3 for value in triple):
            continue
        if has_combination_below(triple, minimum):
            continue
        output.append(triple)
    return tuple(output)


RESIDUAL = (Fraction(0), Fraction(13, 4), Fraction(1), Fraction(13, 4))


def coordinate_minimum(k: int) -> Fraction:
    return Fraction(3 * k * k, 4) + RESIDUAL[k % 4]


@lru_cache(maxsize=None)
def directional_pairs(triple: tuple[int, int, int]) -> tuple[tuple[int, int], ...]:
    s, t, dot = triple
    output = []
    for p in range(-isqrt((4 * s) // 3) - 2, isqrt((4 * s) // 3) + 3):
        s4 = 4 * s - 3 * p * p
        if s4 < 0:
            continue
        for q in range(-isqrt((4 * t) // 3) - 2, isqrt((4 * t) // 3) + 3):
            t4 = 4 * t - 3 * q * q
            cross4 = 4 * dot - 3 * p * q
            if t4 < 0 or s4 * t4 < cross4 * cross4:
                continue
            if any(
                Fraction(i * i * s + j * j * t + 2 * i * j * dot)
                < coordinate_minimum(i * p + j * q)
                for i in (-1, 0, 1)
                for j in (-1, 0, 1)
                if i or j
            ):
                continue
            output.append((p, q))
    return tuple(output)


@dataclass(frozen=True, slots=True)
class DirectionBound:
    control: Control
    type_count: int
    pair_count: int
    maximum: Q2
    maximizing_type: tuple[int, int, int]
    maximizing_pair: tuple[int, int]


def direction_bound(control: Control) -> DirectionBound:
    rows = [
        (control.norm_class.odd_divisor * Q2(p, q), triple, (p, q))
        for triple in gram_types(control)
        for p, q in directional_pairs(triple)
    ]
    if not rows:
        raise AssertionError("empty direction bound")
    maximum, triple, pq = max(rows, key=lambda row: row[0])
    return DirectionBound(control, len(gram_types(control)), len(rows), maximum, triple, pq)


def transformed_gram(triple: tuple[int, int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    s, t, dot = triple
    return (
        (s + t - 2 * dot, 3 * dot - s - 2 * t),
        (3 * dot - s - 2 * t, s + 4 * t - 4 * dot),
    )


@dataclass(frozen=True, slots=True)
class CrossMaximum:
    first_count: int
    second_count: int
    feasible_count: int
    maximum: Q2
    entries: tuple[int, int, int, int]


@lru_cache(maxsize=None)
def cross_maximum(variable: tuple[tuple[int, int], tuple[int, int]]) -> CrossMaximum:
    n1, dot, n2 = variable[0][0], variable[0][1], variable[1][1]
    determinant = 236

    def numerator(left: int, right: int) -> int:
        return 84 * left * left - 92 * left * right + 28 * right * right

    def mixed(a: int, c: int, b: int, d: int) -> int:
        return 84 * a * b - 46 * (a * d + c * b) + 28 * c * d

    def columns(norm: int) -> tuple[tuple[int, int, int], ...]:
        return tuple(
            (left, right, value)
            for left in range(-isqrt(28 * norm) - 2, isqrt(28 * norm) + 3)
            for right in range(-isqrt(84 * norm) - 2, isqrt(84 * norm) + 3)
            if (value := numerator(left, right)) <= norm * determinant
        )

    first, second = columns(n1), columns(n2)
    count = 0
    best: Q2 | None = None
    entries: tuple[int, int, int, int] | None = None
    for a, c, q1 in first:
        r11 = n1 * determinant - q1
        for b, d, q2 in second:
            r22 = n2 * determinant - q2
            r12 = dot * determinant - mixed(a, c, b, d)
            if r11 * r22 < r12 * r12:
                continue
            count += 1
            value = Q2(6 * a + 4 * (b + c) + 3 * d, 4 * a + 3 * (b + c) + 2 * d)
            if best is None or value > best:
                best, entries = value, (a, b, c, d)
    if best is None or entries is None:
        raise AssertionError("empty cross-Gram census")
    return CrossMaximum(len(first), len(second), count, best, entries)


def gp_matrix(matrix: tuple[tuple[int, ...], ...]) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in matrix) + "]"


def pari_replay() -> dict[str, object]:
    executable = shutil.which("gp")
    if executable is None:
        raise RuntimeError("PARI/GP is required for the independent audit")
    gram = gp_matrix(K12_GRAM)
    u = "[" + ",".join(map(str, U_VECTOR)) + "]~"
    v = "[" + ",".join(map(str, V_VECTOR)) + "]~"
    program = (
        f"G={gram};U={u};V={v};"
        "cmp2(a,b,c,d)={my(A=a-c,B=b-d,N);if(A==0,return(sign(B)));"
        "if(B==0,return(sign(A)));if(sign(A)==sign(B),return(sign(A)));"
        "N=A*A-2*B*B;if(N>0,return(sign(A)),return(sign(B)))};\n"
        "R=qfminim(G,10,,2);W=R[3];m1=10^9;m2=10^9;m3=10^9;c1=0;c2=0;c3=0;"
        "for(j=1,matsize(W)[2],forstep(s=-1,1,2,w=s*W[,j];q=w~*G*w;k=w[1];"
        "if(k==1,if(q<m1,m1=q;c1=1,if(q==m1,c1++)));"
        "if(k==2,if(q<m2,m2=q;c2=1,if(q==m2,c2++)));"
        "if(k==3,if(q<m3,m3=q;c3=1,if(q==m3,c3++)))));"
        'print("C|",R[1],"|",matsize(W)[2],"|",m1,"|",c1,"|",m2,"|",c2,"|",m3,"|",c3);'
        "S=qfminim(G,8,,2);W=S[3];a4=-10^9;b4=0;a8=-10^9;b8=0;n8=0;"
        "for(j=1,matsize(W)[2],forstep(s=-1,1,2,w=s*W[,j];q=w~*G*w;"
        "X=U~*G*w;Y=V~*G*w;if(q==4&&cmp2(X,Y,a4,b4)>0,a4=X;b4=Y);"
        "if(q<=8&&cmp2(X,Y,a8,b8)>0,a8=X;b8=Y;n8=q)));"
        'print("M4|",a4,"|",b4);print("M8|",a8,"|",b8,"|",n8);quit()'
    )
    completed = subprocess.run(
        [executable, "-q", "-s", "1500000000"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=240,
        check=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    expected = [
        "C|85680|42840|4|160|4|2|10|160",
        "M4|36|26",
        "M8|40|30|8",
    ]
    if lines != expected:
        raise AssertionError(
            f"independent PARI transcript drifted: {lines!r}; "
            f"stderr={completed.stderr.strip()!r}"
        )
    return {
        "transcript": lines,
        "coset_residue_minima": [0, 4, 4, 10],
        "first_packet_spatial_maxima": [Q2(36, 26).json(), Q2(40, 30).json()],
    }


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stream_hash(rows: Iterable[str]) -> str:
    return hashlib.sha256(("\n".join(rows) + "\n").encode()).hexdigest()


def build_payload() -> dict[str, object]:
    if bareiss(K12_GRAM) != 729:
        raise AssertionError("K12 determinant is not 729")
    gh = tuple(sum(Fraction(K12_GRAM[i][j]) * H[j] for j in range(12)) for i in range(12))
    if gh != (Fraction(1),) + (Fraction(0),) * 11 or pair(H, H) != Fraction(4, 3):
        raise AssertionError("dual direction reconstruction failed")
    if not all((3 * value).denominator == 1 for value in H):
        raise AssertionError("3h is not integral")
    if 3 * H[0] != 4 or pair(tuple(3 * value for value in H), tuple(3 * value for value in H)) != 12:
        raise AssertionError("dual-coordinate translation failed")

    all_controls = controls()
    nonempty = tuple(row for row in all_controls if gram_types(row))
    frontier = tuple((row.height, row.norm_class.norm, len(gram_types(row))) for row in nonempty)
    expected_frontier = (
        (FIRST_HEIGHT, Q2(1), 1),
        (FIRST_HEIGHT, Q2(2), 2),
        (Q2(356, 252), Q2(6, 2), 2),
        (Q2(362, 256), Q2(10, -4), 1),
        (CUT_HEIGHT, Q2(1), 4),
        (CUT_HEIGHT, Q2(2), 11),
    )
    if len(all_controls) != 128 or len({row.height for row in all_controls}) != 87:
        raise AssertionError("independent control chronology differs")
    if frontier != expected_frontier:
        raise AssertionError(f"independent necessary frontier differs: {frontier!r}")

    direction_rows = tuple(direction_bound(row) for row in nonempty)
    expected_pairs = (5, 10, 42, 19, 140, 715)
    expected_maxima = (
        Q2(6, 4), Q2(8, 6), Q2(36, 26), Q2(30, 22), Q2(15, 11), Q2(22, 16)
    )
    if tuple(row.pair_count for row in direction_rows) != expected_pairs:
        raise AssertionError("directional-pair counts differ")
    if tuple(row.maximum for row in direction_rows) != expected_maxima:
        raise AssertionError("directional maxima differ")

    old_margins = tuple(OLD_SCALE * row.maximum - A * row.control.x for row in direction_rows)
    new_margins = tuple(NEW_SCALE * row.maximum - A * row.control.x for row in direction_rows)
    if tuple(value.sign() for value in old_margins[2:]) != (-1, -1, -1, 1):
        raise AssertionError("old-witness cutting chronology differs")
    if any(value.sign() >= 0 for value in new_margins):
        raise AssertionError("7/4 replacement is not strict")

    # Explicit cutter, reconstructed from D,E rather than owner root data.
    if (pair(D_VECTOR, D_VECTOR), pair(E_VECTOR, E_VECTOR), pair(D_VECTOR, E_VECTOR)) != (28, 84, 46):
        raise AssertionError("explicit D,E Gram differs")
    coefficient_gram = (pair(U_VECTOR, U_VECTOR), pair(V_VECTOR, V_VECTOR), pair(U_VECTOR, V_VECTOR))
    if coefficient_gram != (380, 204, 278):
        raise AssertionError("explicit U,V Gram differs")
    root_norm = -A * CUT_CONTROL * CUT_CONTROL + Q2(coefficient_gram[0] + 2 * coefficient_gram[1], 2 * coefficient_gram[2])
    if root_norm != Q2(2) or A * CUT_CONTROL * CUT_CONTROL / 2 != CUT_HEIGHT:
        raise AssertionError("explicit cutter norm/height differs")
    spatial_pairings = tuple(
        sum((Q2(U_VECTOR[i], V_VECTOR[i]) * K12_GRAM[i][j] for i in range(12)), ZERO)
        for j in range(12)
    )
    reflection_integral = (A * CUT_CONTROL).integral and all(value.integral for value in spatial_pairings)
    primitive = any(value % 2 for value in U_VECTOR)
    if not reflection_integral or not primitive:
        raise AssertionError("explicit cutter reflection/primitivity gate failed")
    q2_types = gram_types(nonempty[-1])
    if (380, 204, 278) not in q2_types:
        raise AssertionError("explicit cutter type is absent from the complete q=2 row")

    pari = pari_replay()
    max4, max8 = Q2(36, 26), Q2(40, 30)
    eps2 = Q2(3, 2)
    first_packet_margins = (
        eps2 * (max4 - A * CUT_CONTROL),
        eps2 * (max8 - R * A * CUT_CONTROL),
    )
    if first_packet_margins != (Q2(-18, -13), Q2(-86, -60)):
        raise AssertionError("first-packet cutter margins differ")

    def lookup(height: Q2, norm: Q2) -> Control:
        found = [row for row in all_controls if row.height == height and row.norm_class.norm == norm]
        if len(found) != 1:
            raise AssertionError("ambiguous control lookup")
        return found[0]

    cross_rows = []
    for height, norm in ((Q2(356, 252), Q2(6, 2)), (Q2(362, 256), Q2(10, -4))):
        row = lookup(height, norm)
        for triple in gram_types(row):
            transformed = transformed_gram(triple)
            maximum = cross_maximum(transformed)
            margin = row.norm_class.odd_divisor * maximum.maximum - A * CUT_CONTROL * row.x
            cross_rows.append((row, triple, transformed, maximum, margin))
    if tuple(row[-1] for row in cross_rows) != (Q2(0), Q2(0), Q2(-7, -4)):
        raise AssertionError("intermediate cutter acuteness differs")

    q1 = lookup(CUT_HEIGHT, Q2(1))
    same_q1 = []
    for triple in gram_types(q1):
        maximum = cross_maximum(transformed_gram(triple))
        margin = maximum.maximum - A * CUT_CONTROL * q1.x
        same_q1.append((triple, maximum, margin))
    if tuple(row[-1] for row in same_q1) != (Q2(-2, -1), Q2(-1, -1), Q2(-1, -1), Q2(-1, -1)):
        raise AssertionError("same-height q=1 acuteness differs")

    # Same-height q=2 lemma: enumerate the entire Minkowski box exactly.
    positive_equal_height = tuple(
        Q2(m, n)
        for m in range(-4, 5)
        for n in range(-4, 5)
        if Q2(m, n).sign(0) > 0
        and (Q2(2) - Q2(m, n)).sign(0) >= 0
        and (Q2(2) - Q2(m, n).conjugate()).sign(0) >= 0
        and (Q2(2) + Q2(m, n).conjugate()).sign(0) >= 0
    )
    if positive_equal_height != (R, Q2(1), Q2(2)) and set(positive_equal_height) != {Q2(1), R, Q2(2)}:
        raise AssertionError(f"same-height Minkowski box differs: {positive_equal_height!r}")
    same_height_exclusions = {
        "1": "difference square 2 forces u^2+2v^2=2, impossible at minimum four",
        "sqrt(2)": (
            "difference square 4-2sqrt(2) forces u^2+2v^2=4 and u.v=-1; "
            "both u,v nonzero would contribute at least 12"
        ),
        "2": "difference square zero forces equal spatial vectors",
    }

    h_u, h_v = U_VECTOR[0], V_VECTOR[0]
    old_pairing = -A * CUT_CONTROL + OLD_SCALE * Q2(h_u, h_v)
    new_pairing = -A * CUT_CONTROL + NEW_SCALE * Q2(h_u, h_v)
    if old_pairing != Q2(Fraction(8, 5), Fraction(-1, 5)) or old_pairing.sign() <= 0:
        raise AssertionError("explicit cutter does not cut old witness")
    if new_pairing != Q2(Fraction(1, 2), -1) or new_pairing.sign() >= 0:
        raise AssertionError("explicit cutter does not preserve replacement")
    replacement_norm = Q2(Fraction(4, 3) * NEW_SCALE * NEW_SCALE) - A
    if replacement_norm != Q2(Fraction(205, 12), -12) or replacement_norm.sign() <= 0:
        raise AssertionError("replacement is not spacelike")

    owner = json.loads(OWNER_JSON.read_text())
    owner_summary = {
        "control_count": owner["chronology"]["control_count"],
        "distinct_height_count": owner["chronology"]["distinct_height_count"],
        "direction_counts": [row["finite_directional_pair_count"] for row in owner["chronology"]["nonempty_necessary_rows"]],
        "direction_maxima": [row["maximum_spatial_pairing"]["display"] for row in owner["chronology"]["nonempty_necessary_rows"]],
        "accepted": owner["first_cutting_wall"]["accepted"],
        "replacement_strict": owner["replacement_hyperideal_ray"]["strict_against_every_necessary_row_through_bound"],
    }
    independent_summary = {
        "control_count": len(all_controls),
        "distinct_height_count": len({row.height for row in all_controls}),
        "direction_counts": [row.pair_count for row in direction_rows],
        "direction_maxima": [row.maximum.display() for row in direction_rows],
        "accepted": True,
        "replacement_strict": True,
    }
    if owner_summary != independent_summary:
        raise AssertionError("owner artifact disagrees with clean-room summary")

    control_stream = (
        f"{row.height.display()}|{row.norm_class.norm.display()}|{row.x.display()}|"
        + ";".join(",".join(map(str, triple)) for triple in gram_types(row))
        for row in all_controls
    )
    return {
        "schema_version": "0.1",
        "audit": "independent hostile clean-room replay; no coxeter_polytopes imports",
        "verdict": "PASS",
        "bounded_scope": "complete through H=393+278*sqrt(2), not Vinberg termination",
        "source_dependencies": {
            "owner_module": {"path": str(OWNER_MODULE.relative_to(ROOT)), "sha256": sha(OWNER_MODULE)},
            "owner_freezer": {"path": str(OWNER_SCRIPT.relative_to(ROOT)), "sha256": sha(OWNER_SCRIPT)},
            "owner_json": {"path": str(OWNER_JSON.relative_to(ROOT)), "sha256": sha(OWNER_JSON)},
        },
        "k12": {
            "determinant": 729,
            "dual_vector_square": str(pair(H, H)),
            "G_times_h": [str(value) for value in gh],
            "three_h_integral": True,
            "coordinate_translation_norm": 12,
            "pari": pari,
        },
        "chronology": {
            "norm_class_count": len(NORM_CLASSES),
            "control_count": len(all_controls),
            "distinct_height_count": len({row.height for row in all_controls}),
            "through_first_packet": sum(row.height <= FIRST_HEIGHT for row in all_controls),
            "strictly_between": sum(FIRST_HEIGHT < row.height < CUT_HEIGHT for row in all_controls),
            "at_cutting_height": sum(row.height == CUT_HEIGHT for row in all_controls),
            "control_and_gram_stream_sha256": stream_hash(control_stream),
            "frontier": [
                {
                    "height": row.control.height.json(),
                    "norm": row.control.norm_class.norm.json(),
                    "control": row.control.x.json(),
                    "gram_type_count": row.type_count,
                    "directional_pair_count": row.pair_count,
                    "maximum": row.maximum.json(),
                    "maximizing_type": list(row.maximizing_type),
                    "maximizing_pair": list(row.maximizing_pair),
                    "old_margin": old_margins[index].json(),
                    "replacement_margin": new_margins[index].json(),
                }
                for index, row in enumerate(direction_rows)
            ],
        },
        "explicit_cutter": {
            "D": list(D_VECTOR),
            "E": list(E_VECTOR),
            "U": list(U_VECTOR),
            "V": list(V_VECTOR),
            "difference_gram": [28, 84, 46],
            "coefficient_gram": [int(value) for value in coefficient_gram],
            "h_pairings": [h_u, h_v],
            "norm": root_norm.json(),
            "height": CUT_HEIGHT.json(),
            "reflection_integral": reflection_integral,
            "primitive": primitive,
            "old_witness_pairing": old_pairing.json(),
            "replacement_pairing": new_pairing.json(),
            "first_packet_root_pairing_maxima": [value.json() for value in first_packet_margins],
            "intermediate_cross_gram": [
                {
                    "height": row.height.json(),
                    "norm": row.norm_class.norm.json(),
                    "gram_type": list(triple),
                    "transformed_gram": [list(value) for value in transformed],
                    "column_counts": [maximum.first_count, maximum.second_count],
                    "feasible_count": maximum.feasible_count,
                    "maximum": maximum.maximum.json(),
                    "entries": list(maximum.entries),
                    "root_pairing": margin.json(),
                }
                for row, triple, transformed, maximum, margin in cross_rows
            ],
            "same_height_q1": [
                {
                    "gram_type": list(triple),
                    "feasible_count": maximum.feasible_count,
                    "maximum": maximum.maximum.json(),
                    "root_pairing": margin.json(),
                }
                for triple, maximum, margin in same_q1
            ],
            "same_height_q2": {
                "positive_integral_candidates": [value.display() for value in sorted(set(positive_equal_height))],
                "exclusions": same_height_exclusions,
                "distinct_positive_pairing_impossible": True,
            },
            "accepted": True,
        },
        "replacement": {
            "scale": str(NEW_SCALE),
            "norm": replacement_norm.json(),
            "all_necessary_margins_strict": True,
            "tight_margin": new_margins[-1].json(),
            "prefix_has_open_hyperideal_patch": True,
        },
        "owner_claims_match": True,
        "conclusion": {
            "first_old_witness_cutter_realized": True,
            "cutter_accepted": True,
            "replacement_strict_through_bound": True,
            "finite_volume_example": False,
            "compact_example": False,
            "vinberg_terminated": False,
        },
    }


def render(payload: dict[str, object]) -> str:
    chronology = payload["chronology"]
    cutter = payload["explicit_cutter"]
    replacement = payload["replacement"]
    return f"""# Independent hostile audit: silver K12^2 through H393

## Verdict

**PASS.** A clean-room implementation, importing no project arithmetic or
search backend, reproduces the complete 16-class chronology through
`H=393+278sqrt(2)`: {chronology['control_count']} controls at
{chronology['distinct_height_count']} heights and exactly six nonempty
necessary rows. The explicit norm-two cutter is primitive, reflective, and
acute to the complete earlier prefix. The `7/4` dual-direction replacement
is strictly spacelike and strict against every necessary row through the
bound. This remains a bounded infinite-volume certificate, not Vinberg
termination or a global nonreflectivity result.

## Independent checks

- The literal K12 Gram has determinant 729; `G h=e1`, `h^2=4/3`, and `3h`
  is integral of norm 12. Independent PARI enumeration recovers coordinate
  minima `0,4,4,10` and cutter shell maxima `36+26sqrt(2)` and
  `40+30sqrt(2)`.
- The exact control census has {chronology['through_first_packet']} controls
  through the first packet, {chronology['strictly_between']} strictly between
  it and the cutter, and {chronology['at_cutting_height']} at the cutter.
- Direct coordinate multiplication gives `Gram(D,E)=(28,84,46)`,
  `Gram(U,V)=(380,204,278)`, `r^2=2`, and `(h.U,h.V)=(22,16)`.
  All reflection pairings are integral. An odd rational coefficient of `U`
  excludes divisibility by the unique prime over two, proving primitivity.
- Complete shell maximisation makes the cutter strictly acute to both first
  packet components. Necessary four-vector Schur censuses give intervening
  maxima `0,0,-7-4sqrt(2)`; all four simultaneous norm-one maxima are strict.
- For a distinct simultaneous norm-two root, the two-embedding Cauchy box
  leaves only positive pairings `1,sqrt(2),2`. The first two force forbidden
  coefficient-lattice squares below the K12^2 minimum four; `2` forces equal
  spatial vectors. Thus no distinct positive same-height pairing exists.
- The replacement norm is `{replacement['norm']['display']}` and its tight
  wall margin is `{replacement['tight_margin']['display']}`.

Owner claims and the independently reconstructed summary agree exactly. The
control/Gram stream hash is
`{chronology['control_and_gram_stream_sha256']}`.

## Reproduction

```bash
python3 scripts/independent_hostile_audit_silver_coxeter_todd_squared_h393_replacement.py --verify
PYTHONPATH=src pytest -q tests/test_silver_coxeter_todd_squared_h393_independent_hostile_audit.py
```
"""


def checksum_text() -> str:
    paths = (Path(__file__), AUDIT_JSON, AUDIT_REPORT)
    return "".join(f"{sha(path)}  {path.relative_to(ROOT)}\n" for path in paths)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    report = render(payload)
    if args.verify:
        if AUDIT_JSON.read_text() != serialized:
            raise SystemExit(f"audit JSON drift: {AUDIT_JSON}")
        if AUDIT_REPORT.read_text() != report:
            raise SystemExit(f"audit report drift: {AUDIT_REPORT}")
        expected = checksum_text()
        if AUDIT_SHA.read_text() != expected:
            raise SystemExit(f"audit checksums drift: {AUDIT_SHA}")
        print("PASS independent hostile silver K12^2 H393 audit")
        return
    AUDIT_JSON.write_text(serialized)
    AUDIT_REPORT.write_text(report)
    # The checksum file includes itself only indirectly: source, JSON, report.
    AUDIT_SHA.write_text(checksum_text())
    print(AUDIT_JSON)
    print(AUDIT_REPORT)
    print(AUDIT_SHA)


if __name__ == "__main__":
    main()
