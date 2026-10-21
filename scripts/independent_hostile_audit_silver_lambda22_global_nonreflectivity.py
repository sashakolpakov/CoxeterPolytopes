#!/usr/bin/env python3
"""No-production-import audit of the global silver Lambda22 theorem.

The arithmetic below is deliberately standalone.  The audit reads only the
vendored Lambda22 Gram and the frozen bounded H14/H150 certificates; it does
not import the production global-nonreflectivity module or silver search code.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key
from math import gcd, isqrt, lcm
from pathlib import Path

from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form


ROOT = Path(__file__).resolve().parents[1]
GRAM_PATH = ROOT / "data/lattices/lambda22_nebe_sloane_gram.txt"
H14_PATH = ROOT / "data/enumerations/silver_lambda22_h14.json"
H14_AUDIT_PATH = (
    ROOT / "data/enumerations/silver_lambda22_h14_independent_hostile_audit.json"
)
H150_PATH = ROOT / "data/enumerations/silver_lambda22_h150_frontier.json"
PRODUCTION_PATH = (
    ROOT / "data/enumerations/silver_lambda22_global_nonreflectivity.json"
)
OUTPUT_PATH = (
    ROOT
    / "data/enumerations/silver_lambda22_global_nonreflectivity_independent_hostile_audit.json"
)
REPORT_PATH = (
    ROOT
    / "results/silver_lambda22_global_nonreflectivity_independent_hostile_audit.md"
)

EXPECTED_SHA256 = {
    "gram": "4f3dbdc96e7aead7fae09a834abdbe2d5319620ece24656720f13184dd85b2f1",
    "H14": "5d81e273ae8d826c7b4c26f2cda99bc182ee287860a2fe1b39928ea862fc1a77",
    "H14_audit": "35550b5c534cd00b87dee4490a0a645ae28a6cd9dd68c0790b9181377b76a72c",
    "H150": "df2f0c23321bc1c99ccd4a23157b07bd72ce51b9284399f637f4e30f6943ab00",
}


@dataclass(frozen=True, slots=True)
class K2:
    a: Fraction = Fraction(0)
    b: Fraction = Fraction(0)

    @classmethod
    def make(
        cls, a: int | str | Fraction = 0, b: int | str | Fraction = 0
    ) -> "K2":
        return cls(Fraction(a), Fraction(b))

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
        output = ONE
        base = self
        while exponent:
            if exponent & 1:
                output = output * base
            base = base * base
            exponent >>= 1
        return output

    def sign(self, embedding: int = 0) -> int:
        if embedding not in (0, 1):
            raise ValueError("embedding must be zero or one")
        if embedding == 1:
            return K2(self.a, -self.b).sign()
        sa = (self.a > 0) - (self.a < 0)
        sb = (self.b > 0) - (self.b < 0)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        return sa if self.a * self.a - 2 * self.b * self.b > 0 else sb

    @property
    def integral(self) -> bool:
        return self.a.denominator == self.b.denominator == 1

    def exact(self) -> dict[str, str]:
        return {"a": str(self.a), "b": str(self.b)}

    def display(self) -> str:
        if not self.b:
            return str(self.a)
        if not self.a:
            return f"{self.b}*sqrt(2)"
        operation = "+" if self.b > 0 else "-"
        return f"{self.a} {operation} {abs(self.b)}*sqrt(2)"


ZERO = K2.make()
ONE = K2.make(1)
R = K2.make(0, 1)
COEFFICIENT = R - ONE
A = K2.make(5, 4)
B = K2.make(2, 1)
DOT = K2.make(4, 3)
SHIFT = K2.make(1, 1)
CONTROL_BOUND = K2.make(8, 6)
HEIGHT_BOUND = K2.make(56, 40)
FRONTIER_BOUND = K2.make(150)
H14_BOUND = K2.make(14, 10)
Z_LOWER = K2.make(10, -7) / K2.make(6)
CONTROL_BOX = 45

MIRROR_NORMS = (
    K2.make(1),
    K2.make(2),
    K2.make(4, 2),
    K2.make(4),
    K2.make(3),
    K2.make(6),
    K2.make(12, 6),
    K2.make(12),
)

NORM_GATES = {
    K2.make(1): (1, 1),
    K2.make(2): (1, 1),
    K2.make(4, 2): (2, 1),
    K2.make(4): (2, 2),
    K2.make(3): (3, 3),
    K2.make(6): (3, 3),
    K2.make(12, 6): (6, 3),
    K2.make(12): (6, 6),
}
SUBLATTICE_MINIMUM = {1: 4, 2: 12, 3: 24, 6: 96}
LOCAL_COEFFICIENTS = tuple(
    (first, second, third)
    for first in range(-3, 4)
    for second in range(-3, 4)
    for third in range(-3, 4)
    if first or second or third
)


@dataclass(frozen=True, slots=True)
class Control:
    norm: K2
    x: K2
    target: K2
    height: K2


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_gram() -> Matrix:
    rows: list[list[int]] = []
    header = None
    for raw in GRAM_PATH.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        values = [int(value) for value in line.split()]
        if header is None:
            header = values
        else:
            rows.append(values)
    if header != [22, 0] or len(rows) != 22:
        raise AssertionError("malformed Lambda22 Gram stream")
    return Matrix(22, 22, lambda i, j: rows[max(i, j)][min(i, j)])


def compare_k2(left: K2, right: K2) -> int:
    return (left - right).sign(0)


def compare_controls(left: Control, right: Control) -> int:
    for first, second in (
        (left.height, right.height),
        (left.norm, right.norm),
        (left.x, right.x),
    ):
        comparison = compare_k2(first, second)
        if comparison:
            return comparison
    return 0


def sign_integer_pair(rational: int, radical: int, embedding: int = 0) -> int:
    """Exact sign of ``rational+radical*sqrt(2)`` without floats."""

    if embedding == 1:
        radical = -radical
    elif embedding != 0:
        raise ValueError("embedding must be zero or one")
    sa = (rational > 0) - (rational < 0)
    sb = (radical > 0) - (radical < 0)
    if not sa:
        return sb
    if not sb or sa == sb:
        return sa
    return sa if rational * rational - 2 * radical * radical > 0 else sb


def control_box_proof() -> dict[str, object]:
    """Prove that the symmetric coefficient box contains every H<=150 control."""

    # H<=150 gives a*x_0^2<=150*q_0.  The largest right-hand
    # side occurs for q=12+6sqrt(2), and |x_0|<87 is uniform.
    # Conjugate target positivity is q_1+a_1*x_1^2>=0; the uniform
    # inequality q_1+9*a_1<0 gives |x_1|<3.  Hence
    # |m|=|x_0+x_1|/2<45 and |n|sqrt(2)<45.
    for norm in MIRROR_NORMS:
        if (K2.make(87 * 87) * COEFFICIENT - FRONTIER_BOUND * norm).sign(0) <= 0:
            raise AssertionError("identity control bound 87 is insufficient")
        if (norm + K2.make(9) * COEFFICIENT).sign(1) >= 0:
            raise AssertionError("conjugate control bound 3 is insufficient")
    return {
        "identity_embedding_absolute_bound": "<87",
        "conjugate_embedding_absolute_bound": "<3",
        "derived_rational_coefficient_bound": "<45",
        "derived_radical_coefficient_bound_times_sqrt2": "<45",
        "enumerated_symmetric_coefficient_box": CONTROL_BOX,
    }


def enumerate_controls() -> tuple[Control, ...]:
    """Clean-room enumeration of every reflective arithmetic control H<=150."""

    control_box_proof()
    rows: set[Control] = set()
    for norm in MIRROR_NORMS:
        qa, qb = int(norm.a), int(norm.b)
        qnorm = qa * qa - 2 * qb * qb
        for rational in range(-CONTROL_BOX, CONTROL_BOX + 1):
            for radical in range(-CONTROL_BOX, CONTROL_BOX + 1):
                if sign_integer_pair(rational, radical) <= 0:
                    continue
                # 2*a*x=(-2m+4n)+(2m-2n)sqrt(2).  Division by q is
                # integral precisely when the two conjugate numerators below
                # are divisible by N(q).
                c = -2 * rational + 4 * radical
                d = 2 * rational - 2 * radical
                if (c * qa - 2 * d * qb) % qnorm:
                    continue
                if (d * qa - c * qb) % qnorm:
                    continue
                target_a = (
                    qa
                    - rational * rational
                    - 2 * radical * radical
                    + 4 * rational * radical
                )
                target_b = (
                    qb
                    + rational * rational
                    + 2 * radical * radical
                    - 2 * rational * radical
                )
                if sign_integer_pair(target_a, target_b, 1) < 0:
                    continue
                height = K2.make(target_a - qa, target_b - qb) / norm
                if (height - FRONTIER_BOUND).sign(0) > 0:
                    continue
                rows.add(
                    Control(
                        norm,
                        K2.make(rational, radical),
                        K2.make(target_a, target_b),
                        height,
                    )
                )
    return tuple(sorted(rows, key=cmp_to_key(compare_controls)))


def combination_square(
    gram: tuple[int, int, int], first: int, second: int
) -> int:
    u_square, v_square, dot = gram
    return (
        first * first * u_square
        + second * second * v_square
        + 2 * first * second * dot
    )


def has_forbidden_binary_combination(gram: tuple[int, int, int]) -> bool:
    """Exactly detect all integral binary combinations of norm 1, 2, or 3."""

    u_square, v_square, dot = gram
    determinant = u_square * v_square - dot * dot
    if determinant < 0:
        return True
    if determinant == 0:
        return gcd(gcd(u_square, v_square), abs(dot)) < 4
    first_bound = isqrt((3 * v_square) // determinant) + 1
    second_bound = isqrt((3 * u_square) // determinant) + 1
    return any(
        0 < combination_square(gram, first, second) < 4
        for first in range(-first_bound, first_bound + 1)
        for second in range(-second_bound, second_bound + 1)
        if first or second
    )


def candidate_grams(target: K2) -> tuple[tuple[int, int, int], ...]:
    """Complete even/minimum-four/Cauchy binary necessary census."""

    if not target.integral or target.b.numerator % 2:
        return ()
    rational_target = int(target.a)
    dot = int(target.b / 2)
    rows = []
    for v_square in range(0, rational_target // 2 + 1, 2):
        u_square = rational_target - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if u_square and u_square < 4:
            continue
        if v_square and v_square < 4:
            continue
        gram = (u_square, v_square, dot)
        if dot * dot > u_square * v_square:
            continue
        if has_forbidden_binary_combination(gram):
            continue
        rows.append(gram)
    return tuple(rows)


def extension_determinant(
    gram: tuple[int, int, int], p: int, q: int
) -> int:
    """Determinant of Gram(u,v,d), d^2=4 (one cross term, exactly)."""

    u_square, v_square, dot = gram
    return (
        4 * (u_square * v_square - dot * dot)
        - v_square * p * p
        - u_square * q * q
        + 2 * dot * p * q
    )


def locally_admissible_extension(
    gram: tuple[int, int, int], p: int, q: int
) -> bool:
    """Apply the declared finite [-3,3]^3 minimum-four necessary gate."""

    u_square, v_square, dot = gram
    for first, second, third in LOCAL_COEFFICIENTS:
        square = (
            first * first * u_square
            + second * second * v_square
            + 4 * third * third
            + 2 * first * second * dot
            + 2 * first * third * p
            + 2 * second * third * q
        )
        if 0 < square < 4:
            return False
    return True


def ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def extensions(
    gram: tuple[int, int, int], p_modulus: int, q_modulus: int
) -> tuple[tuple[int, int, int], ...]:
    """Enumerate every gated integral PSD extension by d^2=4."""

    u_square, v_square, dot = gram
    base_determinant = u_square * v_square - dot * dot
    if base_determinant < 0:
        return ()
    output = []
    p_bound = isqrt(4 * u_square)
    q_bound = isqrt(4 * v_square)
    if not u_square:
        if dot:
            return ()
        q_start = ceil_div(-q_bound, q_modulus) * q_modulus
        for q in range(q_start, q_bound + 1, q_modulus):
            if locally_admissible_extension(gram, 0, q):
                output.append((0, q, extension_determinant(gram, 0, q)))
        return tuple(output)
    p_start = ceil_div(-p_bound, p_modulus) * p_modulus
    for p in range(p_start, p_bound + 1, p_modulus):
        radius_square = base_determinant * (4 * u_square - p * p)
        if radius_square < 0:
            continue
        radius = isqrt(radius_square)
        lower = max(-q_bound, ceil_div(dot * p - radius, u_square))
        upper = min(q_bound, (dot * p + radius) // u_square)
        q_start = ceil_div(lower, q_modulus) * q_modulus
        for q in range(q_start, upper + 1, q_modulus):
            determinant = extension_determinant(gram, p, q)
            if determinant < 0:
                raise AssertionError("ellipse interval admitted negative determinant")
            # Direct determinant regression guards against a duplicated cross term.
            direct = int(Matrix([[u_square, dot, p], [dot, v_square, q], [p, q, 4]]).det())
            if determinant != direct:
                raise AssertionError("closed 3x3 determinant formula is wrong")
            if locally_admissible_extension(gram, p, q):
                output.append((p, q, determinant))
    return tuple(output)


def cleanroom_frontier_scan() -> dict[str, object]:
    """Reconstruct the complete H<=150 necessary superset from first principles."""

    controls = enumerate_controls()
    retained: list[dict[str, object]] = []
    dangerous: list[dict[str, object]] = []
    retained_levels: set[tuple[K2, K2, K2]] = set()
    for control in controls:
        u_modulus, v_modulus = NORM_GATES[control.norm]
        cross_modulus = lcm(u_modulus, v_modulus)
        for gram in candidate_grams(control.target):
            u_square, v_square, dot = gram
            # This is deliberately the declared weak necessary gate used by
            # H150: minimum plus divisibility.  S2 norm residues mod 8 enter
            # only later as an exact resolution, never to shrink this census.
            if (
                (u_square and u_square < SUBLATTICE_MINIMUM[u_modulus])
                or (v_square and v_square < SUBLATTICE_MINIMUM[v_modulus])
                or u_square % u_modulus
                or v_square % v_modulus
                or dot % cross_modulus
            ):
                continue
            extension_rows = extensions(gram, u_modulus, v_modulus)
            if not extension_rows:
                continue
            best_p, best_q, best_det = extension_rows[0]
            best_value = K2.make(best_p, best_q)
            for p, q, determinant in extension_rows[1:]:
                value = K2.make(p, q)
                if compare_k2(value, best_value) > 0:
                    best_p, best_q, best_det = p, q, determinant
                    best_value = value
            threshold = K2.make(3) * COEFFICIENT * control.x
            margin = best_value - threshold
            nonnegative = tuple(
                (p, q, determinant)
                for p, q, determinant in extension_rows
                if (K2.make(p, q) - threshold).sign(0) >= 0
            )
            row = {
                "height": control.height.exact(),
                "height_display": control.height.display(),
                "norm": control.norm.exact(),
                "control": control.x.exact(),
                "target": control.target.exact(),
                "gram": list(gram),
                "pairing_gate": {
                    "u_modulus": u_modulus,
                    "v_modulus": v_modulus,
                    "cross_modulus": cross_modulus,
                },
                "extension_count": len(extension_rows),
                "maximizing_extension": [best_p, best_q, best_det],
                "maximum": best_value.exact(),
                "margin": margin.exact(),
                "nonnegative_extensions": [list(item) for item in nonnegative],
            }
            retained.append(row)
            retained_levels.add((control.height, control.norm, control.x))
            if compare_k2(control.height, H14_BOUND) > 0 and margin.sign(0) >= 0:
                dangerous.append(row)

    expected_dangers = (
        (K2.make(21, 15), (60, 36, 45), (15, 12)),
        (K2.make(36, Fraction(51, 2)), (354, 198, 264), (36, 27)),
        (K2.make(41, 29), (88, 40, 58), (18, 12)),
        (K2.make(42, 30), (228, 144, 180), (30, 24)),
        (K2.make(56, 40), (124, 52, 80), (22, 14)),
        (K2.make(56, 40), (100, 64, 80), (20, 16)),
        (K2.make(70, Fraction(99, 2)), (230, 126, 170), (30, 22)),
    )
    actual_dangers = tuple(
        (
            K2.make(row["height"]["a"], row["height"]["b"]),
            tuple(row["gram"]),
            tuple(row["maximizing_extension"][:2]),
        )
        for row in dangerous
    )
    if (len(controls), len(retained_levels), len(retained)) != (58, 18, 72):
        raise AssertionError(
            "clean-room H150 census drifted: "
            f"{len(controls)}, {len(retained_levels)}, {len(retained)}"
        )
    if actual_dangers != expected_dangers:
        raise AssertionError(f"clean-room H150 dangers changed: {actual_dangers!r}")

    # At the two shell-pairing exclusions, checking only the maximizer would
    # be insufficient.  Exhaustion proves the complete nonnegative extension
    # sets are the asserted singletons.
    if dangerous[2]["nonnegative_extensions"] != [[18, 12, 48]]:
        raise AssertionError(
            "H41 nonnegative extension set is not a singleton: "
            f"{dangerous[2]['nonnegative_extensions']!r}"
        )
    if dangerous[4]["nonnegative_extensions"] != [[22, 14, 0]]:
        raise AssertionError("H56 shell threat extension set is not a singleton")

    h36_height = K2.make(36, Fraction(51, 2))
    h36_retained = [
        row for row in retained
        if K2.make(row["height"]["a"], row["height"]["b"]) == h36_height
    ]
    expected_h36_grams = [
        [402, 174, 264],
        [390, 180, 264],
        [378, 186, 264],
        [366, 192, 264],
        [354, 198, 264],
        [342, 204, 264],
    ]
    if [row["gram"] for row in h36_retained] != expected_h36_grams:
        raise AssertionError("complete retained H36 level changed")

    return {
        "coefficient_box_proof": control_box_proof(),
        "control_count": len(controls),
        "retained_level_count": len(retained_levels),
        "retained_necessary_type_count": len(retained),
        "post_H14_nonnegative_type_count": len(dangerous),
        "post_H14_nonnegative_types": dangerous,
        "H41_nonnegative_extension_set": dangerous[2]["nonnegative_extensions"],
        "H56_124_52_80_nonnegative_extension_set": dangerous[4]["nonnegative_extensions"],
        "H36_retained_grams": expected_h36_grams,
        "determinant_formula": "4(UV-M^2)-V*p^2-U*q^2+2*M*p*q",
        "direct_matrix_determinant_regression": True,
    }


def kernel_basis_mod_prime(gram: Matrix, prime: int) -> Matrix:
    """Construct {x:Gx=0 mod p} by clean-room modular row reduction."""

    reduced = [[int(gram[i, j]) % prime for j in range(22)] for i in range(22)]
    pivots: list[int] = []
    row = 0
    for column in range(22):
        pivot = next(
            (index for index in range(row, 22) if reduced[index][column]), None
        )
        if pivot is None:
            continue
        reduced[row], reduced[pivot] = reduced[pivot], reduced[row]
        inverse = pow(reduced[row][column], -1, prime)
        reduced[row] = [(value * inverse) % prime for value in reduced[row]]
        for index in range(22):
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
    for free in (column for column in range(22) if column not in pivots):
        vector = [0] * 22
        vector[free] = 1
        for index, pivot in enumerate(pivots):
            vector[pivot] = (-reduced[index][free]) % prime
        radical.append(Matrix(vector))
    basis = hermite_normal_form(Matrix.hstack(prime * Matrix.eye(22), *radical))
    if basis.shape != (22, 22) or any(int(value) % prime for value in gram * basis):
        raise AssertionError("pairing-radical basis construction failed")
    return basis


def gated_bases(gram: Matrix) -> dict[int, Matrix]:
    s2 = kernel_basis_mod_prime(gram, 2)
    s3 = kernel_basis_mod_prime(gram, 3)
    # Chinese remaindering: S6=S2 intersect S3=3*S2+2*S3.
    s6 = hermite_normal_form(Matrix.hstack(3 * s2, 2 * s3))
    if s6.shape != (22, 22) or any(int(value) % 6 for value in gram * s6):
        raise AssertionError("S6 intersection basis construction failed")
    return {2: s2, 3: s3, 6: s6}


def gp_matrix(matrix: Matrix) -> str:
    return "[" + ";".join(
        ",".join(str(int(matrix[i, j])) for j in range(matrix.cols))
        for i in range(matrix.rows)
    ) + "]"


def independent_gated_shell_audit(gram: Matrix) -> dict[str, object]:
    """Recompute the S2/S3/S6 minima and critical direction pairings in PARI."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required for the clean-room shell audit")
    bases = gated_bases(gram)
    outputs: dict[str, dict[str, object]] = {}
    for modulus, bound in ((2, 12), (3, 24), (6, 96)):
        basis = bases[modulus]
        commands = [
            f"G={gp_matrix(gram)}",
            f"B={gp_matrix(basis)}",
            "H=B~*G*B",
            "d=vector(22,i,i==1)~",
            "ell=d~*G*B",
            f"R=qfminim(H,{bound},,2)",
            "M=R[3]",
            'print("signed_count=",R[1])',
            'print("maximum=",floor(R[2]+1/2))',
            'print("representatives=",matsize(M)[2])',
            "P=List()",
            (
                "for(j=1,matsize(M)[2],"
                f"if(M[,j]~*H*M[,j]=={bound},listput(P,abs(ell*M[,j]))))"
            ),
            f'print("norm_{bound}_representatives=",#P)',
            f'print("norm_{bound}_abs_pairings=",Set(Vec(P)))',
            f"R0=qfminim(H,{bound - 1},1)",
            'print("below_signed_count=",R0[1])',
            'print("below_maximum=",floor(R0[2]+1/2))',
            "quit()",
        ]
        process = subprocess.run(
            [gp, "-q", "-s", "512000000"],
            input=";".join(commands) + "\n",
            text=True,
            capture_output=True,
            timeout=240,
            check=True,
        )
        parsed: dict[str, object] = {}
        for raw in process.stdout.splitlines():
            key, value = raw.strip().split("=", 1)
            parsed[key] = json.loads(value) if key.endswith("abs_pairings") else int(value)
        outputs[str(modulus)] = parsed
    expected = {
        "2": {
            "signed_count": 8_448,
            "maximum": 12,
            "representatives": 4_224,
            "norm_12_representatives": 4_224,
            "norm_12_abs_pairings": [0, 2, 4],
            "below_signed_count": 0,
            "below_maximum": 0,
        },
        "3": {
            "signed_count": 1_782,
            "maximum": 24,
            "representatives": 891,
            "norm_24_representatives": 891,
            "norm_24_abs_pairings": [0, 3, 6],
            "below_signed_count": 0,
            "below_maximum": 0,
        },
        "6": {
            "signed_count": 1_782,
            "maximum": 96,
            "representatives": 891,
            "norm_96_representatives": 891,
            "norm_96_abs_pairings": [0, 6, 12],
            "below_signed_count": 0,
            "below_maximum": 0,
        },
    }
    if outputs != expected:
        raise AssertionError(f"clean-room gated shell census changed: {outputs!r}")
    return {
        "method": "clean-room modular kernels plus PARI qfminim(flag=2)",
        "shells": outputs,
    }


def s2_norm_residues_mod8(gram: Matrix) -> dict[str, object]:
    """Enumerate all four classes of S2/2Lambda22 and their norm residues."""

    basis = gated_bases(gram)[2]
    generators = {
        tuple(int(basis[i, j]) % 2 for i in range(22)) for j in range(22)
    }
    classes = {(0,) * 22}
    for generator in generators:
        classes |= {
            tuple((left[i] + generator[i]) % 2 for i in range(22))
            for left in tuple(classes)
        }
    residues = sorted(
        sum(vector[i] * int(gram[i, j]) * vector[j] for i in range(22) for j in range(22)) % 8
        for vector in classes
    )
    if len(classes) != 4 or residues != [0, 4, 4, 4]:
        raise AssertionError(f"S2/2L norm residues changed: {residues}")
    return {"class_count": len(classes), "residue_multiset": residues}


def resolve_cleanroom_frontier(
    scan: dict[str, object], shells: dict[str, object], residues: dict[str, object]
) -> dict[str, object]:
    """Eliminate every one of the seven independently regenerated threats."""

    dangers = scan["post_H14_nonnegative_types"]
    if len(dangers) != 7:
        raise AssertionError("expected seven clean-room frontier threats")
    shell_rows = shells["shells"]
    if shell_rows["2"]["maximum"] != 12 or shell_rows["2"]["norm_12_abs_pairings"] != [0, 2, 4]:
        raise AssertionError("S2 minimum/pairing exclusion unavailable")
    if shell_rows["3"]["maximum"] != 24:
        raise AssertionError("S3 minimum exclusion unavailable")
    if shell_rows["6"]["maximum"] != 96:
        raise AssertionError("S6 minimum exclusion unavailable")
    if residues != {"class_count": 4, "residue_multiset": [0, 4, 4, 4]}:
        raise AssertionError("S2 residue exclusion unavailable")

    # Each entry is checked numerically rather than accepted as prose.
    g0, g1, g2, g3, g4, g5, g6 = (tuple(row["gram"]) for row in dangers)
    if g0[0] + g0[1] - 2 * g0[2] != 6:
        raise AssertionError("H21 S3 difference calculation changed")
    h36_grams = [tuple(row) for row in scan["H36_retained_grams"]]
    if len(h36_grams) != 6:
        raise AssertionError("H36 retained level was not completely enumerated")
    for u_square, v_square, dot in h36_grams:
        if u_square + 2 * v_square != 750 or dot != 264:
            raise AssertionError("H36 target identities changed")
        if v_square <= 198:
            square = u_square + 4 * v_square - 4 * dot
            if not 0 < square < 96:
                raise AssertionError("H36 S6 branch failed")
        else:
            square = u_square + v_square - 2 * dot
            if not 0 < square < 24:
                raise AssertionError("H36 S3 branch failed")
    if g2[0] + g2[1] - 2 * g2[2] != 12:
        raise AssertionError("H41 S2 difference calculation changed")
    if dangers[2]["nonnegative_extensions"] != [[18, 12, 48]] or 18 - 12 != 6:
        raise AssertionError("H41 complete extension/pairing exclusion changed")
    if g3[0] + g3[1] - 2 * g3[2] != 12:
        raise AssertionError("H42 S6 difference calculation changed")
    if g4[0] + 4 * g4[1] - 4 * g4[2] != 12:
        raise AssertionError("H56 S2 combination calculation changed")
    if dangers[4]["nonnegative_extensions"] != [[22, 14, 0]] or 22 - 2 * 14 != -6:
        raise AssertionError("H56 complete extension/pairing exclusion changed")
    if g5[0] + g5[1] - 2 * g5[2] != 4:
        raise AssertionError("second H56 S2 difference calculation changed")
    if g6[0] % 8 not in (2, 6):
        raise AssertionError("H70 S2 residue obstruction changed")

    resolutions = [
        "u-v is a nonzero S3 vector of norm 6<24",
        "u-2v is a nonzero S6 vector of norm 90<96",
        "the sole nonnegative extension has d.(u-v)=6, absent on the complete S2 norm-12 shell",
        "u-v is a nonzero S6 vector of norm 12<96",
        "the sole nonnegative extension has d.(u-2v)=-6, absent on the complete S2 norm-12 shell",
        "u-v is a nonzero S2 vector of norm 4<12",
        "u lies in S2 but u^2=230 is 6 mod 8, outside the S2/2L residue set {0,4}",
    ]
    return {
        "resolved_type_count": len(resolutions),
        "resolutions": resolutions,
        "all_nonnegative_necessary_extensions_excluded": True,
        "all_seven_resolved": True,
        "strict_against_every_actual_root_through_150": True,
    }


def field_dot(left: tuple[K2, ...], right: tuple[K2, ...], gram: Matrix) -> K2:
    output = ZERO
    for i in range(22):
        for j in range(22):
            if gram[i, j]:
                output = output + left[i] * K2.make(int(gram[i, j])) * right[j]
    return output


def lorentz_dot(left: tuple[K2, ...], right: tuple[K2, ...], gram: Matrix) -> K2:
    return -COEFFICIENT * left[0] * right[0] + field_dot(
        left[1:], right[1:], gram
    )


def direction_dot(vector: tuple[K2, ...], gram: Matrix) -> K2:
    direction = (ONE,) + (ZERO,) * 21
    return field_dot(direction, vector, gram)


def transform(vector: tuple[K2, ...], gram: Matrix, inverse: bool = False) -> tuple[K2, ...]:
    x = vector[0]
    y = vector[1:]
    p = direction_dot(y, gram)
    if inverse:
        new_x = A * x - DOT * p
        shift = -B * x + SHIFT * p
    else:
        new_x = A * x + DOT * p
        shift = B * x + SHIFT * p
    return (new_x, y[0] + shift, *y[1:])


def independent_isometry_audit(gram: Matrix) -> dict[str, object]:
    basis = tuple(
        tuple(ONE if row == column else ZERO for row in range(23))
        for column in range(23)
    )
    images = tuple(transform(vector, gram) for vector in basis)
    inverse_images = tuple(transform(vector, gram, True) for vector in images)
    form_ok = all(
        lorentz_dot(images[i], images[j], gram)
        == lorentz_dot(basis[i], basis[j], gram)
        for i in range(23)
        for j in range(23)
    )
    if inverse_images != basis or not form_ok:
        raise AssertionError("standalone integral-isometry replay failed")
    v0 = basis[0]
    gv0 = transform(v0, gram)
    expected = (A, B) + (ZERO,) * 21
    if gv0 != expected:
        raise AssertionError("basepoint image changed")
    plane_determinant = A**2 - K2.make(4) * DOT * B
    if plane_determinant != ONE or (A - ONE).sign() <= 0:
        raise AssertionError("hyperbolic plane block is not infinite order")
    return {
        "full_basis_form_preservation": True,
        "integral_inverse": True,
        "basepoint_image": [entry.exact() for entry in gv0],
        "plane_determinant": plane_determinant.exact(),
        "plane_trace": (K2.make(2) * A).exact(),
        "trace_identity_embedding_greater_than_two": True,
        "chosen_time_component_preserved": True,
        "infinite_order": True,
    }


def time_zero_audit(h14: dict[str, object]) -> dict[str, object]:
    norms = tuple(
        K2.make(row["exact"]["a"], row["exact"]["b"])
        for row in h14["local_mirror_norms"]["canonical_norms"]
    )
    if norms != MIRROR_NORMS:
        raise AssertionError("complete mirror-norm alphabet changed")
    minima = h14["gated_sublattices"]
    if (
        minima["S2_minimum"],
        minima["S3_minimum"],
        minima["S6_minimum"],
    ) != (12, 24, 96):
        raise AssertionError("gated minima changed")
    # q, Tr(q)/2, minimum for nonzero u, minimum for nonzero v.
    rows = (
        (MIRROR_NORMS[0], 1, 4, 4),
        (MIRROR_NORMS[1], 2, 4, 4),
        (MIRROR_NORMS[2], 4, 12, 4),
        (MIRROR_NORMS[3], 4, 12, 12),
        (MIRROR_NORMS[4], 3, 24, 24),
        (MIRROR_NORMS[5], 6, 24, 24),
        (MIRROR_NORMS[6], 12, 96, 24),
        (MIRROR_NORMS[7], 12, 96, 96),
    )
    # If r=(0,u+sqrt(2)v), then u^2+2v^2=Tr(q)/2.  Every
    # nonzero coefficient already exceeds the right side.
    if not all(q.a == target and umin > target and 2 * vmin > target for q, target, umin, vmin in rows):
        raise AssertionError("a time-zero mirror was not excluded")
    return {
        "complete_norm_class_count": len(rows),
        "identity": "u^2+2v^2=Tr(q)/2",
        "cases": [
            {
                "norm": q.exact(),
                "trace_half_target": target,
                "u_minimum": umin,
                "twice_v_minimum": 2 * vmin,
            }
            for q, target, umin, vmin in rows
        ],
        "no_time_zero_reflective_root": True,
        "basepoint_regular": True,
    }


def frontier_audit(
    h14: dict[str, object],
    h150: dict[str, object],
    gram: Matrix,
) -> dict[str, object]:
    dependency = h150["base_H14_artifact"]
    if dependency["sha256"] != EXPECTED_SHA256["H14"]:
        raise AssertionError("H150 does not pin the audited H14 base")
    scan = h150["necessary_superset_scan"]
    if (scan["control_count"], scan["retained_level_count"], scan["retained_necessary_type_count"]) != (58, 18, 72):
        raise AssertionError("H150 necessary census changed")
    if scan["post_H14_nonnegative_type_count"] != 7:
        raise AssertionError("H150 frontier threat count changed")
    if h150["height_bound"]["display"] != "150":
        raise AssertionError("H150 bound changed")
    if h14["strict_hyperideal_witness"]["vector"] != (
        "z=(1,d/3), d any Lambda22 minimal vector"
    ):
        raise AssertionError("H14 strict witness changed")
    cleanroom = cleanroom_frontier_scan()
    shells = independent_gated_shell_audit(gram)
    residues = s2_norm_residues_mod8(gram)
    resolutions = resolve_cleanroom_frontier(cleanroom, shells, residues)
    production_dangers = tuple(
        (
            row["height"]["exact"],
            row["gram"],
            [
                row["relaxed_maximum"]["d_dot_u"],
                row["relaxed_maximum"]["d_dot_v"],
            ],
        )
        for row in scan["post_H14_nonnegative_types"]
    )
    independent_dangers = tuple(
        (
            row["height"],
            row["gram"],
            row["maximizing_extension"][:2],
        )
        for row in cleanroom["post_H14_nonnegative_types"]
    )
    if independent_dangers != production_dangers:
        raise AssertionError("clean-room dangers disagree with production H150")
    return {
        "H14_sha256": EXPECTED_SHA256["H14"],
        "H14_independent_audit_sha256": EXPECTED_SHA256["H14_audit"],
        "H150_sha256": EXPECTED_SHA256["H150"],
        "H150_control_count": 58,
        "H150_retained_type_count": 72,
        "H150_resolved_nonnegative_types": 7,
        "cleanroom_necessary_superset": cleanroom,
        "cleanroom_gated_shells": shells,
        "cleanroom_S2_norm_residues_mod8": residues,
        "cleanroom_exact_resolutions": resolutions,
        "production_danger_tuple_comparison": True,
        "base_packet_and_post_H14_frontier_complete": True,
        "strict_against_every_actual_root_through_150": True,
    }


def separator_audit() -> dict[str, object]:
    # For an oriented root r=(x,y), x>0.  Separation implies
    # p=(d,y)>=k*x with k=(4-r)/2.  Cauchy p^2<=4(q+a*x^2)
    # then gives x^2/q<=1/(k^2/4-a)=136+96r.
    k = K2.make(4, -1) / K2.make(2)
    denominator = k**2 / K2.make(4) - COEFFICIENT
    control_square = ONE / denominator
    height = COEFFICIENT * control_square
    lower = -COEFFICIENT + k / K2.make(3)
    if control_square != K2.make(136, 96):
        raise AssertionError("separator control square changed")
    if CONTROL_BOUND**2 != control_square:
        raise AssertionError("separator square root changed")
    if height != HEIGHT_BOUND or not (height - FRONTIER_BOUND).sign() < 0:
        raise AssertionError("separator escaped the H150 frontier")
    if lower != Z_LOWER or lower.sign() <= 0:
        raise AssertionError("separator no longer forces positive z pairing")
    # The sign/orientation arithmetic itself is independently pinned:
    # (gv0,r)=-(3+r)x+(2+r)p, because a*A=3+r;
    # dividing by B=2+r gives k=(3+r)/(2+r)=(4-r)/2.
    if COEFFICIENT * A != K2.make(3, 1):
        raise AssertionError("endpoint Lorentz coefficient changed")
    if K2.make(3, 1) / B != k:
        raise AssertionError("separator slope changed")
    return {
        "oriented_root": "x>0, so (v0,r)=-(sqrt(2)-1)x<0",
        "endpoint_pairing": "(gv0,r)=-(3+sqrt(2))x+(2+sqrt(2))(d,y)",
        "separation_slope": k.exact(),
        "cauchy_control_square": control_square.exact(),
        "control_bound": CONTROL_BOUND.exact(),
        "height_bound": height.exact(),
        "height_bound_less_than_150": True,
        "z_pairing_lower_coefficient": lower.exact(),
        "z_pairing_strictly_positive": True,
        "no_separating_reflective_mirror": True,
    }


def chamber_logic_audit() -> dict[str, object]:
    return {
        "mirror_arrangement_invariant_under_g": (
            "An integral isometry conjugates an integral reflection s_r to "
            "the integral reflection s_{g r}; hence it permutes all mirrors."
        ),
        "same_chamber": (
            "The regular endpoints v0 and gv0 have no separating mirror. "
            "Connected components of hyperbolic space minus the locally finite "
            "mirror union are exactly the open reflection chambers."
        ),
        "stabilizer": (
            "Because g permutes mirrors and maps v0 into the same component, "
            "g maps that component to itself."
        ),
        "infinite_volume": (
            "The restriction of g to span(v0,d) is hyperbolic, with expanding "
            "eigenvalue >1 and positive translation length.  A ball of radius "
            "less than half that translation length about any chamber-interior "
            "point has pairwise disjoint g^n-translates, all contained in the "
            "g-invariant chamber.  Therefore its volume is infinite."
        ),
        "index": (
            "The integral arithmetic orthogonal group is a finite-covolume "
            "lattice.  The reflection chamber is a fundamental domain for the "
            "full integral reflection subgroup, so infinite chamber volume is "
            "equivalent to infinite reflection index."
        ),
        "logic_pass": True,
    }


def host_scope_audit(gram: Matrix) -> dict[str, object]:
    if gram.det() != 12:
        raise AssertionError("Lambda22 determinant changed")
    # At the identity place a=sqrt(2)-1>0, hence <-a>+L is Lorentzian.
    # At the other place a'=-sqrt(2)-1<0, so <-a'> is positive and the
    # positive Lambda22 block remains positive.  This is the admissible
    # arithmetic orthogonal construction.  Definite conjugacy also implies
    # K-anisotropy, hence the ambient arithmetic lattice is cocompact.
    if COEFFICIENT.sign(0) <= 0 or COEFFICIENT.sign(1) >= 0:
        raise AssertionError("host is not admissible")
    return {
        "rank": 23,
        "hyperbolic_dimension": 22,
        "spatial_determinant": 12,
        "identity_signature": [22, 1],
        "conjugate_signature": [23, 0],
        "admissible_arithmetic_host": True,
        "anisotropic_and_cocompact": True,
        "scope": (
            "This proves nonreflectivity of this entire integral Lambda22 host, "
            "not nonexistence for every Q(sqrt(2)) host or every chamber in H22."
        ),
    }


def audit() -> dict[str, object]:
    for key, path in (
        ("gram", GRAM_PATH),
        ("H14", H14_PATH),
        ("H14_audit", H14_AUDIT_PATH),
        ("H150", H150_PATH),
    ):
        if sha256(path) != EXPECTED_SHA256[key]:
            raise AssertionError(f"pinned dependency changed: {key}")
    gram = parse_gram()
    h14 = json.loads(H14_PATH.read_text())
    h14_audit = json.loads(H14_AUDIT_PATH.read_text())
    h150 = json.loads(H150_PATH.read_text())
    if h14_audit["verdict"] != "PASS":
        raise AssertionError("the independent H14 base audit is not PASS")
    payload = {
        "schema_version": "0.1.0",
        "verdict": "PASS_GLOBAL_NONREFLECTIVITY",
        "independence": (
            "No import of the production global theorem or silver search. "
            "Standalone Q(sqrt(2)) arithmetic, Gram parsing, 23-basis isometry "
            "checks, time-zero analysis, separator algebra, and theorem logic."
        ),
        "host": host_scope_audit(gram),
        "isometry": independent_isometry_audit(gram),
        "time_zero_mirrors": time_zero_audit(h14),
        "bounded_frontier_dependency": frontier_audit(h14, h150, gram),
        "separator": separator_audit(),
        "chamber_and_volume_logic": chamber_logic_audit(),
        "production_comparison": {
            "path": str(PRODUCTION_PATH.relative_to(ROOT)),
            "sha256": sha256(PRODUCTION_PATH),
            "compared_only_after_independent_derivation": True,
        },
        "conclusion": {
            "full_reflection_chamber_infinite_volume": True,
            "full_integral_reflection_subgroup_infinite_index": True,
            "silver_Lambda22_host_reflective": False,
            "field_wide_nonexistence": False,
        },
    }
    production = json.loads(PRODUCTION_PATH.read_text())
    if not all(
        production["conclusion"][key]
        for key in (
            "reflection_chamber_has_infinite_volume",
            "reflection_subgroup_has_infinite_index",
            "global_nonreflectivity_theorem",
        )
    ):
        raise AssertionError("production theorem conclusion disagrees")
    return payload


def report(payload: dict[str, object]) -> str:
    production = payload["production_comparison"]
    return f"""# Independent hostile audit: global silver Lambda22 nonreflectivity

Terminal verdict: **PASS_GLOBAL_NONREFLECTIVITY**.

For

`M=<-(sqrt(2)-1)> + Lambda22` over `Z[sqrt(2)]`,

the full integral reflection subgroup has infinite index and its complete
reflection chamber has infinite hyperbolic volume.  This is a global theorem
for the displayed host, not a bounded-prefix conclusion.

## Independent arithmetic

The replay imports no production silver module.  It reparses the rank-22
Gram, checks the explicit map on all 23 basis vectors, reconstructs its
integral inverse, and verifies the full form identity.  Its Lorentz-plane
block has determinant one and trace `10+8sqrt(2)>2`, so the map has infinite
order and preserves the chosen time component.

The eight complete mirror-norm classes are replayed separately at time
coordinate zero.  From `y=u+sqrt(2)v` one has
`u^2+2v^2=Tr(q)/2`; the ordinary/S2/S3/S6 minima exclude all eight cases.
Thus the basepoint is regular.

## Separator bound

If an oriented reflective root separates `v0` from `g*v0`, endpoint signs
give

`(d,y) >= ((4-sqrt(2))/2)x`.

Cauchy then gives

`x^2/q <= 136+96sqrt(2)=(8+6sqrt(2))^2`,

and hence `H<=56+40sqrt(2)<150`.  The same endpoint inequality forces

`(z,r)>=(10-7sqrt(2))x/6>0`, for `z=(1,d/3)`.

The H150 frontier is not merely trusted by hash.  This script independently
enumerates all 58 arithmetic controls in a proved Minkowski coefficient box,
all even binary Gram candidates, their exact binary minimum-four test and
pairing gates, and every gated PSD extension by `d`.  It reproduces 18
retained levels, 72 necessary types, and exactly seven post-H14 nonnegative
types.  Clean-room modular kernels and PARI shell enumeration prove minima
`min(S2,S3,S6)=(12,24,96)`, the critical S2 norm-12 pairing spectrum
`{{0,2,4}}`, and the S2 mod-8 residue classes.  These data eliminate all
seven threats; in particular, the complete nonnegative extension sets in the
two shell-pairing cases are singletons `(18,12)` and `(22,14)`.

Thus every actual root through H150 pairs strictly negatively with `z`, and
no mirror separates the regular endpoints.  Since the integral isometry
permutes all integral-reflection mirrors, it stabilizes their common open
chamber.

The restriction is a hyperbolic translation of positive length.  A small
ball in the chamber has infinitely many pairwise disjoint translates inside
that chamber, proving infinite volume.  Since the ambient admissible
arithmetic orthogonal group has finite covolume, the reflection subgroup has
infinite index.  The conjugate form is definite, so the host is anisotropic
and its ambient arithmetic quotient is cocompact.

The production artifact compared after the independent derivation is
`{production['path']}`, SHA-256 `{production['sha256']}`.

This does **not** prove field-wide nonexistence over `Q(sqrt(2))` or exclude
other `H22` hosts.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = audit()
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    rendered = report(payload)
    if args.verify:
        if OUTPUT_PATH.read_text() != serialized or REPORT_PATH.read_text() != rendered:
            raise SystemExit("FAIL: frozen independent hostile audit changed")
        print(payload["verdict"])
        return
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(serialized)
    REPORT_PATH.write_text(rendered)
    print(payload["verdict"])


if __name__ == "__main__":
    main()
