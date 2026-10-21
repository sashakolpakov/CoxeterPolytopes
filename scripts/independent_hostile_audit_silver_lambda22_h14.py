#!/usr/bin/env python3
"""Independent exact audit of the silver Lambda22 H14 freeze.

This replay deliberately does not import ``silver_lambda22``.  It parses the
vendored Gram matrix, reconstructs its Smith data and pairing-radical
sublattices, asks PARI for exact count-only shell censuses, and checks the
closed-form packet and witness inequalities with standalone fraction pairs.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
GRAM_PATH = ROOT / "data/lattices/lambda22_nebe_sloane_gram.txt"
PRODUCTION_PATH = ROOT / "data/enumerations/silver_lambda22_h14.json"
OUTPUT_PATH = (
    ROOT
    / "data/enumerations/silver_lambda22_h14_independent_hostile_audit.json"
)
EXPECTED_GRAM_FILE_SHA256 = (
    "4f3dbdc96e7aead7fae09a834abdbe2d5319620ece24656720f13184dd85b2f1"
)


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

    def sign(self) -> int:
        sa = (self.a > 0) - (self.a < 0)
        sb = (self.b > 0) - (self.b < 0)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        return sa if self.a * self.a - 2 * self.b * self.b > 0 else sb

    def exact(self) -> dict[str, str]:
        return {"a": str(self.a), "b": str(self.b)}


ZERO = K2.make()
ONE = K2.make(1)
SQRT2 = K2.make(0, 1)
A = SQRT2 - ONE
H14 = K2.make(14, 10)


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
    if any(len(row) != index + 1 for index, row in enumerate(rows)):
        raise AssertionError("malformed lower triangle")
    return Matrix(
        22, 22, lambda i, j: rows[max(i, j)][min(i, j)]
    )


def kernel_basis_mod_prime(gram: Matrix, prime: int) -> Matrix:
    """Build ``{x:Gx=0 mod p}`` without project Smith transports."""

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
    generators = [prime * Matrix.eye(22), *radical]
    basis = hermite_normal_form(Matrix.hstack(*generators))
    if any(int(value) % prime for value in gram * basis):
        raise AssertionError("pairing-radical basis violates the gate")
    return basis


def gp_matrix(matrix: Matrix) -> str:
    return "[" + ";".join(
        ",".join(str(int(matrix[i, j])) for j in range(matrix.cols))
        for i in range(matrix.rows)
    ) + "]"


def qfminim_count_only(gram: Matrix, bounds: tuple[int, ...]) -> dict[int, tuple[int, int]]:
    """Exact PARI totals using ``m=1`` and the default algorithm flag 0."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required")
    commands = [
        f"G={gp_matrix(gram)}",
        "U=qflllgram(G)",
        "Q=U~*G*U",
    ]
    for bound in bounds:
        commands.append(
            f"R=qfminim(Q,{bound},1);print(\"ROW {bound} \",R[1],\" \",R[2])"
        )
    commands.append("quit()")
    process = subprocess.run(
        [gp, "-s", "4000000000", "-fq"],
        input=";".join(commands) + "\n",
        text=True,
        capture_output=True,
        check=True,
        timeout=600,
    )
    output: dict[int, tuple[int, int]] = {}
    for line in process.stdout.splitlines():
        if not line.startswith("ROW "):
            continue
        _, bound, count, maximum = line.split()
        output[int(bound)] = (int(count), int(maximum.split(".", 1)[0]))
    if set(output) != set(bounds):
        raise AssertionError(f"truncated PARI output: {process.stdout!r}")
    return output


def audit(production_path: Path = PRODUCTION_PATH) -> dict[str, object]:
    gram = parse_gram()
    if sha256(GRAM_PATH.read_bytes()).hexdigest() != EXPECTED_GRAM_FILE_SHA256:
        raise AssertionError("the vendored Gram file changed")
    smith = tuple(
        abs(int(value))
        for value in smith_normal_form(gram, domain=ZZ).diagonal()
    )
    if gram.det() != 12 or smith != (1,) * 20 + (2, 6):
        raise AssertionError("Gram determinant/Smith form drifted")

    s2_basis = kernel_basis_mod_prime(gram, 2)
    s3_basis = kernel_basis_mod_prime(gram, 3)
    ordinary_counts = qfminim_count_only(gram, (4, 6, 8, 10))
    s2_counts = qfminim_count_only(s2_basis.T * gram * s2_basis, (10, 12, 16, 20))
    s3_counts = qfminim_count_only(s3_basis.T * gram * s3_basis, (23, 24))
    if ordinary_counts != {
        4: (49_896, 4),
        6: (2_871_528, 6),
        8: (53_204_118, 8),
        10: (520_800_918, 10),
    }:
        raise AssertionError(f"ordinary census drifted: {ordinary_counts!r}")
    if s2_counts != {
        10: (0, 0),
        12: (8_448, 12),
        16: (58_344, 16),
        20: (1_426_920, 20),
    }:
        raise AssertionError(f"S2 census drifted: {s2_counts!r}")
    if s3_counts != {23: (0, 0), 24: (1_782, 24)}:
        raise AssertionError(f"S3 census drifted: {s3_counts!r}")

    ordinary_shell = ordinary_counts[10][0] - ordinary_counts[8][0]
    dyadic_shell = s2_counts[20][0] - s2_counts[16][0]
    if (ordinary_shell, dyadic_shell) != (467_596_800, 1_368_576):
        raise AssertionError("H14 shell cardinality drifted")

    # Complete q=4 type resolution: v^2=10,12,14,16,20 and
    # (u-v)^2=20-v^2.  min(S2)=12 leaves only u=v, v^2=20.
    q4_v_squares = (10, 12, 14, 16, 20)
    q4_difference_squares = tuple(20 - value for value in q4_v_squares)
    if q4_difference_squares != (10, 8, 6, 4, 0):
        raise AssertionError("q4 difference identity drifted")

    # The spatial scaling is (1+sqrt(2))^2=3+2sqrt(2).
    spatial_scale = K2.make(3, 2)
    q2_same = -A * K2.make(6, 4) ** 2 + spatial_scale * K2.make(8)
    q4_same = -A * K2.make(8, 6) ** 2 + spatial_scale * K2.make(14)
    cross = -A * K2.make(6, 4) * K2.make(8, 6) + spatial_scale * K2.make(13)
    witness_norm = K2.make(4) / K2.make(9) - A
    witness_q2 = -A * K2.make(6, 4) + (ONE + SQRT2) * K2.make(5) / K2.make(3)
    witness_q4 = -A * K2.make(8, 6) + (ONE + SQRT2) * K2.make(8) / K2.make(3)
    if not all(value.sign() < 0 for value in (q2_same, q4_same, cross)):
        raise AssertionError("same-height acuteness bound failed")
    if witness_norm.sign() <= 0 or not all(
        value.sign() < 0 for value in (witness_q2, witness_q4)
    ):
        raise AssertionError("strict hyperideal witness failed")

    production = json.loads(production_path.read_text())
    packet = production["first_packet"]
    if packet["ordinary_component"]["shell_count"] != ordinary_shell:
        raise AssertionError("production ordinary count disagrees")
    if packet["dyadic_component"]["shell_count"] != dyadic_shell:
        raise AssertionError("production dyadic count disagrees")
    if packet["total_root_count"] != ordinary_shell + dyadic_shell:
        raise AssertionError("production total packet count disagrees")
    if production["gated_sublattices"]["S3_minimum"] != 24:
        raise AssertionError("production retained the obsolete S3 minimum 36")
    if production["conclusion"] != "H14 prefix has infinite hyperbolic volume":
        raise AssertionError("production prefix conclusion drifted")
    next_target = production["next_target"]
    if "H=36" not in next_target or "H=41" in next_target:
        raise AssertionError("production retained the obsolete H41 frontier")

    return {
        "schema_version": "0.1.0",
        "verdict": "PASS",
        "scope": "independent exact audit of the silver Lambda22 H14 prefix",
        "independence": (
            "No import of silver_lambda22; Gram parsing, mod-p kernels, Smith "
            "form, PARI shell counts, and Q(sqrt(2)) inequalities are replayed "
            "with standalone code."
        ),
        "sources": {
            "production_artifact": str(production_path.relative_to(ROOT)),
            "production_sha256": sha256(production_path.read_bytes()).hexdigest(),
            "gram_file": str(GRAM_PATH.relative_to(ROOT)),
            "gram_file_sha256": sha256(GRAM_PATH.read_bytes()).hexdigest(),
        },
        "gram": {
            "determinant": 12,
            "smith_invariants": list(smith),
            "rank_mod_2": 20,
            "rank_mod_3": 21,
            "S2_index": abs(int(s2_basis.det())),
            "S3_index": abs(int(s3_basis.det())),
        },
        "pari_count_only": {
            "method": "qflllgram then qfminim(Q,b,m=1,default flag=0)",
            "ordinary_cumulative": {str(k): list(v) for k, v in ordinary_counts.items()},
            "S2_cumulative": {str(k): list(v) for k, v in s2_counts.items()},
            "S3_cumulative": {str(k): list(v) for k, v in s3_counts.items()},
        },
        "H14_packet": {
            "ordinary_norm10_shell_count": ordinary_shell,
            "dyadic_S2_norm20_shell_count": dyadic_shell,
            "total_root_count": ordinary_shell + dyadic_shell,
            "q4_necessary_v_squares": list(q4_v_squares),
            "q4_difference_squares": list(q4_difference_squares),
            "q4_surviving_description": "u=v in the complete S2 norm-20 shell",
            "same_height_pairings_at_bounds": {
                "ordinary_dot_8": q2_same.exact(),
                "dyadic_dot_14": q4_same.exact(),
                "cross_dot_13": cross.exact(),
            },
        },
        "hyperideal_witness": {
            "description": "z=(1,d/3), d^2=4",
            "norm": witness_norm.exact(),
            "ordinary_max_pairing": witness_q2.exact(),
            "dyadic_max_pairing": witness_q4.exact(),
            "strict": True,
        },
        "conclusion": {
            "H14_prefix_finite_volume": False,
            "H14_prefix_compact": False,
            "eventual_chamber_status": "unresolved",
            "corrected_next_frontier": "H=36 mixed packet; obsolete H41 claim rejected",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--production", type=Path, default=PRODUCTION_PATH)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    payload = audit(args.production)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"verdict": payload["verdict"], **payload["conclusion"]}, indent=2))


if __name__ == "__main__":
    main()
