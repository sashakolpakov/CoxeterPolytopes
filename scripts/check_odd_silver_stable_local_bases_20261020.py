#!/usr/bin/env python3
"""Exact local-basis algebra only. The root's wrapper owns gates and lock.

No search, scientific imports, global integral-isometry matrix, or execution
of the genus/strong-approximation proof. Importing this file runs no check.
"""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse
import json
import signal

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "results/odd_silver_stable_local_bases_20261020.json"
PINS = {
    "results/odd_silver_stable_integral_isometry_and_uniform_obstruction_20261020.md":
        "2ad82b9a52c7b81364a070b985932932eab72438922e7b0d717b392e5e4c5ee3",
    "results/compact_h9_odd_silver_diagonal_two_edges_20261020.md":
        "85095d203aea614961425bf356cfbad21529cd9f6f01d32b77f73cc5b2ce234c",
    "results/fresh_silver_e8_f4_identification_20261016.md":
        "1c9b7df422c598d2cea8d9067b0d9c17047fe93d21b9142757fb1ec8162a4353",
    "results/fresh_silver_e8_chamber_scope_20261016.md":
        "9bb9e4d7f7d02f6bd05f7ecd0589fec2e892ce18718443b4f3299a047b5e2f2e",
}


def q(a=0, b=0):
    return F(a), F(b)


Z, ONE, S, ALPHA = q(), q(1), q(0, 1), q(1, 1)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def neg(a):
    return -a[0], -a[1]


def mul(a, b):
    return a[0]*b[0]+2*a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def inv(a):
    norm = a[0]*a[0]-2*a[1]*a[1]
    assert norm
    return a[0]/norm, -a[1]/norm


def dot(a, b):
    out = Z
    for x, y in zip(a, b, strict=True):
        out = add(out, mul(x, y))
    return out


def transpose(a):
    return [list(row) for row in zip(*a, strict=True)]


def mm(a, b):
    return [[dot(row, col) for col in transpose(b)] for row in a]


def identity(n):
    return [[q(i == j) for j in range(n)] for i in range(n)]


def determinant_inverse(a):
    n = len(a)
    rows = [list(row)+unit for row, unit in zip(a, identity(n), strict=True)]
    determinant = ONE
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j] != Z)
        if pivot != j:
            rows[j], rows[pivot] = rows[pivot], rows[j]
            determinant = neg(determinant)
        value = rows[j][j]
        determinant = mul(determinant, value)
        rows[j] = [mul(x, inv(value)) for x in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [add(x, neg(mul(factor, y)))
                           for x, y in zip(rows[i], rows[j], strict=True)]
    assert [row[:n] for row in rows] == identity(n)
    inverse = [row[n:] for row in rows]
    assert mm(a, inverse) == mm(inverse, a) == identity(n)
    return determinant, inverse


def denominators_supported(matrix, prime):
    for row in matrix:
        for value in row:
            for coefficient in value:
                denominator = coefficient.denominator
                while denominator % prime == 0:
                    denominator //= prime
                assert denominator == 1


def serial_matrix(matrix):
    return [[[str(a), str(b)] for a, b in row] for row in matrix]


def digest(path):
    h = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def check():
    paths = list(PINS)+[str(Path(__file__).resolve().relative_to(ROOT))]
    sources = {path: digest(ROOT/path) for path in paths}
    assert all(sources[path] == pin for path, pin in PINS.items())
    alpha_inverse = inv(ALPHA)
    assert alpha_inverse == q(-1, 1) and mul(S, inv(S)) == ONE
    third = q(F(1, 3))
    columns = [
        [mul(third, v) for v in (ALPHA, q(1, -1), ONE, ONE, ONE)],
        [Z, Z, ONE, q(-1), Z], [Z, Z, Z, ONE, q(-1)],
        [mul(third, v) for v in (alpha_inverse, ALPHA, neg(S), neg(S), q(0, 2))],
        [mul(third, v) for v in (q(1, -2), q(-2, -1), ONE, ONE, ONE)],
    ]
    other_columns = [
        [Z, Z, Z, Z, ONE], [Z, ONE, q(-1), Z, Z],
        [Z, Z, ONE, q(-1), Z], [Z, Z, Z, S, Z],
        [mul(inv(S), v) for v in (ONE, q(-1), q(-1), q(-1), Z)],
    ]
    target = identity(5)
    for i in range(1, 5):
        target[i][i] = q(2)
    for i, value in enumerate((q(-1), neg(S), q(-1)), start=1):
        target[i][i+1] = target[i+1][i] = value
    f4 = [row[1:] for row in target[1:]]
    f4_det, _ = determinant_inverse(f4)
    assert f4_det == ONE
    basis_records = []
    for name, cols, prime in (("away_from_3", columns, 3),
                              ("away_from_sqrt2", other_columns, 2)):
        basis = transpose(cols)
        gram = mm(transpose(basis), basis)
        assert gram == target
        determinant, inverse = determinant_inverse(basis)
        assert determinant in (ONE, q(-1))
        denominators_supported(basis, prime)
        denominators_supported(inverse, prime)
        basis_records.append({"name": name, "rational_denominator_prime": prime,
                              "matrix": serial_matrix(basis), "inverse": serial_matrix(inverse),
                              "determinant": [str(x) for x in determinant]})
    # Actual rank-six ambient norms for the two literal unit-norm vectors.
    def ambient_norm(v):
        return add(neg(mul(ALPHA, mul(v[0], v[0]))), dot(v[1:], v[1:]))
    assert ambient_norm([Z, ONE, q(2), Z, Z, Z]) == q(5)
    assert ambient_norm([ONE, Z, Z, Z, Z, Z]) == neg(ALPHA)
    # Universal identity in O[r]/(r^2+7), not a finite numerical sample.
    minus_one_norm = [q(6), Z, ONE]       # norm of (r,1,1,2) in I4
    assert add(minus_one_norm[0], mul(q(-7), minus_one_norm[2])) == q(-1)
    assert minus_one_norm[1] == Z and (-7) % 8 == 1
    assert all(digest(ROOT/path) == pin for path, pin in sources.items()), "source changed"
    return {"status": "EXACT_LOCAL_BASES_AND_UNIT_WITNESSES_PASS", "owns_lock": False,
            "sources": sources, "field": "Q(sqrt2); exact Fraction pairs",
            "target_gram": serial_matrix(target), "f4_determinant": "1", "bases": basis_records,
            "unit_spinor_witness_norms": ["5", "-1-sqrt2", "-1 conditional on r^2=-7"],
            "conditional_minus_one": {"formal_norm_coefficients": serial_matrix([minus_one_norm])[0],
                                      "relation": "r^2+7", "unit_square_mod8": 1},
            "unverified_written_inputs": ["2-adic square criterion", "local spinor norm theory",
                "PID and sign reduction", "Hasse-Minkowski", "Kneser-Platonov strong approximation",
                "global integral isometry exists; no global matrix constructed", "lower-face geometry",
                "inherited closed silver E8 nonreflectivity"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    for flag in ("write", "verify", "stdout"):
        mode.add_argument("--"+flag, action="store_true")
    args = parser.parse_args()
    if args.write and CERT.exists():
        parser.error("refusing to overwrite existing certificate")
    signal.alarm(45)
    data = check()
    payload = json.dumps(data, sort_keys=True, indent=2)+"\n"
    assert len(payload.encode()) <= 65536
    if args.write:
        with CERT.open("x") as stream:
            stream.write(payload)
    elif args.verify:
        assert CERT.stat().st_size <= 65536
        assert data == json.loads(CERT.read_text()), "saved certificate differs"
    print(payload if args.stdout else json.dumps({"status": "PASS", "bases": 2,
          "spinor_witnesses": 3, "mode": "write" if args.write else "verify"}),
          end="" if args.stdout else "\n")


if __name__ == "__main__":
    main()
