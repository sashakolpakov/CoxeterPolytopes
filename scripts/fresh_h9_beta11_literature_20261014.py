#!/usr/bin/env python3
"""Exact lattice comparison and a sourced McLeod symmetry proposal.

Standalone Fraction arithmetic in Q(phi); no project arithmetic imports.
This does not certify reflectivity, nonreflectivity, or chamber stabilization.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREFIX = "fresh_h9_beta11_literature_20261014"


class K:
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        if isinstance(a, K):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)

    def __add__(self, other):
        o = K(other)
        return K(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-K(other))

    def __rsub__(self, other):
        return K(other) - self

    def __mul__(self, other):
        o = K(other)
        return K(self.a * o.a + self.b * o.b,
                 self.a * o.b + self.b * o.a + self.b * o.b)

    __rmul__ = __mul__

    def conjugate(self):
        return K(self.a + self.b, -self.b)

    def norm(self):
        return self.a * self.a + self.a * self.b - self.b * self.b

    def __truediv__(self, other):
        o = K(other)
        assert o.norm() != 0
        q = self * o.conjugate()
        return K(q.a / o.norm(), q.b / o.norm())

    def __eq__(self, other):
        o = K(other)
        return (self.a, self.b) == (o.a, o.b)

    def __bool__(self):
        return bool(self.a or self.b)

    def integral(self):
        return self.a.denominator == self.b.denominator == 1

    def sign(self):
        # 2(a+b phi) = (2a+b) + b sqrt(5).
        u, v = 2 * self.a + self.b, self.b
        su, sv = (u > 0) - (u < 0), (v > 0) - (v < 0)
        if not su:
            return sv
        if not sv or su == sv:
            return su
        return su if u * u > 5 * v * v else sv

    def encoded(self):
        return [str(self.a), str(self.b)]


Z, O, PHI = K(), K(1), K(0, 1)


def mat(rows):
    return [[K(x) for x in row] for row in rows]


def eye(n):
    return mat([[int(i == j) for j in range(n)] for i in range(n)])


def tr(A):
    return [list(row) for row in zip(*A)]


def mm(A, B):
    return [[sum((x * y for x, y in zip(row, col)), Z)
             for col in zip(*B)] for row in A]


def inv(A):
    n = len(A)
    M = [list(row) + ident for row, ident in zip(A, eye(n))]
    for j in range(n):
        p = next(i for i in range(j, n) if M[i][j])
        M[j], M[p] = M[p], M[j]
        q = M[j][j]
        M[j] = [x / q for x in M[j]]
        for i in range(n):
            if i != j and M[i][j]:
                q = M[i][j]
                M[i] = [x - q * y for x, y in zip(M[i], M[j])]
    assert [row[:n] for row in M] == eye(n)
    return [row[n:] for row in M]


def det(A):
    M = [list(row) for row in A]
    ans = O
    for j in range(len(A)):
        p = next((i for i in range(j, len(A)) if M[i][j]), None)
        if p is None:
            return Z
        if p != j:
            M[j], M[p] = M[p], M[j]
            ans = -ans
        q = M[j][j]
        ans = ans * q
        for i in range(j + 1, len(A)):
            ratio = M[i][j] / q
            M[i] = [x - ratio * y for x, y in zip(M[i], M[j])]
    return ans


def diagonal(values):
    return [[K(x) if i == j else Z for j in range(len(values))]
            for i, x in enumerate(values)]


def enc(x):
    if isinstance(x, K):
        return x.encoded()
    if isinstance(x, list):
        return [enc(y) for y in x]
    return x


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def run():
    host_path = "results/fresh_h9_h8_cap_height_20261014_host_basis.json"
    host = json.loads((ROOT / host_path).read_text())
    for path, digest in host["sources"].items():
        assert sha(path) == digest, path
    def parse(x):
        assert x[1] == ["0", "0"]
        a, b = map(F, x[0])
        return K(a - b, 2 * b)
    G = [[parse(x) for x in row] for row in host["gram"]]
    E8 = mat([[2,0,-1,0,0,0,0,0], [0,2,0,-1,0,0,0,0],
              [-1,0,2,-1,0,0,0,0], [0,-1,-1,2,-1,0,0,0],
              [0,0,0,-1,2,-1,0,0], [0,0,0,0,-1,2,-1,0],
              [0,0,0,0,0,-1,2,-1], [0,0,0,0,0,0,-1,2]])
    expected = diagonal([-2 * PHI] + [0] * 8 + [2])
    for i in range(8):
        for j in range(8):
            expected[i + 1][j + 1] = E8[i][j]
    expected[0][9] = expected[9][0] = K(-1)
    assert G == expected and det(E8) == 1
    beta, d = 3 + PHI, 1 + 4 * PHI
    assert PHI * beta == d and beta.norm() == 11 and det(G) == -d
    assert det(G).norm() == -11
    assert all(x.integral() for row in G for x in row)
    assert all((G[i][i] / 2).integral() for i in range(10))
    assert (8 * 8 - 8 - 1) % 11 == 0 and (3 + 8) % 11 == 0

    # Current coordinates (x, eight E8 coefficients, t); old H11 is (x,-t,E8).
    P = mat([[0] * 10 for _ in range(10)])
    P[0][0], P[9][1] = O, K(-1)
    for j in range(8):
        P[j + 1][j + 2] = O
    old = diagonal([-2 * PHI, 2] + [0] * 8)
    old[0][1] = old[1][0] = O
    for i in range(8):
        for j in range(8):
            old[i + 2][j + 2] = E8[i][j]
    assert mm(tr(P), mm(G, P)) == old and det(P) in (O, K(-1))

    # Standard Euclidean realization of this exact E8 simple-root ordering.
    R = mat([[F(1,2),F(-1,2),F(-1,2),F(-1,2),F(-1,2),F(-1,2),F(-1,2),F(1,2)],
             [1,1,0,0,0,0,0,0], [-1,1,0,0,0,0,0,0],
             [0,-1,1,0,0,0,0,0], [0,0,-1,1,0,0,0,0],
             [0,0,0,-1,1,0,0,0], [0,0,0,0,-1,1,0,0],
             [0,0,0,0,0,-1,1,0]])
    assert mm(R, tr(R)) == E8 and det(R) in (O, K(-1))
    H = mat([[0] * 8 for _ in range(8)])
    for j in range(0, 8, 2):
        H[j][j] = H[j][j + 1] = H[j + 1][j] = O
        H[j + 1][j + 1] = K(-1)
    Q = mm(inv(tr(R)), H)
    assert all(x.integral() for row in Q for x in row)
    T = mat([[0] * 10 for _ in range(10)])
    T[0][0], T[9][0], T[9][1] = K(2), O, O
    for i in range(8):
        for j in range(8):
            T[i + 1][j + 2] = Q[i][j]
    diagonal_form = diagonal([-2 * d] + [2] * 9)
    assert mm(tr(T), mm(G, T)) == diagonal_form
    assert det(T) in (K(32), K(-32))

    euclidean_root = [[K(x)] for x in [1,0,-1,0,0,0,0,0]]
    e8_root = mm(inv(tr(R)), euclidean_root)
    r = [[Z]] + e8_root + [[Z]]
    assert all(row[0].integral() for row in r)
    assert mm(tr(r), mm(G, r)) == [[K(2)]]
    reflection = [[x - y for x, y in zip(row, prod)]
                  for row, prod in zip(eye(10), mm(r, mm(tr(r), G)))]
    assert all(x.integral() for row in reflection for x in row)
    assert mm(tr(reflection), mm(G, reflection)) == G
    m = [[T[i][2]] for i in range(10)]
    reflected = mm(reflection, m)
    reflected_in_M = mm(inv(T), reflected)
    assert any(not row[0].integral() for row in reflected_in_M)

    # McLeod 2013, Appendix B.2.1, printed p.108 / PDF page 116 (index115).
    # Matrix entries are copied literally, then the isometry is independently checked.
    g4 = [[K(463,748),K(-128,-208),K(-78,-126),K(-78,-126)],
          [K(1256,2032),K(-349,-564),K(-212,-342),K(-212,-342)],
          [K(102,166),K(-28,-46),K(-17,-28),K(-18,-28)],
          [K(102,166),K(-28,-46),K(-18,-28),K(-17,-28)]]
    f4 = diagonal([-d, 1, 1, 1])
    assert mm(tr(g4), mm(f4, g4)) == f4
    assert det(g4) in (O, K(-1))
    trace4 = sum((g4[i][i] for i in range(4)), Z)
    assert (trace4 - 4).sign() > 0 and g4[0][0].sign() > 0
    g10 = eye(10)
    for i in range(4):
        for j in range(4):
            g10[i][j] = g4[i][j]
    assert mm(tr(g10), mm(diagonal_form, g10)) == diagonal_form
    transported = mm(T, mm(g10, inv(T)))
    assert mm(tr(transported), mm(G, transported)) == G
    power = eye(10)
    first_integral_power = None
    for exponent in range(1, 61):
        power = mm(power, transported)
        if all(x.integral() for row in power for x in row):
            assert det(power) in (O, K(-1))
            first_integral_power = exponent
            break
    assert first_integral_power is not None
    assert mm(tr(power), mm(G, power)) == G

    sources = [f"scripts/{PREFIX}.py", host_path,
               "src/coxeter_polytopes/search/h11_vinberg.py",
               "results/h11_vinberg_prefix.md",
               "results/h11_h19_finite_index_obstruction.md"]
    return {
        "status": "PASS_HOST_COMPARISON_NO_REFLECTIVITY_TRANSFER",
        "sources": {p: sha(p) for p in sources},
        "field_encoding": "[a,b] denotes a+b*phi, phi^2=phi+1",
        "gram": enc(G), "old_h11_gram": enc(old), "old_h11_isometry": enc(P),
        "beta": enc(beta), "determinant": enc(-d), "determinant_norm": -11,
        "prime_residue": {"prime_norm": 11, "phi_mod_beta": 8},
        "proper_integral_overlattice": False,
        "unit_determinant_rank10_rational_similarity": False,
        "even_lattice": True,
        "diagonal_sublattice": {
            "basis_columns": enc(T), "gram": enc(diagonal_form),
            "basis_determinant": enc(det(T)), "absolute_index": 1024,
            "euclidean_E8_simple_rows": enc(R),
            "reflection_in_L_not_in_M": {
                "root": enc(r), "reflection": enc(reflection),
                "M_generator": enc(m), "image": enc(reflected),
                "image_in_M_coordinates": enc(reflected_in_M)}},
        "mcleod_symmetry_proposal": {
            "source_url": "https://etheses.durham.ac.uk/id/eprint/7743/1/thesis31072013.pdf#page=116",
            "printed_page": 108, "appendix": "B.2.1", "matrix4": enc(g4),
            "trace4": enc(trace4), "determinant4": enc(det(g4)),
            "transported_matrix10": enc(transported),
            "first_integral_power_checked": first_integral_power,
            "integral_power_matrix10": enc(power),
            "preserves_L": True, "infinite_order": True,
            "full_reflection_chamber_stabilization": "NOT_CHECKED"},
        "scope": [
            "The exact host is the existing H11 lattice, not a new host.",
            "Prime determinant excludes integral overlattices and unit-determinant rational similarity.",
            "McLeod diagonal-lattice nonreflectivity is not transferred by the dyadic sublattice.",
            "The transported infinite-order integral isometry is only a proposal for a full mirror-separation test.",
            "No reflectivity or nonreflectivity conclusion for the full norm-2/norm-2beta group."]}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()
    out = ROOT / "results" / f"{PREFIX}.json"
    answer = run()
    if args.verify:
        assert answer == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(answer, indent=2) + "\n")
    print(answer["status"], "McLeod integral power", answer["mcleod_symmetry_proposal"]["first_integral_power_checked"])
