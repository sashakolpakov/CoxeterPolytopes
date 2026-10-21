#!/usr/bin/env python3
"""Exact A8-period lattice algebra; root schedules via parent-lock wrapper.

No scientific libraries, worker imports, subprocesses, inner flock or census.
An isometry PASS is NOT a first-endpoint or nonreflectivity certificate.
Use python -B scripts/run_bounded_research_check_20261018.py --parent-lock
--ledger results/NEW.json -- scripts/check_...py --write (or --verify).
"""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import argparse
import json
import signal

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "results/compact_h9_sqrt3_a8_period_isometry_20261020.json"
CAP = 256 * 1024
PINS = {
    "results/compact_h9_sqrt3_a8_period_isometry_20261020.md":
        "f4d80265de3c167c82e741af40d45cb62957d63ff020211e24b7616f4229137b",
    "results/compact_h9_sqrt3_e8_integral_splitting_20261020.md":
        "0e8f62f9354fb3be65738db2bb6b8ffcbfdd864d535a70c10292f71d4210bacb",
    "results/compact_h9_sqrt3_a7_opposite_edge_20261020.json":
        "226a53992f880c06976bbd9b27ae7dd1ad2ef8b5ca7bbf3312201b6b01d40e76",
    "scripts/run_bounded_research_check_20261018.py":
        "4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132",
    "scripts/build_illustrated_ra_manuscript_20261018.py":
        "27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15",
    "scripts/run_lean_batches_signal_safe_20261018.py":
        "c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d",
    "scripts/run_lean_batches_20261018.py":
        "97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c",
}


class F:
    """Q(sqrt3), with rational coefficients and exact physical sign."""
    def __init__(self, a=0, b=0):
        if isinstance(a, F):
            assert not b
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = Q(a), Q(b)

    def __add__(self, other):
        z = F(other)
        return F(self.a + z.a, self.b + z.b)

    __radd__ = __add__

    def __neg__(self):
        return F(-self.a, -self.b)

    def __sub__(self, other):
        return self + -F(other)

    def __rsub__(self, other):
        return F(other) + -self

    def __mul__(self, other):
        z = F(other)
        return F(self.a*z.a + 3*self.b*z.b,
                 self.a*z.b + self.b*z.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        z = F(other)
        d = z.a*z.a - 3*z.b*z.b
        assert d
        return self * F(z.a/d, -z.b/d)

    def __eq__(self, other):
        z = F(other)
        return self.a == z.a and self.b == z.b

    def __bool__(self):
        return bool(self.a or self.b)

    def sign(self):
        a, b = self.a, self.b
        if not a:
            return (b > 0) - (b < 0)
        if not b or (a > 0) == (b > 0):
            return (a > 0) - (a < 0)
        d = a*a - 3*b*b
        assert d
        return ((a > 0) - (a < 0))*((d > 0) - (d < 0))

    def integral(self):
        return self.a.denominator == self.b.denominator == 1


Z, ONE, S = F(), F(1), F(0, 1)
EPS, ETA, LAMBDA, M = 2+S, 1+S, S-1, 4+S


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(col) for col in zip(*a)]


def mm(a, b):
    assert a and b and len(a[0]) == len(b)
    return [[sum((x*y for x, y in zip(row, col)), Z)
             for col in zip(*b)] for row in a]


def mv(a, v):
    return [row[0] for row in mm(a, [[x] for x in v])]


def va(*vectors):
    return [sum(col, Z) for col in zip(*vectors)]


def vs(k, v):
    return [k*x for x in v]


def inverse(a):
    n = len(a)
    assert all(len(row) == n for row in a)
    work = [list(row) + ident for row, ident in zip(a, eye(n))]
    for col in range(n):
        pivot = next((j for j in range(col, n) if work[j][col]), None)
        assert pivot is not None, "singular exact matrix"
        work[col], work[pivot] = work[pivot], work[col]
        d = work[col][col]
        work[col] = [x/d for x in work[col]]
        for j in range(n):
            if j != col and work[j][col]:
                k = work[j][col]
                work[j] = [x-k*y for x, y in zip(work[j], work[col])]
    assert [row[:n] for row in work] == eye(n)
    return [row[n:] for row in work]


def determinant(a):
    work = [list(row) for row in a]
    n, answer = len(a), ONE
    for col in range(n):
        pivot = next((j for j in range(col, n) if work[j][col]), None)
        if pivot is None:
            return Z
        if pivot != col:
            work[pivot], work[col] = work[col], work[pivot]
            answer = -answer
        d = work[col][col]
        answer = answer*d
        for j in range(col+1, n):
            k = work[j][col]/d
            work[j] = [x-k*y for x, y in zip(work[j], work[col])]
    return answer


def integral_matrix(a):
    return all(x.integral() for row in a for x in row)


def digest(path):
    h = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def encode(x):
    if isinstance(x, F):
        return [str(x.a), str(x.b)]
    if isinstance(x, list):
        return [encode(y) for y in x]
    return x


def check():
    own = str(Path(__file__).resolve().relative_to(ROOT))
    pins = {p: digest(ROOT/p) for p in [*PINS, own]}
    assert all(pins[p] == h for p, h in PINS.items()), "fixed input changed"
    prior_path = ROOT / "results/compact_h9_sqrt3_a7_opposite_edge_20261020.json"
    assert prior_path.stat().st_size <= CAP
    prior = json.loads(prior_path.read_text())
    assert prior["status"] == "PASS_COMPLETE_FIRST_A7_OPPOSITE_EDGE_SQRT3_M"
    assert prior["simple_wall_projective_representatives"] == 1
    assert prior["selected_simple_wall"]["profile"] == "E8"
    assert prior["endpoint_M_homogeneous"] == [[99, 72]]+[[72, 45]]*8+[[36, 27]]
    for p, h in prior["sources"].items():
        target = (ROOT/p).resolve()
        assert target.is_relative_to(ROOT)
        assert digest(target) == h, "A7 certificate source changed"
        if p in pins:
            assert pins[p] == h
        pins[p] = h

    # Reconstruct both lattice Grams without relying on the isometry note's matrix.
    c = [[F(2 if i == j else 0) for j in range(8)] for i in range(8)]
    for i, j in [(i, i+1) for i in range(6)]+[(2, 7)]:
        c[i][j] = c[j][i] = F(-1)
    cinv = inverse(c)
    assert determinant(c) == 1
    assert [cinv[i][i] for i in range(8)] == list(map(F, [4,14,30,20,12,6,2,8]))
    g = [row+[Z] for row in c]+[[Z]*8+[1-S]]
    original = [[17+7*S]+[-M]*8]
    original += [[-M]+[F(1+(i == j)) for j in range(8)] for i in range(8)]
    basis = eye(9)
    old_simple = []
    for i in range(7):
        vector = [Z]*9
        vector[i+1], vector[i+2] = F(-1), ONE
        old_simple.append(vector)
    beta_original = [5+3*S]+[4+2*S]*3+[3+2*S]*5
    v_original = [11+8*S]+[8+5*S]*8
    embedding = transpose(old_simple+[beta_original, v_original])
    assert determinant(embedding) in [ONE, -ONE]
    assert mm(mm(transpose(embedding), original), embedding) == g
    weights = [list(col)+[Z] for col in zip(*cinv)]
    roots, v = basis[:8], basis[8]
    beta = roots[7]
    n8 = va(vs(ETA*EPS, v), vs(-1, weights[6]), vs(-ETA, weights[7]))
    v0 = va(vs(11+8*S, v), vs(-2*EPS, weights[7]))
    r5 = va(vs(23+13*S, v), vs(-(4+3*S), weights[1]),
            vs(-1, weights[0]), vs(-1, weights[7]))
    v1 = va(vs(73+42*S, v), vs(-(16+10*S), weights[1]))
    assert mv(embedding, n8) == [Z]*8+[F(-1)]
    assert mv(embedding, v0) == [F(9)]+[M]*8
    assert mv(embedding, vs(9, v)) == list(map(lambda p: F(*p), prior["endpoint_M_homogeneous"][:9]))
    f = [roots[0], r5, beta, roots[2], roots[3], roots[4], roots[5], roots[6]]
    old = roots[:7]+[n8]
    a8 = [[F(2 if i == j else -1 if abs(i-j) == 1 else 0)
           for j in range(8)] for i in range(8)]
    for frame in (old, f):
        assert mm(mm(frame, g), transpose(frame)) == a8
    for frame, apex in ((old, v0), (f, v1)):
        assert mm(frame, [[x] for x in mv(g, apex)]) == [[Z]]*8
        assert sum((x*y for x, y in zip(apex, mv(g, apex))), Z) == 9*(1-S)
    assert sum((x*y for x, y in zip(v0, mv(g, v1))), Z).sign() < 0

    # Construct independently from the two full rational frames.
    left, right = transpose(old+[v0]), transpose(list(reversed(f))+[v1])
    period = mm(right, inverse(left))
    period_inverse = inverse(period)
    assert integral_matrix(period) and integral_matrix(period_inverse)
    assert determinant(period) in [ONE, -ONE]
    assert mm(period, period_inverse) == mm(period_inverse, period) == eye(9)
    assert mm(mm(transpose(period), g), period) == g
    assert all(mv(period, x) == y for x, y in zip(old, reversed(f)))
    assert mv(period, v0) == v1
    assert period[8][8] == 69+40*S and period[7][7] == -44-25*S
    trace = sum((period[i][i] for i in range(9)), Z)
    assert trace == 5+3*S and (trace-9).sign() > 0
    period_original = mm(mm(embedding, period), inverse(embedding))
    assert integral_matrix(period_original)
    assert mm(mm(transpose(period_original), original), period_original) == original

    # Explicit index-nine glue vector, checked separately from matrix inversion.
    weighted = va(*(vs(i+1, f[i]) for i in range(8)))
    glue = vs(F(Q(1, 9)), va(v1, vs(-M, weighted)))
    expected = [F(34,19),F(68,38),F(96,54),F(76,43),
                F(56,32),F(36,21),F(16,10),F(48,27),F(-21,-12)]
    assert glue == expected and all(x.integral() for x in glue)
    assert pins == {p: digest(ROOT/p) for p in pins}, "input changed during check"
    return {
        "status": "PASS_EXACT_INTEGRAL_INFINITE_ORDER_A8_ISOMETRY",
        "scope": "Exact lattice algebra only; alpha2 first-simple-endpoint and same-chamber interface still required",
        "sources": pins, "owns_lock": False, "rank": 9,
        "matrix_e8_basis": encode(period), "inverse_e8_basis": encode(period_inverse),
        "matrix_original_m_basis": encode(period_original),
        "determinant": encode(determinant(period)), "trace": encode(trace),
        "trace_exceeds_rank": True, "glue_vector": encode(glue),
        "original_apex": encode(v0), "candidate_apex": encode(v1),
        "a7_endpoint_interface": "pinned complete census: original A8 to actual E8",
        "alpha2_endpoint_interface": "PENDING separate complete census with same-chamber simple-cone check",
        "nonreflectivity_claim": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--verify", action="store_true")
    modes.add_argument("--stdout", action="store_true")
    args = parser.parse_args()
    if args.write and CERT.exists():
        parser.error("refusing to overwrite existing certificate")
    signal.alarm(45)
    result = check()
    payload = json.dumps(result, sort_keys=True, indent=2)+"\n"
    assert len(payload.encode()) <= CAP
    if args.write:
        with CERT.open("x") as stream:
            stream.write(payload)
    elif args.verify:
        assert CERT.stat().st_size <= CAP
        assert result == json.loads(CERT.read_text()), "saved certificate differs"
    if args.stdout:
        print(payload, end="")
    else:
        print(json.dumps({"status": "PASS", "rank": 9, "trace": [5,3],
                          "global_nonreflectivity": "conditional; not asserted here"}))


if __name__ == "__main__":
    main()
