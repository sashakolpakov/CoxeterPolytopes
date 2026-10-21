#!/usr/bin/env python3
"""Exact local algebra for the all-finite-side single-top prism obstruction.

This uses only stdlib rational arithmetic, not the cap producer or its arithmetic.
The accompanying report supplies the general field/Schur-complement argument.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CAPS = ROOT / "data/enumerations/fresh_h9_face_bridge_20261013_caps.json"
SOURCE = ROOT / "data/enumerations/fresh_h9_face_bridge_20261013_verify.json"
OUT = ROOT / "results/fresh_h9_prism_lift_20261013.json"
Z, O, PHI = (Q(0), Q(0)), (Q(1), Q(0)), (Q(0), Q(1))


def k(a=0, b=0):
    return Q(a), Q(b)


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    a, b = x
    c, d = y
    return a*c+b*d, a*d+b*c+b*d


def conjugate(x):
    return x[0]+x[1], -x[1]


def norm(x):
    ans = mul(x, conjugate(x))
    assert ans[1] == 0
    return ans[0]


def div(x, y):
    n = norm(y)
    assert n
    z = mul(x, conjugate(y))
    return z[0]/n, z[1]/n


def sign(x):
    """Sign of a+b*phi, reduced to a rational square comparison."""
    u, v = x[0]+x[1]/2, x[1]/2
    if v == 0:
        return (u > 0)-(u < 0)
    if u == 0:
        return (v > 0)-(v < 0)
    if (u > 0) == (v > 0):
        return (u > 0)-(u < 0)
    delta = u*u-5*v*v
    assert delta != 0
    return ((u > 0)-(u < 0))*((delta > 0)-(delta < 0))


def enc(x):
    return {"a": str(x[0]), "b": str(x[1])}


def parse(x):
    return k(x["a"], x["b"])


def matmul(A, B):
    return [[sumk(mul(x, y) for x, y in zip(row, col))
             for col in zip(*B)] for row in A]


def sumk(xs):
    ans = Z
    for x in xs:
        ans = add(ans, x)
    return ans


def inverse(A):
    n = len(A)
    B = [list(row)+[O if i == j else Z for j in range(n)]
         for i, row in enumerate(A)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if B[i][j] != Z)
        B[j], B[pivot] = B[pivot], B[j]
        d = B[j][j]
        B[j] = [div(x, d) for x in B[j]]
        for i in range(n):
            if i != j:
                c = B[i][j]
                B[i] = [sub(x, mul(c, y)) for x, y in zip(B[i], B[j])]
    return [row[n:] for row in B]


def ldl_pivots(A):
    B = [list(row) for row in A]
    pivots = []
    for i in range(len(B)):
        d = B[i][i]
        assert d != Z
        pivots.append(d)
        for j in range(i+1, len(B)):
            for h in range(i+1, len(B)):
                B[j][h] = sub(B[j][h], div(mul(B[j][i], B[i][h]), d))
    return pivots


def square_Q(x):
    if x < 0:
        return False
    return isqrt(x.numerator)**2 == x.numerator and isqrt(x.denominator)**2 == x.denominator


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    source = json.loads(SOURCE.read_text())
    caps = json.loads(CAPS.read_text())
    assert source["source_sha256"] == sha(CAPS)
    assert source["status"] == "PASS_EXACT_96_WALL_COMPACT_H4_PRISM_AND_TARGET_CUT"
    assert source["prism_common_type"] == "5A1"
    assert source["prism_vertices"] == 8
    assert source["prism_old_base"] == "wall17"
    assert source["prism_new_base"] == "facePrism"
    assert caps["caps"][0]["complement_wall"] == "wall17"
    assert caps["caps"][0]["branch"] == -1

    # T=[seed6,crosspair0,nextH4,compact0], a 3--5--3 tetrahedron.
    L = [[k(2), k(-1), Z, Z], [k(-1), k(2), neg(PHI), Z],
         [Z, neg(PHI), k(2), k(-1)], [Z, Z, k(-1), k(2)]]
    Linv = inverse(L)
    tetra_pivots = ldl_pivots(L)
    assert [sign(x) for x in tetra_pivots] == [1, 1, 1, -1]
    assert matmul(L, Linv) == [[O if i == j else Z for j in range(4)] for i in range(4)]
    assert all(sign(x) < 0 for row in Linv for x in row)
    assert Linv[3] == [k('-4/11', '-5/11'), k('-8/11', '-10/11'),
                       k('-3/11', '-12/11'), k('4/11', '-6/11')]
    T = ["seed6", "crosspair0", "nextH4", "compact0"]
    construction = source["new_normals"][0]["construction"]
    coeff = dict(zip(caps["caps"][0]["source_walls"],
                     map(parse, construction["particular_coefficients"])))
    assert [coeff[name] for name in T] == [neg(x) for x in Linv[3]]
    assert all(x == Z for name, x in coeff.items() if name not in T)
    assert all(parse(x) == Z for x in construction["complement_projection_coefficients"])
    a = k('9/11', '3/11')
    assert parse(construction["radicand"]) == a
    assert parse(construction["particular_norm"]) == Linv[3][3]
    assert parse(construction["complement_norm"]) == k(2)
    assert construction["branch"] == -1
    assert add(Linv[3][3], mul(k(2), a)) == k(2)
    assert sign(a) > 0 and sign(conjugate(a)) > 0

    # The polynomial for sqrt(a) is 11X^4-21X^2+9.
    assert add(sub(mul(k(11), mul(a, a)), mul(k(21), a)), k(9)) == Z
    assert add(a, conjugate(a)) == k('21/11')
    assert norm(a) == Q(9, 11)
    assert norm(div(a, k(11))) == Q(9, 1331)
    nonsquares = [norm(a), norm(div(a, k(11))), Q(11), Q(11, 5)]
    assert all(not square_Q(x) for x in nonsquares)

    # The actual existing corner [T1,T2,T3,r] is H4.
    names = [T[1], T[2], T[3], "facePrism"]
    vertex = next(v for v in source["actual_finite_vertices"]
                  if set(names) <= set(v["incident"]))
    indices = [vertex["incident"].index(name) for name in names]
    C = [[parse(vertex["gram"][i][j]) for j in indices] for i in indices]
    expected = [[k(2), neg(PHI), Z, Z], [neg(PHI), k(2), k(-1), Z],
                [Z, k(-1), k(2), k(-1)], [Z, Z, k(-1), k(2)]]
    assert C == expected
    pivots = ldl_pivots(C)
    assert all(sign(x) > 0 for x in pivots)
    det = O
    for p in pivots:
        det = mul(det, p)
    assert det == k(5, -3)
    Cinv = inverse(C)
    assert all(sign(x) > 0 for row in Cinv for x in row)
    diag = [Cinv[i][i] for i in range(4)]
    assert diag == [k(8, 12), k(12, 18), k(6, 8), k(2, 2)]
    assert all(sign(sub(x, k(2))) > 0 for x in diag)
    common = source["prism_common_facets"]
    commonids = [vertex["incident"].index(name) for name in common]
    assert [[parse(vertex["gram"][i][j]) for j in commonids] for i in commonids] == [
        [k(2) if i == j else Z for j in range(5)] for i in range(5)]
    assert all(parse(vertex["gram"][i][j]) == Z for i in indices for j in commonids)
    return {
        "status": "PASS_ALL_LABEL_SINGLE_TOP_H5_PRISM_OBSTRUCTION",
        "source_caps": str(CAPS.relative_to(ROOT)), "source_caps_sha256": sha(CAPS),
        "source_verified_geometry": str(SOURCE.relative_to(ROOT)), "source_geometry_sha256": sha(SOURCE),
        "tetrahedron_order": T, "tetrahedron_inverse": [[enc(x) for x in row] for row in Linv],
        "tetrahedron_ldl_pivots": list(map(enc, tetra_pivots)),
        "a": enc(a), "sqrt_a_polynomial_ascending": [9, 0, -21, 0, 11],
        "nonsquare_rational_checks": list(map(str, nonsquares)),
        "H4_corner": names, "H4_gram": [[enc(x) for x in row] for row in C],
        "H4_ldl_pivots": list(map(enc, pivots)), "H4_determinant": enc(det),
        "H4_inverse": [[enc(x) for x in row] for row in Cinv],
        "H4_inverse_diagonal": list(map(enc, diag)),
        "deductions": [
            "Q(sqrt(a)) is nonnormal quartic, hence sqrt(a) is outside Q^ab.",
            "All six finite Coxeter pairings belong to Q^ab; the side relation forces p_b=0.",
            "At the H4 corner any nonzero finite coupling has magnitude at least 1; the Schur complement is negative.",
            "A finite or ideal endpoint therefore forces p_T1=p_T2=p_T3=p_r=0.",
            "The side relation with nonzero (L^-1)[3,0] forces p_T0=0.",
            "The proposed top is orthogonal to the old H4 ambient space, so its H5 hyperplane repeats the freed A1 base."
        ],
        "scope": "No full-dimensional compact or finite-volume H5 prism obtained by retaining all six side facets of this exact H4 prism and adding a single top with finite Coxeter angles to all six, closing every old vertex edge. This does not exclude dotted side relations, several caps, facet replacements, or other H9 completions.",
        "independence": "Standalone stdlib Fraction algebra. Existing prism incidence/actual-face assertions are inherited from the hashed source geometry; no cap producer or its arithmetic is imported. The all-label argument uses the accompanying mathematical proof, not a finite label census."
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.verify:
        assert result == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])
