#!/usr/bin/env python3
"""Exact universal tensor algebra only; wrapper owns resource gates/lock.

No geometry, code enumeration, optimizer, samples, floats or scientific imports.
Root schedules execution. Merely importing this source runs no check.
"""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import argparse
import json
import signal

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "results/right_angled_sharp_dimension_sos_20261020.json"
PROOF = "results/right_angled_dimension8_descending_height_proof_20261020.md"
PROOF_SHA = "9470fa1fa3a973a12c412fbc1257537655536dba4ffc4e9ac23c320615babdbe"
ZERO = (0,)*12


def constant(value):
    return {ZERO: Q(value)} if value else {}


def variable(index):
    powers = list(ZERO)
    powers[index] = 1
    return {tuple(powers): Q(1)}


def add(*polys):
    out = {}
    for poly in polys:
        for term, value in poly.items():
            out[term] = out.get(term, Q(0))+value
    return {term: value for term, value in out.items() if value}


def scale(value, poly):
    return {term: Q(value)*coefficient for term, coefficient in poly.items() if value}


def mul(a, b):
    out = {}
    for left, x in a.items():
        for right, y in b.items():
            term = tuple(i+j for i, j in zip(left, right))
            out[term] = out.get(term, Q(0))+x*y
    return {term: value for term, value in out.items() if value}


def dot(a, b):
    assert len(a) == len(b)
    return add(*(mul(x, y) for x, y in zip(a, b)))


def tensors(coords):
    rank2, rank3 = [], []
    for a in range(6):
        for b in range(6):
            rank2.append(add(mul(coords[a], coords[b]), constant(-Q(a == b, 6))))
            for c in range(6):
                sym = add(coords[a] if b == c else {}, coords[b] if a == c else {}, coords[c] if a == b else {})
                rank3.append(add(mul(mul(coords[a], coords[b]), coords[c]), scale(-Q(1, 8), sym)))
    return rank2, rank3


def hash_file(path):
    h = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def poly_summary(poly):
    terms = [[list(key), str(value)] for key, value in sorted(poly.items())]
    return {"terms": len(terms), "degree": max((sum(k) for k in poly), default=0), "sha256": sha256(json.dumps(terms, separators=(",", ":")).encode()).hexdigest()}


def check():
    paths = [PROOF, str(Path(__file__).resolve().relative_to(ROOT))]
    pins = {path: hash_file(ROOT/path) for path in paths}
    assert pins[PROOF] == PROOF_SHA, "reviewed proof source changed"
    u, v = [variable(i) for i in range(6)], [variable(i+6) for i in range(6)]
    su, sv, t = dot(u, u), dot(v, v), dot(u, v)
    one, two = constant(1), constant(2)
    au, tu = tensors(u)
    av, tv = tensors(v)
    assert len(au) == len(av) == 36 and len(tu) == len(tv) == 216
    aa, tt = dot(au, av), dot(tu, tv)
    sum_norm = add(su, sv)
    expected_aa = add(mul(t, t), scale(-Q(1, 6), sum_norm), constant(Q(1, 6)))
    expected_tt = add(mul(mul(t, t), t), scale(-Q(3, 8), mul(t, sum_norm)), scale(Q(3, 8), t))
    assert aa == expected_aa and tt == expected_tt
    # Exact tracelessness modulo the two unit-sphere ideal generators.
    for coords, norm, a_tensor, t_tensor in ((u, su, au, tu), (v, sv, av, tv)):
        unit_generator = add(norm, constant(-1))
        assert add(*(a_tensor[6*i+i] for i in range(6))) == unit_generator
        for a in range(6):
            assert add(*(t_tensor[36*a+6*b+b] for b in range(6))) == mul(coords[a], unit_generator)
        assert add(dot(a_tensor, a_tensor), constant(-Q(5, 6))) == mul(unit_generator, add(norm, constant(Q(2, 3))))
        assert add(dot(t_tensor, t_tensor), constant(-Q(5, 8))) == mul(unit_generator, add(mul(norm, norm), scale(Q(1, 4), norm), constant(Q(5, 8))))
    factored = mul(mul(add(t, constant(Q(1, 2))), add(t, constant(Q(1, 2)))), add(t, constant(-Q(1, 4))))
    expanded = add(mul(mul(t, t), t), scale(Q(3, 4), mul(t, t)), constant(-Q(1, 16)))
    assert factored == expanded
    kernel = add(constant(Q(1, 16)), scale(Q(3, 8), t), scale(Q(3, 4), aa), tt)
    remainder = mul(add(two, scale(-1, su), scale(-1, sv)), add(scale(Q(3, 8), t), constant(Q(1, 8))))
    assert add(kernel, scale(-1, expanded)) == remainder
    # The displayed remainder is explicitly in (su-1,sv-1), not sample evidence.
    g = add(scale(Q(3, 8), t), constant(Q(1, 8)))
    assert remainder == add(scale(-1, mul(add(su, constant(-1)), g)), scale(-1, mul(add(sv, constant(-1)), g)))
    p2 = scale(Q(1, 5), add(scale(6, mul(t, t)), constant(-1)))
    p3 = scale(Q(1, 5), add(scale(8, mul(mul(t, t), t)), scale(-3, t)))
    assert add(scale(Q(6, 5), aa), scale(-1, p2)) == scale(Q(1, 5), add(two, scale(-1, sum_norm)))
    assert add(scale(Q(8, 5), tt), scale(-1, p3)) == scale(Q(3, 5), mul(t, add(two, scale(-1, sum_norm))))
    assert expanded == add(constant(Q(1, 16)), scale(Q(3, 8), t), scale(Q(5, 8), p2), scale(Q(5, 8), p3))
    sos_coefficients = [Q(1, 16), Q(3, 8), Q(3, 4), Q(1)]
    assert all(c > 0 for c in sos_coefficients)
    diagonal, lower = Q(27, 16), Q(1, 16)
    assert (1+Q(1, 2))**2*(1-Q(1, 4)) == diagonal and diagonal/lower == 27
    four_norm = 4+12*(-Q(1, 2))
    assert four_norm == -2 < 0
    assert pins == {path: hash_file(ROOT/path) for path in paths}, "input changed during check"
    return {"status": "EXACT_UNIVERSAL_TENSOR_ALGEBRA_PASS", "scope": "Q[u1,...,u6,v1,...,v6]; explicit unit-ideal remainders; no geometry or spherical-code enumeration", "sources": pins, "dimension": 6, "contractions": {"rank2": 36, "rank3": 216}, "polynomials": {"rank2_product": poly_summary(aa), "rank3_product": poly_summary(tt), "scalar_polynomial": poly_summary(expanded), "sos_kernel": poly_summary(kernel), "unit_ideal_remainder": poly_summary(remainder)}, "sos_coefficients": list(map(str, sos_coefficients)), "gegenbauer_coefficients": ["1/16", "3/8", "5/8", "5/8"], "diagonal": str(diagonal), "constant_lower": str(lower), "ratio": 27, "strict_endpoint_four_vector_norm": str(four_norm), "owns_lock": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--stdout", action="store_true")
    args = parser.parse_args()
    if args.write and CERT.exists():
        parser.error("refusing to overwrite existing certificate")
    signal.alarm(45)
    result = check()
    payload = json.dumps(result, sort_keys=True, indent=2)+"\n"
    assert len(payload.encode()) <= 262144
    if args.write:
        with CERT.open("x") as stream:
            stream.write(payload)
    elif args.verify:
        assert CERT.stat().st_size <= 262144
        assert result == json.loads(CERT.read_text()), "saved certificate differs"
    if args.stdout:
        print(payload, end="")
    else:
        print(json.dumps({"status": "PASS", "mode": "write" if args.write else "verify", "rank2": 36, "rank3": 216, "ratio": 27, "strict_endpoint": -2}))


if __name__ == "__main__":
    main()
