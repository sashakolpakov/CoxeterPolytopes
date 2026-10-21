#!/usr/bin/env python3
"""Exact necessary isolated-cap census; see the separately pinned proof.

Run through run_bounded_research_check_20261018.py --parent-lock.
This worker never identifies norm compatibility with geometric existence.
"""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
import hashlib
import json
from math import isqrt
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOMAIN = "results/arithmetic_whole_a18_point52_actual_field_degree_reduction_20261021.md"
DOMAIN_SHA = "b57aa75adbde5c1a1b6f78354c58bba47d03b2047a923438b133be9aee5f2707"
PROOF = "results/arithmetic_whole_a18_isolated_cap_domain_and_norm_screen_20261021.md"
RATIONAL_OBSTRUCTION = "manuscripts/compact_arithmetic_dimension_bound.tex"
RATIONAL_OBSTRUCTION_SHA = "1f09df2976085fd1ad1aa56ed35f1ab182ff95a8332dcfa0b0b80d4849b82203"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            h.update(block)
    return h.hexdigest()


def trim(p):
    p = list(map(Q, p))
    while p and not p[-1]:
        p.pop()
    return tuple(p)


def value(p, x):
    total = Q(0)
    for coefficient in reversed(p):
        total = total * x + coefficient
    return total


def remainder(a, b):
    assert b
    a = list(a)
    while len(a) >= len(b):
        lead = a[-1] / b[-1]
        offset = len(a) - len(b)
        for j, coefficient in enumerate(b):
            a[offset + j] -= lead * coefficient
        a = list(trim(a))
    return tuple(a)


def sturm(p):
    derivative = trim(j * p[j] for j in range(1, len(p)))
    sequence = [p, derivative]
    while sequence[-1]:
        r = trim(-x for x in remainder(sequence[-2], sequence[-1]))
        if not r:
            break
        sequence.append(r)
    assert len(sequence[-1]) == 1, "Irreducibility requires squarefreeness"
    return sequence


def roots(sequence, a, b):
    assert a < b and value(sequence[0], a) and value(sequence[0], b)

    def changes(at):
        signs = [(v > 0) - (v < 0) for p in sequence if (v := value(p, at))]
        return sum(x != y for x, y in zip(signs, signs[1:]))

    answer = changes(a) - changes(b)
    assert 0 <= answer <= len(sequence[0]) - 1
    return answer


def squarefree(n):
    assert n > 0
    original, part, square, prime = n, 1, 1, 2
    while prime * prime <= n:
        exponent = 0
        while n % prime == 0:
            n //= prime
            exponent += 1
        square *= prime ** (exponent // 2)
        if exponent % 2:
            part *= prime
        prime = 3 if prime == 2 else prime + 2
    if n > 1:
        part *= n
    assert original == part * square * square
    return part, square


def reduce_mod(a, p):
    n = len(p) - 1
    r = remainder(trim(a), p)
    return r + (Q(0),) * (n - len(r))


def determinant(matrix):
    a = [row[:] for row in matrix]
    det, n = Q(1), len(a)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return Q(0)
        if pivot != j:
            a[pivot], a[j] = a[j], a[pivot]
            det = -det
        det *= a[j][j]
        for i in range(j + 1, n):
            factor = a[i][j] / a[j][j]
            for k in range(j + 1, n):
                a[i][k] -= factor * a[j][k]
            a[i][j] = Q(0)
    return det


def multiplication_norm(p, k):
    n = len(p) - 1
    f = reduce_mod((0, Q(-38, 19), Q(k, 19)), p)
    columns = [reduce_mod((Q(0),) * j + f, p) for j in range(n)]
    return determinant([[columns[j][i] for j in range(n)] for i in range(n)])


def coefficients(degree):
    if degree == 2:
        for t in range(5, 14):
            for m in range(1, 8):
                yield {"T": t, "M": m}, trim((m, -t, 1))
    else:
        for t in range(5, 15):
            for s in range(1, 15):
                for m in range(1, 4):
                    yield {"T": t, "S": s, "M": m}, trim((-m, s, -t, 1))


def irreducible(p, c, degree):
    t, m = c["T"], c["M"]
    if degree == 2:
        discriminant = t * t - 4 * m
        return discriminant > 0 and isqrt(discriminant) ** 2 != discriminant
    for divisor in range(1, m + 1):
        if m % divisor == 0:
            if not value(p, Q(divisor)) or not value(p, Q(-divisor)):
                return False
    return True


def polynomial_discriminant(c, degree):
    t, m = c["T"], c["M"]
    if degree == 2:
        return t * t - 4 * m
    s = c["S"]
    return t * t * s * s - 4 * s**3 - 4 * t**3 * m - 27 * m * m + 18 * t * s * m


def norm_numerator(c, degree, k):
    t, m = c["T"], c["M"]
    if degree == 2:
        return m * (m * k * k - 38 * t * k + 1444)
    return m * (m * k**3 - 38 * c["S"] * k * k + 1444 * t * k - 54872)


def census(degrees):
    counts = Counter()
    polynomials, cells, grouped = [], [], defaultdict(list)
    for degree in degrees:
        for c, p in coefficients(degree):
            tag = str(degree)
            counts[tag + ":visited"] += 1
            if not irreducible(p, c, degree):
                counts[tag + ":reducible_or_nonreal_quadratic"] += 1
                continue
            sequence = sturm(p)
            if roots(sequence, Q(0), Q(19, 35)) != degree - 1 or roots(sequence, Q(4), Q(40, 3)) != 1:
                counts[tag + ":outside_root_domain"] += 1
                continue
            discriminant = polynomial_discriminant(c, degree)
            assert discriminant > 0
            disc_part, disc_square = squarefree(discriminant)
            polynomial_id = len(polynomials)
            polynomials.append({"id": polynomial_id, "degree": degree,
                "coefficients": c, "ascending": [int(x) for x in p],
                "discriminant": discriminant,
                "discriminant_squarefree_part": disc_part,
                "discriminant_square_factor": disc_square})
            counts[tag + ":retained_polynomials"] += 1
            for i in range(5, 15):
                k = i * (19 - i)
                counts[tag + ":index_tests"] += 1
                if roots(sequence, Q(0), Q(38, k)) != degree - 1:
                    counts[tag + ":failed_index_schur"] += 1
                    continue
                numerator = norm_numerator(c, degree, k)
                assert numerator and (numerator > 0) == (degree == 3)
                assert multiplication_norm(p, k) == Q(numerator, 19**degree)
                norm_part, norm_square = squarefree(abs(numerator))
                cell_id = len(cells)
                cells.append({"id": cell_id, "polynomial": polynomial_id,
                    "degree": degree, "index": i, "k": k,
                    "norm_numerator": numerator, "norm_denominator": 19**degree,
                    "norm_squarefree_part": norm_part,
                    "norm_square_factor": norm_square})
                grouped[(degree, disc_part, norm_part)].append(cell_id)
                counts[tag + ":retained_cells_and_norm_checks"] += 1
    groups = []
    for key, ids in sorted(grouped.items()):
        polynomial_ids = sorted({cells[i]["polynomial"] for i in ids})
        k_values = sorted({cells[i]["k"] for i in ids})
        groups.append({"degree": key[0], "discriminant_squarefree_part": key[1],
            "norm_squarefree_part": key[2], "cells": ids,
            "polynomials": polynomial_ids, "k_values": k_values,
            "one_polynomial_and_one_k": len(polynomial_ids) == len(k_values) == 1})
    for degree in degrees:
        prefix = str(degree) + ":"
        assert counts[prefix + "visited"] == (63 if degree == 2 else 420)
        assert counts[prefix + "visited"] == sum(counts[prefix + k] for k in (
            "reducible_or_nonreal_quadratic", "outside_root_domain", "retained_polynomials"))
        assert counts[prefix + "index_tests"] == 10 * counts[prefix + "retained_polynomials"]
        assert counts[prefix + "index_tests"] == counts[prefix + "failed_index_schur"] + counts[prefix + "retained_cells_and_norm_checks"]
    return {"counts": dict(sorted(counts.items())), "polynomials": polynomials,
        "cells": cells, "necessary_pair_groups": groups,
        "cubic_groups_one_polynomial_and_one_k": (
            all(g["one_polynomial_and_one_k"] for g in groups if g["degree"] == 3)
            if 3 in degrees else None)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", choices=("2", "3", "both"), default="both")
    parser.add_argument("--out", required=True, help="New compact JSON in results/")
    args = parser.parse_args()
    assert os.environ.get("COXETER_CHECK_OWNER"), "Use the bounded parent-lock runner"
    output = (ROOT / args.out).resolve()
    assert output.parent == ROOT / "results" and not output.exists()
    inputs = [ROOT / DOMAIN, ROOT / PROOF, ROOT / RATIONAL_OBSTRUCTION,
              Path(__file__).resolve()]
    pins = [{"path": str(path.relative_to(ROOT)), "sha256": digest(path),
             "bytes": path.stat().st_size} for path in inputs]
    assert pins[0]["sha256"] == DOMAIN_SHA
    assert pins[2]["sha256"] == RATIONAL_OBSTRUCTION_SHA
    degrees = (2, 3) if args.degree == "both" else (int(args.degree),)
    result = census(degrees)
    for path, pin in zip(inputs, pins):
        assert digest(path) == pin["sha256"] and path.stat().st_size == pin["bytes"]
    result.update({"status": "EXACT_NECESSARY_DOMAIN_AND_NORM_SCREEN_COMPLETE",
        "input_pins": pins, "degrees": list(degrees),
        "scope": "Necessary isolated dotted cap domain only; groups do not assert field isomorphism, a K-square, or geometric coexistence."})
    with output.open("x") as stream:
        json.dump(result, stream, separators=(",", ":"), sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "counts": result["counts"],
        "groups": len(result["necessary_pair_groups"]),
        "cubic_groups_one_polynomial_and_one_k": result["cubic_groups_one_polynomial_and_one_k"],
        "output": str(output.relative_to(ROOT)), "output_sha256": digest(output)}))


if __name__ == "__main__":
    main()
