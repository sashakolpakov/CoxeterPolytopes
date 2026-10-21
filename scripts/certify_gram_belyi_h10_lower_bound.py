#!/usr/bin/env python3
"""Exact calibration for the Gram--Belyi Coxeter-dessin complexity.

The theorem report defines a Belyi packet from the scale-invariant squared
Gram entries of an arithmetic Coxeter chamber.  This lightweight producer
checks the nontrivial high-dimensional calibration: the frozen H10 K3
redoubling family contains one packet coordinate whose normalized Belyi
height is at least p on the member N=2(p+1), for every odd prime p.

The only external mathematical input is Rodriguez's theorem that a nonzero
algebraic number with nonzero p-adic valuation has normalized Belyi-polynomial
height at least p.  This script verifies every lattice and valuation
hypothesis exactly; it does not pretend to reprove that theorem.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.json"
FAMILY = ROOT / "data/enumerations/k3_s3_h10_redoubling_family.json"
FAMILY_SCRIPT = ROOT / "scripts/certify_k3_s3_h10_redoubling_family.py"
OUTPUT = ROOT / "data/enumerations/gram_belyi_h10_lower_bound.json"

SCHEMA = "coxeter.gram_belyi.h10_lower_bound.v1"
SOURCE_SHA256 = "93cec35a3c395482f4f99be827b3f5a6ba44764fc7defa22bf23eb6ee951acf0"
FAMILY_SHA256 = "c2cf3c5ee99765ecccc88c027f47666a059eafa521ae45f06a7df46bd47a4df0"
FAMILY_SCRIPT_SHA256 = "a215f2b31a71c1d8617d2cfce699a2fcc2e96d4b2445aadeb1be4049f51293b8"

DOUBLING = (9, 12)  # zero based: a,b
MOVING = (8, 10)  # zero based: c,d
PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pair(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> sp.Integer:
    value = sp.expand((left.T * gram * right)[0])
    if not value.is_Integer:
        raise ArithmeticError(f"nonintegral pairing {value}")
    return sp.Integer(value)


def valuation(integer: int, prime: int) -> int:
    if integer == 0:
        raise ValueError("valuation of zero is not used")
    value = abs(integer)
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def build_certificate() -> dict[str, Any]:
    pinned = {
        "carrier": sha256(SOURCE),
        "redoubling_family": sha256(FAMILY),
        "redoubling_producer": sha256(FAMILY_SCRIPT),
    }
    expected = {
        "carrier": SOURCE_SHA256,
        "redoubling_family": FAMILY_SHA256,
        "redoubling_producer": FAMILY_SCRIPT_SHA256,
    }
    if pinned != expected:
        raise ArithmeticError(f"a pinned H10 source changed: {pinned}")

    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    if not source["all_checks_pass"] or not family["all_checks_pass"]:
        raise ArithmeticError("the frozen H10 family is not green")

    gram = sp.Matrix(source["full_Picard_lattice"]["standard_geometric_gram"])
    roots = tuple(
        sp.Matrix(row["primitive_geometric_standard_F_coordinates"])
        for row in source["transported_H10_13_chamber"]["wall_records"]
    )
    a, b = (roots[index] for index in DOUBLING)
    c, d = (roots[index] for index in MOVING)

    base_pairings = {
        "a2": int(pair(gram, a, a)),
        "b2": int(pair(gram, b, b)),
        "c2": int(pair(gram, c, c)),
        "d2": int(pair(gram, d, d)),
        "ab": int(pair(gram, a, b)),
        "ca": int(pair(gram, c, a)),
        "cb": int(pair(gram, c, b)),
        "cd": int(pair(gram, c, d)),
    }
    if base_pairings != {
        "a2": -2,
        "b2": -2,
        "c2": -4,
        "d2": -2,
        "ab": 2,
        "ca": 2,
        "cb": 0,
        "cd": 0,
    }:
        raise ArithmeticError(f"unexpected redoubling seed: {base_pairings}")

    m = sp.symbols("m", integer=True, nonnegative=True)
    q_even = c + 2 * m**2 * a + 2 * m * (m - 1) * b
    q_square = sp.factor((q_even.T * gram * q_even)[0])
    cq_pair = sp.factor((c.T * gram * q_even)[0])
    normalized_square = sp.factor(cq_pair**2 / (pair(gram, c, c) * q_square))
    expected_symbolics = (-4, 4 * (m**2 - 1), (m**2 - 1) ** 2)
    if any(
        sp.simplify(actual - expected) != 0
        for actual, expected in zip(
            (q_square, cq_pair, normalized_square), expected_symbolics
        )
    ):
        raise ArithmeticError("symbolic moving-wall identity failed")

    rows = []
    for p in PRIMES:
        if not sp.isprime(p) or p == 2:
            raise ArithmeticError(f"bad prime control {p}")
        mm = p + 1
        copies = 2 * mm
        u = (mm * mm - 1) ** 2
        lam = Fraction(u, u + 1)
        vp_num = valuation(lam.numerator, p)
        vp_den = valuation(lam.denominator, p)
        if not (
            copies == 2 * (p + 1)
            and u == (p * (p + 2)) ** 2
            and vp_num == 2
            and vp_den == 0
            and 0 < lam < 1
        ):
            raise ArithmeticError(f"p-adic control failed at p={p}")
        rows.append(
            {
                "prime": p,
                "path_copies_N": copies,
                "facet_count": copies + 12,
                "moving_index_m": mm,
                "normalized_squared_gram_u": u,
                "belyi_coordinate_lambda": f"{lam.numerator}/{lam.denominator}",
                "p_adic_valuation_lambda": vp_num - vp_den,
                "rodriguez_belyi_height_lower_bound": p,
            }
        )

    return {
        "schema": SCHEMA,
        "all_checks_pass": True,
        "pinned_sha256": pinned,
        "seed_pairings": base_pairings,
        "symbolic_identity": {
            "q_2m": "c + 2*m^2*a + 2*m*(m-1)*b",
            "q_2m_square": str(q_square),
            "pair_c_q_2m": str(cq_pair),
            "normalized_squared_gram_u": str(normalized_square),
            "packet_coordinate": "lambda=u/(u+1)",
        },
        "prime_rows": rows,
        "proved_scope": {
            "family": "K3-S3 H10 redoubling P_N",
            "subsequence": "phase A, N=2*(p+1), p odd prime",
            "conclusion": "one normalized Gram-Belyi packet coordinate has degree at least p",
            "resource_claim": "exact arithmetic only; no cone or wall enumeration",
        },
        "external_theorem": {
            "author": "Jose Rodriguez",
            "title": "Bounding the Degree of Belyi Polynomials",
            "arxiv": "1104.2027",
            "statement_used": "nonzero p-adic valuation implies normalized Belyi-polynomial height at least p",
        },
    }


def canonical(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = canonical(payload)
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("frozen Gram-Belyi archive differs from exact replay")
    print(rendered, end="")


if __name__ == "__main__":
    main()
