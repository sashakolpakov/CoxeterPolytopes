#!/usr/bin/env python3
"""Scout Brandhorst--Hofmann rows for the two high-rank K3 folds.

This is a deliberately cheap exact-intersection pass.  It reconstructs the
generic Neron--Severi lattice of every effective action row, filters for the
two complete Kondo chambers of ranks 15 and 16, and then compares every
order-two/order-three Picard element with the archived high fixed actions.
An output hit is a candidate for the stronger integral action/glue audit; it
is not by itself promoted to a Torelli realization.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
import time

import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import _q_histogram
from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    induced_row_action,
)
from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    saturated_row_kernel,
    smith_diagonal,
)


DEFAULT_ROOT = Path("/private/tmp/K3Groups_extract/K3Groups/results")

NS_TARGETS = {
    (
        15,
        8,
        (1,) * 12 + (2, 2, 2),
        (("0", 1), ("1", 3), ("1/2", 3), ("3/2", 1)),
    ): "UE8D4A1",
    (
        16,
        -4,
        (1,) * 14 + (2, 2),
        (("0", 1), ("1", 1), ("1/2", 2)),
    ): "UE8D6",
}

FIXED_TARGETS = {
    (
        "UE8D4A1",
        2,
        14,
        -8,
        (1,) * 12 + (2, 4),
        (("0", 1), ("1", 1), ("1/2", 1), ("3/2", 1), ("3/4", 2), ("5/4", 2)),
    ): "R15-UE8D4A1-C2-H13-19",
    (
        "UE8D4A1",
        3,
        13,
        6,
        (1,) * 12 + (6,),
        (("0", 1), ("3/2", 1), ("4/3", 2), ("5/6", 2)),
    ): "R15-UE8D4A1-C3-H12-17",
    (
        "UE8D6",
        2,
        15,
        4,
        (1,) * 14 + (4,),
        (("0", 1), ("1", 1), ("3/4", 2)),
    ): "R16-UE8D6-C2-H14-17",
}


def parse_row(line: str) -> dict:
    exact = re.sub(r"(?<![\w.])(\d+)/(\d+)", r"Fraction(\1,\2)", line)
    return dict(eval(exact, {"__builtins__": {}, "Fraction": Fraction}))


def rational(value: object) -> sp.Rational:
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    return sp.Rational(value)


def matrix(values: list) -> sp.Matrix:
    return sp.Matrix(22, 22, list(map(rational, values)))


def cyclotomic_evaluation(action: sp.Matrix, order: int) -> sp.Matrix:
    variable = sp.Symbol("x")
    polynomial = sp.Poly(sp.cyclotomic_poly(order, variable), variable)
    result = sp.zeros(22)
    for (exponent,), coefficient in polynomial.terms():
        result += coefficient * action**exponent
    return result


def matrix_order(action: sp.Matrix, bound: int = 256) -> int | None:
    product = sp.eye(action.rows)
    for order in range(1, bound + 1):
        product *= action
        if product == sp.eye(action.rows):
            return order
    return None


def q_key(gram: sp.Matrix) -> tuple[tuple[str, int], ...]:
    return tuple(sorted(_q_histogram(gram).items()))


def generic_ns_rank_fast(row: dict, order: int) -> int:
    symplectic = [matrix(action) for action in row["symplectic"]]
    generator = matrix(row["distinguished_generator"])
    equations = [action - sp.eye(22) for action in symplectic]
    equations.append(cyclotomic_evaluation(generator, order))
    return sp.Matrix.hstack(*equations).rank()


def reconstruct_ns(row: dict, order: int) -> tuple[sp.Matrix, list[sp.Matrix]]:
    basis = matrix(row["basis"])
    ambient_gram = matrix(row["inner_product_matrix"])
    inverse_basis = basis.inv()
    gram = sp.simplify(basis * ambient_gram * basis.T)
    symplectic = [
        sp.simplify(basis * matrix(action) * inverse_basis)
        for action in row["symplectic"]
    ]
    generator = sp.simplify(
        basis * matrix(row["distinguished_generator"]) * inverse_basis
    )
    equations = [action - sp.eye(22) for action in symplectic]
    equations.append(cyclotomic_evaluation(generator, order))
    transcendental = saturated_row_kernel(sp.Matrix.hstack(*equations))
    ns_basis = saturated_row_kernel(gram * transcendental.T)
    ns_gram = sp.simplify(ns_basis * gram * ns_basis.T)
    actions = [induced_row_action(ns_basis, action) for action in symplectic]
    actions.append(induced_row_action(ns_basis, generator))
    return ns_gram, actions


def enumerate_group(generators: list[sp.Matrix], maximum: int = 10000) -> list[sp.Matrix]:
    identity = sp.eye(generators[0].rows)
    seen = {tuple(identity): identity}
    queue = [identity]
    while queue:
        left = queue.pop()
        for right in generators:
            product = left * right
            key = tuple(product)
            if key in seen:
                continue
            seen[key] = product
            queue.append(product)
            if len(seen) > maximum:
                raise RuntimeError(f"Picard image exceeds scout bound {maximum}")
    return list(seen.values())


def scan(paths: list[Path]) -> dict:
    started = time.time()
    total = rank_hits = ns_hits = 0
    found: list[dict] = []
    for path in paths:
        transcendental_order = int(path.stem[-2:])
        for row_number, line in enumerate(path.open(), 1):
            total += 1
            row = parse_row(line)
            if generic_ns_rank_fast(row, transcendental_order) not in (15, 16):
                continue
            ns_gram, generators = reconstruct_ns(row, transcendental_order)
            rank_hits += 1
            ns_key = (
                ns_gram.rows,
                int(ns_gram.det()),
                tuple(smith_diagonal(ns_gram)),
                q_key(ns_gram),
            )
            if ns_key not in NS_TARGETS:
                continue
            ns_hits += 1
            ns_type = NS_TARGETS[ns_key]
            group = enumerate_group(generators)
            for action in group:
                order = matrix_order(action)
                if order not in (2, 3):
                    continue
                fixed_basis = saturated_row_kernel(action - sp.eye(action.rows))
                fixed_gram = sp.simplify(fixed_basis * ns_gram * fixed_basis.T)
                fixed_key = (
                    ns_type,
                    order,
                    fixed_gram.rows,
                    int(fixed_gram.det()),
                    tuple(smith_diagonal(fixed_gram)),
                    q_key(fixed_gram),
                )
                if fixed_key not in FIXED_TARGETS:
                    continue
                hit = {
                    "member": str(path),
                    "row": row_number,
                    "id": row["id"],
                    "transcendental_order": transcendental_order,
                    "ns_type": ns_type,
                    "target": FIXED_TARGETS[fixed_key],
                    "symplectic_generator_count": len(generators) - 1,
                    "picard_image_order": len(group),
                    "element_order": order,
                    "action": [[str(value) for value in action.row(i)] for i in range(action.rows)],
                }
                found.append(hit)
                print("POSITIVE-CANDIDATE", json.dumps(hit, sort_keys=True), flush=True)
        print(
            "DONE",
            path.name,
            "rows",
            total,
            "rank_hits",
            rank_hits,
            "ns_hits",
            ns_hits,
            "found",
            len(found),
            "seconds",
            round(time.time() - started, 1),
            flush=True,
        )
    return {
        "row_count": total,
        "rank_hit_count": rank_hits,
        "ns_hit_count": ns_hits,
        "found": found,
        "elapsed_seconds": round(time.time() - started, 3),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--pure", action="store_true")
    parser.add_argument("--orders", nargs="*", type=int)
    args = parser.parse_args()
    lane = args.root / "purely_ns" if args.pure else args.root
    paths = sorted(lane.glob("order*.txt"))
    if args.orders:
        paths = [path for path in paths if int(path.stem[-2:]) in args.orders]
    print(json.dumps(scan(paths), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
