#!/usr/bin/env python3
"""Exact interface checks for bounded-frame native-host closure."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/ag_field_funnel_bounded_frame_native_host_closure.md"


Matrix = tuple[tuple[Fraction, ...], ...]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_product(left: Matrix, right: Matrix) -> Matrix:
    size = len(left)
    return tuple(
        tuple(
            sum((left[i][k] * right[k][j] for k in range(size)), Fraction(0))
            for j in range(size)
        )
        for i in range(size)
    )


def trace(matrix: Matrix) -> Fraction:
    return sum((matrix[i][i] for i in range(len(matrix))), Fraction(0))


def rank_one_cartan_operator(cartan: Matrix, index: int) -> Matrix:
    """Matrix of T_i=e_i tensor alpha_i in the root basis."""
    size = len(cartan)
    return tuple(
        tuple(cartan[index][j] if i == index else Fraction(0) for j in range(size))
        for i in range(size)
    )


def directed_cycle_product(cartan: Matrix, cycle: Iterable[int]) -> Fraction:
    vertices = tuple(cycle)
    value = Fraction(1)
    for position in range(len(vertices)):
        value *= cartan[vertices[position]][
            vertices[(position + 1) % len(vertices)]
        ]
    return value


def trace_cycle_product(cartan: Matrix, cycle: Iterable[int]) -> Fraction:
    vertices = tuple(cycle)
    operators = [rank_one_cartan_operator(cartan, index) for index in vertices]
    product = operators[0]
    for operator in operators[1:]:
        product = matrix_product(product, operator)
    return trace(product)


def diagonal_rescale(cartan: Matrix, scales: tuple[Fraction, ...]) -> Matrix:
    return tuple(
        tuple(scales[j] * cartan[i][j] / scales[i] for j in range(len(cartan)))
        for i in range(len(cartan))
    )


def recover_tree_scales(
    original: Matrix,
    rescaled: Matrix,
    tree: tuple[tuple[int, int], ...],
    root: int = 0,
) -> tuple[Fraction, ...]:
    size = len(original)
    scales: list[Fraction | None] = [None] * size
    scales[root] = Fraction(1)
    pending = list(tree)
    while pending:
        progressed = False
        next_pending: list[tuple[int, int]] = []
        for i, j in pending:
            if scales[i] is not None and scales[j] is None:
                scales[j] = scales[i] * rescaled[i][j] / original[i][j]
                progressed = True
            elif scales[j] is not None and scales[i] is None:
                scales[i] = scales[j] * rescaled[j][i] / original[j][i]
                progressed = True
            else:
                next_pending.append((i, j))
        if not progressed and next_pending:
            raise ValueError("tree is not connected")
        pending = next_pending
    return tuple(value for value in scales if value is not None)


def matrix_unit_hull_diagonal(diagonal_generators: Iterable[int], rank: int) -> int:
    """Return d for M_rank(Z) L=d Z^rank in a diagonal sample."""
    values = tuple(abs(value) for value in diagonal_generators)
    divisor = 0
    for value in values:
        divisor = gcd(divisor, value)
    if len(values) != rank or divisor == 0:
        raise ValueError("expected a full diagonal lattice")
    return divisor


def report_interface() -> dict[str, object]:
    text = REPORT.read_text(encoding="utf-8")
    required = [
        "simple directed cycle",
        "length-two backtrack",
        "composition is read from right to left",
        "M=BL",
        "matrix units",
        "No native host escapes",
        "AG-DEG",
        "not AG-only",
        "n\\leq995",
        "Nikulin fixed-data finiteness conclusion | none | **not used**",
    ]
    missing = [phrase for phrase in required if phrase not in text]
    return {
        "missing_required_phrases": missing,
        "inline_math_delimiters_balanced": text.count("\\(") == text.count("\\)"),
        "display_math_delimiters_balanced": text.count("\\[") == text.count("\\]"),
        "forbidden_transport_fragments_absent": all(
            fragment not in text for fragment in ["Let (k) be", "(A=62)", "(n\\leq995)"]
        ),
    }


def build_certificate() -> dict[str, object]:
    cartan: Matrix = (
        (Fraction(2), Fraction(-1), Fraction(-2), Fraction(0)),
        (Fraction(-3), Fraction(2), Fraction(-4), Fraction(-1)),
        (Fraction(-5), Fraction(-6), Fraction(2), Fraction(-2)),
        (Fraction(0), Fraction(-7), Fraction(-8), Fraction(2)),
    )
    scales = (Fraction(2), Fraction(3), Fraction(5), Fraction(7))
    changed = diagonal_rescale(cartan, scales)
    cycles = ((0, 1), (1, 2), (2, 3), (0, 1, 2), (1, 2, 3))
    trace_checks = {
        "-".join(map(str, cycle)): {
            "directed_product": str(directed_cycle_product(cartan, cycle)),
            "operator_trace": str(trace_cycle_product(cartan, cycle)),
            "equal": directed_cycle_product(cartan, cycle)
            == trace_cycle_product(cartan, cycle),
        }
        for cycle in cycles
    }
    invariance = {
        "-".join(map(str, cycle)): directed_cycle_product(cartan, cycle)
        == directed_cycle_product(changed, cycle)
        for cycle in cycles
    }
    recovered = recover_tree_scales(cartan, changed, ((0, 1), (1, 2), (2, 3)))
    normalized_expected = tuple(value / scales[0] for value in scales)
    hull_divisor = matrix_unit_hull_diagonal((6, 10, 14), 3)
    interface = report_interface()
    all_checks_pass = all(
        [
            all(row["equal"] for row in trace_checks.values()),
            all(invariance.values()),
            recovered == normalized_expected,
            diagonal_rescale(cartan, recovered) == changed,
            hull_divisor == 2,
            995 - 2 + 1 == 994,
            not interface["missing_required_phrases"],
            interface["inline_math_delimiters_balanced"],
            interface["display_math_delimiters_balanced"],
            interface["forbidden_transport_fragments_absent"],
        ]
    )
    return {
        "schema": "ag-field-funnel-bounded-frame-native-host-closure-v1",
        "all_checks_pass": all_checks_pass,
        "report_sha256": sha256_file(REPORT),
        "cycle_trace_orientation": trace_checks,
        "diagonal_cycle_invariance": invariance,
        "recovered_tree_scales": [str(value) for value in recovered],
        "expected_tree_scales_up_to_common_scalar": [
            str(value) for value in normalized_expected
        ],
        "matrix_unit_hull_sample": {
            "input_diagonal_lattice": [6, 10, 14],
            "rank": 3,
            "full_matrix_order_hull": "2*Z^3",
            "computed_divisor": hull_divisor,
        },
        "dimension_quantifier": {
            "hyperbolic_dimensions": [2, 995],
            "native_ranks": [3, 996],
            "number_of_ranks": 994,
        },
        "report_interface": interface,
        "scope": {
            "direct_AG_DEG_proved": False,
            "fixed_data_arithmetic_core_reconstructed": True,
            "nikulin_fixed_data_conclusion_imported": False,
            "bounded_frame_imported": True,
            "geometric_dimension_ceiling_imported": True,
            "high_dimensional_AG_field_funnel_essential": True,
        },
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))
    if not certificate["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
