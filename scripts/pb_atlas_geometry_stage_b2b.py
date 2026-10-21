#!/usr/bin/env python3
"""PB-Atlas Stage B2b: exact closed-chart section geometry.

This module is deliberately independent of any named Coble surface.  It
consumes one certified Stage-B1 P2 packet, one Stage-B2a chart choice, and an
exact rational closed point of that chart.  It compiles a divisor evaluation
matrix, extracts its unique section when the kernel is one-dimensional, and
replays the strict transform through the complete ordered blowup cluster.

The output certifies exact (not merely lower-bounded) multiplicity at every
centre.  Absolute irreducibility and component-level nonemptiness are separate
finite algebra gates handled by the calling certifier.  No nef, hidden-wall,
same-P, or rank-bound conclusion is made here.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import sympy as sp

from scripts.pb_atlas_cluster_chart_stage_b2a import (
    _build_chart,
    _split_low_high,
    compile_p2_cluster_chart,
    compile_p2_divisor_evaluation_matrix,
)


SCHEMA = "coxeter.pb_atlas.closed_section_geometry_stage_b2b.v1"


def _matrix_from_record(record: Mapping[str, Any]) -> sp.Matrix:
    return sp.Matrix(
        [
            [
                sp.cancel(
                    sp.sympify(entry["numerator"])
                    / sp.sympify(entry["denominator"])
                )
                for entry in row
            ]
            for row in record["entries"]
        ]
    )


def _exact_rational(value: Any, name: str) -> sp.Rational:
    if isinstance(value, bool) or isinstance(value, (float, sp.Float)):
        raise ValueError(f"{name} must be an exact rational")
    try:
        return sp.Rational(value)
    except (TypeError, ValueError) as exception:
        raise ValueError(f"{name} must be an exact rational") from exception


def _closed_substitution(
    chart: Mapping[str, Any], values: Mapping[str | sp.Symbol, Any]
) -> dict[sp.Symbol, sp.Rational]:
    if not isinstance(values, Mapping):
        raise ValueError("closed-chart coordinates must be a mapping")
    parsed: dict[sp.Symbol, sp.Rational] = {}
    for key, value in values.items():
        if not isinstance(key, (str, sp.Symbol)):
            raise ValueError("a closed-chart coordinate key must be a symbol name")
        symbol = sp.Symbol(str(key))
        if symbol in parsed:
            raise ValueError("a closed-chart coordinate occurs more than once")
        parsed[symbol] = _exact_rational(value, f"closed coordinate {symbol}")
    expected = {sp.Symbol(name) for name in chart["coordinate_variables"]}
    if set(parsed) != expected:
        missing = sorted(str(value) for value in expected - set(parsed))
        extra = sorted(str(value) for value in set(parsed) - expected)
        raise ValueError(
            f"closed-chart coordinates do not match the chart; missing={missing}, extra={extra}"
        )
    return parsed


def _specialize_matrix(
    record: Mapping[str, Any], substitution: Mapping[sp.Symbol, sp.Rational]
) -> sp.Matrix:
    matrix = _matrix_from_record(record).applyfunc(
        lambda value: sp.cancel(value.subs(substitution, simultaneous=True))
    )
    if any(value.free_symbols for value in matrix):
        raise ValueError("a closed-chart matrix still contains a free parameter")
    if any(value.has(sp.zoo, sp.nan, sp.oo, -sp.oo) for value in matrix):
        raise ValueError("a closed-chart specialization meets a denominator pole")
    return matrix


def _primitive_kernel(vector: sp.Matrix) -> list[int]:
    denominator = 1
    for value in vector:
        if not value.is_Rational:
            raise ValueError("the closed-chart kernel is not rational")
        denominator = int(sp.ilcm(denominator, int(sp.denom(value))))
    entries = [int(value * denominator) for value in vector]
    divisor = 0
    for value in entries:
        divisor = math.gcd(divisor, abs(value))
    if divisor == 0:
        raise ValueError("the closed-chart kernel vector is zero")
    entries = [value // divisor for value in entries]
    for value in entries:
        if value:
            if value < 0:
                entries = [-entry for entry in entries]
            break
    return entries


def _poly_string(value: sp.Expr) -> str:
    return str(sp.expand(value)).replace("**", "^")


def _exact_cluster_orders(
    *,
    packet: Mapping[str, Any],
    chart_options: Mapping[str, Any] | None,
    substitution: Mapping[sp.Symbol, sp.Rational],
    polynomial: sp.Expr,
    multiplicities: Sequence[int],
) -> list[int]:
    _public, internal = _build_chart(packet, chart_options)
    centers = internal["centers"]
    if len(multiplicities) != len(centers):
        raise ValueError("the multiplicity vector has the wrong length")
    X, Y, Z = sp.symbols("X Y Z")
    current: dict[int, sp.Expr] = {}
    for position, center in enumerate(centers):
        if center["center_kind"] == "proper":
            current[position] = sp.expand(
                polynomial.subs(internal["proper_substitutions"][position])
            ).subs(substitution, simultaneous=True)

    children_by_primary: dict[int, list[int]] = {}
    for child, parent in enumerate(internal["primary_parent"]):
        if parent is not None:
            children_by_primary.setdefault(int(parent), []).append(child)

    orders: list[int] = []
    for position, _center in enumerate(centers):
        if position not in current:
            raise ValueError("the chart did not transport the section to a centre")
        u, v = internal["local_variables"][position]
        expression = sp.cancel(current[position]).subs(substitution, simultaneous=True)
        numerator, denominator = sp.fraction(expression)
        if denominator.has(u) or denominator.has(v):
            raise ValueError("the specialized section has a local Laurent pole")
        parsed = sp.Poly(sp.expand(numerator), u, v)
        nonzero_terms = [
            (exponents, coefficient)
            for exponents, coefficient in parsed.terms()
            if coefficient != 0
        ]
        if not nonzero_terms:
            raise ValueError("the section vanishes identically in a centre chart")
        order = min(sum(exponents) for exponents, _coefficient in nonzero_terms)
        orders.append(order)

        multiplicity = multiplicities[position]
        if isinstance(multiplicity, bool) or not isinstance(multiplicity, int) or multiplicity < 0:
            raise ValueError("section multiplicities must be nonnegative integers")
        low, high = _split_low_high(expression, u, v, multiplicity)
        if any(sp.cancel(value.subs(substitution, simultaneous=True)) != 0 for value in low.values()):
            raise ArithmeticError("the kernel section does not satisfy its declared jets")

        for child in children_by_primary.get(position, []):
            uc, vc = internal["local_variables"][child]
            change = internal["parent_to_child_substitution"][child]
            divisor = uc if internal["chart_mode"][child] == "A" else vc
            transformed = sp.cancel(
                high.subs(change, simultaneous=True) / divisor**multiplicity
            ).subs(substitution, simultaneous=True)
            child_numerator, child_denominator = sp.fraction(transformed)
            if child_denominator.has(uc) or child_denominator.has(vc):
                raise ArithmeticError("strict transform retained a local Laurent denominator")
            current[child] = sp.expand(child_numerator / child_denominator)
    return orders


def compile_closed_section_geometry(
    packet: Mapping[str, Any],
    *,
    degree: int,
    multiplicities: Sequence[int],
    chart_options: Mapping[str, Any] | None,
    closed_coordinates: Mapping[str | sp.Symbol, Any],
) -> dict[str, Any]:
    """Compile one unique section and certify every cluster multiplicity."""

    chart = compile_p2_cluster_chart(packet, chart_options)["chart"]
    substitution = _closed_substitution(chart, closed_coordinates)
    record = compile_p2_divisor_evaluation_matrix(
        packet,
        degree=degree,
        multiplicities=multiplicities,
        chart_options=chart_options,
    )
    matrix = _specialize_matrix(record, substitution)
    kernel = matrix.nullspace()
    if len(kernel) != 1:
        raise ValueError(
            f"the specialized section space is not one-dimensional (nullity={len(kernel)})"
        )
    coefficients = _primitive_kernel(kernel[0])
    X, Y, Z = sp.symbols("X Y Z")
    columns = record["section_column_indices"]
    polynomial = sp.expand(
        sum(
            coefficient * X ** column[0] * Y ** column[1] * Z ** column[2]
            for coefficient, column in zip(coefficients, columns)
        )
    )
    orders = _exact_cluster_orders(
        packet=packet,
        chart_options=chart_options,
        substitution=substitution,
        polynomial=polynomial,
        multiplicities=list(multiplicities),
    )
    return {
        "schema": SCHEMA,
        "degree": degree,
        "matrix_dimensions_rows_columns": [matrix.rows, matrix.cols],
        "matrix_rank": int(matrix.rank()),
        "kernel_dimension": len(kernel),
        "kernel_coefficient_vector": coefficients,
        "homogeneous_polynomial": _poly_string(polynomial),
        "declared_multiplicities_chronological": list(multiplicities),
        "actual_multiplicities_chronological": orders,
        "every_declared_multiplicity_is_exact": orders == list(multiplicities),
        "closed_coordinate_values": {
            str(symbol): str(value) for symbol, value in sorted(substitution.items(), key=lambda item: str(item[0]))
        },
        "absolute_irreducibility_not_asserted_by_this_core": True,
        "all_checks_pass": len(kernel) == 1 and orders == list(multiplicities),
    }


def closed_linear_system_rank(
    packet: Mapping[str, Any],
    *,
    system: str,
    chart_options: Mapping[str, Any] | None,
    closed_coordinates: Mapping[str | sp.Symbol, Any],
) -> dict[str, int]:
    if system not in ("minus_K", "minus_two_K"):
        raise ValueError("system must be minus_K or minus_two_K")
    compiled = compile_p2_cluster_chart(packet, chart_options)
    substitution = _closed_substitution(compiled["chart"], closed_coordinates)
    matrix = _specialize_matrix(compiled["evaluation_matrices"][system], substitution)
    rank = int(matrix.rank())
    return {
        "rows": matrix.rows,
        "columns": matrix.cols,
        "rank": rank,
        "kernel_dimension": matrix.cols - rank,
    }

