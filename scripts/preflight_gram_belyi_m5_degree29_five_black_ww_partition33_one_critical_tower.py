#!/usr/bin/env python3
"""Bounded preflight for the one-critical-endpoint partition-(3^2,1^23) towers.

The endpoint-status convention is ordered as (x=0,x=1):

* CS: K=x(x-r), with 0 and r critical and 1 simple;
* SC: K=(x-1)(x-r), with 1 and r critical and 0 simple.

The four equations D=29*K^2 define a one-dimensional base.  We adjoin the
critical-value equation first, then the remaining endpoint-value equation.
This preflight reports exact Singular dimensions/bases but promotes no
timeout or non-unit computation as a theorem.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_five_black_ww_partition5 import (
    residual_derivative,
)


def equation_strings(weights, derivative, orientation: str) -> list[str]:
    a, b, c, d, e = weights
    g0 = f"(-L)^{a}*(-u)^{b}*(-v)^{c}*(-w)^{d}*(-z)^{e}"
    g1 = f"(1-L)^{a}*(1-u)^{b}*(1-v)^{c}*(1-w)^{d}*(1-z)^{e}"
    gr = f"(r-L)^{a}*(r-u)^{b}*(r-v)^{c}*(r-w)^{d}*(r-z)^{e}"

    if orientation == "CS":
        # 29*x^2*(x-r)^2.
        derivative_equations = [
            str(derivative[0]),
            str(derivative[1]),
            f"({derivative[2]})-29*r^2",
            f"({derivative[3]})+58*r",
        ]
        value_equations = [f"({gr})-({g0})", f"({g1})-({g0})"]
    elif orientation == "SC":
        # 29*(x-1)^2*(x-r)^2.
        derivative_equations = [
            f"({derivative[0]})-29*r^2",
            f"({derivative[1]})+58*r*(1+r)",
            f"({derivative[2]})-29*(1+4*r+r^2)",
            f"({derivative[3]})+58*(1+r)",
        ]
        value_equations = [f"({gr})-({g1})", f"({g0})-({g1})"]
    else:
        raise ValueError(orientation)
    return derivative_equations + value_equations


def singular_tower(prime: int, lam: int, equations, algorithm: str):
    program = "\n".join(
        [
            f"ring R={prime},(r,u,v,w,z),dp;",
            "option(redSB);",
            f"number L={lam};",
            *(f"poly f{i}={value};" for i, value in enumerate(equations)),
            "ideal I=f0,f1,f2,f3;",
            f"ideal B={algorithm}(I);",
            'print("DERIVATIVE_BASE"); print(dim(B));',
            f"ideal C={algorithm}(B+ideal(f4));",
            'print("CRITICAL_VALUE_BASE"); print(dim(C)); print(vdim(C));',
            f"ideal J={algorithm}(C+ideal(f5));",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(dim(J)); print(vdim(J)); }',
            "quit;",
        ]
    )
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    return completed.stdout.strip(), completed.stderr.strip(), time.monotonic() - before


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("orientation", choices=("CS", "SC"))
    parser.add_argument("weights", nargs=5, type=int)
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="std")
    arguments = parser.parse_args()
    weights = tuple(arguments.weights)
    if min(weights) <= 0 or sum(weights) != 29:
        raise SystemExit("positive weights must sum to 29")

    sqrt5, lam = lambda_mod(arguments.prime)
    _context, derivative = residual_derivative(arguments.prime, lam, weights)
    equations = equation_strings(weights, derivative, arguments.orientation)
    stdout, stderr, elapsed = singular_tower(
        arguments.prime, lam, equations, arguments.algorithm
    )
    print(
        json.dumps(
            {
                "orientation": arguments.orientation,
                "weights": weights,
                "prime": arguments.prime,
                "sqrt5": sqrt5,
                "lambda": lam,
                "structural_equation_degrees": [4, 3, 2, 1, 29, 29],
                "stdout": stdout,
                "stderr": stderr,
                "elapsed_seconds": elapsed,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
