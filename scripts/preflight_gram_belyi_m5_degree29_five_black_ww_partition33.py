#!/usr/bin/env python3
"""Exact representative preflight for five-black WW partition (3^2,1^23)."""

from __future__ import annotations

import argparse
import json
import subprocess
import time

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_five_black_ww_partition5 import (
    residual_derivative,
)


def equations(lam: int, weights, derivative) -> list[str]:
    a, b, c_weight, d, e = weights
    desired = (
        "29*t^2*r^2",
        "-58*t*r*(t+r)",
        "29*(t^2+4*t*r+r^2)",
        "-58*(t+r)",
    )
    critical = [f"({derivative[index]})-({desired[index]})" for index in range(4)]
    g0 = f"(-L)^{a}*(-u)^{b}*(-v)^{c_weight}*(-w)^{d}*(-z)^{e}"
    g1 = (
        f"(1-L)^{a}*(1-u)^{b}*(1-v)^{c_weight}*"
        f"(1-w)^{d}*(1-z)^{e}"
    )
    gt = (
        f"(t-L)^{a}*(t-u)^{b}*(t-v)^{c_weight}*"
        f"(t-w)^{d}*(t-z)^{e}"
    )
    gr = (
        f"(r-L)^{a}*(r-u)^{b}*(r-v)^{c_weight}*"
        f"(r-w)^{d}*(r-z)^{e}"
    )
    return critical + [f"({gr})-({gt})", f"({g0})-({gt})", f"({g1})-({gt})"]


def singular_run(prime: int, lam: int, values: list[str], projective: bool):
    variables = "(u,v,w,z,t,r,h)" if projective else "(u,v,w,z,t,r)"
    commands = [
        f"ring R={prime},{variables},dp;",
        "option(redSB);",
        f"number L={lam};",
        *(f"poly f{i}={value};" for i, value in enumerate(values)),
        "ideal I=f0,f1,f2,f3,f4,f5,f6;",
        "ideal J=std(I);",
        'if (size(J)==1 && string(J[1])=="1") { print("AFFINE_UNIT"); } else { print("AFFINE_NONUNIT"); print(dim(J)); print(vdim(J)); }',
    ]
    if projective:
        commands += [*(f"poly H{i}=homog(f{i},h);" for i in range(7))]
        for chart in ("u", "v", "w", "z", "t", "r"):
            chart_equations = ",".join(
                f"subst(subst(H{i},h,0),{chart},1)" for i in range(7)
            )
            commands += [
                f"ideal K={chart_equations}; ideal M=std(K);",
                f'if (size(M)==1 && string(M[1])=="1") {{ print("INF_{chart}_UNIT"); }} else {{ print("INF_{chart}_NONUNIT"); print(dim(M)); }}',
            ]
    commands.append("quit;")
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"],
        input="\n".join(commands),
        text=True,
        capture_output=True,
        check=True,
    )
    return completed.stdout.strip(), completed.stderr.strip(), time.monotonic() - before


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("weights", nargs=5, type=int)
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument("--projective", action="store_true")
    arguments = parser.parse_args()
    weights = tuple(arguments.weights)
    if min(weights) <= 0 or sum(weights) != 29:
        raise SystemExit("positive weights must sum to 29")
    sqrt5, lam = lambda_mod(arguments.prime)
    _context, derivative = residual_derivative(arguments.prime, lam, weights)
    values = equations(lam, weights, derivative)
    stdout, stderr, elapsed = singular_run(
        arguments.prime, lam, values, arguments.projective
    )
    print(
        json.dumps(
            {
                "weights": weights,
                "prime": arguments.prime,
                "sqrt5": sqrt5,
                "lambda": lam,
                "projective": arguments.projective,
                "result": stdout,
                "stderr": stderr,
                "seconds": elapsed,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
