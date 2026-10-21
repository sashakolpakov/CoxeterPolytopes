#!/usr/bin/env python3
"""Exact CC successive-basis preflight for WW partition (3^2,1^23)."""

from __future__ import annotations

import argparse
import json
import subprocess
import time

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_five_black_ww_partition5 import residual_derivative


def equations(lam: int, weights, derivative):
    a, b, c, d, e = weights
    low = [
        str(derivative[0]),
        str(derivative[1]),
        f"({derivative[2]})-29",
        f"({derivative[3]})+58",
    ]
    g0 = f"(-L)^{a}*(-u)^{b}*(-v)^{c}*(-w)^{d}*(-z)^{e}"
    g1 = f"(1-L)^{a}*(1-u)^{b}*(1-v)^{c}*(1-w)^{d}*(1-z)^{e}"
    return low + [f"({g0})-({g1})"]


def singular_run(prime: int, lam: int, values):
    program = "\n".join(
        [
            f"ring R={prime},(u,v,w,z),dp;",
            "option(redSB);",
            f"number L={lam};",
            *(f"poly f{i}={value};" for i, value in enumerate(values)),
            "ideal I=f0,f1,f2,f3; ideal B=std(I);",
            'print("BASE"); print(dim(B)); print(vdim(B));',
            "ideal J=std(B+ideal(f4));",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(dim(J)); print(vdim(J)); }',
            "quit;",
        ]
    )
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    return completed.stdout.strip(), completed.stderr.strip(), time.monotonic() - before


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("weights", nargs=5, type=int)
    parser.add_argument("--prime", type=int, default=101)
    args = parser.parse_args()
    weights = tuple(args.weights)
    sqrt5, lam = lambda_mod(args.prime)
    _context, derivative = residual_derivative(args.prime, lam, weights)
    values = equations(lam, weights, derivative)
    stdout, stderr, elapsed = singular_run(args.prime, lam, values)
    print(json.dumps({"weights": weights, "prime": args.prime, "sqrt5": sqrt5,
        "lambda": lam, "degrees": [4,3,2,1,29-weights[0]], "result": stdout,
        "stderr": stderr, "seconds": elapsed}, sort_keys=True))


if __name__ == "__main__":
    main()
